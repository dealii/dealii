//----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------------

#ifndef __deal2__mg_transfer_templates_h
#define __deal2__mg_transfer_templates_h

#include <lac/sparse_matrix.h>
#include <grid/tria_iterator.h>
#include <fe/fe.h>
#include <dofs/dof_constraints.h>
#include <multigrid/mg_base.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/mg_dof_tools.h>
#include <multigrid/mg_transfer.h>

#include <algorithm>


/* --------------------- MGTransferPrebuilt -------------- */



template <class VECTOR>
template <int dim, class InVector>
void
MGTransferPrebuilt<VECTOR>::copy_to_mg (
  const MGDoFHandler<dim>        &mg_dof_handler,
  MGLevelObject<VECTOR> &dst,
  const InVector                 &src) const
{
				   // forward to the correct
				   // specialization
  copy_to_mg (mg_dof_handler, dst, src, is_1d<(dim==1)>());
}


template <class VECTOR>
template <int dim, class InVector>
void
MGTransferPrebuilt<VECTOR>::copy_to_mg (
  const MGDoFHandler<dim>        &,
  MGLevelObject<VECTOR> &,
  const InVector                 &,
  const is_1d<true>              &) const
{
  Assert(false, ExcNotImplemented());
}


template <class VECTOR>
template <int dim, class InVector>
void
MGTransferPrebuilt<VECTOR>::copy_to_mg (
  const MGDoFHandler<dim>& mg_dof_handler,
  MGLevelObject<VECTOR>& dst,
  const InVector&          src,
  const is_1d<false>&) const
{
				   // Make src a real finite element function
//  InVector src = osrc;
//  constraints->distribute(src);
  const unsigned int dofs_per_cell = mg_dof_handler.get_fe().dofs_per_cell;
  const unsigned int dofs_per_face = mg_dof_handler.get_fe().dofs_per_face;

				   // set the elements of the vectors
				   // on all levels to zero
  unsigned int minlevel = dst.get_minlevel();
  unsigned int maxlevel = dst.get_maxlevel();

  MGTools::reinit_vector(mg_dof_handler, dst);

  Assert(sizes.size()==mg_dof_handler.get_tria().n_levels(),
	 ExcMatricesNotBuilt());
//  for (unsigned int l=minlevel;l<=maxlevel;++l)
//    dst[l].reinit(sizes[l]);
  
  std::vector<unsigned int> global_dof_indices (dofs_per_cell);
  std::vector<unsigned int> level_dof_indices  (dofs_per_cell);
  std::vector<unsigned int> level_face_indices (dofs_per_face);

				   // traverse the grid top-down
				   // (i.e. starting with the most
				   // refined grid). this way, we can
				   // always get that part of one
				   // level of the output vector which
				   // corresponds to a region which is
				   // more refined, by restriction of
				   // the respective vector on the
				   // next finer level, which we then
				   // already have built.
  for (int level=maxlevel; level>=static_cast<signed int>(minlevel); --level)
    {
      typename MGDoFHandler<dim>::active_cell_iterator
	level_cell = mg_dof_handler.begin_active(level);
      const typename MGDoFHandler<dim>::active_cell_iterator
	level_end  = mg_dof_handler.end_active(level);

//TODO:[?] Treat hanging nodes properly
// The single-level vector is not an FE-function, because the values at
// hanging nodes are set to zero. This should be treated before the restriction.

				       // Compute coarse level right hand side
				       // by restricting from fine level.
      for (; level_cell!=level_end; ++level_cell)
	{
	  DoFObjectAccessor<dim, dim>& global_cell = *level_cell;
					   // get the dof numbers of
					   // this cell for the global
					   // and the level-wise
					   // numbering
	  global_cell.get_dof_indices(global_dof_indices);
	  level_cell->get_mg_dof_indices (level_dof_indices);

					   // transfer the global
					   // defect in the vector
					   // into the level-wise one
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    dst[level](level_dof_indices[i]) = src(global_dof_indices[i]);

	  for (unsigned int face_n=0; face_n<GeometryInfo<dim>::faces_per_cell; ++face_n)
	    {
	      const typename MGDoFHandler<dim>::face_iterator
		face = level_cell->face(face_n);
	      if (face->has_children())
		{
		  face->get_mg_dof_indices(level_face_indices);


						   // Delete values on refinement edge,
						   // since restriction will add them again.
		  for (unsigned int i=0; i<dofs_per_face; ++i)
		    dst[level](level_face_indices[i])
		      = 0.;
		};
	    };
	};
				       // for that part of the level
				       // which is further refined:
				       // get the defect by
				       // restriction of the defect on
				       // one level higher
      if (static_cast<unsigned int>(level) < maxlevel)
	{
	  restrict_and_add (level+1, dst[level], dst[level+1]);
	}
    };
}



template <class VECTOR>
template <int dim, class OutVector>
void
MGTransferPrebuilt<VECTOR>::copy_from_mg(
  const MGDoFHandler<dim>&       mg_dof_handler,
  OutVector&                     dst,
  const MGLevelObject<VECTOR>& src) const
{
  const unsigned int dofs_per_cell = mg_dof_handler.get_fe().dofs_per_cell;

  std::vector<unsigned int> global_dof_indices (dofs_per_cell);
  std::vector<unsigned int> level_dof_indices (dofs_per_cell);

  typename MGDoFHandler<dim>::active_cell_iterator
    level_cell = mg_dof_handler.begin_active();
  const typename MGDoFHandler<dim>::active_cell_iterator
    endc = mg_dof_handler.end();

				   // traverse all cells and copy the
				   // data appropriately to the output
				   // vector

				   // Note that the level is
				   // monotonuosly increasing
  for (; level_cell != endc; ++level_cell)
    {
      DoFObjectAccessor<dim, dim>& global_cell = *level_cell;
      const unsigned int level = level_cell->level();
      
				       // get the dof numbers of
				       // this cell for the global
				       // and the level-wise
				       // numbering
      global_cell.get_dof_indices (global_dof_indices);
      level_cell->get_mg_dof_indices(level_dof_indices);

				       // copy level-wise data to
				       // global vector
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	dst(global_dof_indices[i]) = src[level](level_dof_indices[i]);
    };

				   // clear constrained nodes
//TODO:[GK]  constraints->set_zero(dst);
}



template <class VECTOR>
template <int dim, class OutVector>
void
MGTransferPrebuilt<VECTOR>::copy_from_mg_add (
  const MGDoFHandler<dim>              &mg_dof_handler,
  OutVector                            &dst,
  const MGLevelObject<VECTOR> &src) const
{
  const unsigned int dofs_per_cell = mg_dof_handler.get_fe().dofs_per_cell;

  std::vector<unsigned int> global_dof_indices (dofs_per_cell);
  std::vector<unsigned int> level_dof_indices (dofs_per_cell);

  typename MGDoFHandler<dim>::active_cell_iterator
    level_cell = mg_dof_handler.begin_active();
  const typename MGDoFHandler<dim>::active_cell_iterator
    endc = mg_dof_handler.end();

				   // traverse all cells and copy the
				   // data appropriately to the output
				   // vector

				   // Note that the level is
				   // monotonuosly increasing
  for (; level_cell != endc; ++level_cell)
    {
      DoFObjectAccessor<dim, dim>& global_cell = *level_cell;
      const unsigned int level = level_cell->level();
      
				       // get the dof numbers of
				       // this cell for the global
				       // and the level-wise
				       // numbering
      global_cell.get_dof_indices (global_dof_indices);
      level_cell->get_mg_dof_indices(level_dof_indices);

//TODO:[GK] Probably wrong for continuous elements

				       // copy level-wise data to
				       // global vector
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	dst(global_dof_indices[i]) += src[level](level_dof_indices[i]);
    };

				   // clear constrained nodes
//TODO:[GK]  constraints->set_zero(dst);
}


template <class VECTOR>
unsigned int
MGTransferPrebuilt<VECTOR>::memory_consumption () const
{
  unsigned int result = sizeof(*this);
  result += sizeof(unsigned int) * sizes.size();
#ifdef DEAL_PREFER_MATRIX_EZ
  std::vector<boost::shared_ptr<SparseMatrixEZ<double> > >::const_iterator m;
  const std::vector<boost::shared_ptr<SparseMatrixEZ<double> > >::const_iterator end = prolongation_matrices.end();
  for (m = prolongation_matrices.begin(); m != end ; ++m)
    result += *m->memory_consumption();
#else
  for (unsigned int i=0;i<prolongation_matrices.size();++i)
    result += prolongation_matrices[i]->memory_consumption()
	      + prolongation_sparsities[i]->memory_consumption();
#endif
  return result;
}

/* --------------------- MGTransferSelect -------------- */



template <typename number>
template <int dim, typename number2>
void
MGTransferSelect<number>::copy_to_mg (
  const MGDoFHandler<dim>        &mg_dof_handler,
  MGLevelObject<Vector<number> > &dst,
  const BlockVector<number2>     &src) const
{
				   // forward to the correct
				   // specialization
  do_copy_to_mg (mg_dof_handler, dst, src,
		 0, //component_start[selected_component],
		 is_1d<(dim==1)>());
}



template <typename number>
template <int dim, typename number2>
void
MGTransferSelect<number>::copy_to_mg (
  const MGDoFHandler<dim>        &mg_dof_handler,
  MGLevelObject<Vector<number> > &dst,
  const Vector<number2>          &src) const
{
				   // forward to the correct
				   // specialization
  do_copy_to_mg (mg_dof_handler, dst, src,
		 component_start[selected_component],
		 is_1d<(dim==1)>());
}



template <typename number>
template <int dim, class InVector>
void
MGTransferSelect<number>::do_copy_to_mg (
  const MGDoFHandler<dim>&,
  MGLevelObject<Vector<number> >&,
  const InVector&,
  const unsigned int,
  const is_1d<true>&) const
{
  Assert(false, ExcNotImplemented());
}



template <typename number>
template <int dim, class InVector>
void
MGTransferSelect<number>::do_copy_to_mg (
  const MGDoFHandler<dim>&        mg_dof_handler,
  MGLevelObject<Vector<number> >& dst,
  const InVector&                 src,
  const unsigned int              offset,
  const is_1d<false>&) const
{
				   // Make src a real finite element function
//  InVector src = osrc;
//  constraints->distribute(src);

  const FiniteElement<dim>& fe = mg_dof_handler.get_fe();
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int dofs_per_face = fe.dofs_per_face;

				   // set the elements of the vectors
				   // on all levels to zero
  unsigned int minlevel = dst.get_minlevel();
  unsigned int maxlevel = dst.get_maxlevel();
  
  dst.clear();
  
  Assert(sizes.size()==mg_dof_handler.get_tria().n_levels(),
	 ExcMatricesNotBuilt());

  MGTools::reinit_vector(mg_dof_handler, dst,
			 mg_selected, mg_target_component);
  
  std::vector<unsigned int> global_dof_indices (dofs_per_cell);
  std::vector<unsigned int> level_dof_indices  (dofs_per_cell);
  std::vector<unsigned int> level_face_indices (dofs_per_face);

				   // traverse the grid top-down
				   // (i.e. starting with the most
				   // refined grid). this way, we can
				   // always get that part of one
				   // level of the output vector which
				   // corresponds to a region which is
				   // more refined, by restriction of
				   // the respective vector on the
				   // next finer level, which we then
				   // already have built.
  for (int level=maxlevel; level>=static_cast<signed int>(minlevel); --level)
    {
      typename MGDoFHandler<dim>::active_cell_iterator
	level_cell = mg_dof_handler.begin_active(level);
      const typename MGDoFHandler<dim>::active_cell_iterator
	level_end  = mg_dof_handler.end_active(level);

//TODO:[?] Treat hanging nodes properly
// The single-level vector is not an FE-function, because the values at
// hanging nodes are set to zero. This should be treated before the restriction.

				       // Compute coarse level right hand side
				       // by restricting from fine level.
      for (; level_cell!=level_end; ++level_cell)
	{
	  DoFObjectAccessor<dim, dim>& global_cell = *level_cell;
					   // get the dof numbers of
					   // this cell for the global
					   // and the level-wise
					   // numbering
	  global_cell.get_dof_indices(global_dof_indices);
	  level_cell->get_mg_dof_indices (level_dof_indices);

					   // transfer the global
					   // defect in the vector
					   // into the level-wise one
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {
	      const unsigned int component
		= mg_target_component[fe.system_to_component_index(i).first];
	      if (mg_selected[component])
		{
		  const unsigned int level_start
		    = mg_component_start[level][component];
		  
		  dst[level](level_dof_indices[i] - level_start)
		    = src(global_dof_indices[i] - offset);
		}
	    }
	  
// 	  for (unsigned int face_n=0; face_n<GeometryInfo<dim>::faces_per_cell; ++face_n)
// 	    {
// 	      const typename MGDoFHandler<dim>::face_iterator
// 		face = level_cell->face(face_n);
// 	      if (face->has_children())
// 		{
// 		  face->get_mg_dof_indices(level_face_indices);


						   // Delete values on refinement edge,
						   // since restriction will add them again.
//		  for (unsigned int i=0; i<dofs_per_face; ++i)
//		    dst[level](level_face_indices[i])
//		      = 0.;
//		};
//	    };
	}
				       // for that part of the level
				       // which is further refined:
				       // get the defect by
				       // restriction of the defect on
				       // one level higher
      if (static_cast<unsigned int>(level) < maxlevel)
	{
	  restrict_and_add (level+1, dst[level], dst[level+1]);
	}
    };
}



template <typename number>
template <int dim, typename number2>
void
MGTransferSelect<number>::copy_from_mg (
  const MGDoFHandler<dim>&              mg_dof_handler,
  BlockVector<number2>&                 dst,
  const MGLevelObject<Vector<number> >& src) const
{
  do_copy_from_mg (mg_dof_handler, dst, src, 0);
}



template <typename number>
template <int dim, typename number2>
void
MGTransferSelect<number>::copy_from_mg (
  const MGDoFHandler<dim>&              mg_dof_handler,
  Vector<number2>&                      dst,
  const MGLevelObject<Vector<number> >& src) const
{
  do_copy_from_mg (mg_dof_handler, dst, src,
		   component_start[selected_component]);
}



template <typename number>
template <int dim, typename number2>
void
MGTransferSelect<number>::copy_from_mg_add (
  const MGDoFHandler<dim>&              mg_dof_handler,
  BlockVector<number2>&                 dst,
  const MGLevelObject<Vector<number> >& src) const
{
  do_copy_from_mg_add (mg_dof_handler, dst, src, 0);
}



template <typename number>
template <int dim, typename number2>
void
MGTransferSelect<number>::copy_from_mg_add (
  const MGDoFHandler<dim>&              mg_dof_handler,
  Vector<number2>&                      dst,
  const MGLevelObject<Vector<number> >& src) const
{
  do_copy_from_mg_add (mg_dof_handler, dst, src,
		       component_start[selected_component]);
}



template <typename number>
template <int dim, class OutVector>
void
MGTransferSelect<number>::do_copy_from_mg (
  const MGDoFHandler<dim>              &mg_dof_handler,
  OutVector                            &dst,
  const MGLevelObject<Vector<number> > &src,
  const unsigned int offset) const
{

  const FiniteElement<dim>& fe = mg_dof_handler.get_fe();
  const unsigned int dofs_per_cell = fe.dofs_per_cell;

  std::vector<unsigned int> global_dof_indices (dofs_per_cell);
  std::vector<unsigned int> level_dof_indices (dofs_per_cell);

  typename MGDoFHandler<dim>::active_cell_iterator
    level_cell = mg_dof_handler.begin_active();
  const typename MGDoFHandler<dim>::active_cell_iterator
    endc = mg_dof_handler.end();

				   // traverse all cells and copy the
				   // data appropriately to the output
				   // vector

				   // Note that the level is
				   // monotonuosly increasing
  for (; level_cell != endc; ++level_cell)
    {
       DoFObjectAccessor<dim, dim>& global_cell = *level_cell;
      const unsigned int level = level_cell->level();

				       // get the dof numbers of
				       // this cell for the global
				       // and the level-wise
				       // numbering
      global_cell.get_dof_indices (global_dof_indices);
      level_cell->get_mg_dof_indices(level_dof_indices);

				       // copy level-wise data to
				       // global vector
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  const unsigned int component
	    = mg_target_component[fe.system_to_component_index(i).first];
	if (mg_selected[component])
	  {
	    const unsigned int level_start
	      = mg_component_start[level][component];
	    dst(global_dof_indices[i] - offset)
	      = src[level](level_dof_indices[i]-level_start);
	  }
	}
    };

				   // clear constrained nodes
//TODO:[GK]  constraints->set_zero(dst);
}



template <typename number>
template <int dim, class OutVector>
void
MGTransferSelect<number>::do_copy_from_mg_add (
  const MGDoFHandler<dim>              &mg_dof_handler,
  OutVector                            &dst,
  const MGLevelObject<Vector<number> > &src,
  const unsigned int offset) const
{

  const FiniteElement<dim>& fe = mg_dof_handler.get_fe();
  const unsigned int dofs_per_cell = fe.dofs_per_cell;

  std::vector<unsigned int> global_dof_indices (dofs_per_cell);
  std::vector<unsigned int> level_dof_indices (dofs_per_cell);

  typename MGDoFHandler<dim>::active_cell_iterator
    level_cell = mg_dof_handler.begin_active();
  const typename MGDoFHandler<dim>::active_cell_iterator
    endc = mg_dof_handler.end();

				   // traverse all cells and copy the
				   // data appropriately to the output
				   // vector

				   // Note that the level is
				   // monotonuosly increasing
  for (; level_cell != endc; ++level_cell)
    {
       DoFObjectAccessor<dim, dim>& global_cell = *level_cell;
      const unsigned int level = level_cell->level();

				       // get the dof numbers of
				       // this cell for the global
				       // and the level-wise
				       // numbering
      global_cell.get_dof_indices (global_dof_indices);
      level_cell->get_mg_dof_indices(level_dof_indices);

//TODO:[GK+WB] Probably wrong for continuous elements

				       // copy level-wise data to
				       // global vector
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  const unsigned int component
	    = mg_target_component[fe.system_to_component_index(i).first];
	  if (mg_selected[component])
	    {
	      const unsigned int level_start
		= mg_component_start[level][component];
	      dst(global_dof_indices[i] - offset)
		+= src[level](level_dof_indices[i] - level_start);
	    }
	}
    }
				   // clear constrained nodes
//TODO:[GK+WB]  constraints->set_zero(dst);
}


template <typename number>
unsigned int
MGTransferSelect<number>::memory_consumption () const
{
  return sizeof(int) + MGTransferBlockBase::memory_consumption();
}


/* --------------------- MGTransferBlock -------------- */



template <typename number>
template <int dim, class InVector>
void
MGTransferBlock<number>::copy_to_mg (
  const MGDoFHandler<dim>             &mg_dof_handler,
  MGLevelObject<BlockVector<number> > &dst,
  const InVector                      &src) const
{
				   // forward to the correct
				   // specialization
  copy_to_mg (mg_dof_handler, dst, src, is_1d<(dim==1)>());
}



template <typename number>
template <int dim, class InVector>
void
MGTransferBlock<number>::copy_to_mg (
  const MGDoFHandler<dim>             &,
  MGLevelObject<BlockVector<number> > &,
  const InVector                      &,
  const is_1d<true>                   &) const
{
  Assert(false, ExcNotImplemented());
}



template <typename number>
template <int dim, class InVector>
void
MGTransferBlock<number>::copy_to_mg (
  const MGDoFHandler<dim>             &mg_dof_handler,
  MGLevelObject<BlockVector<number> > &dst,
  const InVector                      &src,
  const is_1d<false>                  &) const
{
				   // Make src a real finite element
				   // function
//  InVector src = osrc;
//  constraints->distribute(src);

  const FiniteElement<dim>& fe = mg_dof_handler.get_fe();
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int dofs_per_face = fe.dofs_per_face;

				   // set the elements of the vectors
				   // on all levels to zero
  unsigned int minlevel = dst.get_minlevel();
  unsigned int maxlevel = dst.get_maxlevel();
  
  dst.clear();

  MGTools::reinit_vector(mg_dof_handler, dst, mg_selected, mg_target_component);
  
  std::vector<unsigned int> global_dof_indices (dofs_per_cell);
  std::vector<unsigned int> level_dof_indices  (dofs_per_cell);
  std::vector<unsigned int> level_face_indices (dofs_per_face);

				   // traverse the grid top-down
				   // (i.e. starting with the most
				   // refined grid). this way, we can
				   // always get that part of one
				   // level of the output vector which
				   // corresponds to a region which is
				   // more refined, by restriction of
				   // the respective vector on the
				   // next finer level, which we then
				   // already have built.
  for (int level=maxlevel; level>=static_cast<signed int>(minlevel); --level)
    {
      typename MGDoFHandler<dim>::active_cell_iterator
	level_cell = mg_dof_handler.begin_active(level);
      const typename MGDoFHandler<dim>::active_cell_iterator
	level_end  = mg_dof_handler.end_active(level);

				       // Compute coarse level right hand side
				       // by restricting from fine level.
      for (; level_cell!=level_end; ++level_cell)
	{
	  DoFObjectAccessor<dim, dim>& global_cell = *level_cell;
					   // get the dof numbers of
					   // this cell for the global
					   // and the level-wise
					   // numbering
	  global_cell.get_dof_indices(global_dof_indices);
	  level_cell->get_mg_dof_indices (level_dof_indices);

					   // transfer the global
					   // defect in the vector
					   // into the level-wise one
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {
//TODO:[GK] Not sure if this works if components are selected
	      dst[level](level_dof_indices[i])
		= src(global_dof_indices[i]);
	    }
//TODO: Special treatment of faces? Copy from MGTransferPrebuilt when done there.
	}
				       // for that part of the level
				       // which is further refined:
				       // get the defect by
				       // restriction of the defect on
				       // one level higher
      if (static_cast<unsigned int>(level) < maxlevel)
	{
	  restrict_and_add (level+1, dst[level], dst[level+1]);
	}
    };
}




template <typename number>
template <int dim, class OutVector>
void
MGTransferBlock<number>::copy_from_mg (
  const MGDoFHandler<dim>                   &mg_dof_handler,
  OutVector                                 &dst,
  const MGLevelObject<BlockVector<number> > &src) const
{

  const FiniteElement<dim>& fe = mg_dof_handler.get_fe();
  const unsigned int dofs_per_cell = fe.dofs_per_cell;

  std::vector<unsigned int> global_dof_indices (dofs_per_cell);
  std::vector<unsigned int> level_dof_indices (dofs_per_cell);

  typename MGDoFHandler<dim>::active_cell_iterator
    level_cell = mg_dof_handler.begin_active();
  const typename MGDoFHandler<dim>::active_cell_iterator
    endc = mg_dof_handler.end();
				   // traverse all cells and copy the
				   // data appropriately to the output
				   // vector

				   // Note that level is monotonuosly
				   // increasing
  for (; level_cell != endc; ++level_cell)
    {
       DoFObjectAccessor<dim, dim>& global_cell = *level_cell;
      const unsigned int level = level_cell->level();
				       // get the dof numbers of
				       // this cell for the global
				       // and the level-wise
				       // numbering
      global_cell.get_dof_indices (global_dof_indices);
      level_cell->get_mg_dof_indices(level_dof_indices);

				       // copy level-wise data to
				       // global vector
//TODO:[GK] Does this work with selected components?
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	dst(global_dof_indices[i])
	  = src[level](level_dof_indices[i]);
    }
}



template <typename number>
template <int dim, class OutVector>
void
MGTransferBlock<number>::copy_from_mg_add (
  const MGDoFHandler<dim>                   &mg_dof_handler,
  OutVector                                 &dst,
  const MGLevelObject<BlockVector<number> > &src) const
{
  const FiniteElement<dim>& fe = mg_dof_handler.get_fe();
  const unsigned int dofs_per_cell = fe.dofs_per_cell;

  std::vector<unsigned int> global_dof_indices (dofs_per_cell);
  std::vector<unsigned int> level_dof_indices (dofs_per_cell);

  typename MGDoFHandler<dim>::active_cell_iterator
    level_cell = mg_dof_handler.begin_active();
  const typename MGDoFHandler<dim>::active_cell_iterator
    endc = mg_dof_handler.end();

				   // traverse all cells and copy the
				   // data appropriately to the output
				   // vector

				   // Note that the level monotonuosly
				   // increasing
  for (; level_cell != endc; ++level_cell)
    {
       DoFObjectAccessor<dim, dim>& global_cell = *level_cell;
      const unsigned int level = level_cell->level();

				       // get the dof numbers of
				       // this cell for the global
				       // and the level-wise
				       // numbering
      global_cell.get_dof_indices (global_dof_indices);
      level_cell->get_mg_dof_indices(level_dof_indices);

//TODO:[GK] Probably wrong for continuous elements

				       // copy level-wise data to
				       // global vector
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	dst(global_dof_indices[i])
	  += src[level](level_dof_indices[i]);
    }
}


#endif

//----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------------

#ifndef __deal2__mg_tools_templates_h
#define __deal2__mg_tools_templates_h

#include <dofs/dof_constraints.h>
#include <numerics/data_out.h>
#include <multigrid/multigrid.h>
#include <algorithm>
#include <fstream>

#include <lac/sparse_matrix.h>


//TODO:[?] This function needs to be specially implemented, since in 2d mode we use faces
#if deal_II_dimension == 1

template <int dim, class InVector>
void
MGTransferPrebuilt::copy_to_mg (const MGDoFHandler<dim>&,
				MGLevelObject<Vector<double> >&,
				const InVector&) const
{
  Assert(false, ExcNotImplemented());
}

#else


template <int dim, class InVector>
void
MGTransferPrebuilt::copy_to_mg (const MGDoFHandler<dim>& mg_dof_handler,
				MGLevelObject<Vector<double> >& dst,
				const InVector& src) const
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
  
  dst.clear();

  for (unsigned int l=minlevel;l<=maxlevel;++l)
    dst[l].reinit(sizes[l]);
  
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
  for (int level=maxlevel; level>=static_cast<int>(minlevel); --level)
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
};

#endif


template <int dim, class OutVector>
void
MGTransferPrebuilt::copy_from_mg(const MGDoFHandler<dim>& mg_dof_handler,
				 OutVector &dst,
				 const MGLevelObject<Vector<double> > &src) const
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

				   // Is the level monotonuosly increasing?

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
//TODO:  constraints->set_zero(dst);
}

template <int dim, class OutVector>
void
MGTransferPrebuilt::copy_from_mg_add(const MGDoFHandler<dim>& mg_dof_handler,
				     OutVector &dst,
				     const MGLevelObject<Vector<double> > &src) const
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

				   // Is the level monotonuosly increasing?

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

//TODO: Probably wrong for continuous elements

				       // copy level-wise data to
				       // global vector
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	dst(global_dof_indices[i]) += src[level](level_dof_indices[i]);
    };

				   // clear constrained nodes
//TODO:  constraints->set_zero(dst);
}


//TODO:[?] This function needs to be specially implemented, since in 2d mode we use faces
#if deal_II_dimension == 1

template <int dim, class InVector>
void
MGTransferSelect::copy_to_mg (const MGDoFHandler<dim>&,
				MGLevelObject<Vector<double> >&,
				const InVector&) const
{
  Assert(false, ExcNotImplemented());
}

#else


template <int dim, class InVector>
void
MGTransferSelect::copy_to_mg (const MGDoFHandler<dim>& mg_dof_handler,
				MGLevelObject<Vector<double> >& dst,
				const InVector& osrc) const
{
				   // Make src a real finite element function
  InVector src = osrc;
//  constraints->distribute(src);

  const FiniteElement<dim>& fe = mg_dof_handler.get_fe();
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int dofs_per_face = fe.dofs_per_face;

				   // set the elements of the vectors
				   // on all levels to zero
  unsigned int minlevel = dst.get_minlevel();
  unsigned int maxlevel = dst.get_maxlevel();
  
  dst.clear();
  
  for (unsigned int l=minlevel;l<=maxlevel;++l)
    dst[l].reinit(sizes[l][selected]);
  
  std::vector<unsigned int> global_dof_indices (dofs_per_cell);
  std::vector<unsigned int> level_dof_indices  (dofs_per_cell);
  std::vector<unsigned int> level_face_indices (dofs_per_face);

				   // Start of input vector inside a
				   // block vector.
  const unsigned int start = component_start[selected];

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
  for (int level=maxlevel; level>=static_cast<int>(minlevel); --level)
    {
				       // Start of treated component
				       // for this level
      const unsigned int level_start = mg_component_start[level][selected];

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
	      if (fe.system_to_component_index(i).first == selected)
		dst[level](level_dof_indices[i] - level_start)
		  = src(global_dof_indices[i] - start);
	    }
	  
	  for (unsigned int face_n=0; face_n<GeometryInfo<dim>::faces_per_cell; ++face_n)
	    {
	      const typename MGDoFHandler<dim>::face_iterator
		face = level_cell->face(face_n);
	      if (face->has_children())
		{
		  face->get_mg_dof_indices(level_face_indices);


						   // Delete values on refinement edge,
						   // since restriction will add them again.
//		  for (unsigned int i=0; i<dofs_per_face; ++i)
//		    dst[level](level_face_indices[i])
//		      = 0.;
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
};

#endif


template <int dim, class OutVector>
void
MGTransferSelect::copy_from_mg(const MGDoFHandler<dim>& mg_dof_handler,
				 OutVector &dst,
				 const MGLevelObject<Vector<double> > &src) const
{

  const FiniteElement<dim>& fe = mg_dof_handler.get_fe();
  const unsigned int dofs_per_cell = fe.dofs_per_cell;

  std::vector<unsigned int> global_dof_indices (dofs_per_cell);
  std::vector<unsigned int> level_dof_indices (dofs_per_cell);

  typename MGDoFHandler<dim>::active_cell_iterator
    level_cell = mg_dof_handler.begin_active();
  const typename MGDoFHandler<dim>::active_cell_iterator
    endc = mg_dof_handler.end();

				   // Start of input vector inside a
				   // block vector.
  const unsigned int start = component_start[selected];

				   // traverse all cells and copy the
				   // data appropriately to the output
				   // vector

				   // Is the level monotonuosly increasing?

  for (; level_cell != endc; ++level_cell)
    {
       DoFObjectAccessor<dim, dim>& global_cell = *level_cell;
      const unsigned int level = level_cell->level();

				       // Start of treated component
				       // for this level
      const unsigned int level_start = mg_component_start[level][selected];
      
				       // get the dof numbers of
				       // this cell for the global
				       // and the level-wise
				       // numbering
      global_cell.get_dof_indices (global_dof_indices);
      level_cell->get_mg_dof_indices(level_dof_indices);

				       // copy level-wise data to
				       // global vector
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	if (fe.system_to_component_index(i).first == selected)
	  dst(global_dof_indices[i]-start)
	    = src[level](level_dof_indices[i]-level_start);
    };

				   // clear constrained nodes
//TODO:  constraints->set_zero(dst);
}

template <int dim, class OutVector>
void
MGTransferSelect::copy_from_mg_add(const MGDoFHandler<dim>& mg_dof_handler,
				     OutVector &dst,
				     const MGLevelObject<Vector<double> > &src) const
{

  const FiniteElement<dim>& fe = mg_dof_handler.get_fe();
  const unsigned int dofs_per_cell = fe.dofs_per_cell;

  std::vector<unsigned int> global_dof_indices (dofs_per_cell);
  std::vector<unsigned int> level_dof_indices (dofs_per_cell);

				   // Start of input vector inside a
				   // block vector.
  const unsigned int start = component_start[selected];

  typename MGDoFHandler<dim>::active_cell_iterator
    level_cell = mg_dof_handler.begin_active();
  const typename MGDoFHandler<dim>::active_cell_iterator
    endc = mg_dof_handler.end();

				   // traverse all cells and copy the
				   // data appropriately to the output
				   // vector

				   // Is the level monotonuosly increasing?

  for (; level_cell != endc; ++level_cell)
    {
       DoFObjectAccessor<dim, dim>& global_cell = *level_cell;
      const unsigned int level = level_cell->level();

				       // Start of treated component
				       // for this level
      const unsigned int level_start = mg_component_start[level][selected];
      
				       // get the dof numbers of
				       // this cell for the global
				       // and the level-wise
				       // numbering
      global_cell.get_dof_indices (global_dof_indices);
      level_cell->get_mg_dof_indices(level_dof_indices);

//TODO: Probably wrong for continuous elements

				       // copy level-wise data to
				       // global vector
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	if (fe.system_to_component_index(i).first == selected)
	  dst(global_dof_indices[i]-start)
	    += src[level](level_dof_indices[i]-level_start);
    };

				   // clear constrained nodes
//TODO:  constraints->set_zero(dst);
}


#endif

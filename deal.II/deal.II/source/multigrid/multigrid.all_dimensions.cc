//----------------------------  multigrid.all_dimensions.cc  ---------------------------
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
//----------------------------  multigrid.all_dimensions.cc  ---------------------------



#include <lac/vector.h>
#include <lac/sparse_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <multigrid/mg_transfer.h>
#include <multigrid/mg_smoother.h>
#include <multigrid/mg_transfer.templates.h>
#include <multigrid/multigrid.templates.h>





#ifdef MG_DEBUG
template <>
void
Multigrid<Vector<double> >::print_vector (const unsigned int level,
					  const Vector<double> &v,
					  const char *name) const
{
   if (level!=maxlevel)
     return;
   const unsigned int dim=deal_II_dimension;

//TODO[GK]: How is this supposed to work? the .all_dimension.cc files are supposed to be exactly the same for all space dimensions -- if they aren't you get strange and inconsistent results if you link a program with the 1d, 2d, and 3d libraries at the same time, because you have multiple instances of the exact same function (same name, same template arguments), but they do different things. the linker can't know this, so it may call one or the other, and possible results certainly include crashes
   const DoFHandler<dim> *dof = mg_dof_handler;
  
   Vector<double> out_vector(dof->n_dofs());

   out_vector = -10000;
  
   const unsigned int dofs_per_cell = mg_dof_handler->get_fe().dofs_per_cell;

   std::vector<unsigned int> global_dof_indices (dofs_per_cell);
   std::vector<unsigned int> level_dof_indices (dofs_per_cell);

   DoFHandler<dim>::active_cell_iterator
     global_cell = dof->begin_active(level);
   MGDoFHandler<dim>::active_cell_iterator
     level_cell = mg_dof_handler->begin_active(level);
   const MGDoFHandler<dim>::cell_iterator
     endc = mg_dof_handler->end(level);

 				   // traverse all cells and copy the
 				   // data appropriately to the output
 				   // vector
   for (; level_cell != endc; ++level_cell, ++global_cell)
     {
       global_cell->get_dof_indices (global_dof_indices);
       level_cell->get_mg_dof_indices(level_dof_indices);

 				       // copy level-wise data to
 				       // global vector
       for (unsigned int i=0; i<dofs_per_cell; ++i)
 	out_vector(global_dof_indices[i])
 	  = v(level_dof_indices[i]);
     }

   std::ofstream out_file(name);
   DataOut<dim> out;
   out.attach_dof_handler(*dof);
   out.add_data_vector(out_vector, "v");
   out.build_patches(5);
   out.write_gnuplot(out_file);
}



template <class VECTOR>
void
Multigrid<VECTOR>::print_vector (const unsigned int level,
				 const VECTOR &v,
				 const char *name) const
{
  Assert(false, ExcNotImplemented());
}
#endif




template <class VECTOR>
MGTransferPrebuilt<VECTOR>::~MGTransferPrebuilt () 
{}


template <class VECTOR>
void MGTransferPrebuilt<VECTOR>::prolongate (
  const unsigned int to_level,
  VECTOR&            dst,
  const VECTOR&      src) const 
{
  Assert ((to_level >= 1) && (to_level<=prolongation_matrices.size()),
	  ExcIndexRange (to_level, 1, prolongation_matrices.size()+1));

  prolongation_matrices[to_level-1]->vmult (dst, src);
}


template <class VECTOR>
void MGTransferPrebuilt<VECTOR>::restrict_and_add (
  const unsigned int   from_level,
  VECTOR       &dst,
  const VECTOR &src) const 
{
  Assert ((from_level >= 1) && (from_level<=prolongation_matrices.size()),
	  ExcIndexRange (from_level, 1, prolongation_matrices.size()+1));

  prolongation_matrices[from_level-1]->Tvmult_add (dst, src);
}





template <typename number>
MGTransferBlock<number>::~MGTransferBlock () 
{}


template <typename number>
void MGTransferBlock<number>::prolongate (
  const unsigned int   to_level,
  BlockVector<number>       &dst,
  const BlockVector<number> &src) const 
{
  Assert ((to_level >= 1) && (to_level<=prolongation_matrices.size()),
	  ExcIndexRange (to_level, 1, prolongation_matrices.size()+1));

  unsigned int k=0;
  for (unsigned int b=0; b<src.n_blocks();++b)
    {
      if (!selected[k])
	++k;
      prolongation_matrices[to_level-1]->block(k,k).vmult (dst.block(b), src.block(b));
      ++k;
    }
}


template <typename number>
void MGTransferBlock<number>::restrict_and_add (
  const unsigned int   from_level,
  BlockVector<number>       &dst,
  const BlockVector<number> &src) const 
{
  Assert ((from_level >= 1) && (from_level<=prolongation_matrices.size()),
	  ExcIndexRange (from_level, 1, prolongation_matrices.size()+1));

  unsigned int k=0;
  for (unsigned int b=0; b<src.n_blocks();++b)
    {
      if (!selected[k])
	++k;
      prolongation_matrices[from_level-1]->block(k,k).Tvmult_add (dst.block(b), src.block(b));
      ++k;
    }
}




template <typename number>
MGTransferSelect<number>::~MGTransferSelect () 
{}


template <typename number>
void MGTransferSelect<number>::prolongate (
  const unsigned int   to_level,
  Vector<number>       &dst,
  const Vector<number> &src) const 
{
  Assert ((to_level >= 1) && (to_level<=prolongation_matrices.size()),
	  ExcIndexRange (to_level, 1, prolongation_matrices.size()+1));

      prolongation_matrices[to_level-1]->block(mg_selected_component,
					       mg_selected_component)
	.vmult (dst, src);
}


template <typename number>
void MGTransferSelect<number>::restrict_and_add (
  const unsigned int   from_level,
  Vector<number>       &dst,
  const Vector<number> &src) const
{
  Assert ((from_level >= 1) && (from_level<=prolongation_matrices.size()),
	  ExcIndexRange (from_level, 1, prolongation_matrices.size()+1));

  prolongation_matrices[from_level-1]->block(mg_selected_component,
					     mg_selected_component)
    .Tvmult_add (dst, src);
}


template <class VECTOR>
void
MGSmootherContinuous<VECTOR>::set_zero_interior_boundary (
  const unsigned int level,
  VECTOR&            u) const
{
  if (level==0)
    return;
  else
    for (std::vector<unsigned int>::const_iterator p=interior_boundary_dofs[level-1].begin();
	 p!=interior_boundary_dofs[level-1].end(); ++p)
      u(*p) = 0;
}


// Explicit instantiations

template class Multigrid<Vector<float> >;
template class Multigrid<Vector<double> >;
template class Multigrid<BlockVector<float> >;
template class Multigrid<BlockVector<double> >;

template class MGTransferPrebuilt<Vector<float> >;
template class MGTransferPrebuilt<BlockVector<float> >;
template class MGTransferPrebuilt<Vector<double> >;
template class MGTransferPrebuilt<BlockVector<double> >;
template class MGTransferBlock<float>;
template class MGTransferBlock<double>;
template class MGTransferSelect<float>;
template class MGTransferSelect<double>;

template
void MGSmootherContinuous<Vector<double> >::set_zero_interior_boundary (
  const unsigned int, Vector<double>&) const;
template
void MGSmootherContinuous<Vector<float> >::set_zero_interior_boundary (
  const unsigned int, Vector<float>&) const;
template
void MGSmootherContinuous<BlockVector<double> >::set_zero_interior_boundary (
  const unsigned int, BlockVector<double>&) const;
template
void MGSmootherContinuous<BlockVector<float> >::set_zero_interior_boundary (
  const unsigned int, BlockVector<float>&) const;

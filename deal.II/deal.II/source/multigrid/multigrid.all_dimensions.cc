//----------------------------  multigrid.all_dimensions.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
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




template <typename number>
MGTransferPrebuilt<number>::~MGTransferPrebuilt () 
{}


template <typename number>
void MGTransferPrebuilt<number>::prolongate (
  const unsigned int   to_level,
  Vector<number>       &dst,
  const Vector<number> &src) const 
{
  Assert ((to_level >= 1) && (to_level<=prolongation_matrices.size()),
	  ExcIndexRange (to_level, 1, prolongation_matrices.size()+1));

  prolongation_matrices[to_level-1]->vmult (dst, src);
}


template <typename number>
void MGTransferPrebuilt<number>::restrict_and_add (
  const unsigned int   from_level,
  Vector<number>       &dst,
  const Vector<number> &src) const 
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

      prolongation_matrices[to_level-1]->block(selected, selected)
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

  prolongation_matrices[from_level-1]->block(selected, selected)
    .Tvmult_add (dst, src);
}


// Explicit instantiations

template class Multigrid<Vector<float> >;
template class Multigrid<Vector<double> >;
template class Multigrid<BlockVector<float> >;
template class Multigrid<BlockVector<double> >;

template class MGTransferPrebuilt<float>;
template class MGTransferPrebuilt<double>;
template class MGTransferBlock<float>;
template class MGTransferBlock<double>;
template class MGTransferSelect<float>;
template class MGTransferSelect<double>;

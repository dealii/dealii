//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------



#include <lac/vector.h>
#include <lac/sparse_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <multigrid/mg_transfer.h>
#include <multigrid/mg_transfer_block.h>
#include <multigrid/mg_transfer_component.h>
#include <multigrid/mg_smoother.h>
#include <multigrid/mg_transfer.templates.h>
#include <multigrid/mg_transfer_block.templates.h>
#include <multigrid/mg_transfer_component.templates.h>
#include <multigrid/multigrid.templates.h>

DEAL_II_NAMESPACE_OPEN


// Warning: the following function is for debugging purposes only. It
// will be compiled only, if the additional and undocumented compiler
// flag MG_DEBUG is set. Furthermore, as soon as this function is
// compiled, libraries for different dimensions may NOT be linked
// anymore at the same time.

// If this function is to be used, its declaration must be added to
// the class Multigrid again.

// TODO: deal_II_dimension not set any more, so this will not compile!
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


template<class VECTOR>
MGTransferPrebuilt<VECTOR>::MGTransferPrebuilt ()
{} 


template<class VECTOR>
MGTransferPrebuilt<VECTOR>::MGTransferPrebuilt (const ConstraintMatrix &c, const MGConstrainedDoFs& mg_c)
 :
   constraints(&c),
   mg_constrained_dofs(&mg_c)
{} 

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
MGTransferBlock<number>::MGTransferBlock ()
		:
		memory(0, typeid(*this).name())
{}


template <typename number>
MGTransferBlock<number>::~MGTransferBlock ()
{
  if (memory != 0) memory = 0;
}


template <typename number>
void
MGTransferBlock<number>::initialize (const std::vector<number>& f,
				     VectorMemory<Vector<number> >& mem) 
{
  factors = f;
  memory = &mem;
}


template <typename number>
void MGTransferBlock<number>::prolongate (
  const unsigned int   to_level,
  BlockVector<number>       &dst,
  const BlockVector<number> &src) const 
{
  Assert ((to_level >= 1) && (to_level<=prolongation_matrices.size()),
	  ExcIndexRange (to_level, 1, prolongation_matrices.size()+1));
  Assert (src.n_blocks() == this->n_mg_blocks,
	  ExcDimensionMismatch(src.n_blocks(), this->n_mg_blocks));
  Assert (dst.n_blocks() == this->n_mg_blocks,
	  ExcDimensionMismatch(dst.n_blocks(), this->n_mg_blocks));

				   // Multiplicate with prolongation
				   // matrix, but only those blocks
				   // selected.
  for (unsigned int b=0; b<this->mg_block.size();++b)
    {
      if (this->selected[b])
	prolongation_matrices[to_level-1]->block(b,b).vmult (
	  dst.block(this->mg_block[b]), src.block(this->mg_block[b]));
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
  Assert (src.n_blocks() == this->n_mg_blocks,
	  ExcDimensionMismatch(src.n_blocks(), this->n_mg_blocks));
  Assert (dst.n_blocks() == this->n_mg_blocks,
	  ExcDimensionMismatch(dst.n_blocks(), this->n_mg_blocks));
  
  for (unsigned int b=0; b<this->mg_block.size();++b)
    {
      if (this->selected[b])
	{
	  if (factors.size() != 0)
	    {
	      Assert (memory != 0, ExcNotInitialized());
	      Vector<number>* aux = memory->alloc();
	      aux->reinit(dst.block(this->mg_block[b]));
	      prolongation_matrices[from_level-1]->block(b,b).Tvmult (
		*aux, src.block(this->mg_block[b]));
	      
	      dst.block(this->mg_block[b]).add(factors[b], *aux);
	      memory->free(aux);
	    }
	  else
	    {
	      prolongation_matrices[from_level-1]->block(b,b).Tvmult_add (
		dst.block(this->mg_block[b]), src.block(this->mg_block[b]));
	    }
	}
    }
}



std::size_t
MGTransferComponentBase::memory_consumption () const
{
  std::size_t result = sizeof(*this);
  result += MemoryConsumption::memory_consumption(selected)
	    - sizeof(selected);
  result += MemoryConsumption::memory_consumption(target_component)
	    - sizeof(mg_target_component);
  result += MemoryConsumption::memory_consumption(sizes)
	    - sizeof(sizes);
  result += MemoryConsumption::memory_consumption(component_start)
	    - sizeof(component_start);
  result += MemoryConsumption::memory_consumption(mg_component_start)
	    - sizeof(mg_component_start);
  result += MemoryConsumption::memory_consumption(prolongation_sparsities)
	    - sizeof(prolongation_sparsities);
  result += MemoryConsumption::memory_consumption(prolongation_matrices)
	    - sizeof(prolongation_matrices);
//TODO:[GK] Add this.
//   result += MemoryConsumption::memory_consumption(copy_to_and_from_indices)
// 	    - sizeof(copy_to_and_from_indices);
  return result;
}


//TODO:[GK] Add all those little vectors.
std::size_t
MGTransferBlockBase::memory_consumption () const
{
  std::size_t result = sizeof(*this);
  result += sizeof(unsigned int) * sizes.size();
  result += MemoryConsumption::memory_consumption(selected)
	    - sizeof(selected);
  result += MemoryConsumption::memory_consumption(mg_block)
	    - sizeof(mg_block);
  result += MemoryConsumption::memory_consumption(block_start)
	    - sizeof(block_start);
  result += MemoryConsumption::memory_consumption(mg_block_start)
	    - sizeof(mg_block_start);
  result += MemoryConsumption::memory_consumption(prolongation_sparsities)
	    - sizeof(prolongation_sparsities);
  result += MemoryConsumption::memory_consumption(prolongation_matrices)
	    - sizeof(prolongation_matrices);
//TODO:[GK] Add this.
//   result += MemoryConsumption::memory_consumption(copy_indices)
// 	    - sizeof(copy_indices);
  return result;
}


//----------------------------------------------------------------------//

template<typename number>
MGTransferSelect<number>::MGTransferSelect ()
{} 


template<typename number>
MGTransferSelect<number>::MGTransferSelect (const ConstraintMatrix &c)
 :
   constraints(&c)
{} 

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

      prolongation_matrices[to_level-1]->block(mg_target_component[mg_selected_component],
					       mg_target_component[mg_selected_component])
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

  prolongation_matrices[from_level-1]->block(mg_target_component[mg_selected_component],
					     mg_target_component[mg_selected_component])
    .Tvmult_add (dst, src);
}


//----------------------------------------------------------------------//

template <typename number>
MGTransferBlockSelect<number>::~MGTransferBlockSelect () 
{}


template <typename number>
void MGTransferBlockSelect<number>::prolongate (
  const unsigned int   to_level,
  Vector<number>       &dst,
  const Vector<number> &src) const 
{
  Assert ((to_level >= 1) && (to_level<=prolongation_matrices.size()),
	  ExcIndexRange (to_level, 1, prolongation_matrices.size()+1));

      prolongation_matrices[to_level-1]->block(selected_block,
					       selected_block)
	.vmult (dst, src);
}


template <typename number>
void MGTransferBlockSelect<number>::restrict_and_add (
  const unsigned int   from_level,
  Vector<number>       &dst,
  const Vector<number> &src) const
{
  Assert ((from_level >= 1) && (from_level<=prolongation_matrices.size()),
	  ExcIndexRange (from_level, 1, prolongation_matrices.size()+1));

  prolongation_matrices[from_level-1]->block(selected_block,
					     selected_block)
    .Tvmult_add (dst, src);
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
template class MGTransferBlockSelect<float>;
template class MGTransferBlockSelect<double>;


DEAL_II_NAMESPACE_CLOSE

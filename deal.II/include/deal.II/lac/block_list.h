//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__block_list_h
#define __deal2__block_list_h


#include <base/subscriptor.h>
#include <base/template_constraints.h>
#include <fe/fe.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * @warning This class is still experimental and will most likely be
 * changed in a future release.
 *
 * A vector of index sets listing the indices of small blocks of a
 * linear system. For each block, the indices in that block are
 * listed.
 *
 * The focus of this class is on small blocks of degrees of freedom
 * associated with a single mesh cell or a small patch. These indices
 * may be contiguous or not and we do not optimize for the case they
 * are. For larger sets, the use of a vector of IndexSet objects might
 * be advisable.
 *
 * BlockList objects can conveniently be initialized with iterator
 * ranges from DoFHandler or MGDoFHandler, using either their cell or
 * multigrid indices. Other initializations are possible to be
 * implemented.
 *
 * @author Guido Kanschat
 * @date 2010
 */
class BlockList :
  public Subscriptor
{
  public:
				     /// The container for each index set
    typedef std::vector<unsigned int> block_container;
				     /// The iterator for individual indices
    typedef block_container::const_iterator const_iterator;

				     /**
				      * Add the indices in
				      * <tt>indices</tt> to block
				      * <tt>block</tt>, eliminating
				      * repeated indices.
				      */
    void add(unsigned int block, const std::vector<unsigned int>& indices);
    
				     /**
				      * Add the indices in
				      * <tt>indices</tt> to block
				      * <tt>block</tt>, eliminating
				      * repeated indices. Only add
				      * those indices for which
				      * <tt>selected_indices</tt> is true.
				      */
    void add(unsigned int block,
	     const std::vector<unsigned int>& indices,
	     const std::vector<bool>& selected_indices);
    
				     /**
				      * Set up all index sets using an
				      * DoF iterator range. This
				      * function will call
				      * <tt>begin->get_dof_indices()</tt>
				      * with a signature like
				      * DoFCellAccessor::get_dof_indices().
				      * Typically, the iterators will
				      * loop over active cells of a
				      * triangulation.
				      *
				      * In addition, the function
				      * needs the total number of
				      * blocks as its first argument.
				      */
    template <typename ITERATOR>
    void initialize(unsigned int n_blocks,
		    const ITERATOR begin,
		    const typename identity<ITERATOR>::type end);
    
				     /**
				      * @deprecated This function will
				      * move to DoFTools.
				      *
				      * Set up all index sets using an
				      * DoF iterator range. This
				      * function will call
				      * <tt>begin->get_mg_dof_indices()</tt>
				      * with a signature like
				      * MGDoFCellAccessor::get_mg_dof_indices().
				      * Typically, the iterators loop
				      * over the cells of a single
				      * level or a Triangulation.
				      *
				      * In addition, the function
				      * needs the total number of
				      * blocks as its first argument.
				      */
    template <typename ITERATOR>
    void initialize_mg(unsigned int n_blocks,
		       const ITERATOR begin,
		       const typename identity<ITERATOR>::type end);
    
				     /**
				      * @deprecated This function will
				      * move to DoFTools.
				      *
				      * Set up all index sets using an
				      * DoF iterator range. This
				      * function will call
				      * <tt>begin->get_dof_indices()</tt>
				      * with a signature like
				      * DoFCellAccessor::get_dof_indices().
				      * Typically, the iterators will
				      * loop over active cells of a
				      * triangulation.
				      *
				      * The argument vector
				      * <tt>selected_dofs</tt> should
				      * have the length  of dofs per
				      * cell (thus, this function is
				      * not suitable for hp), and a
				      * true value for each degree of
				      * freedom which should be added
				      * to the index set of this
				      * cell. If you are working on a
				      * single block of a block
				      * system, the <tt>offset</tt> is
				      * the start index of this block.
				      *
				      * In addition, the function
				      * needs the total number of
				      * blocks as its first argument.
				      */
    template <typename ITERATOR>
    void initialize(unsigned int n_blocks,
		    const ITERATOR begin,
		    const typename identity<ITERATOR>::type end,
		    const std::vector<bool>& selected_dofs,
		    unsigned int offset = 0);
				     /**
				      * @deprecated This function will
				      * move to DoFTools.
				      *
				      * Set up all index sets using an
				      * DoF iterator range. This
				      * function will call
				      * <tt>begin->get_mg_dof_indices()</tt>
				      * with a signature like
				      * MGDoFCellAccessor::get_mg_dof_indices().
				      * Typically, the iterators will
				      * loop over cells on a single
				      * level of a triangulation.
				      *
				      * The argument vector
				      * <tt>selected_dofs</tt> should
				      * have the length  of dofs per
				      * cell (thus, this function is
				      * not suitable for hp), and a
				      * true value for each degree of
				      * freedom which should be added
				      * to the index set of this
				      * cell. If you are working on a
				      * single block of a block
				      * system, the <tt>offset</tt> is
				      * the start index of this block.
				      *
				      * In addition, the function
				      * needs the total number of
				      * blocks as its first argument.
				      */
    template <typename ITERATOR>
    void initialize_mg(unsigned int n_blocks,
		       const ITERATOR begin,
		       const typename identity<ITERATOR>::type end,
		       const std::vector<bool>& selected_dofs,
		       unsigned int offset = 0);

				     /**
				      * @deprecated This function will
				      * move to DoFTools.
				      *
				      * Same as initialize_mg(), but
				      * instead of gathering the
				      * degrees of freedom of a single
				      * cell into a block, gather all
				      * degrees of freedom of a patch
				      * around a vertex.
				      */
    template <int dim, typename ITERATOR>
    void initialize_vertex_patches_mg(unsigned int n_blocks,
		       const ITERATOR begin,
		       const typename identity<ITERATOR>::type end,
		       const std::vector<bool>& selected_dofs = std::vector<bool>(),
		       unsigned int offset = 0);

				     /**
				      * @deprecated This function will
				      * move to DoFTools.
				      *
				      * Auxiliary function, counting
				      * the patches around vertices.
				      */
    template <int dim, typename ITERATOR>
    unsigned int count_vertex_patches(
      const ITERATOR begin,
      const typename identity<ITERATOR>::type end,
      bool same_level_only) const;

				     /**
				      * @deprecated This function will
				      * move to DoFTools.
				      *
				      */
    template <int dim, typename ITERATOR>
    bool cell_generates_vertex_patch(const ITERATOR cell,
				     bool same_level_only) const;
    
				     /**
				      * The number of blocks.
				      */
    unsigned int size() const;
    
				     /**
				      * The size of a single block.
				      */
    unsigned int block_size(unsigned int block) const;

				     /**
				      * Iterator to the first index in block.
				      */
    const_iterator begin(unsigned int block) const;
				     /**
				      * End iterator for a single block.
				      */
    const_iterator end(unsigned int block) const;
				     /**
				      * Return the position of
				      * <tt>index</tt> in
				      * <tt>block</tt>, or
				      * numbers::invalid_unsigned_int,
				      * if the index is not in the block.
				      */
    unsigned int local_index(unsigned int block, unsigned int index) const;
    
  private:
				     /**
				      * The container for t he index sets.
				      */
    std::vector<block_container> index_sets;
};


inline
void
BlockList::add(const unsigned int block, const std::vector<unsigned int>& indices)
{
  AssertIndexRange(block, index_sets.size());
  
  for (unsigned int i=0;i<indices.size();++i)
    {
      const unsigned int k = indices[i];
      if (k==numbers::invalid_unsigned_int)
	continue;
      if (std::find(index_sets[block].begin(), index_sets[block].end(), k)
	  == index_sets[block].end())
	index_sets[block].push_back(k);
    }
}


inline
void
BlockList::add(
  const unsigned int block,
  const std::vector<unsigned int>& indices,
  const std::vector<bool>& selected)
{
  AssertIndexRange(block, index_sets.size());
  AssertDimension(indices.size(), selected.size());
  
  for (unsigned int i=0;i<indices.size();++i)
    {
      const unsigned int k = indices[i];
      if (k==numbers::invalid_unsigned_int)
	continue;
      if (selected[i] && std::find(index_sets[block].begin(), index_sets[block].end(), k)
	  == index_sets[block].end())
	index_sets[block].push_back(k);
    }
}


template <typename ITERATOR>
inline
void
BlockList::initialize(unsigned int n_blocks, const ITERATOR begin, const typename identity<ITERATOR>::type end)
{
  index_sets.resize(n_blocks);
  std::vector<unsigned int> indices;
  unsigned int k = 0;
  for (ITERATOR cell = begin; cell != end; ++cell, ++k)
    {
      indices.resize(cell->get_fe().dofs_per_cell);
      cell->get_dof_indices(indices);
      add(k, indices);
    }
}


template <typename ITERATOR>
inline
void
BlockList::initialize_mg(unsigned int n_blocks, const ITERATOR begin, const typename identity<ITERATOR>::type end)
{
  index_sets.resize(n_blocks);
  std::vector<unsigned int> indices;
  unsigned int k = 0;
  for (ITERATOR cell = begin; cell != end; ++cell, ++k)
    {
      indices.resize(cell->get_fe().dofs_per_cell);
      cell->get_mg_dof_indices(indices);
      add(k, indices);
    }
}


template <typename ITERATOR>
inline
void
BlockList::initialize(
  unsigned int n_blocks,
  const ITERATOR begin,
  const typename identity<ITERATOR>::type end,
  const std::vector<bool>& selected_dofs,
  unsigned int offset)
{
  index_sets.resize(n_blocks);
  std::vector<unsigned int> indices;
  unsigned int k = 0;
  for (ITERATOR cell = begin; cell != end; ++cell, ++k)
    {
      indices.resize(cell->get_fe().dofs_per_cell);
      AssertDimension(selected_dofs.size(), cell->get_fe().dofs_per_cell);
      
      cell->get_dof_indices(indices);
      add(k, indices, selected_dofs);
    }
}


template <typename ITERATOR>
inline
void
BlockList::initialize_mg(
  unsigned int n_blocks,
  const ITERATOR begin,
  const typename identity<ITERATOR>::type end,
  const std::vector<bool>& selected_dofs,
  unsigned int offset)
{
  index_sets.resize(n_blocks);
  std::vector<unsigned int> indices;
  unsigned int k = 0;
  for (ITERATOR cell = begin; cell != end; ++cell, ++k)
    {
      indices.resize(cell->get_fe().dofs_per_cell);
      AssertDimension(selected_dofs.size(), cell->get_fe().dofs_per_cell);
      
      cell->get_mg_dof_indices(indices);
      add(k, indices, selected_dofs);
    }
}



template <int dim, typename ITERATOR>
inline
bool
BlockList::cell_generates_vertex_patch(
  const ITERATOR cell,
  bool same_level_only) const
{
  switch(dim)
    {
      case 3:
	    if (cell->at_boundary(4)) break;
	    if (same_level_only && cell->neighbor(4)->level() != cell->level()) break;
	    if (cell->neighbor(4)->at_boundary(0)) break;
	    if (same_level_only && cell->neighbor(4)->neighbor(0)->level() != cell->level()) break;
	    if (cell->neighbor(4)->at_boundary(2)) break;
	    if (same_level_only && cell->neighbor(4)->neighbor(2)->level() != cell->level()) break;
	    if (cell->neighbor(4)->neighbor(0)->at_boundary(2)) break;
	    if (same_level_only && cell->neighbor(4)->neighbor(0)->neighbor(2)->level() != cell->level()) break;
					     // No break here
      case 2:
	    if (cell->at_boundary(2)) break;
	    if (same_level_only && cell->neighbor(2)->level() != cell->level()) break;
	    if (cell->neighbor(2)->at_boundary(0)) break;
	    if (same_level_only && cell->neighbor(2)->neighbor(0)->level() != cell->level()) break;
      case 1:
	    if (cell->at_boundary(0)) break;
	    if (same_level_only && cell->neighbor(0)->level() != cell->level()) break;
	    return true;
	    break;
      default:
	    Assert(false, ExcNotImplemented());
	    break;
    }
  return false;
}


template <int dim, typename ITERATOR>
inline
unsigned int
BlockList::count_vertex_patches(
  const ITERATOR begin,
  const typename identity<ITERATOR>::type end,
  bool same_level_only) const
{
  unsigned int count = 0;
  for (ITERATOR cell = begin; cell != end; ++cell)
    if (cell_generates_vertex_patch<dim>(cell, same_level_only))
      ++count;
  return count;
}


template <int dim, typename ITERATOR>
inline
void
BlockList::initialize_vertex_patches_mg(
  unsigned int n_blocks,
  const ITERATOR begin,
  const typename identity<ITERATOR>::type end,
  const std::vector<bool>& selected_dofs,
  unsigned int offset)
{
  Assert(selected_dofs.size() == 0, ExcNotImplemented());
  Assert(offset==0, ExcNotImplemented());
  const FiniteElement<dim>& fe = begin->get_fe();
  Assert(fe.dofs_per_vertex == 0, ExcNotImplemented());
  
  index_sets.resize(n_blocks);
  std::vector<unsigned int> indices;
  unsigned int k = 0;
  for (ITERATOR cell = begin; cell != end; ++cell)
    {
      if (cell_generates_vertex_patch<dim>(cell, true))
	{
	  indices.resize(cell->get_fe().dofs_per_cell);

	  switch(dim)
	    {
	      case 3:
		    cell->neighbor(4)->get_mg_dof_indices(indices);
		    for (unsigned int i=0;i<fe.dofs_per_face;++i)
		      {
			indices[fe.face_to_cell_index(i,1)] = numbers::invalid_unsigned_int;
			indices[fe.face_to_cell_index(i,3)] = numbers::invalid_unsigned_int;			
			indices[fe.face_to_cell_index(i,4)] = numbers::invalid_unsigned_int;			
		      }		    
		    add(k, indices);
		    cell->neighbor(4)->neighbor(0)->get_mg_dof_indices(indices);
		    for (unsigned int i=0;i<fe.dofs_per_face;++i)
		      {
			indices[fe.face_to_cell_index(i,0)] = numbers::invalid_unsigned_int;
			indices[fe.face_to_cell_index(i,3)] = numbers::invalid_unsigned_int;			
			indices[fe.face_to_cell_index(i,4)] = numbers::invalid_unsigned_int;			
		      }
		    add(k, indices);
		    cell->neighbor(4)->neighbor(2)->get_mg_dof_indices(indices);
		    for (unsigned int i=0;i<fe.dofs_per_face;++i)
		      {
			indices[fe.face_to_cell_index(i,1)] = numbers::invalid_unsigned_int;
			indices[fe.face_to_cell_index(i,2)] = numbers::invalid_unsigned_int;			
			indices[fe.face_to_cell_index(i,4)] = numbers::invalid_unsigned_int;			
		      }
		    add(k, indices);
		    cell->neighbor(4)->neighbor(2)->neighbor(0)->get_mg_dof_indices(indices);
		    for (unsigned int i=0;i<fe.dofs_per_face;++i)
		      {
			indices[fe.face_to_cell_index(i,0)] = numbers::invalid_unsigned_int;
			indices[fe.face_to_cell_index(i,2)] = numbers::invalid_unsigned_int;			
			indices[fe.face_to_cell_index(i,4)] = numbers::invalid_unsigned_int;			
		      }
		    add(k, indices);
	      case 2:
		    cell->neighbor(2)->get_mg_dof_indices(indices);
		    for (unsigned int i=0;i<fe.dofs_per_face;++i)
		      {
			indices[fe.face_to_cell_index(i,1)] = numbers::invalid_unsigned_int;
			indices[fe.face_to_cell_index(i,2)] = numbers::invalid_unsigned_int;			
			if (dim>2)
			  indices[fe.face_to_cell_index(i,5)] = numbers::invalid_unsigned_int;			
		      }
		    add(k, indices);
		    cell->neighbor(2)->neighbor(0)->get_mg_dof_indices(indices);
		    for (unsigned int i=0;i<fe.dofs_per_face;++i)
		      {
			indices[fe.face_to_cell_index(i,0)] = numbers::invalid_unsigned_int;
			indices[fe.face_to_cell_index(i,2)] = numbers::invalid_unsigned_int;			
			if (dim>2)
			  indices[fe.face_to_cell_index(i,5)] = numbers::invalid_unsigned_int;			
		      }
		    add(k, indices);
						     // no break here
	      case 1:
		    cell->get_mg_dof_indices(indices);
		    for (unsigned int i=0;i<fe.dofs_per_face;++i)
		      {
			indices[fe.face_to_cell_index(i,1)] = numbers::invalid_unsigned_int;
			if (dim>1)
			  indices[fe.face_to_cell_index(i,3)] = numbers::invalid_unsigned_int;			
			if (dim>2)
			  indices[fe.face_to_cell_index(i,5)] = numbers::invalid_unsigned_int;			
		      }
		    add(k, indices);
		    cell->neighbor(0)->get_mg_dof_indices(indices);
		    for (unsigned int i=0;i<fe.dofs_per_face;++i)
		      {
			indices[fe.face_to_cell_index(i,0)] = numbers::invalid_unsigned_int;
			if (dim>1)
			  indices[fe.face_to_cell_index(i,3)] = numbers::invalid_unsigned_int;			
			if (dim>2)
			  indices[fe.face_to_cell_index(i,5)] = numbers::invalid_unsigned_int;			
		      }
		    add(k, indices);
		    break;
	      default:
		    Assert(false, ExcNotImplemented());
		    break;
	    }
	  ++k;
	}
    }
}


inline
unsigned int
BlockList::size() const
{
  return index_sets.size();
}


inline
unsigned int
BlockList::block_size(unsigned int block) const
{
  return index_sets[block].size();
}


inline
BlockList::const_iterator
BlockList::begin(unsigned int block) const
{
  AssertIndexRange(block, index_sets.size());
  return index_sets[block].begin();
}


inline
BlockList::const_iterator
BlockList::end(unsigned int block) const
{
  AssertIndexRange(block, index_sets.size());
  return index_sets[block].end();
}


inline
unsigned int
BlockList::local_index(unsigned int block, unsigned int index) const
{
  AssertIndexRange(block, index_sets.size());
  const block_container& b = index_sets[block];
  for (unsigned i=0;i<b.size();++i)
    if (b[i] == index)
      return i;
  return numbers::invalid_unsigned_int;
}


DEAL_II_NAMESPACE_CLOSE

#endif

//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2010 by the deal.II authors
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

#include <vector>
#include <utility>

DEAL_II_NAMESPACE_OPEN

/**
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

      for (std::vector<unsigned int>::const_iterator i=indices.begin();
	   i != indices.end(); ++i)
	index_sets[k].push_back(*i);
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

      for (std::vector<unsigned int>::const_iterator i=indices.begin();
	   i != indices.end(); ++i)
	index_sets[k].push_back(*i);
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
      
      for (unsigned int i=0; i<indices.size(); ++i)
	if (selected_dofs[i])
	  index_sets[k].push_back(indices[i]-offset);
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
      
      for (unsigned int i=0; i<indices.size(); ++i)
	if (selected_dofs[i])
	  index_sets[k].push_back(indices[i]-offset);
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

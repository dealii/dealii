// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__block_list_h
#define __deal2__block_list_h


#include <deal.II/base/subscriptor.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/fe/fe.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * @deprecated This class is experimental and will be
 * removed in a future release.
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
  /**
   * Declare the type for container size.
   */
  typedef types::global_dof_index size_type;

  /// The container for each index set
  typedef std::vector<size_type> block_container;
  /// The iterator for individual indices
  typedef block_container::const_iterator const_iterator;

  /**
   * Since SparsityPattern can
   * handle the tasks of BlockList,
   * this function allows us to
   * create one from an already
   * filled BlockList. A first step
   * to make BlockList obsolete.
   *
   * The additional integer
   * argument is the dimension of
   * the vector space.
   */
  void create_sparsity_pattern(SparsityPattern &sparsity, size_type n) const;

  /**
   * Add the indices in
   * <tt>indices</tt> to block
   * <tt>block</tt>, eliminating
   * repeated indices.
   */
  void add(size_type block, const std::vector<size_type> &indices);

  /**
   * Add the indices in
   * <tt>indices</tt> to block
   * <tt>block</tt>, eliminating
   * repeated indices. Only add
   * those indices for which
   * <tt>selected_indices</tt> is true.
   */
  void add(size_type block,
           const std::vector<size_type> &indices,
           const std::vector<bool> &selected_indices,
           size_type offset = 0);

  /**
   * Just set up the correct size
   * and assign indices to blocks later.
   */
  void initialize(size_type n_blocks);

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
  void initialize(size_type n_blocks,
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
  void initialize_mg(size_type n_blocks,
                     const ITERATOR begin,
                     const typename identity<ITERATOR>::type end) DEAL_II_DEPRECATED;

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
  void initialize(size_type n_blocks,
                  const ITERATOR begin,
                  const typename identity<ITERATOR>::type end,
                  const std::vector<bool> &selected_dofs,
                  size_type offset = 0) DEAL_II_DEPRECATED;
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
  void initialize_mg(size_type n_blocks,
                     const ITERATOR begin,
                     const typename identity<ITERATOR>::type end,
                     const std::vector<bool> &selected_dofs,
                     size_type offset = 0) DEAL_II_DEPRECATED;

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
  void initialize_vertex_patches_mg(size_type n_blocks,
                                    const ITERATOR begin,
                                    const typename identity<ITERATOR>::type end,
                                    const std::vector<bool> &selected_dofs = std::vector<bool>(),
                                    size_type offset = 0) DEAL_II_DEPRECATED;

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
    bool same_level_only) const DEAL_II_DEPRECATED;

  /**
   * @deprecated This function will
   * move to DoFTools.
   *
   */
  template <int dim, typename ITERATOR>
  bool cell_generates_vertex_patch(const ITERATOR cell,
                                   bool same_level_only) const DEAL_II_DEPRECATED;

  /**
   * The number of blocks.
   */
  size_type size() const;

  /**
   * The size of a single block.
   */
  size_type block_size(size_type block) const;

  /**
   * Iterator to the first index in block.
   */
  const_iterator begin(size_type block) const;
  /**
   * End iterator for a single block.
   */
  const_iterator end(size_type block) const;
  /**
   * Return the position of
   * <tt>index</tt> in
   * <tt>block</tt>, or
   * numbers::invalid_size_type,
   * if the index is not in the block.
   */
  size_type local_index(size_type block, size_type index) const;

private:
  /**
   * The container for t he index sets.
   */
  std::vector<block_container> index_sets;
} DEAL_II_DEPRECATED;


inline
void
BlockList::create_sparsity_pattern(SparsityPattern &sparsity, size_type n) const
{
  std::vector<unsigned int> sizes(size());
  for (size_type b=0; b<size(); ++b)
    sizes[b] = block_size(b);

  sparsity.reinit(size(), n, sizes);
  for (size_type b=0; b<size(); ++b)
    {
      for (const_iterator i = begin(b); i != end(b); ++i)
        sparsity.add(b,*i);
    }
  sparsity.compress();
}


inline
void
BlockList::add(const size_type block, const std::vector<size_type> &indices)
{
  AssertIndexRange(block, index_sets.size());

  for (size_type i=0; i<indices.size(); ++i)
    {
      const size_type k = indices[i];
      if (k==numbers::invalid_size_type)
        continue;
      if (std::find(index_sets[block].begin(), index_sets[block].end(), k)
          == index_sets[block].end())
        index_sets[block].push_back(k);
    }
}


inline
void
BlockList::add(
  const size_type block,
  const std::vector<size_type> &indices,
  const std::vector<bool> &selected,
  size_type offset)
{
  AssertIndexRange(block, index_sets.size());
  AssertDimension(indices.size(), selected.size());

  for (size_type i=0; i<indices.size(); ++i)
    {
      const size_type k = indices[i];
      if (k==numbers::invalid_size_type)
        continue;
      if (selected[i] && std::find(index_sets[block].begin(), index_sets[block].end(), k-offset)
          == index_sets[block].end())
        index_sets[block].push_back(k-offset);
    }
}


inline
void
BlockList::initialize(size_type n_blocks)
{
  index_sets.resize(n_blocks);
}


template <typename ITERATOR>
inline
void
BlockList::initialize(size_type n_blocks, const ITERATOR begin, const typename identity<ITERATOR>::type end)
{
  index_sets.resize(n_blocks);
  std::vector<size_type> indices;
  size_type k = 0;
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
BlockList::initialize_mg(size_type n_blocks, const ITERATOR begin, const typename identity<ITERATOR>::type end)
{
  index_sets.resize(n_blocks);
  std::vector<size_type> indices;
  size_type k = 0;
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
  size_type n_blocks,
  const ITERATOR begin,
  const typename identity<ITERATOR>::type end,
  const std::vector<bool> &selected_dofs,
  size_type offset)
{
  index_sets.resize(n_blocks);
  std::vector<size_type> indices;
  size_type k = 0;
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
  size_type n_blocks,
  const ITERATOR begin,
  const typename identity<ITERATOR>::type end,
  const std::vector<bool> &selected_dofs,
  size_type offset)
{
  index_sets.resize(n_blocks);
  std::vector<size_type> indices;
  size_type k = 0;
  for (ITERATOR cell = begin; cell != end; ++cell, ++k)
    {
      indices.resize(cell->get_fe().dofs_per_cell);
      AssertDimension(selected_dofs.size(), cell->get_fe().dofs_per_cell);

      cell->get_mg_dof_indices(indices);
      add(k, indices, selected_dofs, offset);
    }
}



template <int dim, typename ITERATOR>
inline
bool
BlockList::cell_generates_vertex_patch(
  const ITERATOR cell,
  bool same_level_only) const
{
  switch (dim)
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
  size_type n_blocks,
  const ITERATOR begin,
  const typename identity<ITERATOR>::type end,
  const std::vector<bool> &selected_dofs,
  size_type offset)
{
  Assert(selected_dofs.size() == 0, ExcNotImplemented());
  Assert(offset==0, ExcNotImplemented());
  const FiniteElement<dim> &fe = begin->get_fe();
  Assert(fe.dofs_per_vertex == 0, ExcNotImplemented());

  index_sets.resize(n_blocks);
  std::vector<size_type> indices;
  size_type k = 0;
  for (ITERATOR cell = begin; cell != end; ++cell)
    {
      if (cell_generates_vertex_patch<dim>(cell, true))
        {
          indices.resize(cell->get_fe().dofs_per_cell);

          switch (dim)
            {
            case 3:
              cell->neighbor(4)->get_mg_dof_indices(indices);
              for (unsigned int i=0; i<fe.dofs_per_face; ++i)
                {
                  indices[fe.face_to_cell_index(i,1)] = numbers::invalid_dof_index;
                  indices[fe.face_to_cell_index(i,3)] = numbers::invalid_dof_index;
                  indices[fe.face_to_cell_index(i,4)] = numbers::invalid_dof_index;
                }
              add(k, indices);
              cell->neighbor(4)->neighbor(0)->get_mg_dof_indices(indices);
              for (unsigned int i=0; i<fe.dofs_per_face; ++i)
                {
                  indices[fe.face_to_cell_index(i,0)] = numbers::invalid_dof_index;
                  indices[fe.face_to_cell_index(i,3)] = numbers::invalid_dof_index;
                  indices[fe.face_to_cell_index(i,4)] = numbers::invalid_dof_index;
                }
              add(k, indices);
              cell->neighbor(4)->neighbor(2)->get_mg_dof_indices(indices);
              for (unsigned int i=0; i<fe.dofs_per_face; ++i)
                {
                  indices[fe.face_to_cell_index(i,1)] = numbers::invalid_dof_index;
                  indices[fe.face_to_cell_index(i,2)] = numbers::invalid_dof_index;
                  indices[fe.face_to_cell_index(i,4)] = numbers::invalid_dof_index;
                }
              add(k, indices);
              cell->neighbor(4)->neighbor(2)->neighbor(0)->get_mg_dof_indices(indices);
              for (unsigned int i=0; i<fe.dofs_per_face; ++i)
                {
                  indices[fe.face_to_cell_index(i,0)] = numbers::invalid_dof_index;
                  indices[fe.face_to_cell_index(i,2)] = numbers::invalid_dof_index;
                  indices[fe.face_to_cell_index(i,4)] = numbers::invalid_dof_index;
                }
              add(k, indices);
            case 2:
              cell->neighbor(2)->get_mg_dof_indices(indices);
              for (unsigned int i=0; i<fe.dofs_per_face; ++i)
                {
                  indices[fe.face_to_cell_index(i,1)] = numbers::invalid_dof_index;
                  indices[fe.face_to_cell_index(i,2)] = numbers::invalid_dof_index;
                  if (dim>2)
                    indices[fe.face_to_cell_index(i,5)] = numbers::invalid_dof_index;
                }
              add(k, indices);
              cell->neighbor(2)->neighbor(0)->get_mg_dof_indices(indices);
              for (unsigned int i=0; i<fe.dofs_per_face; ++i)
                {
                  indices[fe.face_to_cell_index(i,0)] = numbers::invalid_dof_index;
                  indices[fe.face_to_cell_index(i,2)] = numbers::invalid_dof_index;
                  if (dim>2)
                    indices[fe.face_to_cell_index(i,5)] = numbers::invalid_dof_index;
                }
              add(k, indices);
            // no break here
            case 1:
              cell->get_mg_dof_indices(indices);
              for (unsigned int i=0; i<fe.dofs_per_face; ++i)
                {
                  indices[fe.face_to_cell_index(i,1)] = numbers::invalid_dof_index;
                  if (dim>1)
                    indices[fe.face_to_cell_index(i,3)] = numbers::invalid_dof_index;
                  if (dim>2)
                    indices[fe.face_to_cell_index(i,5)] = numbers::invalid_dof_index;
                }
              add(k, indices);
              cell->neighbor(0)->get_mg_dof_indices(indices);
              for (unsigned int i=0; i<fe.dofs_per_face; ++i)
                {
                  indices[fe.face_to_cell_index(i,0)] = numbers::invalid_dof_index;
                  if (dim>1)
                    indices[fe.face_to_cell_index(i,3)] = numbers::invalid_dof_index;
                  if (dim>2)
                    indices[fe.face_to_cell_index(i,5)] = numbers::invalid_dof_index;
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
BlockList::size_type
BlockList::size() const
{
  return index_sets.size();
}


inline
BlockList::size_type
BlockList::block_size(size_type block) const
{
  return index_sets[block].size();
}


inline
BlockList::const_iterator
BlockList::begin(size_type block) const
{
  AssertIndexRange(block, index_sets.size());
  return index_sets[block].begin();
}


inline
BlockList::const_iterator
BlockList::end(size_type block) const
{
  AssertIndexRange(block, index_sets.size());
  return index_sets[block].end();
}


inline
BlockList::size_type
BlockList::local_index(size_type block, size_type index) const
{
  AssertIndexRange(block, index_sets.size());
  const block_container &b = index_sets[block];
  for (size_type i=0; i<b.size(); ++i)
    if (b[i] == index)
      return i;
  return numbers::invalid_size_type;
}


DEAL_II_NAMESPACE_CLOSE

#endif

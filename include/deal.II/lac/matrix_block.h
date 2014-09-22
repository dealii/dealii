// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2014 by the deal.II authors
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

#ifndef __deal2__matrix_block_h
#define __deal2__matrix_block_h

#include <deal.II/base/config.h>
#include <deal.II/base/named_data.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/std_cxx11/shared_ptr.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/lac/block_indices.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/full_matrix.h>


DEAL_II_NAMESPACE_OPEN

template <class MATRIX> class MatrixBlock;

namespace internal
{
  template <class MATRIX>
  void
  reinit(MatrixBlock<MATRIX> &v,
         const BlockSparsityPattern &p);

  template <typename number>
  void
  reinit(MatrixBlock<dealii::SparseMatrix<number> > &v,
         const BlockSparsityPattern &p);
}

/**
 * A wrapper around a matrix object, storing the coordinates in a
 * block matrix as well.
 *
 * This class is an alternative to BlockMatrixBase, if you only want
 * to generate a single block of the system, not the whole
 * system. Using the add() functions of this class, it is possible to
 * use the standard assembling functions used for block matrices, but
 * only enter in one of the blocks and still avoiding the index
 * computations involved.

 * The reason for this class is, that we may need a different number
 * of matrices for different blocks in a block system. For example, a
 * preconditioner for the Oseen system can be built as a block system,
 * where the pressure block is of the form
 * <b>M</b><sup>-1</sup><b>FA</b><sup>-1</sup> with <b>M</b> the
 * pressure mass matrix, <b>A</b> the pressure Laplacian and <b>F</b>
 * the advection diffusion operator applied to the pressure
 * space. Since only a single matrix is needed for the other blocks,
 * using BlockSparseMatrix or similar would be a waste of memory.
 *
 * While the add() functions make a MatrixBlock appear like a block
 * matrix for assembling, the functions vmult(), Tvmult(),
 * vmult_add(), and Tvmult_add() make it behave like a MATRIX, when it
 * comes to applying it to a vector. This behavior allows us to store
 * MatrixBlock objects in vectors, for instance in MGLevelObject
 * without extracting the #matrix first.
 *
 * MatrixBlock comes handy when using BlockMatrixArray. Once the
 * MatrixBlock has been properly initalized and filled, it can be used
 * in the simplest case as:
 * @code
 * MatrixBlockVector<SparseMatrix<double> > > blocks;
 *
 * ...
 *
 * BlockMatrixArray matrix (n_blocks, n_blocks);
 *
 * for (size_type i=0;i<blocks.size;++i)
 *   matrix.enter(blocks.block(i).row, blocks.block(i).column, blocks.matrix(i));
 * @endcode
 *
 * Here, we have not gained very much, except that we do not need to
 * set up empty blocks in the block system.
 *
 * @note This class expects, that the row and column BlockIndices
 * objects for the system are equal. If they are not, some functions
 * will throw ExcNotImplemented.
 *
 * @todo Example for the product preconditioner of the pressure Schur
 * complement.
 *
 * @ingroup Matrix2
 * @ingroup vector_valued
 * @see @ref GlossBlockLA "Block (linear algebra)"
 * @author Guido Kanschat, 2006
 */
template <class MATRIX>
class MatrixBlock
  : public Subscriptor
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Constructor rendering an
   * uninitialized object.
   */
  MatrixBlock();

  /**
   * Copy constructor.
   */
  MatrixBlock(const MatrixBlock<MATRIX> &M);

  /**
   * Constructor setting block
   * coordinates, but not
   * initializing the matrix.
   */

  MatrixBlock(size_type i, size_type j);

  /**
   * Reinitialize the matrix for a
   * new BlockSparsityPattern. This
   * adjusts the #matrix as well
   * as the #row_indices and
   * #column_indices.
   *
   * @note The row and column block
   * structure of the sparsity
   * pattern must be equal.
   */
  void reinit(const BlockSparsityPattern &sparsity);

  operator MATRIX &();
  operator const MATRIX &() const;

  /**
   * Add <tt>value</tt> to the
   * element (<i>i,j</i>). Throws
   * an error if the entry does not
   * exist or if it is in a
   * different block.
   */
  void add (const size_type i,
            const size_type j,
            const typename MATRIX::value_type value);

  /**
   * Add all elements in a
   * FullMatrix into sparse
   * matrix locations given by
   * <tt>indices</tt>. This function
   * assumes a quadratic sparse
   * matrix and a quadratic
   * full_matrix.  The global
   * locations are translated
   * into locations in this block
   * and ExcBlockIndexMismatch is
   * thrown, if the global index
   * does not point into the
   * block referred to by #row and
   * #column.
   *
   * @todo
   * <tt>elide_zero_values</tt> is
   * currently ignored.
   *
   * The optional parameter
   * <tt>elide_zero_values</tt> can be
   * used to specify whether zero
   * values should be added anyway or
   * these should be filtered away and
   * only non-zero data is added. The
   * default value is <tt>true</tt>,
   * i.e., zero values won't be added
   * into the matrix.
   */
  template <typename number>
  void add (const std::vector<size_type> &indices,
            const FullMatrix<number>     &full_matrix,
            const bool                    elide_zero_values = true);

  /**
   * Add all elements in a
   * FullMatrix into global
   * locations given by
   * <tt>row_indices</tt> and
   * <tt>col_indices</tt>,
   * respectively. The global
   * locations are translated
   * into locations in this block
   * and ExcBlockIndexMismatch is
   * thrown, if the global index
   * does not point into the
   * block referred to by #row and
   * #column.
   *
   * @todo
   * <tt>elide_zero_values</tt> is
   * currently ignored.
   *
   * The optional parameter
   * <tt>elide_zero_values</tt> can be
   * used to specify whether zero
   * values should be added anyway or
   * these should be filtered away and
   * only non-zero data is added. The
   * default value is <tt>true</tt>,
   * i.e., zero values won't be added
   * into the matrix.
   */
  template <typename number>
  void add (const std::vector<size_type> &row_indices,
            const std::vector<size_type> &col_indices,
            const FullMatrix<number>        &full_matrix,
            const bool                       elide_zero_values = true);

  /**
   * Set several elements in the
   * specified row of the matrix
   * with column indices as given
   * by <tt>col_indices</tt> to
   * the respective value. This
   * is the function doing thye
   * actual work for the ones
   * adding full matrices. The
   * global locations
   * <tt>row_index</tt> and
   * <tt>col_indices</tt> are
   * translated into locations in
   * this block and
   * ExcBlockIndexMismatch is
   * thrown, if the global index
   * does not point into the
   * block referred to by #row and
   * #column.
   *
   * @todo
   * <tt>elide_zero_values</tt> is
   * currently ignored.
   *
   * The optional parameter
   * <tt>elide_zero_values</tt> can be
   * used to specify whether zero
   * values should be added anyway or
   * these should be filtered away and
   * only non-zero data is added. The
   * default value is <tt>true</tt>,
   * i.e., zero values won't be added
   * into the matrix.
   */
  template <typename number>
  void add (const size_type               row_index,
            const std::vector<size_type> &col_indices,
            const std::vector<number>    &values,
            const bool                    elide_zero_values = true);

  /**
   * Add an array of values given by
   * <tt>values</tt> in the given
   * global matrix row at columns
   * specified by col_indices in the
   * sparse matrix.
   *
   * The optional parameter
   * <tt>elide_zero_values</tt> can be
   * used to specify whether zero
   * values should be added anyway or
   * these should be filtered away and
   * only non-zero data is added. The
   * default value is <tt>true</tt>,
   * i.e., zero values won't be added
   * into the matrix.
   */
  template <typename number>
  void add (const size_type  row,
            const size_type  n_cols,
            const size_type *col_indices,
            const number    *values,
            const bool       elide_zero_values = true,
            const bool       col_indices_are_sorted = false);

  /**
   * Matrix-vector-multiplication,
   * forwarding to the same
   * function in MATRIX. No index
   * computations are done, thus,
   * the vectors need to have sizes
   * matching #matrix.
   */
  template<class VECTOR>
  void vmult (VECTOR &w, const VECTOR &v) const;

  /**
   * Matrix-vector-multiplication,
   * forwarding to the same
   * function in MATRIX. No index
   * computations are done, thus,
   * the vectors need to have sizes
   * matching #matrix.
   */
  template<class VECTOR>
  void vmult_add (VECTOR &w, const VECTOR &v) const;

  /**
   * Matrix-vector-multiplication,
   * forwarding to the same
   * function in MATRIX. No index
   * computations are done, thus,
   * the vectors need to have sizes
   * matching #matrix.
   */
  template<class VECTOR>
  void Tvmult (VECTOR &w, const VECTOR &v) const;

  /**
   * Matrix-vector-multiplication,
   * forwarding to the same
   * function in MATRIX. No index
   * computations are done, thus,
   * the vectors need to have sizes
   * matching #matrix.
   */
  template<class VECTOR>
  void Tvmult_add (VECTOR &w, const VECTOR &v) const;

  /**
   * The memory used by this object.
   */
  std::size_t memory_consumption () const;

  /**
   * The block number computed from
   * an index by using
   * BlockIndices does not match
   * the block coordinates stored
   * in this object.
   */
  DeclException2(ExcBlockIndexMismatch, size_type, size_type,
                 << "Block index " << arg1 << " does not match " << arg2);

  /**
   * Row coordinate.  This is the
   * position of the data member
   * matrix on the global matrix.
   */
  size_type row;
  /**
   * Column coordinate.  This is
   * the position of the data
   * member matrix on the global
   * matrix.
   */
  size_type column;

  /**
   * The matrix itself
   */
  MATRIX matrix;

private:
  /**
   * The rwo BlockIndices of the
   * whole system. Using row(),
   * this allows us to find the
   * index of the first row degree
   * of freedom for this block.
   */
  BlockIndices row_indices;
  /**
   * The column BlockIndices of the
   * whole system. Using column(),
   * this allows us to find the
   * index of the first column
   * degree of freedom for this
   * block.
   */
  BlockIndices column_indices;

  template <class OTHER_MATRIX>
  friend
  void dealii::internal::reinit(MatrixBlock<OTHER_MATRIX> &,
                                const BlockSparsityPattern &);

  template <typename number>
  friend
  void internal::reinit(MatrixBlock<dealii::SparseMatrix<number> > &v,
                        const BlockSparsityPattern &p);
};


/**
 * A vector of MatrixBlock, which is implemented using shared
 * pointers, in order to allow for copying and rearranging. Each
 * matrix block can be identified by name.
 *
 * @relates MatrixBlock
 * @ingroup vector_valued
 * @author Baerbel Janssen, Guido Kanschat, 2010
 */
template <class MATRIX>
class MatrixBlockVector
  :
  private NamedData<std_cxx11::shared_ptr<MatrixBlock<MATRIX> > >
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * The type of object stored.
   */
  typedef MatrixBlock<MATRIX> value_type;

  /**
   * Add a new matrix block at the
   * position <tt>(row,column)</tt>
   * in the block system.
   */
  void add(size_type row, size_type column, const std::string &name);

  /**
   * For matrices using a
   * SparsityPattern, this function
   * reinitializes each matrix in
   * the vector with the correct
   * pattern from the block system.
   */
  void reinit(const BlockSparsityPattern &sparsity);

  /**
   * Clears the object.
   *
   * Since often only clearing of
   * the individual matrices is
   * desired, but not removing the
   * blocks themselves, there is an
   * optional argument. If the
   * argument is missing or @p
   * false, all matrices will be
   * mepty, but the size of this
   * object and the block positions
   * will not change. If @p
   * really_clean is @p true, then
   * the object will contain no
   * blocks at the end.
   */
  void clear (bool really_clean = false);

  /**
   * The memory used by this object.
   */
  std::size_t memory_consumption () const;

  /**
   * Access a constant reference to
   * the block at position <i>i</i>.
   */
  const value_type &block(size_type i) const;

  /**
   * Access a reference to
   * the block at position <i>i</i>.
   */
  value_type &block(size_type i);

  /**
   * Access the matrix at position
   * <i>i</i> for read and write
   * access.
   */
  MATRIX &matrix(size_type i);

  /**
   * import functions from private base class
   */
  using NamedData<std_cxx11::shared_ptr<value_type> >::subscribe;
  using NamedData<std_cxx11::shared_ptr<value_type> >::unsubscribe;
  using NamedData<std_cxx11::shared_ptr<value_type> >::size;
  using NamedData<std_cxx11::shared_ptr<value_type> >::name;
};


/**
 * A vector of MGLevelObject<MatrixBlock>, which is implemented using shared
 * pointers, in order to allow for copying and rearranging. Each
 * matrix block can be identified by name.
 *
 * @relates MatrixBlock
 * @ingroup vector_valued
 * @author Baerbel Janssen, Guido Kanschat, 2010
 */
template <class MATRIX>
class MGMatrixBlockVector
  : public Subscriptor
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * The type of object stored.
   */
  typedef MGLevelObject<MatrixBlock<MATRIX> > value_type;
  /**
   * Constructor, determining which
   * matrices should be stored.
   *
   * If <tt>edge_matrices</tt> is
   * true, then objects for edge
   * matrices for discretizations
   * with degrees of freedom on
   * faces are allocated.
   *
   * If <tt>edge_flux_matrices</tt>
   * is true, then objects for DG
   * fluxes on the refinement edge
   * are allocated.
   */
  MGMatrixBlockVector(const bool edge_matrices = false,
                      const bool edge_flux_matrices = false);

  /**
   * The number of blocks.
   */
  unsigned int size () const;

  /**
   * Add a new matrix block at the
   * position <tt>(row,column)</tt>
   * in the block system. The third
   * argument allows to give the
   * matrix a name for later
   * identification.
   */
  void add(size_type row, size_type column, const std::string &name);

  /**
   * For matrices using a
   * SparsityPattern, this function
   * reinitializes each matrix in
   * the vector with the correct
   * pattern from the block system.
   *
   * This function reinitializes
   * the level matrices.
   */
  void reinit_matrix(const MGLevelObject<BlockSparsityPattern> &sparsity);
  /**
   * For matrices using a
   * SparsityPattern, this function
   * reinitializes each matrix in
   * the vector with the correct
   * pattern from the block system.
   *
   * This function reinitializes
   * the matrices for degrees of
   * freedom on the refinement edge.
   */
  void reinit_edge(const MGLevelObject<BlockSparsityPattern> &sparsity);
  /**
   * For matrices using a
   * SparsityPattern, this function
   * reinitializes each matrix in
   * the vector with the correct
   * pattern from the block system.
   *
   * This function reinitializes
   * the flux matrices over the
   * refinement edge.
   */
  void reinit_edge_flux(const MGLevelObject<BlockSparsityPattern> &sparsity);

  /**
   * Clears the object.
   *
   * Since often only clearing of
   * the individual matrices is
   * desired, but not removing the
   * blocks themselves, there is an
   * optional argument. If the
   * argument is missing or @p
   * false, all matrices will be
   * mepty, but the size of this
   * object and the block positions
   * will not change. If @p
   * really_clean is @p true, then
   * the object will contain no
   * blocks at the end.
   */
  void clear (bool really_clean = false);

  /**
   * Access a constant reference to
   * the matrix block at position <i>i</i>.
   */
  const value_type &block(size_type i) const;

  /**
   * Access a reference to
   * the matrix block at position <i>i</i>.
   */
  value_type &block(size_type i);

  /**
   * Access a constant reference to
   * the edge matrix block at position <i>i</i>.
   */
  const value_type &block_in(size_type i) const;

  /**
   * Access a reference to
   * the edge matrix block at position <i>i</i>.
   */
  value_type &block_in(size_type i);

  /**
   * Access a constant reference to
   * the edge matrix block at position <i>i</i>.
   */
  const value_type &block_out(size_type i) const;

  /**
   * Access a reference to
   * the edge matrix block at position <i>i</i>.
   */
  value_type &block_out(size_type i);

  /**
   * Access a constant reference to
   * the  edge flux matrix block at position <i>i</i>.
   */
  const value_type &block_up(size_type i) const;

  /**
   * Access a reference to
   * the  edge flux matrix block at position <i>i</i>.
   */
  value_type &block_up(size_type i);

  /**
   * Access a constant reference to
   * the  edge flux matrix block at position <i>i</i>.
   */
  const value_type &block_down(size_type i) const;

  /**
   * Access a reference to
   * the edge flux matrix block at position <i>i</i>.
   */
  value_type &block_down(size_type i);

  /**
   * The memory used by this object.
   */
  std::size_t memory_consumption () const;
private:
  /// Clear one of the matrix objects
  void clear_object(NamedData<MGLevelObject<MatrixBlock<MATRIX> > > &);

  /// Flag for storing matrices_in and matrices_out
  const bool edge_matrices;

  /// Flag for storing flux_matrices_up and flux_matrices_down
  const bool edge_flux_matrices;

  /// The level matrices
  NamedData<MGLevelObject<MatrixBlock<MATRIX> > > matrices;
  /// The matrix from the interior of a level to the refinement edge
  NamedData<MGLevelObject<MatrixBlock<MATRIX> > > matrices_in;
  /// The matrix from the refinement edge to the interior of a level
  NamedData<MGLevelObject<MatrixBlock<MATRIX> > > matrices_out;
  /// The DG flux from a level to the lower level
  NamedData<MGLevelObject<MatrixBlock<MATRIX> > > flux_matrices_down;
  /// The DG flux from the lower level to a level
  NamedData<MGLevelObject<MatrixBlock<MATRIX> > > flux_matrices_up;
};


//----------------------------------------------------------------------//

namespace internal
{
  template <class MATRIX>
  void
  reinit(MatrixBlock<MATRIX> &v,
         const BlockSparsityPattern &p)
  {
    v.row_indices = p.get_row_indices();
    v.column_indices = p.get_column_indices();
  }


  template <typename number>
  void
  reinit(MatrixBlock<dealii::SparseMatrix<number> > &v,
         const BlockSparsityPattern &p)
  {
    v.row_indices = p.get_row_indices();
    v.column_indices = p.get_column_indices();
    v.matrix.reinit(p.block(v.row, v.column));
  }
}


template <class MATRIX>
inline
MatrixBlock<MATRIX>::MatrixBlock()
  :
  row(deal_II_numbers::invalid_size_type),
  column(deal_II_numbers::invalid_size_type)
{}


template <class MATRIX>
inline
MatrixBlock<MATRIX>::MatrixBlock(const MatrixBlock<MATRIX> &M)
  :
  Subscriptor(),
  row(M.row),
  column(M.column),
  matrix(M.matrix),
  row_indices(M.row_indices),
  column_indices(M.column_indices)
{}


template <class MATRIX>
inline
MatrixBlock<MATRIX>::MatrixBlock(size_type i, size_type j)
  :
  row(i), column(j)
{}


template <class MATRIX>
inline
void
MatrixBlock<MATRIX>::reinit(const BlockSparsityPattern &sparsity)
{
  internal::reinit(*this, sparsity);
}


template <class MATRIX>
inline
MatrixBlock<MATRIX>::operator MATRIX &()
{
  return matrix;
}


template <class MATRIX>
inline
MatrixBlock<MATRIX>::operator const MATRIX &() const
{
  return matrix;
}


template <class MATRIX>
inline void
MatrixBlock<MATRIX>::add (
  const size_type gi,
  const size_type gj,
  const typename MATRIX::value_type value)
{
  Assert(row_indices.size() != 0, ExcNotInitialized());
  Assert(column_indices.size() != 0, ExcNotInitialized());

  const std::pair<unsigned int, size_type> bi
    = row_indices.global_to_local(gi);
  const std::pair<unsigned int, size_type> bj
    = column_indices.global_to_local(gj);

  Assert (bi.first == row, ExcBlockIndexMismatch(bi.first, row));
  Assert (bj.first == column, ExcBlockIndexMismatch(bj.first, column));

  matrix.add(bi.second, bj.second, value);
}


template <class MATRIX>
template <typename number>
inline
void
MatrixBlock<MATRIX>::add (const std::vector<size_type> &r_indices,
                          const std::vector<size_type> &c_indices,
                          const FullMatrix<number>     &values,
                          const bool                    elide_zero_values)
{
  Assert(row_indices.size() != 0, ExcNotInitialized());
  Assert(column_indices.size() != 0, ExcNotInitialized());

  AssertDimension (r_indices.size(), values.m());
  AssertDimension (c_indices.size(), values.n());

  for (size_type i=0; i<row_indices.size(); ++i)
    add (r_indices[i], c_indices.size(), &c_indices[0], &values(i,0),
         elide_zero_values);
}


template <class MATRIX>
template <typename number>
inline
void
MatrixBlock<MATRIX>::add (const size_type  b_row,
                          const size_type  n_cols,
                          const size_type *col_indices,
                          const number    *values,
                          const bool,
                          const bool)
{
  Assert(row_indices.size() != 0, ExcNotInitialized());
  Assert(column_indices.size() != 0, ExcNotInitialized());

  const std::pair<unsigned int, size_type> bi
    = row_indices.global_to_local(b_row);

  // In debug mode, we check whether
  // all indices are in the correct
  // block.

  // Actually, for the time being, we
  // leave it at this. While it may
  // not be the most efficient way,
  // it is at least thread safe.
//#ifdef DEBUG
  Assert(bi.first == row, ExcBlockIndexMismatch(bi.first, row));

  for (size_type j=0; j<n_cols; ++j)
    {
      const std::pair<unsigned int, size_type> bj
        = column_indices.global_to_local(col_indices[j]);
      Assert(bj.first == column, ExcBlockIndexMismatch(bj.first, column));

      matrix.add(bi.second, bj.second, values[j]);
    }
//#endif
}


template <class MATRIX>
template <typename number>
inline
void
MatrixBlock<MATRIX>::add (const std::vector<size_type> &indices,
                          const FullMatrix<number>     &values,
                          const bool                    elide_zero_values)
{
  Assert(row_indices.size() != 0, ExcNotInitialized());
  Assert(column_indices.size() != 0, ExcNotInitialized());

  AssertDimension (indices.size(), values.m());
  Assert (values.n() == values.m(), ExcNotQuadratic());

  for (size_type i=0; i<indices.size(); ++i)
    add (indices[i], indices.size(), &indices[0], &values(i,0),
         elide_zero_values);
}



template <class MATRIX>
template <typename number>
inline
void
MatrixBlock<MATRIX>::add (const size_type               row,
                          const std::vector<size_type> &col_indices,
                          const std::vector<number>    &values,
                          const bool                    elide_zero_values)
{
  Assert(row_indices.size() != 0, ExcNotInitialized());
  Assert(column_indices.size() != 0, ExcNotInitialized());

  AssertDimension (col_indices.size(), values.size());
  add (row, col_indices.size(), &col_indices[0], &values[0],
       elide_zero_values);
}


template <class MATRIX>
template <class VECTOR>
inline
void
MatrixBlock<MATRIX>::vmult (VECTOR &w, const VECTOR &v) const
{
  matrix.vmult(w,v);
}


template <class MATRIX>
template <class VECTOR>
inline
void
MatrixBlock<MATRIX>::vmult_add (VECTOR &w, const VECTOR &v) const
{
  matrix.vmult_add(w,v);
}


template <class MATRIX>
template <class VECTOR>
inline
void
MatrixBlock<MATRIX>::Tvmult (VECTOR &w, const VECTOR &v) const
{
  matrix.Tvmult(w,v);
}


template <class MATRIX>
template <class VECTOR>
inline
void
MatrixBlock<MATRIX>::Tvmult_add (VECTOR &w, const VECTOR &v) const
{
  matrix.Tvmult_add(w,v);
}


template <class MATRIX>
inline
std::size_t
MatrixBlock<MATRIX>::memory_consumption () const
{
  return (sizeof(*this)
          + MemoryConsumption::memory_consumption(matrix)
          - sizeof(matrix));
}

//----------------------------------------------------------------------//

template <class MATRIX>
inline void
MatrixBlockVector<MATRIX>::add(
  size_type row, size_type column,
  const std::string &name)
{
  std_cxx11::shared_ptr<value_type> p(new value_type(row, column));
  NamedData<std_cxx11::shared_ptr<value_type> >::add(p, name);
}


template <class MATRIX>
inline void
MatrixBlockVector<MATRIX>::reinit(const BlockSparsityPattern &sparsity)
{
  for (size_type i=0; i<this->size(); ++i)
    {
      block(i).reinit(sparsity);
    }
}


template <class MATRIX>
inline void
MatrixBlockVector<MATRIX>::clear(bool really_clean)
{
  if (really_clean)
    {
      Assert(false, ExcNotImplemented());
    }
  else
    {
      for (size_type i=0; i<this->size(); ++i)
        matrix(i).clear();
    }
}



template <class MATRIX>
inline const MatrixBlock<MATRIX> &
MatrixBlockVector<MATRIX>::block(size_type i) const
{
  return *this->read(i);
}


template <class MATRIX>
inline MatrixBlock<MATRIX> &
MatrixBlockVector<MATRIX>::block(size_type i)
{
  return *(*this)(i);
}


template <class MATRIX>
inline MATRIX &
MatrixBlockVector<MATRIX>::matrix(size_type i)
{
  return (*this)(i)->matrix;
}



//----------------------------------------------------------------------//

template <class MATRIX>
inline
MGMatrixBlockVector<MATRIX>::MGMatrixBlockVector(
  const bool e, const bool f)
  :
  edge_matrices(e),
  edge_flux_matrices(f)
{}


template <class MATRIX>
inline
unsigned int
MGMatrixBlockVector<MATRIX>::size () const
{
  return matrices.size();
}


template <class MATRIX>
inline void
MGMatrixBlockVector<MATRIX>::add(
  size_type row, size_type column,
  const std::string &name)
{
  MGLevelObject<MatrixBlock<MATRIX> > p(0, 1);
  p[0].row = row;
  p[0].column = column;

  matrices.add(p, name);
  if (edge_matrices)
    {
      matrices_in.add(p, name);
      matrices_out.add(p, name);
    }
  if (edge_flux_matrices)
    {
      flux_matrices_up.add(p, name);
      flux_matrices_down.add(p, name);
    }
}


template <class MATRIX>
inline const MGLevelObject<MatrixBlock<MATRIX> > &
MGMatrixBlockVector<MATRIX>::block(size_type i) const
{
  return matrices.read(i);
}


template <class MATRIX>
inline MGLevelObject<MatrixBlock<MATRIX> > &
MGMatrixBlockVector<MATRIX>::block(size_type i)
{
  return matrices(i);
}


template <class MATRIX>
inline const MGLevelObject<MatrixBlock<MATRIX> > &
MGMatrixBlockVector<MATRIX>::block_in(size_type i) const
{
  return matrices_in.read(i);
}


template <class MATRIX>
inline MGLevelObject<MatrixBlock<MATRIX> > &
MGMatrixBlockVector<MATRIX>::block_in(size_type i)
{
  return matrices_in(i);
}


template <class MATRIX>
inline const MGLevelObject<MatrixBlock<MATRIX> > &
MGMatrixBlockVector<MATRIX>::block_out(size_type i) const
{
  return matrices_out.read(i);
}


template <class MATRIX>
inline MGLevelObject<MatrixBlock<MATRIX> > &
MGMatrixBlockVector<MATRIX>::block_out(size_type i)
{
  return matrices_out(i);
}


template <class MATRIX>
inline const MGLevelObject<MatrixBlock<MATRIX> > &
MGMatrixBlockVector<MATRIX>::block_up(size_type i) const
{
  return flux_matrices_up.read(i);
}


template <class MATRIX>
inline MGLevelObject<MatrixBlock<MATRIX> > &
MGMatrixBlockVector<MATRIX>::block_up(size_type i)
{
  return flux_matrices_up(i);
}


template <class MATRIX>
inline const MGLevelObject<MatrixBlock<MATRIX> > &
MGMatrixBlockVector<MATRIX>::block_down(size_type i) const
{
  return flux_matrices_down.read(i);
}


template <class MATRIX>
inline MGLevelObject<MatrixBlock<MATRIX> > &
MGMatrixBlockVector<MATRIX>::block_down(size_type i)
{
  return flux_matrices_down(i);
}


template <class MATRIX>
inline void
MGMatrixBlockVector<MATRIX>::reinit_matrix(const MGLevelObject<BlockSparsityPattern> &sparsity)
{
  for (size_type i=0; i<this->size(); ++i)
    {
      MGLevelObject<MatrixBlock<MATRIX> > &o = block(i);
      const size_type row = o[o.min_level()].row;
      const size_type col = o[o.min_level()].column;

      o.resize(sparsity.min_level(), sparsity.max_level());
      for (size_type level = o.min_level(); level <= o.max_level(); ++level)
        {
          o[level].row = row;
          o[level].column = col;
          internal::reinit(o[level], sparsity[level]);
        }
    }
}


template <class MATRIX>
inline void
MGMatrixBlockVector<MATRIX>::reinit_edge(const MGLevelObject<BlockSparsityPattern> &sparsity)
{
  for (size_type i=0; i<this->size(); ++i)
    {
      MGLevelObject<MatrixBlock<MATRIX> > &o = block(i);
      const size_type row = o[o.min_level()].row;
      const size_type col = o[o.min_level()].column;

      block_in(i).resize(sparsity.min_level(), sparsity.max_level());
      block_out(i).resize(sparsity.min_level(), sparsity.max_level());
      for (size_type level = o.min_level(); level <= o.max_level(); ++level)
        {
          block_in(i)[level].row = row;
          block_in(i)[level].column = col;
          internal::reinit(block_in(i)[level], sparsity[level]);
          block_out(i)[level].row = row;
          block_out(i)[level].column = col;
          internal::reinit(block_out(i)[level], sparsity[level]);
        }
    }
}


template <class MATRIX>
inline void
MGMatrixBlockVector<MATRIX>::reinit_edge_flux(const MGLevelObject<BlockSparsityPattern> &sparsity)
{
  for (size_type i=0; i<this->size(); ++i)
    {
      MGLevelObject<MatrixBlock<MATRIX> > &o = block(i);
      const size_type row = o[o.min_level()].row;
      const size_type col = o[o.min_level()].column;

      block_up(i).resize(sparsity.min_level(), sparsity.max_level());
      block_down(i).resize(sparsity.min_level(), sparsity.max_level());
      for (size_type level = o.min_level(); level <= o.max_level(); ++level)
        {
          block_up(i)[level].row = row;
          block_up(i)[level].column = col;
          internal::reinit(block_up(i)[level], sparsity[level]);
          block_down(i)[level].row = row;
          block_down(i)[level].column = col;
          internal::reinit(block_down(i)[level], sparsity[level]);
        }

    }
}


template <class MATRIX>
inline void
MGMatrixBlockVector<MATRIX>::clear_object(NamedData<MGLevelObject<MatrixBlock<MATRIX> > > &mo)
{
  for (size_type i=0; i<mo.size(); ++i)
    {
      MGLevelObject<MatrixBlock<MATRIX> > &o = mo(i);
      for (size_type level = o.min_level(); level <= o.max_level(); ++level)
        o[level].matrix.clear();
    }
}


template <class MATRIX>
inline void
MGMatrixBlockVector<MATRIX>::clear(bool really_clean)
{
  if (really_clean)
    {
      Assert(false, ExcNotImplemented());
    }
  else
    {
      clear_object(matrices);
      clear_object(matrices_in);
      clear_object(matrices_out);
      clear_object(flux_matrices_up);
      clear_object(flux_matrices_down);
    }
}





DEAL_II_NAMESPACE_CLOSE

#endif

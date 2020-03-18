// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_matrix_block_h
#define dealii_matrix_block_h

#include <deal.II/base/config.h>

#include <deal.II/algorithms/any_data.h>

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/smartpointer.h>

#include <deal.II/lac/block_indices.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>

#include <memory>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <typename MatrixType>
class MatrixBlock;
#endif

namespace internal
{
  template <typename MatrixType>
  void
  reinit(MatrixBlock<MatrixType> &v, const BlockSparsityPattern &p);

  template <typename number>
  void
  reinit(MatrixBlock<dealii::SparseMatrix<number>> &v,
         const BlockSparsityPattern &               p);
} // namespace internal

/**
 * A wrapper around a matrix object, storing the coordinates in a block matrix
 * as well.
 *
 * This class is an alternative to BlockMatrixBase, if you only want to
 * generate a single block of the system, not the whole system. Using the
 * add() functions of this class, it is possible to use the standard
 * assembling functions used for block matrices, but only enter in one of the
 * blocks and still avoiding the index computations involved.  The reason for
 * this class is, that we may need a different number of matrices for
 * different blocks in a block system. For example, a preconditioner for the
 * Oseen system can be built as a block system, where the pressure block is of
 * the form <b>M</b><sup>-1</sup><b>FA</b><sup>-1</sup> with <b>M</b> the
 * pressure mass matrix, <b>A</b> the pressure Laplacian and <b>F</b> the
 * advection diffusion operator applied to the pressure space. Since only a
 * single matrix is needed for the other blocks, using BlockSparseMatrix or
 * similar would be a waste of memory.
 *
 * While the add() functions make a MatrixBlock appear like a block matrix for
 * assembling, the functions vmult(), Tvmult(), vmult_add(), and Tvmult_add()
 * make it behave like a MatrixType, when it comes to applying it to a vector.
 * This behavior allows us to store MatrixBlock objects in vectors, for
 * instance in MGLevelObject without extracting the #matrix first.
 *
 * MatrixBlock comes handy when using BlockMatrixArray. Once the MatrixBlock
 * has been properly initialized and filled, it can be used in the simplest
 * case as:
 * @code
 * MatrixBlockVector<SparseMatrix<double> > > blocks;
 *
 * ...
 *
 * BlockMatrixArray matrix (n_blocks, n_blocks);
 *
 * for (size_type i=0;i<blocks.size;++i)
 *   matrix.enter(blocks.block(i).row, blocks.block(i).column,
 * blocks.matrix(i));
 * @endcode
 *
 * Here, we have not gained very much, except that we do not need to set up
 * empty blocks in the block system.
 *
 * @note This class expects, that the row and column BlockIndices objects for
 * the system are equal. If they are not, some functions will throw
 * ExcNotImplemented.
 *
 * @todo Example for the product preconditioner of the pressure Schur
 * complement.
 *
 * @ingroup Matrix2
 * @ingroup vector_valued
 *
 * @see
 * @ref GlossBlockLA "Block (linear algebra)"
 * @author Guido Kanschat, 2006
 */
template <typename MatrixType>
class MatrixBlock : public Subscriptor
{
public:
  /**
   * Declare type for container size.
   */
  using size_type = types::global_dof_index;

  /**
   * Declare a type for matrix entries.
   */
  using value_type = typename MatrixType::value_type;

  /**
   * Constructor rendering an uninitialized object.
   */
  MatrixBlock();

  /**
   * Copy constructor.
   */
  MatrixBlock(const MatrixBlock<MatrixType> &M) = default;

  /**
   * Assignment operator.
   */
  MatrixBlock<MatrixType> &
  operator=(const MatrixBlock<MatrixType> &) = default;

  /**
   * Constructor setting block coordinates, but not initializing the matrix.
   */

  MatrixBlock(size_type i, size_type j);

  /**
   * Reinitialize the matrix for a new BlockSparsityPattern. This adjusts the
   * #matrix as well as the #row_indices and #column_indices.
   *
   * @note The row and column block structure of the sparsity pattern must be
   * equal.
   */
  void
  reinit(const BlockSparsityPattern &sparsity);

  operator MatrixType &();
  operator const MatrixType &() const;

  /**
   * Add <tt>value</tt> to the element (<i>i,j</i>). Throws an error if the
   * entry does not exist or if it is in a different block.
   */
  void
  add(const size_type                       i,
      const size_type                       j,
      const typename MatrixType::value_type value);

  /**
   * Add all elements in a FullMatrix into sparse matrix locations given by
   * <tt>indices</tt>. This function assumes a quadratic sparse matrix and a
   * quadratic full_matrix.  The global locations are translated into
   * locations in this block and ExcBlockIndexMismatch is thrown, if the
   * global index does not point into the block referred to by #row and
   * #column.
   *
   * @todo <tt>elide_zero_values</tt> is currently ignored.
   *
   * The optional parameter <tt>elide_zero_values</tt> can be used to specify
   * whether zero values should be added anyway or these should be filtered
   * away and only non-zero data is added. The default value is <tt>true</tt>,
   * i.e., zero values won't be added into the matrix.
   */
  template <typename number>
  void
  add(const std::vector<size_type> &indices,
      const FullMatrix<number> &    full_matrix,
      const bool                    elide_zero_values = true);

  /**
   * Add all elements in a FullMatrix into global locations given by
   * <tt>row_indices</tt> and <tt>col_indices</tt>, respectively. The global
   * locations are translated into locations in this block and
   * ExcBlockIndexMismatch is thrown, if the global index does not point into
   * the block referred to by #row and #column.
   *
   * @todo <tt>elide_zero_values</tt> is currently ignored.
   *
   * The optional parameter <tt>elide_zero_values</tt> can be used to specify
   * whether zero values should be added anyway or these should be filtered
   * away and only non-zero data is added. The default value is <tt>true</tt>,
   * i.e., zero values won't be added into the matrix.
   */
  template <typename number>
  void
  add(const std::vector<size_type> &row_indices,
      const std::vector<size_type> &col_indices,
      const FullMatrix<number> &    full_matrix,
      const bool                    elide_zero_values = true);

  /**
   * Set several elements in the specified row of the matrix with column
   * indices as given by <tt>col_indices</tt> to the respective value. This is
   * the function doing the actual work for the ones adding full matrices. The
   * global locations <tt>row_index</tt> and <tt>col_indices</tt> are
   * translated into locations in this block and ExcBlockIndexMismatch is
   * thrown, if the global index does not point into the block referred to by
   * #row and #column.
   *
   * @todo <tt>elide_zero_values</tt> is currently ignored.
   *
   * The optional parameter <tt>elide_zero_values</tt> can be used to specify
   * whether zero values should be added anyway or these should be filtered
   * away and only non-zero data is added. The default value is <tt>true</tt>,
   * i.e., zero values won't be added into the matrix.
   */
  template <typename number>
  void
  add(const size_type               row_index,
      const std::vector<size_type> &col_indices,
      const std::vector<number> &   values,
      const bool                    elide_zero_values = true);

  /**
   * Add an array of values given by <tt>values</tt> in the given global
   * matrix row at columns specified by col_indices in the sparse matrix.
   *
   * The optional parameter <tt>elide_zero_values</tt> can be used to specify
   * whether zero values should be added anyway or these should be filtered
   * away and only non-zero data is added. The default value is <tt>true</tt>,
   * i.e., zero values won't be added into the matrix.
   */
  template <typename number>
  void
  add(const size_type  row,
      const size_type  n_cols,
      const size_type *col_indices,
      const number *   values,
      const bool       elide_zero_values      = true,
      const bool       col_indices_are_sorted = false);

  /**
   * Matrix-vector-multiplication, forwarding to the same function in
   * MatrixType. No index computations are done, thus, the vectors need to
   * have sizes matching #matrix.
   */
  template <class VectorType>
  void
  vmult(VectorType &w, const VectorType &v) const;

  /**
   * Matrix-vector-multiplication, forwarding to the same function in
   * MatrixType. No index computations are done, thus, the vectors need to
   * have sizes matching #matrix.
   */
  template <class VectorType>
  void
  vmult_add(VectorType &w, const VectorType &v) const;

  /**
   * Matrix-vector-multiplication, forwarding to the same function in
   * MatrixType. No index computations are done, thus, the vectors need to
   * have sizes matching #matrix.
   */
  template <class VectorType>
  void
  Tvmult(VectorType &w, const VectorType &v) const;

  /**
   * Matrix-vector-multiplication, forwarding to the same function in
   * MatrixType. No index computations are done, thus, the vectors need to
   * have sizes matching #matrix.
   */
  template <class VectorType>
  void
  Tvmult_add(VectorType &w, const VectorType &v) const;

  /**
   * The memory used by this object.
   */
  std::size_t
  memory_consumption() const;

  /**
   * The block number computed from an index by using BlockIndices does not
   * match the block coordinates stored in this object.
   */
  DeclException2(ExcBlockIndexMismatch,
                 size_type,
                 size_type,
                 << "Block index " << arg1 << " does not match " << arg2);

  /**
   * Row coordinate.  This is the position of the data member matrix on the
   * global matrix.
   */
  size_type row;
  /**
   * Column coordinate.  This is the position of the data member matrix on the
   * global matrix.
   */
  size_type column;

  /**
   * The matrix itself
   */
  MatrixType matrix;

private:
  /**
   * The row BlockIndices of the whole system. Using row(), this allows us to
   * find the index of the first row degree of freedom for this block.
   */
  BlockIndices row_indices;
  /**
   * The column BlockIndices of the whole system. Using column(), this allows
   * us to find the index of the first column degree of freedom for this
   * block.
   */
  BlockIndices column_indices;

  template <class OTHER_MatrixType>
  friend void
  dealii::internal::reinit(MatrixBlock<OTHER_MatrixType> &,
                           const BlockSparsityPattern &);

  template <typename number>
  friend void
  internal::reinit(MatrixBlock<dealii::SparseMatrix<number>> &v,
                   const BlockSparsityPattern &               p);
};


/**
 * A vector of MatrixBlock, which is implemented using shared pointers, in
 * order to allow for copying and rearranging. Each matrix block can be
 * identified by name.
 *
 * @relatesalso MatrixBlock
 * @ingroup vector_valued
 * @author Baerbel Janssen, Guido Kanschat, 2010
 */
template <typename MatrixType>
class MatrixBlockVector : private AnyData
{
public:
  /**
   * Declare type for container size.
   */
  using size_type = types::global_dof_index;

  /**
   * The type of object stored.
   */
  using value_type = MatrixBlock<MatrixType>;

  /**
   * The pointer type used for storing the objects. We use a shard pointer,
   * such that they get deleted automatically when not used anymore.
   */
  using ptr_type = std::shared_ptr<value_type>;

  /**
   * Add a new matrix block at the position <tt>(row,column)</tt> in the block
   * system.
   */
  void
  add(size_type row, size_type column, const std::string &name);

  /**
   * For matrices using a SparsityPattern, this function reinitializes each
   * matrix in the vector with the correct pattern from the block system.
   */
  void
  reinit(const BlockSparsityPattern &sparsity);

  /**
   * Clear the object.
   *
   * Since often only clearing of the individual matrices is desired, but not
   * removing the blocks themselves, there is an optional argument. If the
   * argument is missing or @p false, all matrices will be empty, but the size
   * of this object and the block positions will not change. If @p
   * really_clean is @p true, then the object will contain no blocks at the
   * end.
   */
  void
  clear(bool really_clean = false);

  /**
   * The memory used by this object.
   */
  std::size_t
  memory_consumption() const;

  /**
   * Access a constant reference to the block at position <i>i</i>.
   */
  const value_type &
  block(size_type i) const;

  /**
   * Access a reference to the block at position <i>i</i>.
   */
  value_type &
  block(size_type i);

  /**
   * Access the matrix at position <i>i</i> for read and write access.
   */
  MatrixType &
  matrix(size_type i);

  /**
   * import functions from private base class
   */
  using AnyData::name;
  using AnyData::size;
  using AnyData::subscribe;
  using AnyData::unsubscribe;
};


/**
 * A vector of MGLevelObject<MatrixBlock>, which is implemented using shared
 * pointers, in order to allow for copying and rearranging. Each matrix block
 * can be identified by name.
 *
 * @relatesalso MatrixBlock
 * @ingroup vector_valued
 * @author Baerbel Janssen, Guido Kanschat, 2010
 */
template <typename MatrixType>
class MGMatrixBlockVector : public Subscriptor
{
public:
  /**
   * Declare type for container size.
   */
  using size_type = types::global_dof_index;

  /**
   * The type of object stored.
   */
  using value_type = MGLevelObject<MatrixBlock<MatrixType>>;
  /**
   * Constructor, determining which matrices should be stored.
   *
   * If <tt>edge_matrices</tt> is true, then objects for edge matrices for
   * discretizations with degrees of freedom on faces are allocated.
   *
   * If <tt>edge_flux_matrices</tt> is true, then objects for DG fluxes on the
   * refinement edge are allocated.
   */
  MGMatrixBlockVector(const bool edge_matrices      = false,
                      const bool edge_flux_matrices = false);

  /**
   * The number of blocks.
   */
  unsigned int
  size() const;

  /**
   * Add a new matrix block at the position <tt>(row,column)</tt> in the block
   * system. The third argument allows to give the matrix a name for later
   * identification.
   */
  void
  add(size_type row, size_type column, const std::string &name);

  /**
   * For matrices using a SparsityPattern, this function reinitializes each
   * matrix in the vector with the correct pattern from the block system.
   *
   * This function reinitializes the level matrices.
   */
  void
  reinit_matrix(const MGLevelObject<BlockSparsityPattern> &sparsity);
  /**
   * For matrices using a SparsityPattern, this function reinitializes each
   * matrix in the vector with the correct pattern from the block system.
   *
   * This function reinitializes the matrices for degrees of freedom on the
   * refinement edge.
   */
  void
  reinit_edge(const MGLevelObject<BlockSparsityPattern> &sparsity);
  /**
   * For matrices using a SparsityPattern, this function reinitializes each
   * matrix in the vector with the correct pattern from the block system.
   *
   * This function reinitializes the flux matrices over the refinement edge.
   */
  void
  reinit_edge_flux(const MGLevelObject<BlockSparsityPattern> &sparsity);

  /**
   * Clear the object.
   *
   * Since often only clearing of the individual matrices is desired, but not
   * removing the blocks themselves, there is an optional argument. If the
   * argument is missing or @p false, all matrices will be empty, but the size
   * of this object and the block positions will not change. If @p
   * really_clean is @p true, then the object will contain no blocks at the
   * end.
   */
  void
  clear(bool really_clean = false);

  /**
   * Access a constant reference to the matrix block at position <i>i</i>.
   */
  const value_type &
  block(size_type i) const;

  /**
   * Access a reference to the matrix block at position <i>i</i>.
   */
  value_type &
  block(size_type i);

  /**
   * Access a constant reference to the edge matrix block at position
   * <i>i</i>.
   */
  const value_type &
  block_in(size_type i) const;

  /**
   * Access a reference to the edge matrix block at position <i>i</i>.
   */
  value_type &
  block_in(size_type i);

  /**
   * Access a constant reference to the edge matrix block at position
   * <i>i</i>.
   */
  const value_type &
  block_out(size_type i) const;

  /**
   * Access a reference to the edge matrix block at position <i>i</i>.
   */
  value_type &
  block_out(size_type i);

  /**
   * Access a constant reference to the  edge flux matrix block at position
   * <i>i</i>.
   */
  const value_type &
  block_up(size_type i) const;

  /**
   * Access a reference to the  edge flux matrix block at position <i>i</i>.
   */
  value_type &
  block_up(size_type i);

  /**
   * Access a constant reference to the  edge flux matrix block at position
   * <i>i</i>.
   */
  const value_type &
  block_down(size_type i) const;

  /**
   * Access a reference to the edge flux matrix block at position <i>i</i>.
   */
  value_type &
  block_down(size_type i);

  /**
   * The memory used by this object.
   */
  std::size_t
  memory_consumption() const;

private:
  /// Clear one of the matrix objects
  void
  clear_object(AnyData &);

  /// Flag for storing matrices_in and matrices_out
  const bool edge_matrices;

  /// Flag for storing flux_matrices_up and flux_matrices_down
  const bool edge_flux_matrices;

  /// The level matrices
  AnyData matrices;
  /// The matrix from the interior of a level to the refinement edge
  AnyData matrices_in;
  /// The matrix from the refinement edge to the interior of a level
  AnyData matrices_out;
  /// The DG flux from a level to the lower level
  AnyData flux_matrices_down;
  /// The DG flux from the lower level to a level
  AnyData flux_matrices_up;
};


//----------------------------------------------------------------------//

namespace internal
{
  template <typename MatrixType>
  void
  reinit(MatrixBlock<MatrixType> &v, const BlockSparsityPattern &p)
  {
    v.row_indices    = p.get_row_indices();
    v.column_indices = p.get_column_indices();
  }


  template <typename number>
  void
  reinit(MatrixBlock<dealii::SparseMatrix<number>> &v,
         const BlockSparsityPattern &               p)
  {
    v.row_indices    = p.get_row_indices();
    v.column_indices = p.get_column_indices();
    v.matrix.reinit(p.block(v.row, v.column));
  }
} // namespace internal


template <typename MatrixType>
inline MatrixBlock<MatrixType>::MatrixBlock()
  : row(numbers::invalid_size_type)
  , column(numbers::invalid_size_type)
{}


template <typename MatrixType>
inline MatrixBlock<MatrixType>::MatrixBlock(size_type i, size_type j)
  : row(i)
  , column(j)
{}


template <typename MatrixType>
inline void
MatrixBlock<MatrixType>::reinit(const BlockSparsityPattern &sparsity)
{
  internal::reinit(*this, sparsity);
}


template <typename MatrixType>
inline MatrixBlock<MatrixType>::operator MatrixType &()
{
  return matrix;
}


template <typename MatrixType>
inline MatrixBlock<MatrixType>::operator const MatrixType &() const
{
  return matrix;
}


template <typename MatrixType>
inline void
MatrixBlock<MatrixType>::add(const size_type                       gi,
                             const size_type                       gj,
                             const typename MatrixType::value_type value)
{
  Assert(row_indices.size() != 0, ExcNotInitialized());
  Assert(column_indices.size() != 0, ExcNotInitialized());

  const std::pair<unsigned int, size_type> bi = row_indices.global_to_local(gi);
  const std::pair<unsigned int, size_type> bj =
    column_indices.global_to_local(gj);

  Assert(bi.first == row, ExcBlockIndexMismatch(bi.first, row));
  Assert(bj.first == column, ExcBlockIndexMismatch(bj.first, column));

  matrix.add(bi.second, bj.second, value);
}


template <typename MatrixType>
template <typename number>
inline void
MatrixBlock<MatrixType>::add(const std::vector<size_type> &r_indices,
                             const std::vector<size_type> &c_indices,
                             const FullMatrix<number> &    values,
                             const bool                    elide_zero_values)
{
  Assert(row_indices.size() != 0, ExcNotInitialized());
  Assert(column_indices.size() != 0, ExcNotInitialized());

  AssertDimension(r_indices.size(), values.m());
  AssertDimension(c_indices.size(), values.n());

  for (size_type i = 0; i < row_indices.size(); ++i)
    add(r_indices[i],
        c_indices.size(),
        c_indices.data(),
        &values(i, 0),
        elide_zero_values);
}


template <typename MatrixType>
template <typename number>
inline void
MatrixBlock<MatrixType>::add(const size_type  b_row,
                             const size_type  n_cols,
                             const size_type *col_indices,
                             const number *   values,
                             const bool,
                             const bool)
{
  Assert(row_indices.size() != 0, ExcNotInitialized());
  Assert(column_indices.size() != 0, ExcNotInitialized());

  const std::pair<unsigned int, size_type> bi =
    row_indices.global_to_local(b_row);

  // In debug mode, we check whether
  // all indices are in the correct
  // block.

  // Actually, for the time being, we
  // leave it at this. While it may
  // not be the most efficient way,
  // it is at least thread safe.
  //#ifdef DEBUG
  Assert(bi.first == row, ExcBlockIndexMismatch(bi.first, row));

  for (size_type j = 0; j < n_cols; ++j)
    {
      const std::pair<unsigned int, size_type> bj =
        column_indices.global_to_local(col_indices[j]);
      Assert(bj.first == column, ExcBlockIndexMismatch(bj.first, column));

      matrix.add(bi.second, bj.second, values[j]);
    }
  //#endif
}


template <typename MatrixType>
template <typename number>
inline void
MatrixBlock<MatrixType>::add(const std::vector<size_type> &indices,
                             const FullMatrix<number> &    values,
                             const bool                    elide_zero_values)
{
  Assert(row_indices.size() != 0, ExcNotInitialized());
  Assert(column_indices.size() != 0, ExcNotInitialized());

  AssertDimension(indices.size(), values.m());
  Assert(values.n() == values.m(), ExcNotQuadratic());

  for (size_type i = 0; i < indices.size(); ++i)
    add(indices[i],
        indices.size(),
        indices.data(),
        &values(i, 0),
        elide_zero_values);
}



template <typename MatrixType>
template <typename number>
inline void
MatrixBlock<MatrixType>::add(const size_type               row,
                             const std::vector<size_type> &col_indices,
                             const std::vector<number> &   values,
                             const bool                    elide_zero_values)
{
  Assert(row_indices.size() != 0, ExcNotInitialized());
  Assert(column_indices.size() != 0, ExcNotInitialized());

  AssertDimension(col_indices.size(), values.size());
  add(row,
      col_indices.size(),
      col_indices.data(),
      values.data(),
      elide_zero_values);
}


template <typename MatrixType>
template <class VectorType>
inline void
MatrixBlock<MatrixType>::vmult(VectorType &w, const VectorType &v) const
{
  matrix.vmult(w, v);
}


template <typename MatrixType>
template <class VectorType>
inline void
MatrixBlock<MatrixType>::vmult_add(VectorType &w, const VectorType &v) const
{
  matrix.vmult_add(w, v);
}


template <typename MatrixType>
template <class VectorType>
inline void
MatrixBlock<MatrixType>::Tvmult(VectorType &w, const VectorType &v) const
{
  matrix.Tvmult(w, v);
}


template <typename MatrixType>
template <class VectorType>
inline void
MatrixBlock<MatrixType>::Tvmult_add(VectorType &w, const VectorType &v) const
{
  matrix.Tvmult_add(w, v);
}


template <typename MatrixType>
inline std::size_t
MatrixBlock<MatrixType>::memory_consumption() const
{
  return (sizeof(*this) + MemoryConsumption::memory_consumption(matrix) -
          sizeof(matrix));
}

//----------------------------------------------------------------------//

template <typename MatrixType>
inline void
MatrixBlockVector<MatrixType>::add(size_type          row,
                                   size_type          column,
                                   const std::string &name)
{
  ptr_type p(new value_type(row, column));
  AnyData::add(p, name);
}


template <typename MatrixType>
inline void
MatrixBlockVector<MatrixType>::reinit(const BlockSparsityPattern &sparsity)
{
  for (size_type i = 0; i < this->size(); ++i)
    {
      block(i).reinit(sparsity);
    }
}


template <typename MatrixType>
inline void
MatrixBlockVector<MatrixType>::clear(bool really_clean)
{
  if (really_clean)
    {
      Assert(false, ExcNotImplemented());
    }
  else
    {
      for (size_type i = 0; i < this->size(); ++i)
        matrix(i).clear();
    }
}



template <typename MatrixType>
inline const MatrixBlock<MatrixType> &
MatrixBlockVector<MatrixType>::block(size_type i) const
{
  return *this->read<ptr_type>(i);
}


template <typename MatrixType>
inline MatrixBlock<MatrixType> &
MatrixBlockVector<MatrixType>::block(size_type i)
{
  return *this->entry<ptr_type>(i);
}


template <typename MatrixType>
inline MatrixType &
MatrixBlockVector<MatrixType>::matrix(size_type i)
{
  return this->entry<ptr_type>(i)->matrix;
}



//----------------------------------------------------------------------//

template <typename MatrixType>
inline MGMatrixBlockVector<MatrixType>::MGMatrixBlockVector(const bool e,
                                                            const bool f)
  : edge_matrices(e)
  , edge_flux_matrices(f)
{}


template <typename MatrixType>
inline unsigned int
MGMatrixBlockVector<MatrixType>::size() const
{
  return matrices.size();
}


template <typename MatrixType>
inline void
MGMatrixBlockVector<MatrixType>::add(size_type          row,
                                     size_type          column,
                                     const std::string &name)
{
  MGLevelObject<MatrixBlock<MatrixType>> p(0, 1);
  p[0].row    = row;
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


template <typename MatrixType>
inline const MGLevelObject<MatrixBlock<MatrixType>> &
MGMatrixBlockVector<MatrixType>::block(size_type i) const
{
  return *matrices.read<const MGLevelObject<MatrixType> *>(i);
}


template <typename MatrixType>
inline MGLevelObject<MatrixBlock<MatrixType>> &
MGMatrixBlockVector<MatrixType>::block(size_type i)
{
  return *matrices.entry<MGLevelObject<MatrixType> *>(i);
}


template <typename MatrixType>
inline const MGLevelObject<MatrixBlock<MatrixType>> &
MGMatrixBlockVector<MatrixType>::block_in(size_type i) const
{
  return *matrices_in.read<const MGLevelObject<MatrixType> *>(i);
}


template <typename MatrixType>
inline MGLevelObject<MatrixBlock<MatrixType>> &
MGMatrixBlockVector<MatrixType>::block_in(size_type i)
{
  return *matrices_in.entry<MGLevelObject<MatrixType> *>(i);
}


template <typename MatrixType>
inline const MGLevelObject<MatrixBlock<MatrixType>> &
MGMatrixBlockVector<MatrixType>::block_out(size_type i) const
{
  return *matrices_out.read<const MGLevelObject<MatrixType> *>(i);
}


template <typename MatrixType>
inline MGLevelObject<MatrixBlock<MatrixType>> &
MGMatrixBlockVector<MatrixType>::block_out(size_type i)
{
  return *matrices_out.entry<MGLevelObject<MatrixType> *>(i);
}


template <typename MatrixType>
inline const MGLevelObject<MatrixBlock<MatrixType>> &
MGMatrixBlockVector<MatrixType>::block_up(size_type i) const
{
  return *flux_matrices_up.read<const MGLevelObject<MatrixType> *>(i);
}


template <typename MatrixType>
inline MGLevelObject<MatrixBlock<MatrixType>> &
MGMatrixBlockVector<MatrixType>::block_up(size_type i)
{
  return *flux_matrices_up.entry<MGLevelObject<MatrixType> *>(i);
}


template <typename MatrixType>
inline const MGLevelObject<MatrixBlock<MatrixType>> &
MGMatrixBlockVector<MatrixType>::block_down(size_type i) const
{
  return *flux_matrices_down.read<const MGLevelObject<MatrixType> *>(i);
}


template <typename MatrixType>
inline MGLevelObject<MatrixBlock<MatrixType>> &
MGMatrixBlockVector<MatrixType>::block_down(size_type i)
{
  return *flux_matrices_down.entry<MGLevelObject<MatrixType> *>(i);
}


template <typename MatrixType>
inline void
MGMatrixBlockVector<MatrixType>::reinit_matrix(
  const MGLevelObject<BlockSparsityPattern> &sparsity)
{
  for (size_type i = 0; i < this->size(); ++i)
    {
      MGLevelObject<MatrixBlock<MatrixType>> &o   = block(i);
      const size_type                         row = o[o.min_level()].row;
      const size_type                         col = o[o.min_level()].column;

      o.resize(sparsity.min_level(), sparsity.max_level());
      for (size_type level = o.min_level(); level <= o.max_level(); ++level)
        {
          o[level].row    = row;
          o[level].column = col;
          internal::reinit(o[level], sparsity[level]);
        }
    }
}


template <typename MatrixType>
inline void
MGMatrixBlockVector<MatrixType>::reinit_edge(
  const MGLevelObject<BlockSparsityPattern> &sparsity)
{
  for (size_type i = 0; i < this->size(); ++i)
    {
      MGLevelObject<MatrixBlock<MatrixType>> &o   = block(i);
      const size_type                         row = o[o.min_level()].row;
      const size_type                         col = o[o.min_level()].column;

      block_in(i).resize(sparsity.min_level(), sparsity.max_level());
      block_out(i).resize(sparsity.min_level(), sparsity.max_level());
      for (size_type level = o.min_level(); level <= o.max_level(); ++level)
        {
          block_in(i)[level].row    = row;
          block_in(i)[level].column = col;
          internal::reinit(block_in(i)[level], sparsity[level]);
          block_out(i)[level].row    = row;
          block_out(i)[level].column = col;
          internal::reinit(block_out(i)[level], sparsity[level]);
        }
    }
}


template <typename MatrixType>
inline void
MGMatrixBlockVector<MatrixType>::reinit_edge_flux(
  const MGLevelObject<BlockSparsityPattern> &sparsity)
{
  for (size_type i = 0; i < this->size(); ++i)
    {
      MGLevelObject<MatrixBlock<MatrixType>> &o   = block(i);
      const size_type                         row = o[o.min_level()].row;
      const size_type                         col = o[o.min_level()].column;

      block_up(i).resize(sparsity.min_level(), sparsity.max_level());
      block_down(i).resize(sparsity.min_level(), sparsity.max_level());
      for (size_type level = o.min_level(); level <= o.max_level(); ++level)
        {
          block_up(i)[level].row    = row;
          block_up(i)[level].column = col;
          internal::reinit(block_up(i)[level], sparsity[level]);
          block_down(i)[level].row    = row;
          block_down(i)[level].column = col;
          internal::reinit(block_down(i)[level], sparsity[level]);
        }
    }
}


template <typename MatrixType>
inline void
MGMatrixBlockVector<MatrixType>::clear_object(AnyData &mo)
{
  for (size_type i = 0; i < mo.size(); ++i)
    {
      MGLevelObject<MatrixBlock<MatrixType>> &o =
        mo.entry<MGLevelObject<MatrixType> *>(i);
      for (size_type level = o.min_level(); level <= o.max_level(); ++level)
        o[level].matrix.clear();
    }
}


template <typename MatrixType>
inline void
MGMatrixBlockVector<MatrixType>::clear(bool really_clean)
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

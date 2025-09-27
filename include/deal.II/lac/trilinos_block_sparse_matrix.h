// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_trilinos_block_sparse_matrix_h
#define dealii_trilinos_block_sparse_matrix_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/template_constraints.h>

#  include <deal.II/lac/block_matrix_base.h>
#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/full_matrix.h>
#  include <deal.II/lac/trilinos_parallel_block_vector.h>
#  include <deal.II/lac/trilinos_sparse_matrix.h>

#  include <cmath>

DEAL_II_NAMESPACE_OPEN

// forward declarations
#  ifndef DOXYGEN
class BlockSparsityPattern;
template <typename number>
class BlockSparseMatrix;
#  endif

namespace TrilinosWrappers
{
  /**
   * @addtogroup TrilinosWrappers
   * @{
   */

  /**
   * Blocked sparse matrix based on the TrilinosWrappers::SparseMatrix class.
   * This class implements the functions that are specific to the Trilinos
   * SparseMatrix base objects for a blocked sparse matrix, and leaves the
   * actual work relaying most of the calls to the individual blocks to the
   * functions implemented in the base class. See there also for a description
   * of when this class is useful.
   *
   * In contrast to the deal.II-type SparseMatrix class, the Trilinos matrices
   * do not have external objects for the sparsity patterns. Thus, one does
   * not determine the size of the individual blocks of a block matrix of this
   * type by attaching a block sparsity pattern, but by calling reinit() to
   * set the number of blocks and then by setting the size of each block
   * separately. In order to fix the data structures of the block matrix, it
   * is then necessary to let it know that we have changed the sizes of the
   * underlying matrices. For this, one has to call the collect_sizes()
   * function, for much the same reason as is documented with the
   * BlockSparsityPattern class.
   *
   * @ingroup Matrix1
   * @see @ref GlossBlockLA "Block (linear algebra)"
   */
  class BlockSparseMatrix : public BlockMatrixBase<SparseMatrix>
  {
  public:
    /**
     * Typedef the base class for simpler access to its own alias.
     */
    using BaseClass = BlockMatrixBase<SparseMatrix>;

    /**
     * Typedef the type of the underlying matrix.
     */
    using BlockType = BaseClass::BlockType;

    /**
     * Import the alias from the base class.
     */
    using value_type      = BaseClass::value_type;
    using pointer         = BaseClass::pointer;
    using const_pointer   = BaseClass::const_pointer;
    using reference       = BaseClass::reference;
    using const_reference = BaseClass::const_reference;
    using size_type       = BaseClass::size_type;
    using iterator        = BaseClass::iterator;
    using const_iterator  = BaseClass::const_iterator;

    /**
     * Constructor; initializes the matrix to be empty, without any structure,
     * i.e.  the matrix is not usable at all. This constructor is therefore
     * only useful for matrices which are members of a class. All other
     * matrices should be created at a point in the data flow where all
     * necessary information is available.
     *
     * You have to initialize the matrix before usage with
     * reinit(BlockSparsityPattern). The number of blocks per row and column
     * are then determined by that function.
     */
    BlockSparseMatrix() = default;

    /**
     * Destructor.
     */
    ~BlockSparseMatrix() override;

    /**
     * Pseudo copy operator only copying empty objects. The sizes of the block
     * matrices need to be the same.
     */
    BlockSparseMatrix &
    operator=(const BlockSparseMatrix &) = default;

    /**
     * This operator assigns a scalar to a matrix. Since this does usually not
     * make much sense (should we set all matrix entries to this value? Only
     * the nonzero entries of the sparsity pattern?), this operation is only
     * allowed if the actual value to be assigned is zero. This operator only
     * exists to allow for the obvious notation <tt>matrix=0</tt>, which sets
     * all elements of the matrix to zero, but keep the sparsity pattern
     * previously used.
     */
    BlockSparseMatrix &
    operator=(const double d);

    /**
     * Resize the matrix, by setting the number of block rows and columns.
     * This deletes all blocks and replaces them with uninitialized ones, i.e.
     * ones for which also the sizes are not yet set. You have to do that by
     * calling the @p reinit functions of the blocks themselves. Do not forget
     * to call collect_sizes() after that on this object.
     *
     * The reason that you have to set sizes of the blocks yourself is that
     * the sizes may be varying, the maximum number of elements per row may be
     * varying, etc. It is simpler not to reproduce the interface of the @p
     * SparsityPattern class here but rather let the user call whatever
     * function they desire.
     */
    void
    reinit(const size_type n_block_rows, const size_type n_block_columns);

    /**
     * Resize the matrix, by using an array of index sets to determine the
     * %parallel distribution of the individual matrices. This function
     * assumes that a quadratic block matrix is generated.
     */
    template <typename BlockSparsityPatternType>
    void
    reinit(const std::vector<IndexSet>    &input_maps,
           const BlockSparsityPatternType &block_sparsity_pattern,
           const MPI_Comm                  communicator  = MPI_COMM_WORLD,
           const bool                      exchange_data = false);

    /**
     * Resize the matrix and initialize it by the given sparsity pattern.
     * Since no distribution map is given, the result is a block matrix for
     * which all elements are stored locally.
     */
    template <typename BlockSparsityPatternType>
    void
    reinit(const BlockSparsityPatternType &block_sparsity_pattern);

    /**
     * This function initializes the Trilinos matrix using the deal.II sparse
     * matrix and the entries stored therein. It uses a threshold to copy only
     * elements whose modulus is larger than the threshold (so zeros in the
     * deal.II matrix can be filtered away).
     */
    void
    reinit(
      const std::vector<IndexSet>               &parallel_partitioning,
      const ::dealii::BlockSparseMatrix<double> &dealii_block_sparse_matrix,
      const MPI_Comm                             communicator = MPI_COMM_WORLD,
      const double                               drop_tolerance = 1e-13);

    /**
     * This function initializes the Trilinos matrix using the deal.II sparse
     * matrix and the entries stored therein. It uses a threshold to copy only
     * elements whose modulus is larger than the threshold (so zeros in the
     * deal.II matrix can be filtered away). Since no Epetra_Map is given, all
     * the elements will be locally stored.
     */
    void
    reinit(const ::dealii::BlockSparseMatrix<double> &deal_ii_sparse_matrix,
           const double                               drop_tolerance = 1e-13);

    /**
     * Return the state of the matrix, i.e., whether compress() needs to be
     * called after an operation requiring data exchange. Does only return
     * non-true values when used in <tt>debug</tt> mode, since it is quite
     * expensive to keep track of all operations that lead to the need for
     * compress().
     */
    bool
    is_compressed() const;

    /**
     * This function collects the sizes of the sub-objects and stores them in
     * internal arrays, in order to be able to relay global indices into the
     * matrix to indices into the subobjects. You *must* call this function
     * each time after you have changed the size of the sub-objects. Note that
     * this is a @ref GlossCollectiveOperation "collective operation", i.e., it needs to be called on all MPI
     * processes. This command internally calls the method
     * <tt>compress()</tt>, so you don't need to call that function in case
     * you use <tt>collect_sizes()</tt>.
     */
    void
    collect_sizes();

    /**
     * Return the total number of nonzero elements of this matrix (summed
     * over all MPI processes).
     */
    std::uint64_t
    n_nonzero_elements() const;

    /**
     * Return the underlying MPI communicator.
     */
    MPI_Comm
    get_mpi_communicator() const;

    /**
     * Return the partitioning of the domain space for the individual blocks of
     * this matrix, i.e., the partitioning of the block vectors this matrix has
     * to be multiplied with.
     */
    std::vector<IndexSet>
    locally_owned_domain_indices() const;

    /**
     * Return the partitioning of the range space for the individual blocks of
     * this matrix, i.e., the partitioning of the block vectors that result
     * from matrix-vector products.
     */
    std::vector<IndexSet>
    locally_owned_range_indices() const;

    /**
     * Matrix-vector multiplication: let $dst = M*src$ with $M$ being this
     * matrix. The vector types can be block vectors or non-block vectors
     * (only if the matrix has only one row or column, respectively), and need
     * to define TrilinosWrappers::SparseMatrix::vmult.
     */
    template <typename VectorType1, typename VectorType2>
    void
    vmult(VectorType1 &dst, const VectorType2 &src) const;

    /**
     * Matrix-vector multiplication: let $dst = M^T*src$ with $M$ being this
     * matrix. This function does the same as vmult() but takes the transposed
     * matrix.
     */
    template <typename VectorType1, typename VectorType2>
    void
    Tvmult(VectorType1 &dst, const VectorType2 &src) const;

    /**
     * Compute the residual of an equation <i>Mx=b</i>, where the residual is
     * defined to be <i>r=b-Mx</i>. Write the residual into @p dst. The
     * <i>l<sub>2</sub></i> norm of the residual vector is returned.
     *
     * Source <i>x</i> and destination <i>dst</i> must not be the same vector.
     *
     * Note that both vectors have to be distributed vectors generated using
     * the same Map as was used for the matrix.
     *
     * This function only applicable if the matrix only has one block row.
     */
    TrilinosScalar
    residual(MPI::BlockVector       &dst,
             const MPI::BlockVector &x,
             const MPI::BlockVector &b) const;

    /**
     * Compute the residual of an equation <i>Mx=b</i>, where the residual is
     * defined to be <i>r=b-Mx</i>. Write the residual into @p dst. The
     * <i>l<sub>2</sub></i> norm of the residual vector is returned.
     *
     * This function is only applicable if the matrix only has one block row.
     */
    TrilinosScalar
    residual(MPI::BlockVector       &dst,
             const MPI::Vector      &x,
             const MPI::BlockVector &b) const;

    /**
     * Compute the residual of an equation <i>Mx=b</i>, where the residual is
     * defined to be <i>r=b-Mx</i>. Write the residual into @p dst. The
     * <i>l<sub>2</sub></i> norm of the residual vector is returned.
     *
     * This function is only applicable if the matrix only has one block column.
     */
    TrilinosScalar
    residual(MPI::Vector            &dst,
             const MPI::BlockVector &x,
             const MPI::Vector      &b) const;

    /**
     * Compute the residual of an equation <i>Mx=b</i>, where the residual is
     * defined to be <i>r=b-Mx</i>. Write the residual into @p dst. The
     * <i>l<sub>2</sub></i> norm of the residual vector is returned.
     *
     * This function is only applicable if the matrix only has one block.
     */
    TrilinosScalar
    residual(MPI::Vector       &dst,
             const MPI::Vector &x,
             const MPI::Vector &b) const;

    /**
     * Make the clear() function in the base class visible, though it is
     * protected.
     */
    using BlockMatrixBase<SparseMatrix>::clear;

    /**
     * @addtogroup Exceptions
     * @{
     */

    /**
     * Exception
     */
    DeclException4(ExcIncompatibleRowNumbers,
                   int,
                   int,
                   int,
                   int,
                   << "The blocks [" << arg1 << ',' << arg2 << "] and [" << arg3
                   << ',' << arg4 << "] have differing row numbers.");

    /**
     * Exception
     */
    DeclException4(ExcIncompatibleColNumbers,
                   int,
                   int,
                   int,
                   int,
                   << "The blocks [" << arg1 << ',' << arg2 << "] and [" << arg3
                   << ',' << arg4 << "] have differing column numbers.");
    /** @} */

  private:
    /**
     * Internal version of (T)vmult with two block vectors
     */
    template <typename VectorType1, typename VectorType2>
    void
    vmult(VectorType1       &dst,
          const VectorType2 &src,
          const bool         transpose,
          const std::bool_constant<true>,
          const std::bool_constant<true>) const;

    /**
     * Internal version of (T)vmult where the source vector is a block vector
     * but the destination vector is a non-block vector
     */
    template <typename VectorType1, typename VectorType2>
    void
    vmult(VectorType1       &dst,
          const VectorType2 &src,
          const bool         transpose,
          const std::bool_constant<false>,
          const std::bool_constant<true>) const;

    /**
     * Internal version of (T)vmult where the source vector is a non-block
     * vector but the destination vector is a block vector
     */
    template <typename VectorType1, typename VectorType2>
    void
    vmult(VectorType1       &dst,
          const VectorType2 &src,
          const bool         transpose,
          const std::bool_constant<true>,
          const std::bool_constant<false>) const;

    /**
     * Internal version of (T)vmult where both source vector and the
     * destination vector are non-block vectors (only defined if the matrix
     * consists of only one block)
     */
    template <typename VectorType1, typename VectorType2>
    void
    vmult(VectorType1       &dst,
          const VectorType2 &src,
          const bool         transpose,
          const std::bool_constant<false>,
          const std::bool_constant<false>) const;
  };



  /** @} */

  // ------------- inline and template functions -----------------



  inline BlockSparseMatrix &
  BlockSparseMatrix::operator=(const double d)
  {
    Assert(d == 0, ExcScalarAssignmentOnlyForZeroValue());

    for (size_type r = 0; r < this->n_block_rows(); ++r)
      for (size_type c = 0; c < this->n_block_cols(); ++c)
        this->block(r, c) = d;

    return *this;
  }



  inline bool
  BlockSparseMatrix::is_compressed() const
  {
    bool compressed = true;
    for (size_type row = 0; row < n_block_rows(); ++row)
      for (size_type col = 0; col < n_block_cols(); ++col)
        if (block(row, col).is_compressed() == false)
          {
            compressed = false;
            break;
          }

    return compressed;
  }



  template <typename VectorType1, typename VectorType2>
  inline void
  BlockSparseMatrix::vmult(VectorType1 &dst, const VectorType2 &src) const
  {
    vmult(dst,
          src,
          false,
          std::bool_constant<IsBlockVector<VectorType1>::value>(),
          std::bool_constant<IsBlockVector<VectorType2>::value>());
  }



  template <typename VectorType1, typename VectorType2>
  inline void
  BlockSparseMatrix::Tvmult(VectorType1 &dst, const VectorType2 &src) const
  {
    vmult(dst,
          src,
          true,
          std::bool_constant<IsBlockVector<VectorType1>::value>(),
          std::bool_constant<IsBlockVector<VectorType2>::value>());
  }



  template <typename VectorType1, typename VectorType2>
  inline void
  BlockSparseMatrix::vmult(VectorType1       &dst,
                           const VectorType2 &src,
                           const bool         transpose,
                           std::bool_constant<true>,
                           std::bool_constant<true>) const
  {
    if (transpose == true)
      BaseClass::Tvmult_block_block(dst, src);
    else
      BaseClass::vmult_block_block(dst, src);
  }



  template <typename VectorType1, typename VectorType2>
  inline void
  BlockSparseMatrix::vmult(VectorType1       &dst,
                           const VectorType2 &src,
                           const bool         transpose,
                           std::bool_constant<false>,
                           std::bool_constant<true>) const
  {
    if (transpose == true)
      BaseClass::Tvmult_nonblock_block(dst, src);
    else
      BaseClass::vmult_nonblock_block(dst, src);
  }



  template <typename VectorType1, typename VectorType2>
  inline void
  BlockSparseMatrix::vmult(VectorType1       &dst,
                           const VectorType2 &src,
                           const bool         transpose,
                           std::bool_constant<true>,
                           std::bool_constant<false>) const
  {
    if (transpose == true)
      BaseClass::Tvmult_block_nonblock(dst, src);
    else
      BaseClass::vmult_block_nonblock(dst, src);
  }



  template <typename VectorType1, typename VectorType2>
  inline void
  BlockSparseMatrix::vmult(VectorType1       &dst,
                           const VectorType2 &src,
                           const bool         transpose,
                           std::bool_constant<false>,
                           std::bool_constant<false>) const
  {
    if (transpose == true)
      BaseClass::Tvmult_nonblock_nonblock(dst, src);
    else
      BaseClass::vmult_nonblock_nonblock(dst, src);
  }



  inline std::vector<IndexSet>
  BlockSparseMatrix::locally_owned_domain_indices() const
  {
    Assert(this->n_block_cols() != 0, ExcNotInitialized());
    Assert(this->n_block_rows() != 0, ExcNotInitialized());

    std::vector<IndexSet> domain_indices;
    for (size_type c = 0; c < this->n_block_cols(); ++c)
      domain_indices.push_back(
        this->sub_objects[0][c]->locally_owned_domain_indices());

    return domain_indices;
  }



  inline std::vector<IndexSet>
  BlockSparseMatrix::locally_owned_range_indices() const
  {
    Assert(this->n_block_cols() != 0, ExcNotInitialized());
    Assert(this->n_block_rows() != 0, ExcNotInitialized());

    std::vector<IndexSet> range_indices;
    for (size_type r = 0; r < this->n_block_rows(); ++r)
      range_indices.push_back(
        this->sub_objects[r][0]->locally_owned_range_indices());

    return range_indices;
  }



  namespace internal
  {
    namespace BlockLinearOperatorImplementation
    {
      /**
       * This is an extension class to BlockLinearOperators for Trilinos block
       * sparse matrices.
       *
       * @note This class does very little at the moment other than to check
       * that the correct Payload type for each subblock has been chosen
       * correctly. Further extensions to the class may be necessary in the
       * future in order to add further functionality to BlockLinearOperators
       * while retaining compatibility with the Trilinos sparse matrix and
       * preconditioner classes.
       *
       *
       * @ingroup TrilinosWrappers
       */
      template <typename PayloadBlockType>
      class TrilinosBlockPayload
      {
      public:
        /**
         * Type of payload held by each subblock
         */
        using BlockType = PayloadBlockType;

        /**
         * Default constructor
         *
         * This simply checks that the payload for each block has been chosen
         * correctly (i.e. is of type TrilinosPayload). Apart from this, this
         * class does not do anything in particular and needs no special
         * configuration, we have only one generic constructor that can be
         * called under any conditions.
         */
        template <typename... Args>
        TrilinosBlockPayload(const Args &...)
        {
          static_assert(
            std::is_same_v<
              PayloadBlockType,
              internal::LinearOperatorImplementation::TrilinosPayload>,
            "TrilinosBlockPayload can only accept a payload of type TrilinosPayload.");
        }
      };

    } // namespace BlockLinearOperatorImplementation
  }   /* namespace internal */


} /* namespace TrilinosWrappers */


DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS

#endif // dealii_trilinos_block_sparse_matrix_h

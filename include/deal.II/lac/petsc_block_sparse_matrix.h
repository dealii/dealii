// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_petsc_block_sparse_matrix_h
#define dealii_petsc_block_sparse_matrix_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/block_matrix_base.h>
#  include <deal.II/lac/block_sparsity_pattern.h>
#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/petsc_block_vector.h>
#  include <deal.II/lac/petsc_sparse_matrix.h>

#  include <cmath>
#  include <cstddef>

DEAL_II_NAMESPACE_OPEN



namespace PETScWrappers
{
  namespace MPI
  {
    /**
     * @addtogroup PETScWrappers
     * @{
     */

    /**
     * Blocked sparse matrix based on the PETScWrappers::MPI::SparseMatrix
     * class. This class implements the functions that are specific to the
     * PETSc SparseMatrix base objects for a blocked sparse matrix, and leaves
     * the actual work relaying most of the calls to the individual blocks to
     * the functions implemented in the base class. See there also for a
     * description of when this class is useful.
     *
     * In contrast to the deal.II-type SparseMatrix class, the PETSc matrices
     * do not have external objects for the sparsity patterns. Thus, one does
     * not determine the size of the individual blocks of a block matrix of
     * this type by attaching a block sparsity pattern, but by calling
     * reinit() to set the number of blocks and then by setting the size of
     * each block separately. In order to fix the data structures of the block
     * matrix, it is then necessary to let it know that we have changed the
     * sizes of the underlying matrices. For this, one has to call the
     * collect_sizes() function, for much the same reason as is documented
     * with the BlockSparsityPattern class.
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
       * Constructor; initializes the matrix to be empty, without any
       * structure, i.e.  the matrix is not usable at all. This constructor is
       * therefore only useful for matrices which are members of a class. All
       * other matrices should be created at a point in the data flow where
       * all necessary information is available.
       *
       * You have to initialize the matrix before usage with
       * reinit(BlockSparsityPattern). The number of blocks per row and column
       * are then determined by that function.
       */
      BlockSparseMatrix();

      /**
       * Create a BlockSparseMatrix with a PETSc Mat that describes the entire
       * block matrix.
       * It infers the number of blocks from the Mat if it is of type MATNEST,
       * otherwise the block operator will only have a single block.
       * Internally, we always store a MATNEST matrix.
       */
      explicit BlockSparseMatrix(const Mat &);

      /**
       * Create a BlockSparseMatrix with an array of PETSc matrices.
       */
      template <std::size_t block_rows, std::size_t block_columns>
      explicit BlockSparseMatrix(
        const std::array<std::array<Mat, block_columns>, block_rows> &);

      /**
       * Destructor.
       */
      ~BlockSparseMatrix() override;

      /**
       * Pseudo copy operator only copying empty objects. The sizes of the
       * block matrices need to be the same.
       */
      BlockSparseMatrix &
      operator=(const BlockSparseMatrix &);

      /**
       * This operator assigns a scalar to a matrix. Since this does usually
       * not make much sense (should we set all matrix entries to this value?
       * Only the nonzero entries of the sparsity pattern?), this operation is
       * only allowed if the actual value to be assigned is zero. This
       * operator only exists to allow for the obvious notation
       * <tt>matrix=0</tt>, which sets all elements of the matrix to zero, but
       * keep the sparsity pattern previously used.
       */
      BlockSparseMatrix &
      operator=(const double d);

      /**
       * Resize the matrix, by setting the number of block rows and columns.
       * This deletes all blocks and replaces them with uninitialized ones,
       * i.e.  ones for which also the sizes are not yet set. You have to do
       * that by calling the @p reinit functions of the blocks themselves. Do
       * not forget to call collect_sizes() after that on this object.
       *
       * The reason that you have to set sizes of the blocks yourself is that
       * the sizes may be varying, the maximum number of elements per row may
       * be varying, etc. It is simpler not to reproduce the interface of the
       * SparsityPattern class here but rather let the user call whatever
       * function they desire.
       */
      void
      reinit(const size_type n_block_rows, const size_type n_block_columns);


      /**
       * Efficiently reinit the block matrix for a parallel computation. Only
       * the BlockSparsityPattern of the Simple type can efficiently store
       * large sparsity patterns in parallel, so this is the only supported
       * argument. The IndexSets describe the locally owned range of DoFs for
       * each block. Note that the IndexSets needs to be ascending and 1:1.
       * For a symmetric structure hand in the same vector for the first two
       * arguments.
       */
      void
      reinit(const std::vector<IndexSet>       &rows,
             const std::vector<IndexSet>       &cols,
             const BlockDynamicSparsityPattern &bdsp,
             const MPI_Comm                     com);


      /**
       * Same as above but for a symmetric structure only.
       */
      void
      reinit(const std::vector<IndexSet>       &sizes,
             const BlockDynamicSparsityPattern &bdsp,
             const MPI_Comm                     com);


      /**
       * This method associates the PETSc Mat to the instance of the class.
       * Infers the number of blocks from A if it is of type MATNEST, otherwise
       * the block operator will only have a single block.
       */
      void
      reinit(Mat A);


      /**
       * Matrix-vector multiplication: let $dst = M*src$ with $M$ being this
       * matrix.
       */
      void
      vmult(BlockVector &dst, const BlockVector &src) const;

      /**
       * Matrix-vector multiplication. Just like the previous function, but
       * only applicable if the matrix has only one block column.
       */
      void
      vmult(BlockVector &dst, const Vector &src) const;

      /**
       * Matrix-vector multiplication. Just like the previous function, but
       * only applicable if the matrix has only one block row.
       */
      void
      vmult(Vector &dst, const BlockVector &src) const;

      /**
       * Matrix-vector multiplication. Just like the previous function, but
       * only applicable if the matrix has only one block.
       */
      void
      vmult(Vector &dst, const Vector &src) const;

      /**
       * Matrix-vector multiplication: let $dst = M^T*src$ with $M$ being this
       * matrix. This function does the same as vmult() but takes the
       * transposed matrix.
       */
      void
      Tvmult(BlockVector &dst, const BlockVector &src) const;

      /**
       * Matrix-vector multiplication. Just like the previous function, but
       * only applicable if the matrix has only one block row.
       */
      void
      Tvmult(BlockVector &dst, const Vector &src) const;

      /**
       * Matrix-vector multiplication. Just like the previous function, but
       * only applicable if the matrix has only one block column.
       */
      void
      Tvmult(Vector &dst, const BlockVector &src) const;

      /**
       * Matrix-vector multiplication. Just like the previous function, but
       * only applicable if the matrix has only one block.
       */
      void
      Tvmult(Vector &dst, const Vector &src) const;

      /**
       * This function collects the sizes of the sub-objects and stores them
       * in internal arrays, in order to be able to relay global indices into
       * the matrix to indices into the subobjects. You *must* call this
       * function each time after you have changed the size of the sub-objects.
       */
      void
      collect_sizes();

      /**
       * Call the compress() function on all the subblocks of the matrix
       * and update the internal state of the PETSc object.
       *
       * See
       * @ref GlossCompress "Compressing distributed objects"
       * for more information.
       */
      void
      compress(VectorOperation::values operation);

      /**
       * Return the partitioning of the domain space of this matrix, i.e., the
       * partitioning of the vectors this matrix has to be multiplied with.
       */
      std::vector<IndexSet>
      locally_owned_domain_indices() const;

      /**
       * Return the partitioning of the range space of this matrix, i.e., the
       * partitioning of the vectors that are result from matrix-vector
       * products.
       */
      std::vector<IndexSet>
      locally_owned_range_indices() const;

      /**
       * Return the number of nonzero elements of this matrix. Actually, it
       * returns the number of entries in the sparsity pattern; if any of the
       * entries should happen to be zero, it is counted anyway.
       */
      std::uint64_t
      n_nonzero_elements() const;

      /**
       * Return the underlying MPI communicator.
       */
      MPI_Comm
      get_mpi_communicator() const;

      /**
       * Make the clear() function in the base class visible, though it is
       * protected.
       */
      using BlockMatrixBase<SparseMatrix>::clear;

      /**
       * Conversion operator to gain access to the underlying PETSc type. If you
       * do this, you cut this class off some information it may need, so this
       * conversion operator should only be used if you know what you do. In
       * particular, it should only be used for read-only operations into the
       * matrix.
       */
      operator const Mat &() const;

      /**
       * Return a reference to the underlying PETSc type. It can be used to
       * modify the underlying data, so use it only when you know what you
       * are doing.
       *
       * The PETSc Mat object returned here describes the *entire* matrix,
       * not just one of its blocks. Internally, this is done using
       * a "nested" matrix using PETSc's MATNEST object whose individual
       * blocks are the blocks of this matrix.
       */
      Mat &
      petsc_matrix();

    private:
      /**
       * A PETSc Mat object that describes the entire block matrix.
       * Internally, this is done by creating
       * a "nested" matrix using PETSc's MATNEST object whose individual
       * blocks are the blocks of this matrix.
       */
      Mat petsc_nest_matrix;

      /**
       * Utility to set up the MATNEST object.
       */
      void
      setup_nest_mat();

      /**
       * An utility method to populate empty blocks with actual objects.
       * This is needed because MATNEST supports nullptr as a block,
       * while the BlockMatrixBase class does not.
       */
      void
      create_empty_matrices_if_needed();
    };



    /** @} */

    // ------------- inline and template functions -----------------
    inline BlockSparseMatrix::BlockSparseMatrix()
      : BaseClass()
      , petsc_nest_matrix(nullptr)
    {}



    inline BlockSparseMatrix::BlockSparseMatrix(const Mat &A)
      : BlockSparseMatrix()
    {
      this->reinit(A);
    }



    template <std::size_t block_rows, std::size_t block_columns>
    inline BlockSparseMatrix::BlockSparseMatrix(
      const std::array<std::array<Mat, block_columns>, block_rows> &arrayA)
      : BlockSparseMatrix()
    {
      this->reinit(block_rows, block_columns);
      this->sub_objects.reinit(block_rows, block_columns);
      for (auto r = 0; r < block_rows; ++r)
        for (auto c = 0; c < block_columns; ++c)
          {
            if (arrayA[r][c])
              this->sub_objects[r][c] = new BlockType(arrayA[r][c]);
            else
              this->sub_objects[r][c] = nullptr;
          }
      this->collect_sizes();
    }



    inline BlockSparseMatrix &
    BlockSparseMatrix::operator=(const double d)
    {
      Assert(d == 0, ExcScalarAssignmentOnlyForZeroValue());

      for (size_type r = 0; r < this->n_block_rows(); ++r)
        for (size_type c = 0; c < this->n_block_cols(); ++c)
          this->block(r, c) = d;

      return *this;
    }



    inline void
    BlockSparseMatrix::vmult(BlockVector &dst, const BlockVector &src) const
    {
      BaseClass::vmult_block_block(dst, src);
    }



    inline void
    BlockSparseMatrix::vmult(BlockVector &dst, const Vector &src) const
    {
      BaseClass::vmult_block_nonblock(dst, src);
    }



    inline void
    BlockSparseMatrix::vmult(Vector &dst, const BlockVector &src) const
    {
      BaseClass::vmult_nonblock_block(dst, src);
    }



    inline void
    BlockSparseMatrix::vmult(Vector &dst, const Vector &src) const
    {
      BaseClass::vmult_nonblock_nonblock(dst, src);
    }


    inline void
    BlockSparseMatrix::Tvmult(BlockVector &dst, const BlockVector &src) const
    {
      BaseClass::Tvmult_block_block(dst, src);
    }



    inline void
    BlockSparseMatrix::Tvmult(BlockVector &dst, const Vector &src) const
    {
      BaseClass::Tvmult_block_nonblock(dst, src);
    }



    inline void
    BlockSparseMatrix::Tvmult(Vector &dst, const BlockVector &src) const
    {
      BaseClass::Tvmult_nonblock_block(dst, src);
    }



    inline void
    BlockSparseMatrix::Tvmult(Vector &dst, const Vector &src) const
    {
      BaseClass::Tvmult_nonblock_nonblock(dst, src);
    }

  } // namespace MPI

} // namespace PETScWrappers


DEAL_II_NAMESPACE_CLOSE


#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC

#endif // dealii_petsc_block_sparse_matrix_h

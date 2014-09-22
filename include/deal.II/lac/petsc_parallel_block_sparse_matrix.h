// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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

#ifndef __deal2__petsc_parallel_block_sparse_matrix_h
#define __deal2__petsc_parallel_block_sparse_matrix_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/base/table.h>
#  include <deal.II/lac/block_matrix_base.h>
#  include <deal.II/lac/block_sparsity_pattern.h>
#  include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#  include <deal.II/lac/petsc_parallel_block_vector.h>
#  include <deal.II/lac/exceptions.h>
#  include <cmath>

DEAL_II_NAMESPACE_OPEN



namespace PETScWrappers
{
  namespace MPI
  {

    /*! @addtogroup PETScWrappers
     *@{
     */

    /**
     * Blocked sparse matrix based on the PETScWrappers::SparseMatrix class. This
     * class implements the functions that are specific to the PETSc SparseMatrix
     * base objects for a blocked sparse matrix, and leaves the actual work
     * relaying most of the calls to the individual blocks to the functions
     * implemented in the base class. See there also for a description of when
     * this class is useful.
     *
     * In contrast to the deal.II-type SparseMatrix class, the PETSc matrices do
     * not have external objects for the sparsity patterns. Thus, one does not
     * determine the size of the individual blocks of a block matrix of this type
     * by attaching a block sparsity pattern, but by calling reinit() to set the
     * number of blocks and then by setting the size of each block separately. In
     * order to fix the data structures of the block matrix, it is then necessary
     * to let it know that we have changed the sizes of the underlying
     * matrices. For this, one has to call the collect_sizes() function, for much
     * the same reason as is documented with the BlockSparsityPattern class.
     *
     * @ingroup Matrix1
     * @see @ref GlossBlockLA "Block (linear algebra)"
     * @author Wolfgang Bangerth, 2004
     */
    class BlockSparseMatrix : public BlockMatrixBase<SparseMatrix>
    {
    public:
      /**
       * Typedef the base class for simpler
       * access to its own typedefs.
       */
      typedef BlockMatrixBase<SparseMatrix> BaseClass;

      /**
       * Typedef the type of the underlying
       * matrix.
       */
      typedef BaseClass::BlockType  BlockType;

      /**
       * Import the typedefs from the base
       * class.
       */
      typedef BaseClass::value_type      value_type;
      typedef BaseClass::pointer         pointer;
      typedef BaseClass::const_pointer   const_pointer;
      typedef BaseClass::reference       reference;
      typedef BaseClass::const_reference const_reference;
      typedef BaseClass::size_type       size_type;
      typedef BaseClass::iterator        iterator;
      typedef BaseClass::const_iterator  const_iterator;

      /**
       * Constructor; initializes the
       * matrix to be empty, without
       * any structure, i.e.  the
       * matrix is not usable at
       * all. This constructor is
       * therefore only useful for
       * matrices which are members of
       * a class. All other matrices
       * should be created at a point
       * in the data flow where all
       * necessary information is
       * available.
       *
       * You have to initialize the
       * matrix before usage with
       * reinit(BlockSparsityPattern). The
       * number of blocks per row and
       * column are then determined by
       * that function.
       */
      BlockSparseMatrix ();

      /**
       * Destructor.
       */
      ~BlockSparseMatrix ();

      /**
       * Pseudo copy operator only copying
       * empty objects. The sizes of the
       * block matrices need to be the
       * same.
       */
      BlockSparseMatrix &
      operator = (const BlockSparseMatrix &);

      /**
       * This operator assigns a scalar to
       * a matrix. Since this does usually
       * not make much sense (should we set
       * all matrix entries to this value?
       * Only the nonzero entries of the
       * sparsity pattern?), this operation
       * is only allowed if the actual
       * value to be assigned is zero. This
       * operator only exists to allow for
       * the obvious notation
       * <tt>matrix=0</tt>, which sets all
       * elements of the matrix to zero,
       * but keep the sparsity pattern
       * previously used.
       */
      BlockSparseMatrix &
      operator = (const double d);

      /**
       * Resize the matrix, by setting
       * the number of block rows and
       * columns. This deletes all
       * blocks and replaces them by
       * unitialized ones, i.e. ones
       * for which also the sizes are
       * not yet set. You have to do
       * that by calling the @p reinit
       * functions of the blocks
       * themselves. Do not forget to
       * call collect_sizes() after
       * that on this object.
       *
       * The reason that you have to
       * set sizes of the blocks
       * yourself is that the sizes may
       * be varying, the maximum number
       * of elements per row may be
       * varying, etc. It is simpler
       * not to reproduce the interface
       * of the SparsityPattern
       * class here but rather let the
       * user call whatever function
       * she desires.
       */
      void reinit (const size_type n_block_rows,
                   const size_type n_block_columns);


      /**
       * Efficiently reinit the block matrix for a parallel computation.
       * Only the BlockSparsityPattern of the Simple type can efficiently
       * store large sparsity patterns in parallel, so this is the only
       * supported argument.
       * The IndexSets describe the locally owned range of DoFs for each block.
       * Note that each IndexSet needs to be contiguous. For a symmetric structure
       * hand in the same vector for the first two arguments.
       */
      void reinit(const std::vector<IndexSet> &rows,
                  const std::vector<IndexSet> &cols,
                  const BlockCompressedSimpleSparsityPattern &bcsp,
                  const MPI_Comm &com);


      /**
       * Same as above but for a symmetric structure only.
       */
      void reinit(const std::vector<IndexSet> &sizes,
                  const BlockCompressedSimpleSparsityPattern &bcsp,
                  const MPI_Comm &com);



      /**
       * Matrix-vector multiplication:
       * let $dst = M*src$ with $M$
       * being this matrix.
       */
      void vmult (BlockVector       &dst,
                  const BlockVector &src) const;

      /**
       * Matrix-vector
       * multiplication. Just like the
       * previous function, but only
       * applicable if the matrix has
       * only one block column.
       */
      void vmult (BlockVector          &dst,
                  const Vector &src) const;

      /**
       * Matrix-vector
       * multiplication. Just like the
       * previous function, but only
       * applicable if the matrix has
       * only one block row.
       */
      void vmult (Vector    &dst,
                  const BlockVector &src) const;

      /**
       * Matrix-vector
       * multiplication. Just like the
       * previous function, but only
       * applicable if the matrix has
       * only one block.
       */
      void vmult (Vector       &dst,
                  const Vector &src) const;

      /**
       * Matrix-vector multiplication:
       * let $dst = M^T*src$ with $M$
       * being this matrix. This
       * function does the same as
       * vmult() but takes the
       * transposed matrix.
       */
      void Tvmult (BlockVector       &dst,
                   const BlockVector &src) const;

      /**
       * Matrix-vector
       * multiplication. Just like the
       * previous function, but only
       * applicable if the matrix has
       * only one block row.
       */
      void Tvmult (BlockVector  &dst,
                   const Vector &src) const;

      /**
       * Matrix-vector
       * multiplication. Just like the
       * previous function, but only
       * applicable if the matrix has
       * only one block column.
       */
      void Tvmult (Vector    &dst,
                   const BlockVector &src) const;

      /**
       * Matrix-vector
       * multiplication. Just like the
       * previous function, but only
       * applicable if the matrix has
       * only one block.
       */
      void Tvmult (Vector       &dst,
                   const Vector &src) const;

      /**
       * This function collects the
       * sizes of the sub-objects and
       * stores them in internal
       * arrays, in order to be able to
       * relay global indices into the
       * matrix to indices into the
       * subobjects. You *must* call
       * this function each time after
       * you have changed the size of
       * the sub-objects.
       */
      void collect_sizes ();

      /**
       * Return a reference to the MPI
       * communicator object in use with
       * this matrix.
       */
      const MPI_Comm &get_mpi_communicator () const;

      /**
       * Make the clear() function in the
       * base class visible, though it is
       * protected.
       */
      using BlockMatrixBase<SparseMatrix>::clear;
    };



    /*@}*/

// ------------- inline and template functions -----------------

    inline
    BlockSparseMatrix &
    BlockSparseMatrix::operator = (const double d)
    {
      Assert (d==0, ExcScalarAssignmentOnlyForZeroValue());

      for (size_type r=0; r<this->n_block_rows(); ++r)
        for (size_type c=0; c<this->n_block_cols(); ++c)
          this->block(r,c) = d;

      return *this;
    }



    inline
    void
    BlockSparseMatrix::vmult (BlockVector       &dst,
                              const BlockVector &src) const
    {
      BaseClass::vmult_block_block (dst, src);
    }



    inline
    void
    BlockSparseMatrix::vmult (BlockVector  &dst,
                              const Vector &src) const
    {
      BaseClass::vmult_block_nonblock (dst, src);
    }



    inline
    void
    BlockSparseMatrix::vmult (Vector            &dst,
                              const BlockVector &src) const
    {
      BaseClass::vmult_nonblock_block (dst, src);
    }



    inline
    void
    BlockSparseMatrix::vmult (Vector       &dst,
                              const Vector &src) const
    {
      BaseClass::vmult_nonblock_nonblock (dst, src);
    }


    inline
    void
    BlockSparseMatrix::Tvmult (BlockVector       &dst,
                               const BlockVector &src) const
    {
      BaseClass::Tvmult_block_block (dst, src);
    }



    inline
    void
    BlockSparseMatrix::Tvmult (BlockVector  &dst,
                               const Vector &src) const
    {
      BaseClass::Tvmult_block_nonblock (dst, src);
    }



    inline
    void
    BlockSparseMatrix::Tvmult (Vector            &dst,
                               const BlockVector &src) const
    {
      BaseClass::Tvmult_nonblock_block (dst, src);
    }



    inline
    void
    BlockSparseMatrix::Tvmult (Vector       &dst,
                               const Vector &src) const
    {
      BaseClass::Tvmult_nonblock_nonblock (dst, src);
    }

  }

}


DEAL_II_NAMESPACE_CLOSE


#endif    // DEAL_II_WITH_PETSC

#endif    // __deal2__petsc_parallel_block_sparse_matrix_h

//----------------------------  petsc_block_sparse_matrix.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  petsc_block_sparse_matrix.h  ---------------------------
#ifndef __deal2__petsc_block_sparse_matrix_h
#define __deal2__petsc_block_sparse_matrix_h


#include <base/config.h>
#include <base/table.h>
#include <lac/block_matrix_base.h>
#include <lac/petsc_sparse_matrix.h>
#include <cmath>


#ifdef DEAL_II_USE_PETSC


/*! @addtogroup PETSc
 *@{
 */

namespace PETScWrappers
{
  

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
 * @author Wolfgang Bangerth, 2004
 */
  class BlockSparseMatrix : public BlockMatrixBase<PETScWrappers::SparseMatrix>
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
                                        * empty objects. The sizes of the block
                                        * matrices need to be the same.
                                        */
      BlockSparseMatrix &
      operator = (const BlockSparseMatrix &);

                                       /**
                                        * This operator assigns a scalar to a
                                        * matrix. Since this does usually not
                                        * make much sense (should we set all
                                        * matrix entries to this value? Only
                                        * the nonzero entries of the sparsity
                                        * pattern?), this operation is only
                                        * allowed if the actual value to be
                                        * assigned is zero. This operator only
                                        * exists to allow for the obvious
                                        * notation <tt>matrix=0</tt>, which
                                        * sets all elements of the matrix to
                                        * zero, but keep the sparsity pattern
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
                                        * call @p collect_sizes after
                                        * that on this object.
                                        *
                                        * The reason that you have to
                                        * set sizes of the blocks
                                        * yourself is that the sizes may
                                        * be varying, the maximum number
                                        * of elements per row may be
                                        * varying, etc. It is simpler
                                        * not to reproduce the interface
                                        * of the @p SparsityPattern
                                        * class here but rather let the
                                        * user call whatever function
                                        * she desires.
                                        */
      void reinit (const unsigned int n_block_rows,
                   const unsigned int n_block_columns);
      

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
                                        * Make the clear() function in the
                                        * base class visible, though it is
                                        * protected.
                                        */
      BlockMatrixBase<SparseMatrix>::clear;
      
                                       /**
                                        * Exception
                                        */
      DeclException4 (ExcIncompatibleRowNumbers,
                      int, int, int, int,
                      << "The blocks [" << arg1 << ',' << arg2 << "] and ["
                      << arg3 << ',' << arg4 << "] have differing row numbers.");
                                       /**
                                        * Exception
                                        */
      DeclException4 (ExcIncompatibleColNumbers,
                      int, int, int, int,
                      << "The blocks [" << arg1 << ',' << arg2 << "] and ["
                      << arg3 << ',' << arg4 << "] have differing column numbers.");
  };



/*@}*/
}


#endif    // DEAL_II_USE_PETSC

#endif    // __deal2__petsc_block_sparse_matrix_h

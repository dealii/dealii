//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__trilinos_block_sparse_matrix_h
#define __deal2__trilinos_block_sparse_matrix_h


#include <base/config.h>

#ifdef DEAL_II_USE_TRILINOS

#  include <base/table.h>
#  include <lac/block_matrix_base.h>
#  include <lac/trilinos_sparse_matrix.h>
#  include <lac/trilinos_block_vector.h>
#  include <lac/full_matrix.h>
#  include <lac/exceptions.h>

#  include <cmath>

#  define TrilinosScalar double

DEAL_II_NAMESPACE_OPEN

                                   // forward declarations
class BlockSparsityPattern;
class BlockCompressedSparsityPattern;
class BlockCompressedSetSparsityPattern;
class BlockCompressedSimpleSparsityPattern;
template <typename number> class BlockSparseMatrix;


namespace TrilinosWrappers
{

/*! @addtogroup TrilinosWrappers
 *@{
 */

/**
 * Blocked sparse matrix based on the TrilinosWrappers::SparseMatrix class. This
 * class implements the functions that are specific to the Trilinos SparseMatrix
 * base objects for a blocked sparse matrix, and leaves the actual work
 * relaying most of the calls to the individual blocks to the functions
 * implemented in the base class. See there also for a description of when
 * this class is useful.
 *
 * In contrast to the deal.II-type SparseMatrix class, the Trilinos matrices do
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
 * @author Martin Kronbichler, Wolfgang Bangerth, 2008
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
                                        * of the @p SparsityPattern
                                        * class here but rather let the
                                        * user call whatever function
                                        * she desires.
                                        */
      void reinit (const unsigned int n_block_rows,
                   const unsigned int n_block_columns);

                                       /**
                                        * Resize the matrix, by using an
				        * array of Epetra maps to determine
				        * the %parallel distribution of the
				        * individual matrices. This function
				        * assumes that a quadratic block
				        * matrix is generated.
                                        */
      template <typename BlockSparsityType>
      void reinit (const std::vector<Epetra_Map> &input_maps,
		   const BlockSparsityType       &block_sparsity_pattern);

                                       /**
                                        * Resize the matrix, by using an
				        * array of index sets to determine
				        * the %parallel distribution of the
				        * individual matrices. This function
				        * assumes that a quadratic block
				        * matrix is generated.
                                        */
      template <typename BlockSparsityType>
      void reinit (const std::vector<IndexSet> &input_maps,
		   const BlockSparsityType     &block_sparsity_pattern,
		   const MPI_Comm              &communicator = MPI_COMM_WORLD);

                                       /**
                                        * Resize the matrix and initialize it
                                        * by the given sparsity pattern. Since
                                        * no distribution map is given, the
                                        * result is a block matrix for which
                                        * all elements are stored locally.
                                        */
      template <typename BlockSparsityType>
      void reinit (const BlockSparsityType &block_sparsity_pattern);

                                       /**
                                        * This function initializes the
				        * Trilinos matrix using the deal.II
				        * sparse matrix and the entries stored
				        * therein. It uses a threshold
				        * to copy only elements whose
				        * modulus is larger than the
				        * threshold (so zeros in the
				        * deal.II matrix can be filtered
				        * away).
                                        */
      void reinit (const std::vector<Epetra_Map>             &input_maps,
		   const ::dealii::BlockSparseMatrix<double> &deal_ii_sparse_matrix,
		   const double                               drop_tolerance=1e-13);

                                       /**
                                        * This function initializes
				        * the Trilinos matrix using
				        * the deal.II sparse matrix
				        * and the entries stored
				        * therein. It uses a threshold
				        * to copy only elements whose
				        * modulus is larger than the
				        * threshold (so zeros in the
				        * deal.II matrix can be
				        * filtered away). Since no
				        * Epetra_Map is given, all the
				        * elements will be locally
				        * stored.
                                        */
      void reinit (const ::dealii::BlockSparseMatrix<double> &deal_ii_sparse_matrix,
		   const double                               drop_tolerance=1e-13);

				       /**
					* This function calls the compress()
				        * command of all matrices after
				        * the assembly is
				        * completed. Note that all MPI
				        * processes need to call this
				        * command (whereas the individual
				        * assembly routines will most probably
				        * only be called on each processor
				        * individually) before any
				        * can complete it.
				        */
      void compress ();

				       /**
					* Returns the state of the
					* matrix, i.e., whether
					* compress() needs to be called
					* after an operation requiring
					* data exchange. Does only
					* return non-true values when
					* used in <tt>debug</tt> mode,
					* since it is quite expensive to
					* keep track of all operations
					* that lead to the need for
					* compress().
					*/
      bool is_compressed () const;

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
                                        * the sub-objects. Note that
                                        * this is a collective
                                        * operation, i.e., it needs to
                                        * be called on all MPI
                                        * processes. This command
                                        * internally calls the method
                                        * <tt>compress()</tt>, so you
                                        * don't need to call that
                                        * function in case you use
                                        * <tt>collect_sizes()</tt>.
                                        */
      void collect_sizes ();

                                       /**
                                        * Return the number of nonzero
                                        * elements of this
                                        * matrix.
                                        */
      unsigned int n_nonzero_elements () const;

                                       /**
                                        * Matrix-vector multiplication:
                                        * let $dst = M*src$ with $M$
                                        * being this matrix.
                                        */
      void vmult (MPI::BlockVector       &dst,
                  const MPI::BlockVector &src) const;


                                       /**
                                        * Matrix-vector multiplication:
                                        * let $dst = M*src$ with $M$
                                        * being this matrix, now applied
                                        * to localized block vectors
                                        * (works only when run on one
                                        * processor).
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
      void vmult (MPI::BlockVector  &dst,
                  const MPI::Vector &src) const;

                                       /**
                                        * Matrix-vector
                                        * multiplication. Just like the
                                        * previous function, but only
                                        * applicable if the matrix has
                                        * only one block column, now
                                        * applied to localized vectors
                                        * (works only when run on one
                                        * processor).
                                        */
      void vmult (BlockVector  &dst,
		  const Vector &src) const;

                                       /**
                                        * Matrix-vector
                                        * multiplication. Just like the
                                        * previous function, but only
                                        * applicable if the matrix has
                                        * only one block row.
                                        */
      void vmult (MPI::Vector            &dst,
		  const MPI::BlockVector &src) const;

                                       /**
                                        * Matrix-vector
                                        * multiplication. Just like the
                                        * previous function, but only
                                        * applicable if the matrix has
                                        * only one block row, now
                                        * applied to localized vectors
                                        * (works only when run on one
                                        * processor).
                                        */
      void vmult (Vector            &dst,
		  const BlockVector &src) const;

                                       /**
                                        * Matrix-vector
                                        * multiplication. Just like the
                                        * previous function, but only
                                        * applicable if the matrix has
                                        * only one block.
                                        */
      void vmult (VectorBase       &dst,
                  const VectorBase &src) const;

                                       /**
                                        * Matrix-vector multiplication:
                                        * let $dst = M^T*src$ with $M$
                                        * being this matrix. This
                                        * function does the same as
                                        * vmult() but takes the
                                        * transposed matrix.
                                        */
      void Tvmult (MPI::BlockVector       &dst,
                   const MPI::BlockVector &src) const;

                                       /**
                                        * Matrix-vector multiplication:
                                        * let $dst = M^T*src$ with $M$
                                        * being this matrix. This
                                        * function does the same as
                                        * vmult() but takes the
                                        * transposed matrix, now applied
                                        * to localized Trilinos vectors
                                        * (works only when run on one
                                        * processor).
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
      void Tvmult (MPI::BlockVector  &dst,
                   const MPI::Vector &src) const;

                                       /**
                                        * Matrix-vector
                                        * multiplication. Just like the
                                        * previous function, but only
                                        * applicable if the matrix has
                                        * only one block row, now
                                        * applied to localized Trilinos
                                        * vectors (works only when run
                                        * on one processor).
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
      void Tvmult (MPI::Vector    &dst,
                   const MPI::BlockVector &src) const;

                                       /**
                                        * Matrix-vector
                                        * multiplication. Just like the
                                        * previous function, but only
                                        * applicable if the matrix has
                                        * only one block column, now
                                        * applied to localized Trilinos
                                        * vectors (works only when run
                                        * on one processor).
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
      void Tvmult (VectorBase       &dst,
                   const VectorBase &src) const;

                                       /**
                                        * Compute the residual of an
                                        * equation <i>Mx=b</i>, where
                                        * the residual is defined to
                                        * be <i>r=b-Mx</i>. Write the
                                        * residual into @p dst. The
                                        * <i>l<sub>2</sub></i> norm of
                                        * the residual vector is
                                        * returned.
                                        *
                                        * Source <i>x</i> and
                                        * destination <i>dst</i> must
                                        * not be the same vector.
					*
					* Note that both vectors have
					* to be distributed vectors
					* generated using the same Map
					* as was used for the matrix
					* in case you work on a
					* distributed memory
					* architecture, using the
					* interface in the
					* TrilinosWrappers::MPI::BlockVector
					* class.
                                        */
      TrilinosScalar residual (MPI::BlockVector       &dst,
			       const MPI::BlockVector &x,
			       const MPI::BlockVector &b) const;

                                       /**
                                        * Compute the residual of an
                                        * equation <i>Mx=b</i>, where
                                        * the residual is defined to
                                        * be <i>r=b-Mx</i>. Write the
                                        * residual into @p dst. The
                                        * <i>l<sub>2</sub></i> norm of
                                        * the residual vector is
                                        * returned.
                                        *
                                        * Source <i>x</i> and
                                        * destination <i>dst</i> must
                                        * not be the same vector.
					*
					* Note that both vectors have
					* to be distributed vectors
					* generated using the same Map
					* as was used for the matrix
					* in case you work on a
					* distributed memory
					* architecture, using the
					* interface in the
					* TrilinosWrappers::BlockVector
					* class. Since the block
					* matrix is in general
					* distributed among processes,
					* this function only works
					* when running the program on
					* one processor.
                                        */
      TrilinosScalar residual (BlockVector       &dst,
			       const BlockVector &x,
			       const BlockVector &b) const;

                                       /**
                                        * Compute the residual of an
                                        * equation <i>Mx=b</i>, where
                                        * the residual is defined to
                                        * be <i>r=b-Mx</i>. Write the
                                        * residual into @p dst. The
                                        * <i>l<sub>2</sub></i> norm of
                                        * the residual vector is
                                        * returned. Just like the
                                        * previous function, but only
                                        * applicable if the matrix
                                        * only has one block row.
                                        */
      TrilinosScalar residual (MPI::BlockVector       &dst,
			       const MPI::Vector      &x,
			       const MPI::BlockVector &b) const;

                                       /**
                                        * Compute the residual of an
                                        * equation <i>Mx=b</i>, where
                                        * the residual is defined to
                                        * be <i>r=b-Mx</i>. Write the
                                        * residual into @p dst. The
                                        * <i>l<sub>2</sub></i> norm of
                                        * the residual vector is
                                        * returned. Just like the
                                        * previous function, but only
                                        * applicable if the matrix
                                        * only has one block row.
                                        */
      TrilinosScalar residual (BlockVector       &dst,
			       const Vector      &x,
			       const BlockVector &b) const;

                                       /**
                                        * Compute the residual of an
                                        * equation <i>Mx=b</i>, where
                                        * the residual is defined to
                                        * be <i>r=b-Mx</i>. Write the
                                        * residual into @p dst. The
                                        * <i>l<sub>2</sub></i> norm of
                                        * the residual vector is
                                        * returned. Just like the
                                        * previous function, but only
                                        * applicable if the matrix
                                        * only has one block column.
                                        */
      TrilinosScalar residual (MPI::Vector            &dst,
			       const MPI::BlockVector &x,
			       const MPI::Vector      &b) const;

                                       /**
                                        * Compute the residual of an
                                        * equation <i>Mx=b</i>, where
                                        * the residual is defined to
                                        * be <i>r=b-Mx</i>. Write the
                                        * residual into @p dst. The
                                        * <i>l<sub>2</sub></i> norm of
                                        * the residual vector is
                                        * returned. Just like the
                                        * previous function, but only
                                        * applicable if the matrix
                                        * only has one block column.
                                        */
      TrilinosScalar residual (Vector            &dst,
			       const BlockVector &x,
			       const Vector      &b) const;

                                       /**
                                        * Compute the residual of an
                                        * equation <i>Mx=b</i>, where
                                        * the residual is defined to
                                        * be <i>r=b-Mx</i>. Write the
                                        * residual into @p dst. The
                                        * <i>l<sub>2</sub></i> norm of
                                        * the residual vector is
                                        * returned. Just like the
                                        * previous function, but only
                                        * applicable if the matrix
                                        * only has one block.
                                        */
      TrilinosScalar residual (VectorBase       &dst,
			       const VectorBase &x,
			       const VectorBase &b) const;

                                       /**
                                        * Make the clear() function in the
                                        * base class visible, though it is
                                        * protected.
                                        */
      using BlockMatrixBase<SparseMatrix>::clear;

				       /** @addtogroup Exceptions
					* @{
					*/

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
				       ///@}
  };



/*@}*/

// ------------- inline and template functions -----------------

  inline
  bool
  BlockSparseMatrix::is_compressed () const
  {
    bool compressed = true;
    for (unsigned int row=0; row<n_block_rows(); ++row)
      for (unsigned int col=0; col<n_block_cols(); ++col)
	if (block(row, col).is_compressed() == false)
	  {
	    compressed = false;
	    break;
	  }

    return compressed;
  }



  inline
  void
  BlockSparseMatrix::vmult (MPI::BlockVector       &dst,
			    const MPI::BlockVector &src) const
  {
    BaseClass::vmult_block_block (dst, src);
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
  BlockSparseMatrix::vmult (MPI::BlockVector  &dst,
			    const MPI::Vector &src) const
  {
    BaseClass::vmult_block_nonblock (dst, src);
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
  BlockSparseMatrix::vmult (MPI::Vector            &dst,
                            const MPI::BlockVector &src) const
  {
    BaseClass::vmult_nonblock_block (dst, src);
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
  BlockSparseMatrix::vmult (VectorBase       &dst,
			    const VectorBase &src) const
  {
    BaseClass::vmult_nonblock_nonblock (dst, src);
  }



  inline
  void
  BlockSparseMatrix::Tvmult (MPI::BlockVector       &dst,
			     const MPI::BlockVector &src) const
  {
    BaseClass::Tvmult_block_block (dst, src);
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
  BlockSparseMatrix::Tvmult (MPI::BlockVector  &dst,
			     const MPI::Vector &src) const
  {
    BaseClass::Tvmult_block_nonblock (dst, src);
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
  BlockSparseMatrix::Tvmult (MPI::Vector            &dst,
			     const MPI::BlockVector &src) const
  {
    BaseClass::Tvmult_nonblock_block (dst, src);
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
  BlockSparseMatrix::Tvmult (VectorBase       &dst,
			     const VectorBase &src) const
  {
    BaseClass::Tvmult_nonblock_nonblock (dst, src);
  }

}

DEAL_II_NAMESPACE_CLOSE

#endif    // DEAL_II_USE_TRILINOS

#endif    // __deal2__trilinos_block_sparse_matrix_h

//----------------------------  petsc_parallel_sparse_matrix.h  ---------------------------
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
//----------------------------  petsc_parallel_sparse_matrix.h  ---------------------------
#ifndef __deal2__petsc_parallel_sparse_matrix_h
#define __deal2__petsc_parallel_sparse_matrix_h


#include <base/config.h>
#include <lac/exceptions.h>

#ifdef DEAL_II_USE_PETSC

#include <lac/petsc_matrix_base.h>
#include <lac/block_sparsity_pattern.h>

#include <vector>



namespace PETScWrappers
{
  namespace MPI
  {
    
/**
 * Implementation of a parallel sparse matrix class based on PETSC, with rows
 * of the matrix distributed across an MPI network. All the functionality is
 * actually in the base class, except for the calls to generate a parallel
 * sparse matrix. This is possible since PETSc only works on an abstract
 * matrix type and internally distributes to functions that do the actual work
 * depending on the actual matrix type (much like using virtual
 * functions). Only the functions creating a matrix of specific type differ,
 * and are implemented in this particular class.
 *
 * There are a number of comments on the communication model as well as access
 * to individual elements in the documentation to the parallel vector
 * class. These comments apply here as well.
 *
 *
 * @section Partitioning Partitioning of matrices
 *
 * PETSc partitions parallel matrices so that each MPI process "owns" a
 * certain number of rows (i.e. only this process stores the respective
 * entries in these rows). The number of rows each process owns has to be
 * passed to the constructors and reinit() functions via the argument @p
 * local_rows. The individual values passed as @p local_rows on all the MPI
 * processes of course have to add up to the global number of rows of the
 * matrix.
 *
 * In addition to this, PETSc also partitions the rectangular chunk of the
 * matrix it owns (i.e. the @p local_rows times n() elements in the matrix),
 * so that matrix vector multiplications can be performed efficiently. This
 * column-partitioning therefore has to match the partitioning of the vectors
 * with which the matrix is multiplied, just as the row-partitioning has to
 * match the partitioning of destination vectors. This partitioning is passed
 * to the constructors and reinit() functions through the @p local_columns
 * variable, which again has to add up to the global number of columns in the
 * matrix. The name @p local_columns may be named inappropriately since it
 * does not reflect that only these columns are stored locally, but it
 * reflects the fact that these are the columns for which the elements of
 * incoming vectors are stored locally.
 *
 * To make things even more complicated, PETSc needs a very good estimate of
 * the number of elements to be stored in each row to be efficient. Otherwise
 * it spends most of the time with allocating small chunks of memory, a
 * process that can slow down programs to a crawl if it happens to often. As
 * if a good estimate of the number of entries per row isn't even, it even
 * needs to split this as follows: for each row it owns, it needs an estimate
 * for the number of elements in this row that fall into the columns that are
 * set apart for this process (see above), and the number of elements that are
 * in the rest of the columns.
 *
 * Since in general this information is not readily available, most of the
 * initializing functions of this class assume that all of the number of
 * elements you give as an argument to @p n_nonzero_per_row or by @p
 * row_lengths fall into the columns "owned" by this process, and none into
 * the other ones. This is a fair guess for most of the rows, since in a good
 * domain partitioning, nodes only interact with nodes that are within the
 * same subdomain. It does not hold for nodes on the interfaces of subdomain,
 * however, and for the rows corresponding to these nodes, PETSc will have to
 * allocate additional memory, a costly process.
 *
 * The only way to avoid this is to tell PETSc where the actual entries of the
 * matrix will be. For this, there are constructors and reinit() functions of
 * this class that take a CompressedSparsityPattern object containing all this
 * information. While in the general case it is sufficient if the constructors
 * and reinit() functions know the number of local rows and columns, the
 * functions getting a sparsity pattern also need to know the number of local
 * rows (@p local_rows_per_process) and columns (@p local_columns_per_process)
 * for all other processes, in order to compute which parts of the matrix are
 * which. Thus, it is not sufficient to just count the number of degrees of
 * freedom that belong to a particular process, but you have to have the
 * numbers for all processes available at all processes.
 * 
 * @ingroup PETScWrappers
 * @author Wolfgang Bangerth, 2004
 */
    class SparseMatrix : public MatrixBase
    {
      public:
                                         /**
                                          * A structure that describes some of
                                          * the traits of this class in terms
                                          * of its run-time behavior. Some
                                          * other classes (such as the block
                                          * matrix classes) that take one or
                                          * other of the matrix classes as its
                                          * template parameters can tune their
                                          * behavior based on the variables in
                                          * this class.
                                          */
        struct Traits
        {
                                             /**
                                              * It is not safe to elide
                                              * additions of zeros to
                                              * individual elements of this
                                              * matrix. The reason is that
                                              * additions to the matrix may
                                              * trigger collective operations
                                              * synchronising buffers on
                                              * multiple processes. If an
                                              * addition is elided on one
                                              * process, this may lead to
                                              * other processes hanging in an
                                              * infinite waiting loop.
                                              */
            static const bool zero_addition_can_be_elided = false;
        };

                                         /**
                                          * Default constructor. Create an
                                          * empty matrix.
                                          */
        SparseMatrix ();
      
                                         /**
                                          * Create a sparse matrix of
                                          * dimensions @p m times @p n, with
                                          * an initial guess of @p
                                          * n_nonzero_per_row nonzero elements
                                          * per row. PETSc is able to cope
                                          * with the situation that more than
                                          * this number of elements are later
                                          * allocated for a row, but this
                                          * involves copying data, and is thus
                                          * expensive.
                                          *
                                          * For the meaning of the @p
                                          * local_row and @p local_columns
                                          * parameters, see the class
                                          * documentation.
                                          * 
                                          * The @p is_symmetric flag determines
                                          * whether we should tell PETSc that
                                          * the matrix is going to be symmetric
                                          * (as indicated by the call
                                          * <tt>MatSetOption(mat,
                                          * MAT_SYMMETRIC)</tt>. Note that the
                                          * PETSc documentation states that one
                                          * cannot form an ILU decomposition of
                                          * a matrix for which this flag has
                                          * been set to @p true, only an
                                          * ICC. The default value of this flag
                                          * is @p false.
                                          */
        SparseMatrix (const MPI_Comm     &communicator,
                      const unsigned int  m,
                      const unsigned int  n,
                      const unsigned int  local_rows,
                      const unsigned int  local_columns,
                      const unsigned int  n_nonzero_per_row,
                      const bool          is_symmetric = false);

                                         /**
                                          * Initialize a rectangular matrix
                                          * with @p m rows and @p n columns.
                                          * The maximal number of nonzero
                                          * entries for each row separately is
                                          * given by the @p row_lengths array.
                                          *
                                          * For the meaning of the @p
                                          * local_row and @p local_columns
                                          * parameters, see the class
                                          * documentation.
                                          * 
                                          * Just as for the other
                                          * constructors: PETSc is able to
                                          * cope with the situation that more
                                          * than this number of elements are
                                          * later allocated for a row, but
                                          * this involves copying data, and is
                                          * thus expensive.
                                          *
                                          * The @p is_symmetric flag
                                          * determines whether we should tell
                                          * PETSc that the matrix is going to
                                          * be symmetric (as indicated by the
                                          * call <tt>MatSetOption(mat,
                                          * MAT_SYMMETRIC)</tt>. Note that the
                                          * PETSc documentation states that
                                          * one cannot form an ILU
                                          * decomposition of a matrix for
                                          * which this flag has been set to @p
                                          * true, only an ICC. The default
                                          * value of this flag is @p false.
                                          */
        SparseMatrix (const MPI_Comm                  &communicator,
                      const unsigned int               m,
                      const unsigned int               n,
                      const unsigned int               local_rows,
                      const unsigned int               local_columns,
                      const std::vector<unsigned int> &row_lengths,
                      const bool                       is_symmetric = false);
        
                                         /**
                                          * Initialize using the given
                                          * sparsity pattern with
                                          * communication happening over the
                                          * provided @p communicator.
                                          *
                                          * For the meaning of the @p
                                          * local_rows_per_process and @p
                                          * local_columns_per_process
                                          * parameters, see the class
                                          * documentation.
                                          * 
                                          * Note that PETSc can be very slow
                                          * if you do not provide it with a
                                          * good estimate of the lengths of
                                          * rows. Using the present function
                                          * is a very efficient way to do
                                          * this, as it uses the exact number
                                          * of nonzero entries for each row of
                                          * the matrix by using the given
                                          * sparsity pattern argument. If the
                                          * @p preset_nonzero_locations flag
                                          * is @p true, this function in
                                          * addition not only sets the correct
                                          * row sizes up front, but also
                                          * pre-allocated the correct nonzero
                                          * entries in the matrix.
                                          *
                                          * PETsc allows to later add
                                          * additional nonzero entries to a
                                          * matrix, by simply writing to these
                                          * elements. However, this will then
                                          * lead to additional memory
                                          * allocations which are very
                                          * inefficient and will greatly slow
                                          * down your program. It is therefore
                                          * significantly more efficient to
                                          * get memory allocation right from
                                          * the start.
                                          *
                                          * Despite the fact that it would
                                          * seem to be an obvious win, setting
                                          * the @p preset_nonzero_locations
                                          * flag to @p true doesn't seem to
                                          * accelerate program. Rather on the
                                          * contrary, it seems to be able to
                                          * slow down entire programs
                                          * somewhat. This is suprising, since
                                          * we can use efficient function
                                          * calls into PETSc that allow to
                                          * create multiple entries at once;
                                          * nevertheless, given the fact that
                                          * it is inefficient, the respective
                                          * flag has a default value equal to
                                          * @p false.
                                          */
        SparseMatrix (const MPI_Comm                  &communicator,
                      const CompressedSparsityPattern &sparsity_pattern,
                      const std::vector<unsigned int> &local_rows_per_process,
                      const std::vector<unsigned int> &local_columns_per_process,
                      const unsigned int               this_process,
                      const bool                       preset_nonzero_locations = false);

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
        SparseMatrix & operator = (const double d);

                                         /**
                                          * Throw away the present matrix and
                                          * generate one that has the same
                                          * properties as if it were created by
                                          * the constructor of this class with
                                          * the same argument list as the
                                          * present function.
                                          */
        void reinit (const MPI_Comm     &communicator,
                     const unsigned int  m,
                     const unsigned int  n,
                     const unsigned int  local_rows,
                     const unsigned int  local_columns,
                     const unsigned int  n_nonzero_per_row,
                     const bool          is_symmetric = false);

                                         /**
                                          * Throw away the present matrix and
                                          * generate one that has the same
                                          * properties as if it were created by
                                          * the constructor of this class with
                                          * the same argument list as the
                                          * present function.
                                          */
        void reinit (const MPI_Comm                  &communicator,
                     const unsigned int               m,
                     const unsigned int               n,
                     const unsigned int               local_rows,
                     const unsigned int               local_columns,
                     const std::vector<unsigned int> &row_lengths,
                     const bool                       is_symmetric = false);

                                         /**
                                          * Initialize using the given
                                          * sparsity pattern with
                                          * communication happening over the
                                          * provided @p communicator.
                                          *
                                          * Note that PETSc can be very slow
                                          * if you do not provide it with a
                                          * good estimate of the lengths of
                                          * rows. Using the present function
                                          * is a very efficient way to do
                                          * this, as it uses the exact number
                                          * of nonzero entries for each row of
                                          * the matrix by using the given
                                          * sparsity pattern argument. If the
                                          * @p preset_nonzero_locations flag
                                          * is @p true, this function in
                                          * addition not only sets the correct
                                          * row sizes up front, but also
                                          * pre-allocated the correct nonzero
                                          * entries in the matrix.
                                          *
                                          * PETsc allows to later add
                                          * additional nonzero entries to a
                                          * matrix, by simply writing to these
                                          * elements. However, this will then
                                          * lead to additional memory
                                          * allocations which are very
                                          * inefficient and will greatly slow
                                          * down your program. It is therefore
                                          * significantly more efficient to
                                          * get memory allocation right from
                                          * the start.
                                          *
                                          * Despite the fact that it would
                                          * seem to be an obvious win, setting
                                          * the @p preset_nonzero_locations
                                          * flag to @p true doesn't seem to
                                          * accelerate program. Rather on the
                                          * contrary, it seems to be able to
                                          * slow down entire programs
                                          * somewhat. This is suprising, since
                                          * we can use efficient function
                                          * calls into PETSc that allow to
                                          * create multiple entries at once;
                                          * nevertheless, given the fact that
                                          * it is inefficient, the respective
                                          * flag has a default value equal to
                                          * @p false.
                                          */
        void reinit (const MPI_Comm                  &communicator,
                     const CompressedSparsityPattern &sparsity_pattern,
                     const std::vector<unsigned int> &local_rows_per_process,
                     const std::vector<unsigned int> &local_columns_per_process,
                     const unsigned int               this_process,
                     const bool                       preset_nonzero_locations = false);
        
                                         /**
                                          * Return a reference to the MPI
                                          * communicator object in use with
                                          * this matrix.
                                          */
        const MPI_Comm & get_mpi_communicator () const;
        
				     /** @addtogroup Exceptions
				      * @{ */
					 /**
					  * Exception
					  */
	DeclException2 (ExcLocalRowsTooLarge,
			int, int,
			<< "The number of local rows " << arg1
			<< " must be larger than the total number of rows " << arg2);
					 //@}        
      private:

                                         /**
                                          * Copy of the communicator object to
                                          * be used for this parallel vector.
                                          */
        MPI_Comm communicator;

                                         /**
                                          * Do the actual work for the
                                          * respective reinit() function and
                                          * the matching constructor,
                                          * i.e. create a matrix. Getting rid
                                          * of the previous matrix is left to
                                          * the caller.
                                          */
        void do_reinit (const unsigned int m,
                        const unsigned int n,
                        const unsigned int local_rows,
                        const unsigned int local_columns,
                        const unsigned int n_nonzero_per_row,
                        const bool         is_symmetric = false);

                                         /**
                                          * Same as previous function.
                                          */
        void do_reinit (const unsigned int               m,
                        const unsigned int               n,
                        const unsigned int               local_rows,
                        const unsigned int               local_columns,
                        const std::vector<unsigned int> &row_lengths,
                        const bool                       is_symmetric = false);

                                         /**
                                          * Same as previous functions.
                                          */
        void do_reinit (const CompressedSparsityPattern &sparsity_pattern,
                        const std::vector<unsigned int> &local_rows_per_process,
                        const std::vector<unsigned int> &local_columns_per_process,
                        const unsigned int               this_process,
                        const bool                       preset_nonzero_locations);
    };



// -------- template and inline functions ----------

    inline
    const MPI_Comm &    
    SparseMatrix::get_mpi_communicator () const
    {
      return communicator;
    }
  }
}

#endif // DEAL_II_USE_PETSC

/*----------------------------   petsc_parallel_sparse_matrix.h     ---------------------------*/

#endif
/*----------------------------   petsc_parallel_sparse_matrix.h     ---------------------------*/

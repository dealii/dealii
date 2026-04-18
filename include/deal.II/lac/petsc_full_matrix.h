// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2004 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_petsc_full_matrix_h
#define dealii_petsc_full_matrix_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/petsc_matrix_base.h>

DEAL_II_NAMESPACE_OPEN



namespace PETScWrappers
{
  /**
   * @addtogroup PETScWrappers
   * @{
   */

  /**
   * Implementation of a sequential dense matrix class based on PETSc. All the
   * functionality is actually in the base class, except for the calls to
   * generate a sequential dense matrix. This is possible since PETSc only
   * works on an abstract matrix type and internally distributes to functions
   * that do the actual work depending on the actual matrix type (much like
   * using virtual functions). Only the functions creating a matrix of
   * specific type differ, and are implemented in this particular class.
   *
   * @ingroup Matrix1
   */
  class FullMatrix : public MatrixBase
  {
  public:
    /**
     * Declare type for container size.
     */
    using size_type = types::global_dof_index;


    /**
     * Default constructor. Create an empty matrix.
     */
    FullMatrix();


    /**
     * Create a full matrix of dimensions @p m times @p n.
     */
    FullMatrix(const size_type m, const size_type n);


    /**
     * Throw away the present matrix and generate one that has the same
     * properties as if it were created by the constructor of this class with
     * the same argument list as the present function.
     */
    void
    reinit(const size_type m, const size_type n);


  private:
    /**
     * Do the actual work for the respective reinit() function and the
     * matching constructor, i.e. create a matrix. Getting rid of the previous
     * matrix is left to the caller.
     */
    void
    do_reinit(const size_type m, const size_type n);
  };

  namespace MPI
  {
    /* Implementation of a parallel full matrix class based on PETSc,
    with rows of the matrix distributed across an MPI network.
    *
    * @ingroup PETScWrappers
    * @ingroup Matrix1
    */
    class FullMatrix : public MatrixBase
    {
    public:
      /**
       * Declare type for container size.
       */
      using size_type = types::global_dof_index;


      /**
       * Default constructor. Create an empty matrix.
       */
      FullMatrix();

      /* Initialize with communication
       * happening over the provided @p communicator.
       * For the meaning of the @p local_rows_per_process and @p
       * local_columns_per_process parameters, see the class documentation. */
      FullMatrix(const MPI_Comm                communicator,
                 const size_type               m,
                 const size_type               n,
                 const std::vector<size_type> &local_rows_per_process,
                 const std::vector<size_type> &local_columns_per_process,
                 const unsigned int            this_process);

      /**
       * Initialize a FullMatrix from a PETSc Mat object. Note that we do not
       * copy the matrix. The Mat object is referenced by the newly created
       * instance of the class using PetscObjectReference. This is in line
       * with the PETSc approach to object ownership, which mimics
       * std::shared_ptr.
       */
      explicit FullMatrix(const Mat &);

      /**
       * Destructor to free the PETSc object.
       */
      ~FullMatrix() override;

      void
      reinit(const MPI_Comm                communicator,
             const size_type               m,
             const size_type               n,
             const std::vector<size_type> &local_rows_per_process,
             const std::vector<size_type> &local_columns_per_process,
             const unsigned int            this_process);

    private:
      /**
       * Do the actual work for the respective reinit() function and the
       * matching constructor, i.e. create a matrix. Getting rid of the
       * previous matrix is left to the caller.
       */
      void
      do_reinit(const MPI_Comm                communicator,
                const size_type               m,
                const size_type               n,
                const std::vector<size_type> &local_rows_per_process,
                const std::vector<size_type> &local_columns_per_process,
                const unsigned int            this_process);
    };
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

#endif

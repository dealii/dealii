// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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

  /** @} */
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

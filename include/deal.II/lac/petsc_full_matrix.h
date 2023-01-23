// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
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

#endif // DEAL_II_WITH_PETSC

#endif

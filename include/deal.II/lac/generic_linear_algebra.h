// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2017 by the deal.II authors
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

#ifndef dealii_generic_linear_algebra_h
#define dealii_generic_linear_algebra_h

#include <deal.II/base/config.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>


DEAL_II_NAMESPACE_OPEN


/**
 * A namespace in which the deal.II linear algebra classes are typedef'ed to
 * generic names. There are similar namespaces LinearAlgebraPETSc and
 * LinearAlgebraTrilinos for typedefs to classes that interface with the PETSc
 * and Trilinos libraries.
 */
namespace LinearAlgebraDealII
{
  typedef Vector<double>      Vector;
  typedef BlockVector<double> BlockVector;

  typedef SparseMatrix<double> SparseMatrix;

  typedef PreconditionSSOR<SparseMatrix> PreconditionSSOR;
} // namespace LinearAlgebraDealII


DEAL_II_NAMESPACE_CLOSE


#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/block_sparsity_pattern.h>
#  include <deal.II/lac/petsc_parallel_block_sparse_matrix.h>
#  include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#  include <deal.II/lac/petsc_precondition.h>
#  include <deal.II/lac/petsc_solver.h>

DEAL_II_NAMESPACE_OPEN

/**
 * A namespace in which the wrappers to the PETSc linear algebra classes are
 * typedef'ed to generic names. There are similar namespaces
 * LinearAlgebraDealII and LinearAlgebraTrilinos for typedefs to deal.II's own
 * classes and classes that interface with Trilinos.
 */
namespace LinearAlgebraPETSc
{
  using namespace dealii;

  typedef PETScWrappers::SolverCG    SolverCG;
  typedef PETScWrappers::SolverGMRES SolverGMRES;

  /**
   * A namespace with typedefs to generic names for parallel PETSc linear
   * algebra objects.
   */
  namespace MPI
  {
    /**
     * Typedef for the vector type used.
     */
    typedef PETScWrappers::MPI::Vector Vector;

    /**
     * Typedef for the type used to describe vectors that consist of multiple
     * blocks.
     */
    typedef PETScWrappers::MPI::BlockVector BlockVector;

    /**
     * Typedef for the sparse matrix type used.
     */
    typedef PETScWrappers::MPI::SparseMatrix SparseMatrix;

    /**
     * Typedef for the type used to describe sparse matrices that consist of
     * multiple blocks.
     */
    typedef PETScWrappers::MPI::BlockSparseMatrix BlockSparseMatrix;

    typedef dealii::BlockDynamicSparsityPattern BlockCompressedSparsityPattern;

    /**
     * Typedef for the AMG preconditioner type.
     */
    typedef PETScWrappers::PreconditionBoomerAMG PreconditionAMG;

    /**
     * Typedef for the Incomplete Cholesky preconditioner.
     */
    typedef PETScWrappers::PreconditionICC PreconditionIC;

    /**
     * Typedef for the Incomplete LU decomposition preconditioner.
     */
    typedef PETScWrappers::PreconditionILU PreconditionILU;

    /**
     * Typedef for the Incomplete Jacobi decomposition preconditioner.
     */
    typedef PETScWrappers::PreconditionJacobi PreconditionJacobi;

    /**
     * Typedef for the SSOR preconditioner.
     */
    typedef PETScWrappers::PreconditionSSOR PreconditionSSOR;

  } // namespace MPI

} // namespace LinearAlgebraPETSc
DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_WITH_PETSC

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/lac/block_sparsity_pattern.h>
#  include <deal.II/lac/trilinos_block_sparse_matrix.h>
#  include <deal.II/lac/trilinos_precondition.h>
#  include <deal.II/lac/trilinos_solver.h>
#  include <deal.II/lac/trilinos_sparse_matrix.h>

DEAL_II_NAMESPACE_OPEN

/**
 * A namespace in which the wrappers to the Trilinos linear algebra classes
 * are typedef'ed to generic names. There are similar namespaces
 * LinearAlgebraDealII and LinearAlgebraPETSc for typedefs to deal.II's own
 * classes and classes that interface with PETSc.
 */
namespace LinearAlgebraTrilinos
{
  using namespace dealii;

  typedef TrilinosWrappers::SolverCG    SolverCG;
  typedef TrilinosWrappers::SolverGMRES SolverGMRES;

  /**
   * A namespace with typedefs to generic names for parallel Trilinos linear
   * algebra objects.
   */
  namespace MPI
  {
    /**
     * Typedef for the vector type used.
     */
    typedef TrilinosWrappers::MPI::Vector Vector;

    /**
     * Typedef for the type used to describe vectors that consist of multiple
     * blocks.
     */
    typedef TrilinosWrappers::MPI::BlockVector BlockVector;

    /**
     * Typedef for the sparse matrix type used.
     */
    typedef TrilinosWrappers::SparseMatrix SparseMatrix;

    /**
     * Typedef for the type used to describe sparse matrices that consist of
     * multiple blocks.
     */
    typedef TrilinosWrappers::BlockSparseMatrix BlockSparseMatrix;

    typedef TrilinosWrappers::BlockSparsityPattern
      BlockCompressedSparsityPattern;

    /**
     * Typedef for the AMG preconditioner type.
     */
    typedef TrilinosWrappers::PreconditionAMG PreconditionAMG;

    /**
     * Typedef for the Incomplete Cholesky preconditioner.
     */
    typedef TrilinosWrappers::PreconditionIC PreconditionIC;

    /**
     * Typedef for the Incomplete LU decomposition preconditioner.
     */
    typedef TrilinosWrappers::PreconditionILU PreconditionILU;

    /**
     * Typedef for the Incomplete Jacobi decomposition preconditioner.
     */
    typedef TrilinosWrappers::PreconditionJacobi PreconditionJacobi;

    /**
     * Typedef for the SSOR preconditioner
     */
    typedef TrilinosWrappers::PreconditionSSOR PreconditionSSOR;


  } // namespace MPI

} // namespace LinearAlgebraTrilinos

DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_WITH_TRILINOS



#endif

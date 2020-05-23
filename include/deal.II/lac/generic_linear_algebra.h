// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2020 by the deal.II authors
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

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>


DEAL_II_NAMESPACE_OPEN


/**
 * A namespace in which the deal.II linear algebra classes are aliased to
 * generic names. There are similar namespaces LinearAlgebraPETSc and
 * LinearAlgebraTrilinos for alias to classes that interface with the PETSc
 * and Trilinos libraries.
 */
namespace LinearAlgebraDealII
{
  /**
   * Typedef for the vector type used
   */
  using Vector = Vector<double>;

  /**
   * Typedef for the block-vector type used
   */
  using BlockVector = BlockVector<double>;

  /**
   * Typedef for sparse matrix type used
   */
  using SparseMatrix = SparseMatrix<double>;

  /**
   * Typedef describing sparse matrices that consist of multiple blocks.
   */
  using BlockSparseMatrix = BlockSparseMatrix<double>;

  /**
   * Typedef for the SSOR preconditioner used
   */
  using PreconditionSSOR = PreconditionSSOR<SparseMatrix>;
} // namespace LinearAlgebraDealII


DEAL_II_NAMESPACE_CLOSE


#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/petsc_block_sparse_matrix.h>
#  include <deal.II/lac/petsc_precondition.h>
#  include <deal.II/lac/petsc_solver.h>
#  include <deal.II/lac/petsc_sparse_matrix.h>

DEAL_II_NAMESPACE_OPEN

/**
 * A namespace in which the wrappers to the PETSc linear algebra classes are
 * aliased to generic names. There are similar namespaces
 * LinearAlgebraDealII and LinearAlgebraTrilinos for alias to deal.II's own
 * classes and classes that interface with Trilinos.
 */
namespace LinearAlgebraPETSc
{
  /**
   * Typedef for the CG solver type used.
   */
  using SolverCG = PETScWrappers::SolverCG;

  /**
   * Typedef for the GMRES solver type used.
   */
  using SolverGMRES = PETScWrappers::SolverGMRES;

  /**
   * A namespace with alias to generic names for parallel PETSc linear
   * algebra objects.
   */
  namespace MPI
  {
    /**
     * Typedef for the vector type used.
     */
    using Vector = PETScWrappers::MPI::Vector;

    /**
     * Typedef for the type used to describe vectors that consist of multiple
     * blocks.
     */
    using BlockVector = PETScWrappers::MPI::BlockVector;

    /**
     * Typedef for the sparse matrix type used.
     */
    using SparseMatrix = PETScWrappers::MPI::SparseMatrix;

    /**
     * Typedef for the type used to describe sparse matrices that consist of
     * multiple blocks.
     */
    using BlockSparseMatrix = PETScWrappers::MPI::BlockSparseMatrix;

    /**
     * Typedef for the compressed block sparsity pattern used.
     */
    using BlockCompressedSparsityPattern = dealii::BlockDynamicSparsityPattern;

    /**
     * Typedef for the AMG preconditioner type.
     */
    using PreconditionAMG = PETScWrappers::PreconditionBoomerAMG;

    /**
     * Typedef for the Incomplete Cholesky preconditioner.
     */
    using PreconditionIC = PETScWrappers::PreconditionICC;

    /**
     * Typedef for the Incomplete LU decomposition preconditioner.
     */
    using PreconditionILU = PETScWrappers::PreconditionILU;

    /**
     * Typedef for the Incomplete Jacobi decomposition preconditioner.
     */
    using PreconditionJacobi = PETScWrappers::PreconditionJacobi;

    /**
     * Typedef for the SSOR preconditioner.
     */
    using PreconditionSSOR = PETScWrappers::PreconditionSSOR;

  } // namespace MPI

} // namespace LinearAlgebraPETSc
DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_WITH_PETSC

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/lac/trilinos_block_sparse_matrix.h>
#  include <deal.II/lac/trilinos_precondition.h>
#  include <deal.II/lac/trilinos_solver.h>
#  include <deal.II/lac/trilinos_sparse_matrix.h>

DEAL_II_NAMESPACE_OPEN

/**
 * A namespace in which the wrappers to the Trilinos linear algebra classes
 * are aliased to generic names. There are similar namespaces
 * LinearAlgebraDealII and LinearAlgebraPETSc for alias to deal.II's own
 * classes and classes that interface with PETSc.
 */
namespace LinearAlgebraTrilinos
{
  /**
   * Typedef for the CG solver type used.
   */
  using SolverCG = TrilinosWrappers::SolverCG;

  /**
   * Typdef for the GMRES solver type used.
   */
  using SolverGMRES = TrilinosWrappers::SolverGMRES;

  /**
   * A namespace with alias to generic names for parallel Trilinos linear
   * algebra objects.
   */
  namespace MPI
  {
    /**
     * Typedef for the vector type used.
     */
    using Vector = TrilinosWrappers::MPI::Vector;

    /**
     * Typedef for the type used to describe vectors that consist of multiple
     * blocks.
     */
    using BlockVector = TrilinosWrappers::MPI::BlockVector;

    /**
     * Typedef for the sparse matrix type used.
     */
    using SparseMatrix = TrilinosWrappers::SparseMatrix;

    /**
     * Typedef for the type used to describe sparse matrices that consist of
     * multiple blocks.
     */
    using BlockSparseMatrix = TrilinosWrappers::BlockSparseMatrix;

    /**
     * Typedef for the type used for compressed block sparsity pattern.
     */
    using BlockCompressedSparsityPattern =
      TrilinosWrappers::BlockSparsityPattern;

    /**
     * Typedef for the AMG preconditioner type.
     */
    using PreconditionAMG = TrilinosWrappers::PreconditionAMG;

    /**
     * Typedef for the Incomplete Cholesky preconditioner.
     */
    using PreconditionIC = TrilinosWrappers::PreconditionIC;

    /**
     * Typedef for the Incomplete LU decomposition preconditioner.
     */
    using PreconditionILU = TrilinosWrappers::PreconditionILU;

    /**
     * Typedef for the Incomplete Jacobi decomposition preconditioner.
     */
    using PreconditionJacobi = TrilinosWrappers::PreconditionJacobi;

    /**
     * Typedef for the SSOR preconditioner
     */
    using PreconditionSSOR = TrilinosWrappers::PreconditionSSOR;


  } // namespace MPI

} // namespace LinearAlgebraTrilinos

DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_WITH_TRILINOS



#endif

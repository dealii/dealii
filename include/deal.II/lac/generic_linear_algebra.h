// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2016 by the deal.II authors
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

#ifndef dealii__generic_linear_algebra_h
#define dealii__generic_linear_algebra_h

#include <deal.II/base/config.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_selector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/precondition_selector.h>
#include <deal.II/lac/vector.h>

#include <complex>


DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup GenericLinearAlgebra
 * @{
 */

/**
 * A namespace in which the real valued deal.II linear algebra classes are
 * typedef'ed to generic names. There are similar namespaces LinearAlgebraPETSc
 * and LinearAlgebraTrilinos for typedefs to classes that interface with the
 * PETSc and Trilinos libraries.
 *
 * @ingroup GenericLinearAlgebra
 * @author Timo Heister 2008, Jean-Paul Pelteret 2016
 */
namespace LinearAlgebraDealII
{
  /**
   * Typedef for the vector type used.
   */
  typedef Vector<double> Vector;

  /**
  * Typedef for the type used to describe vectors that consist of multiple
  * blocks.
  */
  typedef BlockVector<double> BlockVector;

  /**
   * Typedef for the sparse matrix type used.
   */
  typedef SparseMatrix<double> SparseMatrix;

  /**
   * Typedef for the type used to describe sparse matrices that consist of
   * multiple blocks.
   */
  typedef BlockSparseMatrix<double> BlockSparseMatrix;

  /**
   * Typedef for the Jacobi preconditioner type.
   */
  typedef PreconditionJacobi<SparseMatrix> PreconditionJacobi;

  /**
  * Typedef for the generic preconditioner selector.
  */
  typedef PreconditionSelector<SparseMatrix> PreconditionSelector;

  /**
  * Typedef for the SOR preconditioner type.
  */
  typedef PreconditionSOR<SparseMatrix> PreconditionSOR;
  /**
  * Typedef for the SSOR preconditioner type.
  */
  typedef PreconditionSSOR<SparseMatrix> PreconditionSSOR;

  /**
   * Typedef for the base class for linear solvers.
   *
   * This typedef is useful when wanting to create a generic solver interface
   * wherein one can switch between linear solvers. For example:
   * @code
   * std_cxx11::unique_ptr<LA::SolverBase> solver;
   * if (parameters.solver_type == "CG")
   * {
   *   solver.reset(new LA::SolverCG(solver_control));
   * }
   * @endcode
   */
  typedef Solver<Vector> SolverBase;

  /**
   * Typedef for the BiCGStab linear solver.
   */
  typedef SolverBicgstab<Vector> SolverBicgstab;

  /**
   * Typedef for the CG linear solver.
   */
  typedef SolverCG<Vector> SolverCG;

  /**
   * Typedef for the GMRES linear solver.
   */
  typedef SolverGMRES<Vector> SolverGMRES;

  /**
   * Typedef for the generic linear solver selector.
   */
  typedef SolverSelector<Vector> SolverSelector;

  // ----------------------------------------------------
  // Below are the common elements between this namespace
  // and any namespaces nested within this one
  // ----------------------------------------------------

  /**
   * Typedef for dynamic sparsity patterns.
   */
  typedef DynamicSparsityPattern DynamicSparsityPattern;

  /**
   * Typedef for block dynamic sparsity patterns.
   */
  typedef BlockDynamicSparsityPattern BlockDynamicSparsityPattern;

#ifdef DEAL_II_WITH_UMFPACK

  /**
   * Typedef for the UMFPACK direct linear solver.
   */
  typedef SparseDirectUMFPACK SparseDirectUMFPACK;

#endif

  /**
  * Typedef for the identity preconditioner.
  */
  typedef PreconditionIdentity PreconditionIdentity;

  /**
   * A namespace in which the complex valued deal.II linear algebra classes are
   * typedef'ed to generic names.
   *
   * @ingroup GenericLinearAlgebra
   * @author Jean-Paul Pelteret 2016
   */
  namespace Complex
  {
    /**
     * Typedef for the vector type used.
     */
    typedef dealii::Vector< std::complex<double> > Vector;

    /**
    * Typedef for the type used to describe vectors that consist of multiple
    * blocks.
    */
    typedef dealii::BlockVector< std::complex<double> > BlockVector;

    /**
     * Typedef for the sparse matrix type used.
     */
    typedef dealii::SparseMatrix< std::complex<double> > SparseMatrix;

    /**
     * Typedef for the type used to describe sparse matrices that consist of
     * multiple blocks.
     */
    typedef dealii::BlockSparseMatrix< std::complex<double> > BlockSparseMatrix;

    /**
     * Typedef for the Jacobi preconditioner type.
     */
    typedef dealii::PreconditionJacobi<SparseMatrix> PreconditionJacobi;

    /**
    * Typedef for the generic preconditioner selector.
    */
    typedef dealii::PreconditionSelector<SparseMatrix> PreconditionSelector;

    /**
    * Typedef for the SOR preconditioner type.
    */
    typedef dealii::PreconditionSOR<SparseMatrix> PreconditionSOR;
    /**
    * Typedef for the SSOR preconditioner type.
    */
    typedef dealii::PreconditionSSOR<SparseMatrix> PreconditionSSOR;

    /**
     * Typedef for the base class for linear solvers.
     *
     * This typedef is useful when wanting to create a generic solver interface
     * wherein one can switch between linear solvers. For example:
     * @code
     * std_cxx11::unique_ptr<LA::SolverBase> solver;
     * if (parameters.solver_type == "CG")
     * {
     *   solver.reset(new LA::SolverCG(solver_control));
     * }
     * @endcode
     */
    typedef dealii::Solver<Vector> SolverBase;

    /**
     * Typedef for the BiCGStab linear solver.
     */
    typedef dealii::SolverBicgstab<Vector> SolverBicgstab;

    /**
     * Typedef for the CG linear solver.
     */
    typedef dealii::SolverCG<Vector> SolverCG;

    /**
     * Typedef for the GMRES linear solver.
     */
    typedef dealii::SolverGMRES<Vector> SolverGMRES;

    /**
     * Typedef for the generic linear solver selector.
     */
    typedef dealii::SolverSelector<Vector> SolverSelector;

    // ----------------------------------------------------
    // Below are the common elements between this namespace
    // and the namespaces that it is nested within
    // ----------------------------------------------------

    /**
     * Typedef for dynamic sparsity patterns.
     */
    typedef DynamicSparsityPattern DynamicSparsityPattern;

    /**
     * Typedef for block dynamic sparsity patterns.
     */
    typedef BlockDynamicSparsityPattern BlockDynamicSparsityPattern;

#ifdef DEAL_II_WITH_UMFPACK

    /**
     * Typedef for the UMFPACK direct linear solver.
     */
    typedef SparseDirectUMFPACK SparseDirectUMFPACK;

#endif

    /**
    * Typedef for the identity preconditioner.
    */
    typedef PreconditionIdentity PreconditionIdentity;
  }
}


DEAL_II_NAMESPACE_CLOSE


#ifdef DEAL_II_WITH_PETSC

#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/petsc_block_sparse_matrix.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_block_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

DEAL_II_NAMESPACE_OPEN

/**
 * A namespace in which the wrappers to the PETSc linear algebra classes are
 * typedef'ed to generic names. There are similar namespaces
 * LinearAlgebraDealII and LinearAlgebraTrilinos for typedefs to deal.II's own
 * classes and classes that interface with Trilinos.
 *
 * @ingroup GenericLinearAlgebra
 * @author Timo Heister 2008, Jean-Paul Pelteret 2016
 */
namespace LinearAlgebraPETSc
{
  /**
   * Typedef for the vector type used.
   */
  typedef dealii::PETScWrappers::Vector Vector;

  /**
   * Typedef for the type used to describe vectors that consist of multiple
   * blocks.
   */
  typedef dealii::PETScWrappers::BlockVector BlockVector;

  /**
   * Typedef for the sparse matrix type used.
   */
  typedef dealii::PETScWrappers::SparseMatrix SparseMatrix;

  /**
   * Typedef for the type used to describe sparse matrices that consist of
   * multiple blocks.
   */
  typedef dealii::PETScWrappers::BlockSparseMatrix BlockSparseMatrix;

  // ----------------------------------------------------
  // Below are the common elements between this namespace
  // and any namespaces nested within this one
  // ----------------------------------------------------

  /**
   * Typedef for dynamic sparsity patterns.
   */
  typedef dealii::DynamicSparsityPattern DynamicSparsityPattern;

  /**
   * Typedef for block dynamic sparsity patterns.
   */
  typedef dealii::BlockDynamicSparsityPattern BlockDynamicSparsityPattern;

  /**
   * Typedef for the base class for linear solvers.
   *
   * This typedef is useful when wanting to create a generic solver interface
   * wherein one can switch between linear solvers. For example:
   * @code
   * std_cxx11::unique_ptr<LA::SolverBase> solver;
   * if (parameters.solver_type == "CG")
   * {
   *   solver.reset(new LA::SolverCG(solver_control));
   * }
   * @endcode
   */
  typedef dealii::PETScWrappers::SolverBase SolverBase;

  /**
   * Typedef for the BiCGStab linear solver.
   */
  typedef dealii::PETScWrappers::SolverBicgstab SolverBicgstab;

  /**
   * Typedef for the CG linear solver.
   */
  typedef dealii::PETScWrappers::SolverCG SolverCG;

  /**
   * Typedef for the CG Squared linear solver.
   */
  typedef dealii::PETScWrappers::SolverCGS SolverCGS;

#ifdef DEAL_II_PETSC_WITH_MUMPS

  /**
   * Typedef for the MUMPS direct linear solver.
   */
  typedef dealii::PETScWrappers::SparseDirectMUMPS SparseDirectMUMPS;

#endif

  /**
   * Typedef for the GMRES linear solver.
   */
  typedef dealii::PETScWrappers::SolverGMRES SolverGMRES;

  /**
   * Typedef for the TFQMR linear solver.
   */
  typedef dealii::PETScWrappers::SolverTFQMR SolverTFQMR;

  /**
   * Typedef for the base class for preconditioners.
   *
   * This typedef is useful when wanting to create a generic solver interface
   * wherein one can switch between preconditioners. For example:
   * @code
   * std_cxx11::unique_ptr<LA::PreconditionBase> preconditioner;
   * if (parameters.preconditioner_type == "Jacobi")
   * {
   *   LA::MPI::PreconditionJacobi* ptr_prec
   *     = new LA::MPI::PreconditionJacobi ();
   *
   *   LA::MPI::PreconditionJacobi::AdditionalData
   *     additional_data (parameters.preconditioner_relaxation);
   *
   *   ptr_prec->initialize(tangent_matrix.block(u_block,u_block),
   *                        additional_data);
   *   preconditioner.reset(ptr_prec);
   * }
   * @endcode
   */
  typedef dealii::PETScWrappers::PreconditionerBase PreconditionBase;

  /**
   * Typedef for the AMG preconditioner type.
   */
  typedef dealii::PETScWrappers::PreconditionBoomerAMG PreconditionAMG;

  /**
   * Typedef for the Incomplete Cholesky preconditioner.
   */
  typedef dealii::PETScWrappers::PreconditionICC PreconditionIC;

  /**
   * Typedef for the identity preconditioner.
   */
  typedef dealii::PETScWrappers::PreconditionNone PreconditionIdentity;

  /**
   * Typedef for the Incomplete LU decomposition preconditioner.
   */
  typedef dealii::PETScWrappers::PreconditionILU PreconditionILU;

  /**
   * Typedef for the Incomplete Jacobi decomposition preconditioner.
   */
  typedef dealii::PETScWrappers::PreconditionJacobi PreconditionJacobi;

  /**
   * Typedef for the SOR preconditioner.
   */
  typedef dealii::PETScWrappers::PreconditionSOR PreconditionSOR;

  /**
   * Typedef for the SSOR preconditioner.
   */
  typedef dealii::PETScWrappers::PreconditionSSOR PreconditionSSOR;

  /**
   * A namespace with typedefs to generic names for parallel PETSc linear
   * algebra objects.
   *
   * @ingroup GenericLinearAlgebra
   * @author Timo Heister 2008, Jean-Paul Pelteret 2016
   */
  namespace MPI
  {
    /**
     * Typedef for the vector type used.
     */
    typedef dealii::PETScWrappers::MPI::Vector Vector;

    /**
     * Typedef for the type used to describe vectors that consist of multiple
     * blocks.
     */
    typedef dealii::PETScWrappers::MPI::BlockVector BlockVector;

    /**
     * Typedef for the sparse matrix type used.
     */
    typedef dealii::PETScWrappers::MPI::SparseMatrix SparseMatrix;

    /**
     * Typedef for the type used to describe sparse matrices that consist of
     * multiple blocks.
     */
    typedef dealii::PETScWrappers::MPI::BlockSparseMatrix BlockSparseMatrix;

    // ----------------------------------------------------
    // Below are the common elements between this namespace
    // and the namespaces that it is nested within
    // ----------------------------------------------------

    /**
     * Typedef for dynamic sparsity patterns.
     */
    typedef dealii::DynamicSparsityPattern DynamicSparsityPattern;

    /**
     * Typedef for block dynamic sparsity patterns.
     */
    typedef dealii::BlockDynamicSparsityPattern BlockDynamicSparsityPattern;

    /**
     * Typedef for the base class for linear solvers.
     *
     * This typedef is useful when wanting to create a generic solver interface
     * wherein one can switch between linear solvers. For example:
     * @code
     * std_cxx11::unique_ptr<LA::SolverBase> solver;
     * if (parameters.solver_type == "CG")
     * {
     *   solver.reset(new LA::SolverCG(solver_control));
     * }
     * @endcode
     */
    typedef dealii::PETScWrappers::SolverBase SolverBase;

    /**
     * Typedef for the BiCGStab linear solver.
     */
    typedef dealii::PETScWrappers::SolverBicgstab SolverBicgstab;

    /**
     * Typedef for the CG linear solver.
     */
    typedef dealii::PETScWrappers::SolverCG SolverCG;

    /**
     * Typedef for the CG Squared linear solver.
     */
    typedef dealii::PETScWrappers::SolverCGS SolverCGS;

#ifdef DEAL_II_PETSC_WITH_MUMPS

    /**
     * Typedef for the MUMPS direct linear solver.
     */
    typedef dealii::PETScWrappers::SparseDirectMUMPS SparseDirectMUMPS;

#endif

    /**
     * Typedef for the GMRES linear solver.
     */
    typedef dealii::PETScWrappers::SolverGMRES SolverGMRES;

    /**
     * Typedef for the TFQMR linear solver.
     */
    typedef dealii::PETScWrappers::SolverTFQMR SolverTFQMR;

    /**
     * Typedef for the base class for preconditioners.
     *
     * This typedef is useful when wanting to create a generic solver interface
     * wherein one can switch between preconditioners. For example:
     * @code
     * std_cxx11::unique_ptr<LA::PreconditionBase> preconditioner;
     * if (parameters.preconditioner_type == "Jacobi")
     * {
     *   LA::MPI::PreconditionJacobi* ptr_prec
     *     = new LA::MPI::PreconditionJacobi ();
     *
     *   LA::MPI::PreconditionJacobi::AdditionalData
     *     additional_data (parameters.preconditioner_relaxation);
     *
     *   ptr_prec->initialize(tangent_matrix.block(u_block,u_block),
     *                        additional_data);
     *   preconditioner.reset(ptr_prec);
     * }
     * @endcode
     */
    typedef dealii::PETScWrappers::PreconditionerBase PreconditionBase;

    /**
     * Typedef for the AMG preconditioner type.
     */
    typedef dealii::PETScWrappers::PreconditionBoomerAMG PreconditionAMG;

    /**
     * Typedef for the Incomplete Cholesky preconditioner.
     */
    typedef dealii::PETScWrappers::PreconditionICC PreconditionIC;

    /**
     * Typedef for the identity preconditioner.
     */
    typedef dealii::PETScWrappers::PreconditionNone PreconditionIdentity;

    /**
     * Typedef for the Incomplete LU decomposition preconditioner.
     */
    typedef dealii::PETScWrappers::PreconditionILU PreconditionILU;

    /**
     * Typedef for the Incomplete Jacobi decomposition preconditioner.
     */
    typedef dealii::PETScWrappers::PreconditionJacobi PreconditionJacobi;

    /**
     * Typedef for the SOR preconditioner.
     */
    typedef dealii::PETScWrappers::PreconditionSOR PreconditionSOR;

    /**
     * Typedef for the SSOR preconditioner.
     */
    typedef dealii::PETScWrappers::PreconditionSSOR PreconditionSSOR;
  }

}
DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_WITH_PETSC

#ifdef DEAL_II_WITH_TRILINOS

#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_vector.h>

DEAL_II_NAMESPACE_OPEN

/**
 * A namespace in which the wrappers to the Trilinos linear algebra classes
 * are typedef'ed to generic names. There are similar namespaces
 * LinearAlgebraDealII and LinearAlgebraPETSc for typedefs to deal.II's own
 * classes and classes that interface with PETSc.
 *
 * @ingroup GenericLinearAlgebra
 * @author Timo Heister 2008, Jean-Paul Pelteret 2016
 */
namespace LinearAlgebraTrilinos
{
  /**
   * Typedef for the vector type used.
   */
  typedef dealii::TrilinosWrappers::Vector Vector;

  /**
   * Typedef for the type used to describe vectors that consist of multiple
   * blocks.
   */
  typedef dealii::TrilinosWrappers::BlockVector BlockVector;

  // ----------------------------------------------------
  // Below are the common elements between this namespace
  // and any namespaces nested within this one
  // ----------------------------------------------------

  /**
   * Typedef for the sparse matrix type used.
   */
  typedef dealii::TrilinosWrappers::SparseMatrix SparseMatrix;

  /**
   * Typedef for dynamic sparsity patterns.
   */
  typedef dealii::DynamicSparsityPattern DynamicSparsityPattern;

  /**
   * Typedef for sparsity patterns.
   */
  typedef dealii::TrilinosWrappers::SparsityPattern SparsityPattern;

  /**
   * Typedef for the type used to describe sparse matrices that consist of
   * multiple blocks.
   */
  typedef dealii::TrilinosWrappers::BlockSparseMatrix BlockSparseMatrix;

  /**
   * Typedef for block dynamic sparsity patterns.
   */
  typedef dealii::BlockDynamicSparsityPattern BlockDynamicSparsityPattern;

  /**
   * Typedef for block sparsity patterns.
   */
  typedef dealii::TrilinosWrappers::BlockSparsityPattern BlockSparsityPattern;

  /**
   * Typedef for the base class for linear solvers.
   *
   * This typedef is useful when wanting to create a generic solver interface
   * wherein one can switch between linear solvers. For example:
   * @code
   * std_cxx11::unique_ptr<LA::SolverBase> solver;
   * if (parameters.solver_type == "CG")
   * {
   *   solver.reset(new LA::SolverCG(solver_control));
   * }
   * @endcode
   */
  typedef dealii::TrilinosWrappers::SolverBase SolverBase;

  /**
   * Typedef for the BiCGStab linear solver.
   */
  typedef dealii::TrilinosWrappers::SolverBicgstab SolverBicgstab;

  /**
   * Typedef for the CG linear solver.
   */
  typedef dealii::TrilinosWrappers::SolverCG SolverCG;

  /**
   * Typedef for the CG Squared linear solver.
   */
  typedef dealii::TrilinosWrappers::SolverCGS SolverCGS;

  /**
   * Typedef for the direct linear solver.
   */
  typedef dealii::TrilinosWrappers::SolverDirect SolverDirect;

  /**
   * Typedef for the GMRES linear solver.
   */
  typedef dealii::TrilinosWrappers::SolverGMRES SolverGMRES;

  /**
   * Typedef for the TFQMR linear solver.
   */
  typedef dealii::TrilinosWrappers::SolverTFQMR SolverTFQMR;

  /**
   * Typedef for the base class for preconditioners.
   *
   * This typedef is useful when wanting to create a generic solver interface
   * wherein one can switch between preconditioners. For example:
   * @code
   * std_cxx11::unique_ptr<LA::PreconditionBase> preconditioner;
   * if (parameters.preconditioner_type == "Jacobi")
   * {
   *   LA::MPI::PreconditionJacobi* ptr_prec
   *     = new LA::MPI::PreconditionJacobi ();
   *
   *   LA::MPI::PreconditionJacobi::AdditionalData
   *     additional_data (parameters.preconditioner_relaxation);
   *
   *   ptr_prec->initialize(tangent_matrix.block(u_block,u_block),
   *                        additional_data);
   *   preconditioner.reset(ptr_prec);
   * }
   * @endcode
   */
  typedef dealii::TrilinosWrappers::PreconditionBase PreconditionBase;

  /**
   * Typedef for the AMG preconditioner type.
   */
  typedef dealii::TrilinosWrappers::PreconditionAMG PreconditionAMG;

  /**
   * Typedef for the Incomplete Cholesky preconditioner.
   */
  typedef dealii::TrilinosWrappers::PreconditionIC PreconditionIC;

  /**
   * Typedef for the identity preconditioner.
   */
  typedef dealii::TrilinosWrappers::PreconditionIdentity PreconditionIdentity;

  /**
   * Typedef for the Incomplete LU decomposition preconditioner.
   */
  typedef dealii::TrilinosWrappers::PreconditionILU PreconditionILU;

  /**
   * Typedef for the Incomplete Jacobi decomposition preconditioner.
   */
  typedef dealii::TrilinosWrappers::PreconditionJacobi PreconditionJacobi;

  /**
   * Typedef for the SOR preconditioner
   */
  typedef dealii::TrilinosWrappers::PreconditionSOR PreconditionSOR;

  /**
   * Typedef for the SSOR preconditioner
   */
  typedef dealii::TrilinosWrappers::PreconditionSSOR PreconditionSSOR;

  /**
   * A namespace with typedefs to generic names for parallel Trilinos linear
   * algebra objects.
   *
   * @ingroup GenericLinearAlgebra
   * @author Timo Heister 2008, Jean-Paul Pelteret 2016
   */
  namespace MPI
  {
    /**
     * Typedef for the vector type used.
     */
    typedef dealii::TrilinosWrappers::MPI::Vector Vector;

    /**
     * Typedef for the type used to describe vectors that consist of multiple
     * blocks.
     */
    typedef dealii::TrilinosWrappers::MPI::BlockVector BlockVector;

    // ----------------------------------------------------
    // Below are the common elements between this namespace
    // and the namespaces that it is nested within
    // ----------------------------------------------------

    /**
     * Typedef for the sparse matrix type used.
     */
    typedef dealii::TrilinosWrappers::SparseMatrix SparseMatrix;

    /**
     * Typedef for dynamic sparsity patterns.
     */
    typedef dealii::DynamicSparsityPattern DynamicSparsityPattern;

    /**
     * Typedef for sparsity patterns.
     */
    typedef dealii::TrilinosWrappers::SparsityPattern SparsityPattern;

    /**
     * Typedef for the type used to describe sparse matrices that consist of
     * multiple blocks.
     */
    typedef dealii::TrilinosWrappers::BlockSparseMatrix BlockSparseMatrix;

    /**
     * Typedef for block dynamic sparsity patterns.
     */
    typedef dealii::BlockDynamicSparsityPattern BlockDynamicSparsityPattern;

    /**
     * Typedef for block sparsity patterns.
     */
    typedef dealii::TrilinosWrappers::BlockSparsityPattern BlockSparsityPattern;

    /**
     * Typedef for the base class for linear solvers.
     *
     * This typedef is useful when wanting to create a generic solver interface
     * wherein one can switch between linear solvers. For example:
     * @code
     * std_cxx11::unique_ptr<LA::SolverBase> solver;
     * if (parameters.solver_type == "CG")
     * {
     *   solver.reset(new LA::SolverCG(solver_control));
     * }
     * @endcode
     */
    typedef dealii::TrilinosWrappers::SolverBase SolverBase;

    /**
     * Typedef for the BiCGStab linear solver.
     */
    typedef dealii::TrilinosWrappers::SolverBicgstab SolverBicgstab;

    /**
     * Typedef for the CG linear solver.
     */
    typedef dealii::TrilinosWrappers::SolverCG SolverCG;

    /**
     * Typedef for the CG Squared linear solver.
     */
    typedef dealii::TrilinosWrappers::SolverCGS SolverCGS;

    /**
     * Typedef for the direct linear solver.
     */
    typedef dealii::TrilinosWrappers::SolverDirect SolverDirect;

    /**
     * Typedef for the GMRES linear solver.
     */
    typedef dealii::TrilinosWrappers::SolverGMRES SolverGMRES;

    /**
     * Typedef for the TFQMR linear solver.
     */
    typedef dealii::TrilinosWrappers::SolverTFQMR SolverTFQMR;

    /**
     * Typedef for the base class for preconditioners.
     *
     * This typedef is useful when wanting to create a generic solver interface
     * wherein one can switch between preconditioners. For example:
     * @code
     * std_cxx11::unique_ptr<LA::PreconditionBase> preconditioner;
     * if (parameters.preconditioner_type == "Jacobi")
     * {
     *   LA::MPI::PreconditionJacobi* ptr_prec
     *     = new LA::MPI::PreconditionJacobi ();
     *
     *   LA::MPI::PreconditionJacobi::AdditionalData
     *     additional_data (parameters.preconditioner_relaxation);
     *
     *   ptr_prec->initialize(tangent_matrix.block(u_block,u_block),
     *                        additional_data);
     *   preconditioner.reset(ptr_prec);
     * }
     * @endcode
     */
    typedef dealii::TrilinosWrappers::PreconditionBase PreconditionBase;

    /**
     * Typedef for the AMG preconditioner type.
     */
    typedef dealii::TrilinosWrappers::PreconditionAMG PreconditionAMG;

    /**
     * Typedef for the Incomplete Cholesky preconditioner.
     */
    typedef dealii::TrilinosWrappers::PreconditionIC PreconditionIC;

    /**
     * Typedef for the identity preconditioner.
     */
    typedef dealii::TrilinosWrappers::PreconditionIdentity PreconditionIdentity;

    /**
     * Typedef for the Incomplete LU decomposition preconditioner.
     */
    typedef dealii::TrilinosWrappers::PreconditionILU PreconditionILU;

    /**
     * Typedef for the Incomplete Jacobi decomposition preconditioner.
     */
    typedef dealii::TrilinosWrappers::PreconditionJacobi PreconditionJacobi;

    /**
     * Typedef for the SOR preconditioner
     */
    typedef dealii::TrilinosWrappers::PreconditionSOR PreconditionSOR;

    /**
     * Typedef for the SSOR preconditioner
     */
    typedef dealii::TrilinosWrappers::PreconditionSSOR PreconditionSSOR;
  }

}

/*@}*/

DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_WITH_TRILINOS



#endif

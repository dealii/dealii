// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_trilinos_tpetra_solver_direct_h
#define dealii_trilinos_tpetra_solver_direct_h

#include <deal.II/base/config.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/trilinos_tpetra_types.h>

#include <string>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA
#  ifdef DEAL_II_TRILINOS_WITH_AMESOS2

#    include <deal.II/base/types.h>

#    include <deal.II/lac/la_parallel_vector.h>
#    include <deal.II/lac/solver_control.h>
#    include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#    include <Amesos2.hpp>
#    include <Teuchos_ConfigDefs.hpp>
#    include <Teuchos_ParameterList.hpp>
#    include <Teuchos_RCPDecl.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
    // forward declarations
#    ifndef DOXYGEN
    template <typename Number, typename MemorySpace>
    class SparseMatrix;
#    endif

    /// The chosen Solver is not supported or configured with Amesos2.
    DeclException1(
      ExcTrilinosAmesos2SolverUnsupported,
      std::string,
      << "You tried to select the solver type <" << arg1 << ">\n"
      << "but this solver is not supported by Trilinos/Amesos2\n"
      << "due to one of the following reasons:\n"
      << "* This solver does not exist\n"
      << "* This solver is not (yet) supported by Trilinos/Amesos2\n"
      << "* Trilinos/Amesos2 was not configured for its use.");


    /**
     * The base class for all direct solvers based on the Amesos2 package
     * of Trilinos.
     *
     * @ingroup TpetraWrappers
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class SolverDirectBase
    {
    public:
      /**
       * Destructor.
       */
      virtual ~SolverDirectBase() = default;

      /**
       * Initializes the direct solver for the matrix <tt>A</tt> and creates a
       * factorization for it with the package chosen from the additional
       * data structure. Note that there is no need for a preconditioner
       * here and solve() is not called.
       */
      void
      initialize(const SparseMatrix<Number, MemorySpace> &A);

      /**
       * Solve the linear system <tt>Ax=b</tt> based on the
       * package set in initialize(). Note the matrix is not refactorized during
       * this call.
       */
      void
      solve(Vector<Number, MemorySpace>       &x,
            const Vector<Number, MemorySpace> &b);

      /**
       * Solve the linear system <tt>Ax=b</tt>. Creates a factorization of the
       * matrix with the package chosen from the additional data structure and
       * performs the solve. Note that there is no need for a preconditioner
       * here.
       */
      void
      solve(const SparseMatrix<Number, MemorySpace> &A,
            Vector<Number, MemorySpace>             &x,
            const Vector<Number, MemorySpace>       &b);


      /**
       * Access to object that controls convergence.
       */
      SolverControl &
      control() const;

      /**
       * Exception
       */
      DeclException1(ExcTrilinosError,
                     int,
                     << "An error with error number " << arg1
                     << " occurred while calling a Trilinos function");

    protected:
      /**
       * Constructor. Takes the solver control object and name
       * and creates the solver.
       */
      SolverDirectBase(SolverControl &cn, const std::string &solver_type);
      /**
       * Actually performs the operations for solving the linear system,
       * including the factorization and forward and backward substitution.
       */
      void
      do_solve();

      /**
       * Reference to the object that controls convergence of the iterative
       * solver. In fact, for these Trilinos wrappers, Trilinos does so itself,
       * but we copy the data from this object before starting the solution
       * process, and copy the data back into it afterwards.
       */
      SolverControl &solver_control;

      /**
       * A structure that contains the Trilinos solver object.
       */
      Teuchos::RCP<
        Amesos2::Solver<TpetraTypes::MatrixType<Number, MemorySpace>,
                        TpetraTypes::MultiVectorType<Number, MemorySpace>>>
        solver;

      /*
       * The set solver type to be handed to the solver factory of Amesos2.
       */
      std::string solver_type;

      /**
       * An optional Teuchos::ParameterList for fine tuning the solver.
       * Please refer to the Amesos2 manual to see which parameters
       * to set for each individual solver.
       */
      Teuchos::ParameterList parameter_list;
    }; // Base



    /**
     * A general purpose class to support any solver that the Amesos2 package
     * of Trilinos provides.
     *
     * Notes for users switching from TrilinosWrappers to TpetraWrappers:
     * The general interface of this class is kept identical to the Amesos(1)
     * wrapper with the notable exception of being templated.
     *
     * A further addition is the option to fine tune each solver with a
     * Teuchos::ParameterList, to change default parameters.
     *
     * @ingroup TpetraWrappers
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class SolverDirect : public SolverDirectBase<Number, MemorySpace>
    {
    public:
      /**
       * A structure whose member variables describe details of the algorithm
       * to be employed by the surrounding class. An object of this type can be
       * used to initialize an object of the surrounding class.
       */
      struct AdditionalData
      {
        AdditionalData(const std::string &solver_name);


        /**
         * Set the solver type (for third party solver support of Trilinos
         * Amesos package). Current possibilities are:
         * <ul>
         * <li>  "Basker" </li>
         * <li>  "Cholmod" </li>
         * <li>  "cuSOLVER" </li>
         * <li>  "KLU2" </li>
         * <li>  "LAPACK" </li>
         * <li>  "MUMPS" </li>
         * <li>  "PardisoMKL" </li>
         * <li>  "ShyLUBasker" </li>
         * <li>  "SuperLU" </li>
         * <li>  "SuperLU_DIST" </li>
         * <li>  "SuperLU_MT" </li>
         * <li>  "Tacho" </li>
         * <li>  "UMFPACK" </li>
         * </ul>
         * Note that the availability of these solvers in deal.II depends on
         * which solvers were set when configuring Trilinos.
         * Additionally, Amesos2 may add support for solvers not listed here
         * that can nevertheless be used if configured.
         */
        std::string solver_name;
      };

      /**
       * Constructor. Takes the solver control object and creates the solver.
       */
      SolverDirect(SolverControl        &cn,
                   const AdditionalData &additional_data = AdditionalData());

      /**
       * Set a parameter list to fine tune the solver.
       * The valid parameters depend on the solver you use
       * and the possible parameters can be found in the
       * documentation of Amesos2.
       */
      void
      set_pararameter_list(Teuchos::ParameterList &parameter_list);
    }; // SolverDirect


    /**
     * A wrapper class for the solver KLU2 that works in serial and parallel.
     * This solver is part of Amesos2 and enabled by default.
     *
     * The AdditionalData structure allows to pass options specific to this
     * solver and the default will result in the same solver as constructing
     * SolverDirect with solver_name KLU2.
     *
     * @ingroup TpetraWrappers
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class SolverDirectKLU2 : public SolverDirectBase<Number, MemorySpace>
    {
    public:
      /**
       * A structure whose member variables describe details of the algorithm
       * to be employed by the surrounding class. An object of this type can be
       * used to initialize an object of the surrounding class.
       */
      struct AdditionalData
      {
        AdditionalData(const std::string &transpose_mode       = "NOTRANS",
                       const bool         symmetric_mode       = false,
                       const bool         equilibrate_matrix   = true,
                       const std::string &column_permutation   = "COLAMD",
                       const std::string &iterative_refinement = "NO");
        /**
         * Decide which system to solve
         * "NOTRANS": Ax=b (default)
         * "TRANS": A^Tx=b
         * "CONJ": A^*x=b
         */
        std::string transpose_mode;

        /**
         * Is A symmetric or not?
         */
        bool symmetric_mode;

        /**
         * Equilibrate matrix before solving? (default: true)
         */
        bool equilibrate_matrix;

        /**
         * Choose an ordering strategy
         * "COLAMD": approximate minimum degree column ordering (default).
         * "NATURAL": natural ordering.
         * "MMD_ATA": minimum degree ordering on the structure of A^TA.
         * "MMD_AT_PLUS_A": minimum degree ordering on the structure of A^T+A.
         */
        std::string column_permutation;

        /**
         * Perform iterative refinement?
         * "NO": no (default)
         * "SINGLE": yes, use single precision for the residual
         * "DOUBLE": yes, use double precision for the residual
         * "EXTRA": ??
         */
        std::string iterative_refinement;
      };

      /**
       * Constructor. Takes the solver control object and creates the solver.
       */
      SolverDirectKLU2(
        SolverControl        &cn,
        const AdditionalData &additional_data = AdditionalData());
    }; // KLU2

  } // namespace TpetraWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_TRILINOS_WITH_AMESOS2
#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_TPETRA

#endif

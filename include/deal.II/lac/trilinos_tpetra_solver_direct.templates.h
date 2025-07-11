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

#ifndef dealii_trilinos_tpetra_solver_templates_h
#define dealii_trilinos_tpetra_solver_templates_h

#include <deal.II/base/config.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/trilinos_tpetra_types.h>

#include <string>


#ifdef DEAL_II_TRILINOS_WITH_TPETRA
#  ifdef DEAL_II_TRILINOS_WITH_AMESOS2

#    include <deal.II/base/conditional_ostream.h>
#    include <deal.II/base/template_constraints.h>

#    include <deal.II/lac/trilinos_tpetra_solver_direct.h>


DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
    /* ---------------------- SolverDirectBase ------------------------ */

    template <typename Number, typename MemorySpace>
    SolverDirectBase<Number, MemorySpace>::SolverDirectBase(
      SolverControl     &cn,
      const std::string &solver_type)
      : solver_control(cn)
      , solver_type(solver_type)
    {
      AssertThrow(Amesos2::query(solver_type),
                  ExcTrilinosAmesos2SolverUnsupported(solver_type));
    }



    template <typename Number, typename MemorySpace>
    SolverControl &
    SolverDirectBase<Number, MemorySpace>::control() const
    {
      return solver_control;
    }



    template <typename Number, typename MemorySpace>
    void
    SolverDirectBase<Number, MemorySpace>::initialize(
      const SparseMatrix<Number, MemorySpace> &A)
    {
      // Allocate the Amesos2 solver with the concrete solver, if possible.
      solver =
        Amesos2::create<TpetraTypes::MatrixType<Number, MemorySpace>,
                        TpetraTypes::MultiVectorType<Number, MemorySpace>>(
          solver_type, A.trilinos_rcp());

      solver->setParameters(Teuchos::rcpFromRef(parameter_list));
      // Now do the actual factorization, which is a two step procedure.
      // The symbolic factorization determines the structure of the inverse,
      // while the numeric factorization does to actual computation of L and U
      solver->symbolicFactorization();

      solver->numericFactorization();
    }



    template <typename Number, typename MemorySpace>
    void
    SolverDirectBase<Number, MemorySpace>::solve(
      Vector<Number, MemorySpace>       &x,
      const Vector<Number, MemorySpace> &b)
    {
      // Assign the empty solution vector
      solver->setX(x.trilinos_rcp());

      // Assign the RHS vector
      solver->setB(b.trilinos_rcp());

      solver->solve();

      // Finally, force the SolverControl object to report convergence
      solver_control.check(0, 0);
      if (solver_control.last_check() != SolverControl::success)
        AssertThrow(false,
                    SolverControl::NoConvergence(solver_control.last_step(),
                                                 solver_control.last_value()));
    }



    template <typename Number, typename MemorySpace>
    void
    SolverDirectBase<Number, MemorySpace>::do_solve()
    {
      // Set the parameter list. If it is empty defaults will be used.
      solver->setParameters(Teuchos::rcpFromRef(parameter_list));
      // Now do the actual factorization, which is a two step procedure.
      // The symbolic factorization determines the structure of the inverse,
      // while the numeric factorization does to actual computation of L and U
      solver->symbolicFactorization();

      solver->numericFactorization();

      solver->solve();

      // Finally, force the SolverControl object to report convergence
      solver_control.check(0, 0);

      if (solver_control.last_check() != SolverControl::success)
        AssertThrow(false,
                    SolverControl::NoConvergence(solver_control.last_step(),
                                                 solver_control.last_value()));
    }



    template <typename Number, typename MemorySpace>
    void
    SolverDirectBase<Number, MemorySpace>::solve(
      const SparseMatrix<Number, MemorySpace> &A,
      Vector<Number, MemorySpace>             &x,
      const Vector<Number, MemorySpace>       &b)
    {
      solver =
        Amesos2::create<TpetraTypes::MatrixType<Number, MemorySpace>,
                        TpetraTypes::MultiVectorType<Number, MemorySpace>>(
          solver_type, A.trilinos_rcp(), x.trilinos_rcp(), b.trilinos_rcp());
      do_solve();
    }



    template <typename Number, typename MemorySpace>
    SolverDirect<Number, MemorySpace>::AdditionalData::AdditionalData(
      const std::string &solver_name)
      : solver_name(solver_name)
    {}



    template <typename Number, typename MemorySpace>
    SolverDirect<Number, MemorySpace>::SolverDirect(SolverControl        &cn,
                                                    const AdditionalData &ad)
      : SolverDirectBase<Number, MemorySpace>(cn, ad.solver_name)
    {}



    template <typename Number, typename MemorySpace>
    void
    SolverDirect<Number, MemorySpace>::set_pararameter_list(
      Teuchos::ParameterList &parameter_list)
    {
      this->parameter_list.setParameters(parameter_list);
    }



    template <typename Number, typename MemorySpace>
    SolverDirectKLU2<Number, MemorySpace>::AdditionalData::AdditionalData(
      const std::string &transpose_mode,
      const bool         symmetric_mode,
      const bool         equilibrate_matrix,
      const std::string &column_permutation,
      const std::string &iterative_refinement)
      : transpose_mode(transpose_mode)
      , symmetric_mode(symmetric_mode)
      , equilibrate_matrix(equilibrate_matrix)
      , column_permutation(column_permutation)
      , iterative_refinement(iterative_refinement)
    {}



    template <typename Number, typename MemorySpace>
    SolverDirectKLU2<Number, MemorySpace>::SolverDirectKLU2(
      SolverControl        &cn,
      const AdditionalData &ad)
      : SolverDirectBase<Number, MemorySpace>(cn, "KLU2")
    {
      this->parameter_list               = Teuchos::ParameterList("Amesos2");
      Teuchos::ParameterList klu2_params = this->parameter_list.sublist("KLU2");
      klu2_params.set("Trans", ad.transpose_mode);
      klu2_params.set("Equil", ad.equilibrate_matrix);
      klu2_params.set("IterRefine", ad.iterative_refinement);
      klu2_params.set("SymmetricMode", ad.symmetric_mode);
      klu2_params.set("ColPerm", ad.column_permutation);
    }

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

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check TrilinosWrappers::NOXSolver by solving f(x) = x^2 = 0 with initial
// condition x=2.
//
// This test runs the same solver twice: Once via the deal.II interface to NOX
// in the TrilinosWrappers::NOXSolver class and once using the native,
// Epetra-based interface to NOX. The output should of course be the same.

#include <deal.II/base/mpi.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/trilinos/nox.h>

// as reference solution
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <NOX_Epetra_Group.H>
#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_LinearSystem_AztecOO.H>
#include <NOX_Solver_Factory.H>
#include <NOX_Solver_Generic.H>
#include <NOX_StatusTest_Combo.H>
#include <NOX_StatusTest_MaxIters.H>
#include <NOX_StatusTest_NormF.H>
#include <NOX_StatusTest_RelativeNormF.H>

#include "../tests.h"


int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  using Number     = double;
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  // set up solver control
  const unsigned int n_max_iterations  = 100;
  const double       abs_tolerance     = 1e-9;
  const double       rel_tolerance     = 1e-5;
  const double       lin_rel_tolerance = 1e-3;

  TrilinosWrappers::NOXSolver<VectorType>::AdditionalData additional_data(
    n_max_iterations, abs_tolerance, rel_tolerance);

  // set up parameters
  Teuchos::RCP<Teuchos::ParameterList> non_linear_parameters =
    Teuchos::rcp(new Teuchos::ParameterList);

  non_linear_parameters->set("Nonlinear Solver", "Line Search Based");
  non_linear_parameters->sublist("Printing").set("Output Information", 15);
  non_linear_parameters->sublist("Direction").set("Method", "Newton");
  non_linear_parameters->sublist("Direction")
    .sublist("Newton")
    .sublist("Linear Solver")
    .set("Tolerance", lin_rel_tolerance);
  non_linear_parameters->sublist("Line Search").set("Method", "Polynomial");

  /*
   * First check: Use the deal.II-based wrappers for the test.
   */
  {
    // set up solver
    TrilinosWrappers::NOXSolver<VectorType> solver(additional_data,
                                                   non_linear_parameters);

    // ... helper functions
    double J = 0.0;

    solver.residual = [](const auto &src, auto &dst) {
      // compute residual
      dst[0] = src[0] * src[0];
    };

    solver.setup_jacobian = [&](const auto &src) {
      // compute Jacobian
      J = 2.0 * src[0];
    };

    solver.apply_jacobian = [&](const auto &src, auto &dst) {
      // solve with Jacobian
      dst[0] = src[0] * J;
    };

    solver.solve_with_jacobian = [&](const auto &src, auto &dst, const auto) {
      // solve with Jacobian
      dst[0] = src[0] / J;
    };

    // initial guess
    VectorType solution(1);
    solution[0] = 2.0;

    // solve with the given initial guess
    solver.solve(solution);

    deallog << "The solution is: " << solution[0] << std::endl;
  }

  /*
   * Second check: Run the same test through the native NOX interfaces.
   */
  {
    class NoxInterface : public NOX::Epetra::Interface::Required,
                         public NOX::Epetra::Interface::Jacobian
    {
    public:
      bool
      computeF(const Epetra_Vector                       &x,
               Epetra_Vector                             &f,
               NOX::Epetra::Interface::Required::FillType F) override
      {
        (void)F;

        f[0] = x[0] * x[0];

        return true;
      }

      bool
      computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac) override
      {
        auto jac = dynamic_cast<Epetra_CrsMatrix *>(&Jac);

        AssertThrow(jac, ExcNotImplemented());

        jac->PutScalar(2.0 * x[0]);

        return true;
      }
    };


    // convert data structures to Epetra structures
    IndexSet is(1);
    is.add_index(0);

    TrilinosWrappers::MPI::Vector solution(is, MPI_COMM_WORLD);
    solution[0] = 2.0;

    TrilinosWrappers::SparsityPattern dsp(is, MPI_COMM_WORLD);
    dsp.add(0, 0);
    dsp.compress();

    TrilinosWrappers::SparseMatrix system_matrix;
    system_matrix.reinit(dsp);


    // setup linear system and group

    Epetra_Vector InitialGuess(Copy, solution.trilinos_vector(), 0);

    auto A =
      Teuchos::rcp(new Epetra_CrsMatrix(system_matrix.trilinos_matrix()));

    Teuchos::ParameterList &printParams =
      non_linear_parameters->sublist("Printing");
    Teuchos::ParameterList &lsParams =
      non_linear_parameters->sublist("Newton").sublist("Linear Solver");

    Teuchos::RCP<NoxInterface> interface = Teuchos::rcp(new NoxInterface());

    Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;

    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(
        printParams, lsParams, iReq, iJac, A, InitialGuess));

    NOX::Epetra::Vector              noxInitGuess(InitialGuess, NOX::DeepCopy);
    Teuchos::RCP<NOX::Epetra::Group> group = Teuchos::rcp(
      new NOX::Epetra::Group(printParams, iReq, noxInitGuess, linSys));

    // setup solver control
    auto check =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

    if (additional_data.abs_tol > 0.0)
      {
        const auto additional_data_norm_f_abs =
          Teuchos::rcp(new NOX::StatusTest::NormF(additional_data.abs_tol));
        check->addStatusTest(additional_data_norm_f_abs);
      }

    if (additional_data.rel_tol > 0.0)
      {
        const auto additional_data_norm_f_rel = Teuchos::rcp(
          new NOX::StatusTest::RelativeNormF(additional_data.rel_tol));
        check->addStatusTest(additional_data_norm_f_rel);
      }

    if (additional_data.max_iter > 0)
      {
        const auto additional_data_max_iterations =
          Teuchos::rcp(new NOX::StatusTest::MaxIters(additional_data.max_iter));
        check->addStatusTest(additional_data_max_iterations);
      }

    Teuchos::RCP<NOX::Solver::Generic> solver =
      NOX::Solver::buildSolver(group, check, non_linear_parameters);
    auto status = solver->solve();

    AssertThrow(status == NOX::StatusTest::Converged, ExcInternalError());

    deallog << "The solution is: " << group->getX().norm() << std::endl;
  }
}

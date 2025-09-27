// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector.h>

#include <deal.II/sundials/arkode.h>

#include "../tests.h"


// Test implicit-explicit time stepper. solve_jacobian + custom linear solver +
// custom mass solver

/**
 * Test problem inspired by linear 1D FE problem with two unknowns u = [u1 u2]:
 *
 * M*u' = K*u + f
 *
 * where the mass matrix M = [2/3 1/3; 1/3 2/3]
 * and the stiffness matrix K = [1/2 -1/2; 1/2 -1/2]
 * and the force vector f = [1,2].
 * The initial values is chosen as u0 = [1,0].
 *
 * The implicit part is taken as K*u while f is explicit.
 *
 * Note that this is not the classical second-order mechanical system and thus
 * the term "mass matrix" may be misleading. It is still used here to be in line
 * with SUNDIALS nomenclature.
 */
int
main()
{
  initlog();

  using VectorType = Vector<double>;

  ParameterHandler                             prm;
  SUNDIALS::ARKode<VectorType>::AdditionalData data;
  data.add_parameters(prm);

  if (false)
    {
      std::ofstream ofile(SOURCE_DIR "/arkode_08_in.prm");
      prm.print_parameters(ofile, ParameterHandler::ShortPRM);
      ofile.close();
    }

  std::ifstream ifile(SOURCE_DIR "/arkode_08_in.prm");
  prm.parse_input(ifile);

  SUNDIALS::ARKode<VectorType> ode(data);

  // Explicit jacobian = stiffness matrix
  FullMatrix<double> K(2, 2);
  K(0, 0) = K(1, 1) = 0.5;
  K(1, 0) = K(0, 1) = -0.5;

  // mass matrix
  FullMatrix<double> M(2, 2);
  M(0, 0) = M(1, 1) = 2.0 / 3;
  M(1, 0) = M(0, 1) = 1.0 / 3;

  ode.implicit_function = [&](double, const VectorType &y, VectorType &ydot) {
    K.vmult(ydot, y);
  };


  ode.explicit_function = [&](double, const VectorType &y, VectorType &ydot) {
    ydot[0] = 1;
    ydot[1] = 2;
  };

  ode.jacobian_times_vector = [&](const VectorType &v,
                                  VectorType       &Jv,
                                  double            t,
                                  const VectorType &y,
                                  const VectorType &fy) { K.vmult(Jv, v); };

  const auto solve_function =
    [&](SUNDIALS::SundialsOperator<VectorType>       &op,
        SUNDIALS::SundialsPreconditioner<VectorType> &prec,
        VectorType                                   &x,
        const VectorType                             &b,
        double                                        tol) {
      ReductionControl     control;
      SolverCG<VectorType> solver_cg(control);
      solver_cg.solve(op, x, b, prec);
    };

  ode.solve_linearized_system = solve_function;

  ode.solve_mass = solve_function;

  FullMatrix<double> M_inv(2, 2);

  ode.mass_preconditioner_solve =
    [&](double t, const VectorType &r, VectorType &z, double gamma, int lr) {
      LogStream::Prefix prefix("mass_preconditioner_solve");
      deallog << "applied" << std::endl;
      M_inv.vmult(z, r);
    };

  ode.mass_preconditioner_setup = [&](double t) {
    LogStream::Prefix prefix("mass_preconditioner_setup");
    deallog << "applied" << std::endl;
    M_inv.invert(M);
  };

  ode.mass_times_vector = [&](const double      t,
                              const VectorType &v,
                              VectorType       &Mv) { M.vmult(Mv, v); };


  ode.output_step =
    [&](const double t, const VectorType &sol, const unsigned int step_number) {
      deallog << t << ' ' << sol[0] << ' ' << sol[1] << std::endl;
    };

  Vector<double> y(2);
  y[0] = 1;
  y[1] = 0;
  ode.solve_ode(y);
}

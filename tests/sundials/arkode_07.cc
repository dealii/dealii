//-----------------------------------------------------------
//
//    Copyright (C) 2020 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//-----------------------------------------------------------

#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector.h>

#include <deal.II/sundials/arkode.h>

#include <arkode/arkode_arkstep.h>

#include "../tests.h"


// Test implicit-explicit time stepper. jac_times_vector/setup +
// custom linear solver + custom preconditioner supplied through SUNDIALS

/**
 * This test problem is called "brusselator", and is a typical benchmark for
 * ODE solvers. This problem has 3 dependent variables u, v and w, that depend
 * on the independent variable t via the IVP system
 *
 * du/dt = a − (w + 1)u + v u^2
 * dv/dt = w u − v u^2
 * dw/dt = (b − w)/eps -w u
 *
 * We integrate over the interval 0 ≤ t ≤ 10, with the initial conditions
 *
 * u(0) = 3.9, v(0) = 1.1, w(0) = 2.8,
 *
 * and parameters
 *
 * a = 1.2, b = 2.5, and eps = 10−5
 *
 * The implicit part only contains the stiff part of the problem (the part with
 * eps in right hand side of the third equation).
 */
int
main(int argc, char **argv)
{
  initlog();
  // restrict output to highest level
  deallog.depth_file(1);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  using VectorType = Vector<double>;

  ParameterHandler                             prm;
  SUNDIALS::ARKode<VectorType>::AdditionalData data;
  data.add_parameters(prm);

  if (false)
    {
      std::ofstream ofile(SOURCE_DIR "/arkode_07.prm");
      prm.print_parameters(ofile, ParameterHandler::ShortText);
      ofile.close();
    }

  std::ifstream ifile(SOURCE_DIR "/arkode_07.prm");
  prm.parse_input(ifile);

  SUNDIALS::ARKode<VectorType> ode(data);

  // Parameters
  double u0 = 3.9, v0 = 1.1, w0 = 2.8, a = 1.2, b = 2.5, eps = 1e-5;
  // Explicit jacobian.
  FullMatrix<double> J(3, 3);

  ode.implicit_function =
    [&](double, const VectorType &y, VectorType &ydot) -> int {
    ydot[0] = 0;
    ydot[1] = 0;
    ydot[2] = (b - y[2]) / eps;
    return 0;
  };


  ode.explicit_function =
    [&](double, const VectorType &y, VectorType &ydot) -> int {
    ydot[0] = a - (y[2] + 1) * y[0] + y[1] * y[0] * y[0];
    ydot[1] = y[2] * y[0] - y[1] * y[0] * y[0];
    ydot[2] = -y[2] * y[0];
    return 0;
  };


  ode.jacobian_times_setup =
    [&](realtype t, const VectorType &y, const VectorType &fy) -> int {
    J       = 0;
    J(2, 2) = -1.0 / eps;
    return 0;
  };

  ode.jacobian_times_vector = [&](const VectorType &v,
                                  VectorType &      Jv,
                                  double            t,
                                  const VectorType &y,
                                  const VectorType &fy) -> int {
    J.vmult(Jv, v);
    return 0;
  };

  ode.solve_linearized_system =
    [&](SUNDIALS::SundialsOperator<VectorType> &      op,
        SUNDIALS::SundialsPreconditioner<VectorType> &prec,
        VectorType &                                  x,
        const VectorType &                            b,
        double                                        tol) -> int {
    ReductionControl     control;
    SolverCG<VectorType> solver_cg(control);
    solver_cg.solve(op, x, b, prec);
    return 0;
  };

  ode.jacobian_preconditioner_setup = [&](double            t,
                                          const VectorType &y,
                                          const VectorType &fy,
                                          int               jok,
                                          int &             jcur,
                                          double            gamma) -> int {
    deallog << "jacobian_preconditioner_setup called\n";
    return 0;
  };

  ode.jacobian_preconditioner_solve = [&](double            t,
                                          const VectorType &y,
                                          const VectorType &fy,
                                          const VectorType &r,
                                          VectorType &      z,
                                          double            gamma,
                                          double            delta,
                                          int               lr) -> int {
    deallog << "jacobian_preconditioner_solve called\n";
    z = r;
    return 0;
  };


  ode.output_step = [&](const double       t,
                        const VectorType & sol,
                        const unsigned int step_number) -> int {
    deallog << t << " " << sol[0] << " " << sol[1] << " " << sol[2]
            << std::endl;
    return 0;
  };

  // after 5.2.0 a special interpolation mode should be used for stiff problems
#if DEAL_II_SUNDIALS_VERSION_GTE(5, 2, 0)
  ode.custom_setup = [&](void *arkode_mem) {
    ARKStepSetInterpolantType(arkode_mem, ARK_INTERP_LAGRANGE);
  };
#endif

  Vector<double> y(3);
  y[0] = u0;
  y[1] = v0;
  y[2] = w0;
  ode.solve_ode(y);
  return 0;
}

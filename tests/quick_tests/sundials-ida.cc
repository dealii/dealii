//-----------------------------------------------------------
//
//    Copyright (C) 2015 - 2020 by the deal2lkit authors
//
//    This file is part of the deal2lkit library.
//
//    The deal2lkit library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level of the deal2lkit distribution.
//
//-----------------------------------------------------------

#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/sundials/ida.h>

using namespace dealii;


/**
 * Solve the Harmonic oscillator problem.
 *
 * u'' = -k^2 u
 * u (0) = 0
 * u'(0) = k
 *
 * write in terms of a first order ode:
 *
 * y[0]' -     y[1]  = 0
 * y[1]' + k^2 y[0]  = 0
 *
 * That is
 *
 * F(y', y, t) = y' + A y = 0
 *
 * A = [ 0 , -1; k^2, 0 ]
 *
 * y_0  = 0, k
 * y_0' = 0, 0
 *
 * The exact solution is
 *
 * y[0](t) = sin(k t)
 * y[1](t) = k cos(k t)
 *
 * The Jacobian to assemble is the following:
 *
 * J = alpha I + A
 */
class HarmonicOscillator
{
public:
  HarmonicOscillator(double _kappa = 1.0)
    : y(2)
    , y_dot(2)
    , diff(2)
    , J(2, 2)
    , A(2, 2)
    , Jinv(2, 2)
    , kappa(_kappa)
  {
    diff[0] = 1.0;
    diff[1] = 1.0;

    time_stepper.reinit_vector = [&](Vector<double> &v) { v.reinit(2); };


    using VectorType = Vector<double>;

    time_stepper.residual = [&](const double /*t*/,
                                const VectorType &y,
                                const VectorType &y_dot,
                                VectorType &      res) -> int {
      res = y_dot;
      A.vmult_add(res, y);
      return 0;
    };

    time_stepper.setup_jacobian = [&](const double,
                                      const VectorType &,
                                      const VectorType &,
                                      const double alpha) -> int {
      A(0, 1) = -1.0;
      A(1, 0) = kappa * kappa;

      J = A;

      J(0, 0) = alpha;
      J(1, 1) = alpha;

      Jinv.invert(J);
      return 0;
    };

    time_stepper.solve_jacobian_system = [&](const VectorType &src,
                                             VectorType &      dst) -> int {
      Jinv.vmult(dst, src);
      return 0;
    };

    time_stepper.output_step = [&](const double /*t*/,
                                   const VectorType & /*sol*/,
                                   const VectorType & /*sol_dot*/,
                                   const unsigned int /*step_number*/) -> int {
      // In this test, don't output anything.
      return 0;
    };

    time_stepper.solver_should_restart =
      [](const double, VectorType &, VectorType &) -> bool { return false; };


    time_stepper.differential_components = [&]() -> IndexSet {
      return complete_index_set(2);
    };
  }

  void
  run()
  {
    y[1] = kappa;
    time_stepper.solve_dae(y, y_dot);
  }
  SUNDIALS::IDA<Vector<double>> time_stepper;

private:
  Vector<double>     y;
  Vector<double>     y_dot;
  Vector<double>     diff;
  FullMatrix<double> J;
  FullMatrix<double> A;
  FullMatrix<double> Jinv;
  double             kappa;
};


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  HarmonicOscillator ode(2 * numbers::PI);
  ode.run();
  return 0;
}

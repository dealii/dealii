//-----------------------------------------------------------
//
//    Copyright (C) 2017 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE at
//    the top level of the deal.II distribution.
//
//-----------------------------------------------------------

#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/sundials/ida.h>

#include "../tests.h"


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
 * y_0' = k, 0
 *
 * The exact solution is
 *
 * y[0](t) = sin(k t)
 * y[1](t) = k cos(k t)
 *
 * The Jacobian to assemble is the following:
 *
 * J = alpha I + A
 *
 */
class HarmonicOscillator
{
public:
  HarmonicOscillator(
    double                                                        _kappa,
    const typename SUNDIALS::IDA<Vector<double>>::AdditionalData &data) :
    time_stepper(data),
    y(2),
    y_dot(2),
    J(2, 2),
    A(2, 2),
    Jinv(2, 2),
    kappa(_kappa),
    out("output")
  {
    typedef Vector<double> VectorType;

    time_stepper.reinit_vector = [&](VectorType &v) { v.reinit(2); };


    time_stepper.residual = [&](const double      t,
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

    time_stepper.output_step = [&](const double       t,
                                   const VectorType & sol,
                                   const VectorType & sol_dot,
                                   const unsigned int step_number) -> int {
      out << t << " " << sol[0] << " " << sol[1] << " " << sol_dot[0] << " "
          << sol_dot[1] << std::endl;
      return 0;
    };
  }

  void
  run()
  {
    y[1]     = kappa;
    y_dot[0] = kappa;
    time_stepper.solve_dae(y, y_dot);
  }
  SUNDIALS::IDA<Vector<double>> time_stepper;

private:
  Vector<double>     y;
  Vector<double>     y_dot;
  FullMatrix<double> J;
  FullMatrix<double> A;
  FullMatrix<double> Jinv;
  double             kappa;

  std::ofstream out;
};


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  SUNDIALS::IDA<Vector<double>>::AdditionalData data;
  ParameterHandler                              prm;
  data.add_parameters(prm);

  // std::ofstream ofile(SOURCE_DIR "/harmonic_oscillator_01.prm");
  // prm.print_parameters(ofile, ParameterHandler::ShortText);
  // ofile.close();

  std::ifstream ifile(SOURCE_DIR "/harmonic_oscillator_01.prm");
  prm.parse_input(ifile);


  HarmonicOscillator ode(1.0, data);
  ode.run();
  return 0;
}

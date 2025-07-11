// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/vector.h>

#include <deal.II/sundials/kinsol.h>

#include "../tests.h"


// Like the kinsol_06 test, but with a case where KINSOL, after a few
// tries with the given residual/Jacobian, gives up and terminates
// because it can't find a solution. We need to make sure that we
// throw a catchable exception.
//
// This testcase is a variation of the kinsol_06 test, modified by
// Simon Wiesheier, and posted on the mailing list. Then further
// adapted by Wolfgang Bangerth.


int
main()
{
  initlog();

  using VectorType = Vector<double>;

  // Size of the problem
  const unsigned int N = 1;

  SUNDIALS::KINSOL<VectorType>::AdditionalData data;
  ParameterHandler                             prm;
  data.add_parameters(prm);

  std::ifstream ifile(SOURCE_DIR "/kinsol_06_v3_input.prm");
  prm.parse_input(ifile);

  SUNDIALS::KINSOL<VectorType> kinsol(data);

  kinsol.reinit_vector = [N](VectorType &v) { v.reinit(N); };

  kinsol.residual = [&](const VectorType &u, VectorType &F) {
    // Count how many times this function has been called. Let it call
    // *every* time after the first attempt:
    static int count = 0;
    deallog << "Computing residual for the " << count + 1
            << "th time, at u=" << u[0] << std::endl;
    if ((u[0] < -10) || (u[0] > 20))
      {
        deallog << "Reporting recoverable failure." << std::endl;
        throw RecoverableUserCallbackError();
      }
    count++;


    F.reinit(u);
    F[0] = std::atan(u[0]) - 0.5;
  };

  double J_inverse;

  kinsol.setup_jacobian = [&J_inverse](const VectorType &u,
                                       const VectorType &F) {
    deallog << "Setting up Jacobian system at u=" << u[0] << std::endl;

    const double J = 1. / (1 + u[0] * u[0]);
    J_inverse      = 1. / J;
  };


  kinsol.solve_with_jacobian = [&](const VectorType &rhs,
                                   VectorType       &dst,
                                   double) { dst[0] = J_inverse * rhs[0]; };

  VectorType v(N);
  v[0] = 10;

  try
    {
      auto niter = kinsol.solve(v);
    }
  catch (const ExceptionBase &e)
    {
      deallog << "KINSOL threw an exception with the following message:"
              << std::endl;
      e.print_info(deallog.get_file_stream());
    }
}

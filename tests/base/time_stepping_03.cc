// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test the coefficients of the low-storage Runge-Kutta methods in TimeStepping
#include <deal.II/base/time_stepping.h>

#include "../tests.h"

int
main()
{
  initlog();

  const std::vector<double> bi = {{1153189308089. / 22510343858157.,
                                   1772645290293. / 4653164025191.,
                                   -1672844663538. / 4480602732383.,
                                   2114624349019. / 3568978502595.,
                                   5198255086312. / 14908931495163.}};
  std::vector<double>       ai;
  ai = {{970286171893. / 4311952581923.,
         6584761158862. / 12103376702013.,
         2251764453980. / 15575788980749.,
         26877169314380. / 34165994151039.}};

  deallog << "Check low-storage Runge-Kutta coefficients" << std::endl;
  TimeStepping::LowStorageRungeKutta<Vector<double>> lsrk54(
    TimeStepping::LOW_STORAGE_RK_STAGE5_ORDER4);
  std::vector<double> a, b, c;
  lsrk54.get_coefficients(a, b, c);

  double sum_previous_bi = 0.;
  for (unsigned int i = 0; i < b.size(); ++i)
    {
      if (i != b.size() - 1)
        Assert(std::fabs(a[i] - ai[i]) < 1e-10, ExcInternalError());
      Assert(std::fabs(b[i] - bi[i]) < 1e-10, ExcInternalError());

      if (i > 0)
        {
          const double ci = sum_previous_bi + ai[i - 1];
          deallog << " c: " << c[i] - ci << std::endl;
          sum_previous_bi += bi[i - 1];
        }
      else
        deallog << " c: " << c[0] << std::endl;
    }

  return 0;
}

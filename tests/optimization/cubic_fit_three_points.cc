// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// check minimization of the cubic fit based on f(x1), f(x2) and f'(x1)
// and f'(x2)

#include <deal.II/base/logstream.h>

#include <deal.II/optimization/line_minimization.h>

#include <fstream>
#include <functional>
#include <iostream>

#include "../tests.h"


void
test()
{
  // test 1:
  {
    auto f = [](double x) {
      return std::pow(x, 4) - 20. * std::pow(x, 3) + 0.1 * x;
    };
    auto g = [](double x) {
      return 4. * std::pow(x, 3) - 60. * std::pow(x, 2) + 0.1;
    };

    const double x1 = 17;
    const double x2 = 10;
    const double x3 = 5;
    const double f1 = f(x1);
    const double f2 = f(x2);
    const double f3 = f(x3);
    const double g1 = g(x1);
    const double res =
      *LineMinimization::cubic_fit_three_points(x1, f1, g1, x2, f2, x3, f3);
    deallog << x1 << ' ' << f1 << ' ' << g1 << std::endl
            << x2 << ' ' << f2 << std::endl
            << x3 << ' ' << f3 << std::endl
            << res << std::endl;
  }
}


int
main(int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile, /*do not print job id*/ false);
  deallog.depth_console(0);

  test();
}

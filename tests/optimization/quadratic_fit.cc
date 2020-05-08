//-----------------------------------------------------------
//
//    Copyright (C) 2018 - 2020 by the deal.II authors
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
//---------------------------------------------------------------

// check minimization of the quadratic fit based on f(x1), f(x2) and f'(x1)

#include <deal.II/base/logstream.h>

#include <deal.II/optimization/line_minimization.h>

#include <fstream>
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

    const double x1  = 10;
    const double x2  = 17;
    const double f1  = f(x1);
    const double f2  = f(x2);
    const double g1  = g(x1);
    const double g2  = g(x2);
    const double res = *LineMinimization::quadratic_fit(x1, f1, g1, x2, f2);
    deallog << x1 << " " << x2 << std::endl
            << f1 << " " << f2 << std::endl
            << g1 << " " << g2 << std::endl
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

//-----------------------------------------------------------
//
//    Copyright (C) 2018 by the deal.II authors
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

// check line minimization with strong Wolfe conditions

#include <deal.II/base/logstream.h>

#include <deal.II/optimization/line_minimization.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

/*
 * MWE in Maxima

func(x):=100 * x^4 + (1-x)^2;
gfunc(x):=''(diff(func(x),x));
eta : 0.1;
mu : 0.01;
w1(x) := func(x) - func(0) - x * mu * gfunc(0);
w2(x) := abs(gfunc(x)) - eta * abs(gfunc(0));
plot2d([func(x),gfunc(x),w1(x), w2(x)], [x,0.1,0.2]);
w1(0.159668);
w2(0.159668);
bfloat(solve(gfunc(x)=0)[3]);

 */


void
test()
{
  // test 1:
  {
    const double min_x = 0.161262023139589;
    auto         func  = [](const double x) {
      const double f = 100. * std::pow(x, 4) + std::pow(1. - x, 2);
      const double g = 400. * std::pow(x, 3) - 2. * (1. - x);
      return std::make_pair(f, g);
    };

    const auto fg0 = func(0);
    const auto res =
      LineMinimization::line_search<double>(func,
                                            fg0.first,
                                            fg0.second,
                                            LineMinimization::poly_fit<double>,
                                            0.1,
                                            0.1,
                                            0.01,
                                            100,
                                            20,
                                            true);
    deallog << "Solution: " << res.first << std::endl
            << "Distance: " << std::fabs(res.first - min_x) << std::endl;
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

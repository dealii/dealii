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
// similar to line_minimization.cc but a different function and
// different initial steps

#include <deal.II/base/logstream.h>

#include <deal.II/optimization/line_minimization.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

/*
 * MWE in Maxima

func(x):=-3*x/(x^2+2)-0.03*x;
gfunc(x):=''(diff(func(x),x));
mu:0.025;
eta : 0.1;
a_max : 30;
f_min : func(0) + a_max * mu * gfunc(0);
w1(x) := func(x) - func(0) - x * mu * gfunc(0);
w2(x) := abs(gfunc(x)) - eta * abs(gfunc(0));
w(x) := signum( signum(w1(x)) + signum(w2(x)) + 1 );
pline(x):=func(0)+x*mu*gfunc(0);
plot2d([func(x),pline(x),w(x)], [x,0,20]);
bfloat(solve(gfunc(x)=0)[4]);

 */


void
test()
{
  // test 1:
  {
    const double min_x = 1.474531468108294;
    auto         func  = [](const double x) {
      const double f = (-3. * x) / (x * x + 2.) - 0.03 * x;
      const double g =
        -3. / (x * x + 2.) + (6. * x * x) / std::pow(x * x + 2., 2) - 0.03;
      return std::make_pair(f, g);
    };

    const auto fg0 = func(0);

    {
      deallog << "Case 1:" << std::endl;
      // First, overshoot and get to solution immediately
      const auto res = LineMinimization::line_search<double>(
        func,
        fg0.first,
        fg0.second,
        LineMinimization::poly_fit<double>,
        13,
        0.1,
        0.025,
        30,
        20,
        true);
      deallog << "Solution: " << res.first << std::endl
              << "Distance: " << std::fabs(res.first - min_x) << std::endl;
    }

    {
      deallog << "Case 2:" << std::endl;
      // Now a small step to converge where needed:
      const auto res = LineMinimization::line_search<double>(
        func,
        fg0.first,
        fg0.second,
        LineMinimization::poly_fit<double>,
        0.1,
        0.1,
        0.025,
        30,
        20,
        true);
      deallog << "Solution: " << res.first << std::endl
              << "Distance: " << std::fabs(res.first - min_x) << std::endl;
    }

    {
      deallog << "Case 3:" << std::endl;
      // Now do a big step so that next one in bracketing satisfies both Wolf:
      // at the termination point the derivative is alos negative!
      // Also that point satisfies both Wolfe conditions as well, but
      // we are interested in another segment, which contains local
      // minimizer
      const auto res = LineMinimization::line_search<double>(
        func,
        fg0.first,
        fg0.second,
        LineMinimization::poly_fit<double>,
        1,
        0.1,
        0.025,
        30,
        20,
        true);
      deallog << "Solution: " << res.first << std::endl
              << "Distance: " << std::fabs(res.first - min_x) << std::endl;
    }
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

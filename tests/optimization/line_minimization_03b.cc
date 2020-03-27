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

// same as line_minimization_03 but use three points cubic fit.

#include <deal.II/base/logstream.h>

#include <deal.II/optimization/line_minimization.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

/*
 * MWE in Maxima


Case 1 (compare to Table I) Function 5.1 Figure 3:

b:2;
func(x):=-x/(x^2+b);
gfunc(x):=''(diff(func(x),x));
mu:0.001;
eta : 0.1;
w1(x) := func(x) - func(0) - x * mu * gfunc(0);
w2(x) := abs(gfunc(x)) - eta * abs(gfunc(0));
w(x) := signum( signum(w1(x)) + signum(w2(x)) + 1 );
pline(x):=func(0)+x*mu*gfunc(0);
plot2d([func(x),pline(x),w(x)], [x,0,16]);

Case 2, Function 5.2 Figure 4 (compared to Table II we have more iterations, up
to 8 vs 16):

b:0.004;
func(x):=(x+b)^5-2*(x+b)^4;
gfunc(x):=''(diff(func(x),x));
mu:0.1;
eta : 0.1;
w1(x) := func(x) - func(0) - x * mu * gfunc(0);
w2(x) := abs(gfunc(x)) - eta * abs(gfunc(0));
w(x) := signum( signum(w1(x)) + signum(w2(x)) + 1 );
pline(x):=func(0)+x*mu*gfunc(0);
plot2d([func(x),pline(x),w(x)], [x,0,2]);

Case 3-5 Function 5.4 Figure 6 (Table IV-V-VI):

b1:0.001;
b2:0.001;

b1:0.01;
b2:0.001;

b1:0.001;
b2:0.01;

g(x):=(1+x^2)^(1/2)-x;
func(x):=g(b1)*((1-x)^2+b2^2)^(1/2) + g(b2)*(x^2+b1^2)^(1/2);
gfunc(x):=''(diff(func(x),x));
plot2d([func(x)], [x,0,1]);

 */


void
test()
{
  const std::vector<double> values = {{1e-3, 1e-1, 1e+1, 1e+3}};
  double                    f0, g0, fi, gi;

  {
    deallog << "Table 1:" << std::endl;
    const double b    = 2;
    auto         func = [&](const double x) {
      const double f = -x / (x * x + b);
      const double g = 2. * x * x / std::pow(x * x + 2., 2) - 1. / (x * x + 2.);
      return std::make_pair(f, g);
    };

    const auto fg0 = func(0);

    for (auto a1 : values)
      {
        const auto res = LineMinimization::line_search<double>(
          func,
          fg0.first,
          fg0.second,
          LineMinimization::poly_fit_three_points<double>,
          a1,
          0.1,
          0.001);

        const auto fgi = func(res.first);
        deallog << res.second << " " << res.first << " " << fgi.second
                << std::endl;
      }
  }

  {
    deallog << "Table 2:" << std::endl;
    const double b    = 0.004;
    auto         func = [&](const double x) {
      const double f = std::pow(x + b, 5) - 2. * std::pow(x + b, 4);
      const double g = 5. * std::pow(x + b, 4) - 8. * std::pow(x + b, 3);
      return std::make_pair(f, g);
    };

    const auto fg0 = func(0);

    for (auto a1 : values)
      {
        const auto res = LineMinimization::line_search<double>(
          func,
          fg0.first,
          fg0.second,
          LineMinimization::poly_fit_three_points<double>,
          a1,
          0.100001,
          0.1);

        const auto fgi = func(res.first);
        deallog << res.second << " " << res.first << " " << fgi.second
                << std::endl;
      }
  }

  {
    const std::vector<std::pair<double, double>> params = {
      {{0.001, 0.001}, {0.01, 0.001}, {0.001, 0.01}}};

    unsigned int ind = 4;
    for (auto p : params)
      {
        deallog << "Table " << ind++ << ":" << std::endl;
        const double b1 = p.first;
        const double b2 = p.second;

        const double gb1 = std::sqrt(1. + b1 * b1) - b1;
        const double gb2 = std::sqrt(1. + b2 * b2) - b2;

        auto func = [&](const double x) {
          const double f = gb1 * std::sqrt(std::pow(1. - x, 2) + b2 * b2) +
                           gb2 * std::sqrt(x * x + b1 * b1);
          const double g =
            gb2 * x / sqrt(x * x + b1 * b1) -
            gb1 * (1. - x) / std::sqrt(std::pow(1 - x, 2) + b2 * b2);
          return std::make_pair(f, g);
        };

        const auto fg0 = func(0);

        for (auto a1 : values)
          {
            const auto res = LineMinimization::line_search<double>(
              func,
              fg0.first,
              fg0.second,
              LineMinimization::poly_fit_three_points<double>,
              a1,
              0.00100001,
              0.001);

            const auto fgi = func(res.first);
            deallog << res.second << " " << res.first << " " << fgi.second
                    << std::endl;
          }
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

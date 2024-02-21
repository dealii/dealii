// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test InterpolatedUniformGridData::gradient

#include <deal.II/base/function_lib.h>

#include "../tests.h"

// now interpolate the function x*y*z onto points. note that this function is
// (bi/tri)linear and so we can later know what the correct value is that the
// function should provide
Table<1, double>
fill(const std::array<std::vector<double>, 1> &coordinates)
{
  Table<1, double> data(coordinates[0].size());
  for (unsigned int i = 0; i < coordinates[0].size(); ++i)
    data[i] = coordinates[0][i];
  return data;
}

Table<2, double>
fill(const std::array<std::vector<double>, 2> &coordinates)
{
  Table<2, double> data(coordinates[0].size(), coordinates[1].size());
  for (unsigned int i = 0; i < coordinates[0].size(); ++i)
    for (unsigned int j = 0; j < coordinates[1].size(); ++j)
      data[i][j] = coordinates[0][i] * coordinates[1][j];
  return data;
}

Table<3, double>
fill(const std::array<std::vector<double>, 3> &coordinates)
{
  Table<3, double> data(coordinates[0].size(),
                        coordinates[1].size(),
                        coordinates[2].size());
  for (unsigned int i = 0; i < coordinates[0].size(); ++i)
    for (unsigned int j = 0; j < coordinates[1].size(); ++j)
      for (unsigned int k = 0; k < coordinates[2].size(); ++k)
        data[i][j][k] =
          coordinates[0][i] * coordinates[1][j] * coordinates[2][k];
  return data;
}


template <int dim>
void
check()
{
  // have coordinate arrays that span an interval starting at d+1
  // d+5 nonuniform intervals
  std::array<std::pair<double, double>, dim> intervals;
  std::array<unsigned int, dim>              n_subintervals;
  for (unsigned int d = 0; d < dim; ++d)
    {
      intervals[d]      = std::make_pair(d + 2., 2 * d + 5.);
      n_subintervals[d] = d + 1 + d * d;
    }

  std::array<std::vector<double>, dim> coordinates;
  for (unsigned int d = 0; d < dim; ++d)
    {
      const double x = intervals[d].first;
      const double dx =
        (intervals[d].second - intervals[d].first) / n_subintervals[d];

      for (unsigned int i = 0; i < n_subintervals[d] + 1; ++i)
        coordinates[d].push_back(x + dx * i);
    }

  const Table<dim, double> data = fill(coordinates);

  Functions::InterpolatedUniformGridData<dim> f(intervals,
                                                n_subintervals,
                                                data);

  // now choose a number of randomly chosen points inside the box and
  // verify that the functions returned are correct
  for (unsigned int i = 0; i < 10; ++i)
    {
      Point<dim> p;
      for (unsigned int d = 0; d < dim; ++d)
        p[d] =
          coordinates[d][0] + (random_value<double>()) *
                                (coordinates[d].back() - coordinates[d][0]);

      // The function is x*y*z, so we can compute the gradient pretty easily
      Tensor<1, dim> exact_gradient;
      for (unsigned int d = 0; d < dim; ++d)
        {
          exact_gradient[d] = 1;
          for (unsigned int e = 0; e < dim; ++e)
            if (e != d)
              exact_gradient[d] *= p[e];
        }

      AssertThrow((exact_gradient - f.gradient(p)).norm() < 1e-12,
                  ExcInternalError());
    }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  check<1>();
  check<2>();
  check<3>();
}

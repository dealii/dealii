// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2020 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Test the support points of FE_SimplexP for consistency with equidistant
// suppoort points.

#include <deal.II/fe/fe_simplex_p.h>

#include "../tests.h"


template <int dim>
std::vector<Point<dim>>
equidistant_support_points(const unsigned int degree)
{
  std::vector<Point<dim>> unit_points;
  if constexpr (dim == 1)
    for (unsigned int i = 0; i < degree + 1; ++i)
      unit_points.emplace_back(double(i) / degree);
  else if constexpr (dim == 2)
    for (unsigned int i = 0; i < degree + 1; ++i)
      for (unsigned int j = 0; j < degree + 1 - i; ++j)
        unit_points.emplace_back(double(i) / degree, double(j) / degree);
  else if constexpr (dim == 3)
    for (unsigned int i = 0; i < degree + 1; ++i)
      for (unsigned int j = 0; j < degree + 1 - i; ++j)
        for (unsigned int k = 0; k < degree + 1 - i - j; ++k)
          unit_points.emplace_back(double(i) / degree,
                                   double(j) / degree,
                                   double(k) / degree);
  return unit_points;
}

template <int dim>
void
test(const unsigned int degree)
{
  deallog << "Support points of degree " << degree << std::endl;
  const auto fe             = FE_SimplexP<dim>(degree);
  const auto support_points = fe.get_unit_support_points();

  const unsigned int n = dim == 2 ?
                           (degree + 1) * (degree + 2) / 2 :
                           (degree + 1) * (degree + 2) * (degree + 3) / 6;
  if (support_points.size() - n == 0)
    deallog << "Correct size" << std::endl;

  const auto equidistant_points = equidistant_support_points<dim>(degree);

  bool consistent = true;
  for (const auto &p : support_points)
    {
      bool found_node = false;
      for (const auto &eqi : equidistant_points)
        if (p.distance(eqi) < 1e-12)
          found_node = true;

      if (found_node == false)
        consistent = false;
    }

  if (consistent)
    deallog << "All nodes correct" << std::endl;
  else
    deallog << "Error in the nodes" << std::endl;

  deallog << std::endl;
}


int
main()
{
  initlog();

  deallog.push("2D");
  for (unsigned int i = 1; i < 4; ++i)
    test<2>(i);
  deallog.pop();

  deallog.push("3D");
  for (unsigned int i = 1; i < 4; ++i)
    test<3>(i);
  deallog.pop();
}

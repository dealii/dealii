// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Extract vertex rtree from the cache, and try to use it

#include <deal.II/base/patterns.h>

#include <deal.II/boost_adaptors/point.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria.h>

#include <algorithm>

#include "../tests.h"

using Patterns::Tools::to_string;

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

template <int dim, int spacedim>
void
test(const unsigned int ref = 2, const unsigned int n_points = 10)
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(ref);
  GridTools::Cache<dim, spacedim> cache(tria);

  const auto &v_tree = cache.get_used_vertices_rtree();

  std::vector<Point<spacedim>> points(n_points);
  std::generate(points.begin(), points.end(), []() {
    return random_point<spacedim>();
  });

  deallog << "Testing dim = " << dim << ", spacedim = " << spacedim
          << std::endl;

  for (const auto &p : points)
    {
      std::vector<std::pair<Point<spacedim>, unsigned int>> res;
      v_tree.query(bgi::nearest(p, 1), std::back_inserter(res));
      deallog << "Nearest vertex to " << p << ": v[" << res[0].second
              << "] = " << res[0].first << std::endl;
    }
}

int
main(int argc, char **argv)
{
  initlog();

  test<1, 1>();
  test<1, 2>();
  test<1, 3>();
  test<2, 2>();
  test<2, 3>();
  test<3, 3>();
}

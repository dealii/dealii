// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Read a Surface_mesh from an .off file, then create a coarse mesh filled with
// tets

#include <deal.II/base/config.h>

#include <deal.II/cgal/utilities.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <CGAL/IO/File_medit.h>
#include <CGAL/IO/io.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include "../tests.h"

// Compute volumes by splitting polyhedra into simplices.

using K                 = CGAL::Exact_predicates_inexact_constructions_kernel;
using CGALPoint         = CGAL::Point_3<K>;
using CGALTriangulation = CGAL::Triangulation_3<K>;
using namespace CGALWrappers;



void
test()
{
  const std::vector<std::string> fnames{SOURCE_DIR "/input_grids/cube.off",
                                        SOURCE_DIR
                                        "/input_grids/tetrahedron.off",
                                        SOURCE_DIR "/input_grids/hedra.off",
                                        SOURCE_DIR
                                        "/input_grids/octahedron.off"};
  CGAL::Surface_mesh<CGALPoint>  sm;
  CGALTriangulation              tria;
  constexpr int                  degree = 3;
  for (const auto &name : fnames)
    {
      std::ifstream input(name);
      input >> sm;
      tria.insert(sm.points().begin(), sm.points().end());
      auto b = compute_quadrature(tria, degree);
      deallog << "Volume of poly with Quadrature: " << std::setprecision(12)
              << std::accumulate(b.get_weights().begin(),
                                 b.get_weights().end(),
                                 0.)
              << "\t Expected:" << std::setprecision(12)
              << CGAL::to_double(CGAL::Polygon_mesh_processing::volume(sm))
              << std::endl;
      sm.clear(); // reset surface
      tria.clear();
    }
}

int
main()
{
  initlog();
  test();
}

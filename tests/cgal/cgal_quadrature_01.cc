// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// Read a Surface_mesh from an .off file, then create a coarse mesh filled with
// tets

#include <deal.II/base/config.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <CGAL/IO/File_medit.h>
#include <CGAL/IO/io.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <deal.II/cgal/utilities.h>

#include "../tests.h"

// Compute volumes by splitting polyhedra into simplices.

using K         = CGAL::Exact_predicates_inexact_constructions_kernel;
using CGALPoint = CGAL::Point_3<K>;
using namespace CGALWrappers;
using Mesh_domain =
  CGAL::Polyhedral_mesh_domain_with_features_3<K,
                                               CGAL::Surface_mesh<CGALPoint>>;
using Tr = typename CGAL::
  Mesh_triangulation_3<Mesh_domain, CGAL::Default, ConcurrencyTag>::type;

using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;
using C3t3          = CGAL::Mesh_complex_3_in_triangulation_3<Tr,
                                                     Mesh_domain::Corner_index,
                                                     Mesh_domain::Curve_index>;


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
  C3t3                           tria;
  constexpr int                  degree = 3;
  for (const auto &name : fnames)
    {
      std::ifstream input(name);
      input >> sm;
      cgal_surface_mesh_to_cgal_triangulation(sm, tria);
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

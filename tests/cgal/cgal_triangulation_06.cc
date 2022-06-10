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

// Read a surface mesh, make a coarse CGAL triangulation out of it, and
// translate the result into a deal.II Triangulation.

#include <deal.II/base/config.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/io.h>
#include <deal.II/cgal/surface_mesh.h>
#include <deal.II/cgal/triangulation.h>
#include <deal.II/cgal/utilities.h>

#include "../tests.h"

using namespace CGALWrappers;

using K         = CGAL::Exact_predicates_inexact_constructions_kernel;
using CGALPoint = K::Point_3;

using Mesh_domain =
  CGAL::Polyhedral_mesh_domain_with_features_3<K,
                                               CGAL::Surface_mesh<CGALPoint>>;
using Tr = typename CGAL::
  Mesh_triangulation_3<Mesh_domain, CGAL::Default, ConcurrencyTag>::type;
using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;
using C3t3          = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;

int
main()
{
  initlog();
  // Build a deal.II triangulation
  CGAL::Surface_mesh<CGALPoint> cgal_surface_mesh;
  {
    std::ifstream input_mesh(SOURCE_DIR "/input_grids/tetrahedron.off");
    input_mesh >> cgal_surface_mesh;
  }

  C3t3 cgal_triangulation;
  cgal_surface_mesh_to_cgal_triangulation(cgal_surface_mesh,
                                          cgal_triangulation);

  Triangulation<3> tria;
  cgal_triangulation_to_dealii_triangulation(cgal_triangulation, tria);

  deallog << "OK" << std::endl;
}

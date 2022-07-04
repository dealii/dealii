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

// Remesh the union of two deal.II hyper_spheres in order to avoid bad
// triangles.

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/cgal/surface_mesh.h>
#include <deal.II/cgal/triangulation.h>
#include <string.h>

#include "../tests.h"

using namespace CGALWrappers;
using K         = CGAL::Exact_predicates_inexact_constructions_kernel;
using CGALPoint = CGAL::Point_3<K>;
using CGALMesh  = CGAL::Surface_mesh<CGALPoint>;

template <int dim, int spacedim>
void
test()
{
  deallog << "dim= " << dim << ",\t spacedim= " << spacedim << std::endl;

  Triangulation<spacedim> tria0, tria1;
  Triangulation<2, 3>     tria_out;
  GridOut                 go;
  CGALMesh                surface_mesh0, surface_mesh1, out_mesh;

  GridGenerator::hyper_ball(tria0, {0., 0., 0.}, 0.4);
  GridGenerator::hyper_ball(tria1, {0.3, 0.3, 0.3}, 0.4);
  tria0.refine_global(3);
  tria1.refine_global(3);

  // Move to CGAL surfaces
  dealii_tria_to_cgal_surface_mesh(tria0, surface_mesh0);
  dealii_tria_to_cgal_surface_mesh(tria1, surface_mesh1);

  // close the surfaces
  CGAL::Polygon_mesh_processing::stitch_borders(surface_mesh0);
  CGAL::Polygon_mesh_processing::stitch_borders(surface_mesh1);

  CGAL::Polygon_mesh_processing::triangulate_faces(surface_mesh0);
  CGAL::Polygon_mesh_processing::triangulate_faces(surface_mesh1);

  compute_boolean_operation(surface_mesh0,
                            surface_mesh1,
                            BooleanOperation::compute_union,
                            out_mesh);
  AdditionalData<3> data;
  data.edge_size      = .02;
  data.facet_angle    = 25;
  data.facet_size     = 0.05;
  data.facet_distance = 0.05;

  remesh_surface<CGALPoint>(out_mesh, data);

  //  Now back to deal.II
  cgal_surface_mesh_to_dealii_triangulation(out_mesh, tria_out);
  std::ofstream out_name_spheres("remeshed_union_spheres.vtk");
  go.write_vtk(tria_out, out_name_spheres);
  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();
  test<3, 3>();
}

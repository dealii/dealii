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

// Performs the following: deal.II Triangulations -> CGAL::Surface_mesh(es) ->
// Perform boolean operation -> deal.II tria.

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

#include <deal.II/cgal/surface_mesh.h>
#include <deal.II/cgal/triangulation.h>
#include <deal.II/cgal/utilities.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

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

  GridGenerator::hyper_ball(tria0, {0., 0., 0.}, 0.5);
  GridGenerator::hyper_ball(tria1, {0.2, 0.2, 0.2}, 0.5);
  tria0.refine_global(3);
  tria1.refine_global(3);

  // Move to CGAL surfaces
  dealii_tria_to_cgal_surface_mesh(tria0, surface_mesh0);
  dealii_tria_to_cgal_surface_mesh(tria1, surface_mesh1);

  // Ensure the meshes are closed
  Assert(CGAL::is_closed(surface_mesh0),
         ExcMessage("The CGAL mesh 0 is not closed"));
  Assert(CGAL::is_closed(surface_mesh1),
         ExcMessage("The CGAL mesh 1 is not closed"));

  // Surfaces automatically closed but still need to be triangulated
  // before using compute_boolean_operation
  CGAL::Polygon_mesh_processing::triangulate_faces(surface_mesh0);
  CGAL::Polygon_mesh_processing::triangulate_faces(surface_mesh1);

  compute_boolean_operation(surface_mesh0,
                            surface_mesh1,
                            BooleanOperation::compute_intersection,
                            out_mesh);
  // Now back to deal.II
  cgal_surface_mesh_to_dealii_triangulation(out_mesh, tria_out);
  std::ofstream out_name_spheres("boolean_intersection_hyper_spheres.vtk");
  go.write_vtk(tria_out, out_name_spheres);
  deallog << "OK" << std::endl;
  remove("boolean_intersection_hyper_spheres.vtk");

  // Clear everything
  tria0.clear();
  tria1.clear();
  tria_out.clear();
  surface_mesh0.clear();
  surface_mesh1.clear();
  out_mesh.clear();

  GridGenerator::hyper_cube(tria0);
  GridGenerator::hyper_cube(tria1, 0.5, 1.5);
  tria0.refine_global(3);
  tria1.refine_global(3);
  GridTools::rotate(Tensor<1, 3>{{0., 0., 1.}}, numbers::PI_4, tria1);

  // Move to CGAL surfaces
  dealii_tria_to_cgal_surface_mesh(tria0, surface_mesh0);
  dealii_tria_to_cgal_surface_mesh(tria1, surface_mesh1);

  // Ensure the meshes are closed
  Assert(CGAL::is_closed(surface_mesh0),
         ExcMessage("The CGAL mesh 0 is not closed"));
  Assert(CGAL::is_closed(surface_mesh1),
         ExcMessage("The CGAL mesh 1 is not closed"));

  // Surfaces automatically closed but still need to be triangulated
  // before using compute_boolean_operation
  CGAL::Polygon_mesh_processing::triangulate_faces(surface_mesh0);
  CGAL::Polygon_mesh_processing::triangulate_faces(surface_mesh1);

  compute_boolean_operation(surface_mesh0,
                            surface_mesh1,
                            BooleanOperation::compute_intersection,
                            out_mesh);
  // Now back to deal.II
  cgal_surface_mesh_to_dealii_triangulation(out_mesh, tria_out);
  std::ofstream out_name_cubes("boolean_intersection_cubes.vtk");
  go.write_vtk(tria_out, out_name_cubes);
  deallog << "OK" << std::endl;
  remove("boolean_intersection_cubes.vtk");
}

int
main()
{
  initlog();
  test<3, 3>();
}

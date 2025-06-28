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

// Remesh the union of two deal.II hyper_spheres in order to avoid bad
// triangles.

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

#include <deal.II/cgal/surface_mesh.h>
#include <deal.II/cgal/triangulation.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <CGAL/Polygon_mesh_processing/orientation.h>
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
  // This test might trigger spurious floating point exception despite
  // functioning properly. Simply disable floating point exceptions again
  // (after they had been enabled int tests.h)
#if defined(DEBUG) && defined(DEAL_II_HAVE_FP_EXCEPTIONS)
  {
    const int current_fe_except = fegetexcept();
    fedisableexcept(current_fe_except);
  }
#endif

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

  Assert(surface_mesh0.is_valid() && surface_mesh1.is_valid(),
         ExcMessage("The CGAL surface mesh is not valid."));
  if (dim == 3)
    {
      Assert(CGAL::is_closed(surface_mesh0) && CGAL::is_closed(surface_mesh1),
             dealii::ExcMessage("The CGAL mesh is not closed"));
      Assert(
        CGAL::Polygon_mesh_processing::is_outward_oriented(surface_mesh0) &&
          CGAL::Polygon_mesh_processing::is_outward_oriented(surface_mesh1),
        dealii::ExcMessage(
          "The normal vectors of the CGAL mesh are not oriented outwards"));
    }

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

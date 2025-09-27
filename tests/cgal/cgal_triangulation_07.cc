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

// Create a Tria<3> out of a Tria<2,3>.

#include <deal.II/base/config.h>

#include <deal.II/cgal/triangulation.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

using namespace CGALWrappers;
int
main()
{
  initlog();
  Triangulation<2, 3> surface_tria;
  GridGenerator::hyper_sphere(surface_tria, {0., 1., 0.}, 1.);
  surface_tria.refine_global(3);

  Triangulation<3>  out_tria;
  AdditionalData<3> data;
  data.facet_size             = 0.1;
  data.facet_distance         = 0.;
  data.cell_radius_edge_ratio = 2.;
  data.cell_size              = 0.1;
  GridGenerator::surface_mesh_to_volumetric_mesh(surface_tria, out_tria, data);

  GridOut       go;
  std::ofstream out_tria_name("tria_surface_to_volumetric.vtk");
  go.write_vtk(out_tria, out_tria_name);

  remove("tria_surface_to_volumetric.vtk");
  deallog << "OK" << std::endl;
}

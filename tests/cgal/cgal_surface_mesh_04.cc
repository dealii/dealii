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

// Convert a deal.II triangulation to a CGAL surface mesh. In the 2D case,
// the whole triangulation is a 2D surface mesh. In 3D, the surface mesh
// describes the boundary of the deal.II Triangulation.

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

#include <deal.II/cgal/surface_mesh.h>
#include <deal.II/cgal/triangulation.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
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
  deallog << "dim= " << dim << ",\t spacedim= " << spacedim << std::endl;

  Triangulation<spacedim> tria_in;
  Triangulation<2, 3>     tria_out;
  GridOut                 go;
  CGALMesh                surface_mesh;

  std::vector<std::pair<std::string, std::string>> names_and_args;
  if constexpr (spacedim == 2)
    {
      names_and_args = {{"hyper_cube", "0.0 : 1.0 : false"},
                        {"hyper_ball", "0.,0. : 1. : false"},
                        {"hyper_L", "0.0 : 1.0 : false"},
                        {"channel_with_cylinder", "0.03 : 2 : 2.0 : false"}};
    }
  else
    {
      names_and_args = {{"hyper_cube", "0.0 : 1.0 : false"},
                        {"hyper_ball", "0.,0.,0. : 1. : false"},
                        {"hyper_L", "0.0 : 1.0 : false"},
                        {"channel_with_cylinder", "0.03 : 2 : 2.0 : false"}};
    }


  for (const auto &info_pair : names_and_args)
    {
      auto name = info_pair.first;
      auto args = info_pair.second;
      deallog << "dim = " << dim << ", spacedim = " << spacedim
              << " name: " << name << std::endl;
      GridGenerator::generate_from_name_and_arguments(tria_in, name, args);
      tria_in.refine_global(2);
      dealii_tria_to_cgal_surface_mesh(tria_in, surface_mesh);
      Assert(surface_mesh.is_valid(),
             ExcMessage("The CGAL surface mesh is not valid."));

      if (dim == 3)
        {
          Assert(CGAL::is_closed(surface_mesh),
                 ExcMessage("The CGAL mesh is not closed"));
          Assert(
            CGAL::Polygon_mesh_processing::is_outward_oriented(surface_mesh),
            ExcMessage(
              "The normal vectors of the CGAL mesh are not oriented outwards"));
        }

      // Now back to the original dealii tria.
      cgal_surface_mesh_to_dealii_triangulation(surface_mesh, tria_out);
      std::ofstream out_name(name + std::to_string(spacedim) + ".vtk");
      go.write_vtk(tria_out, out_name);
      // If we got here, everything was ok, including writing the grid.
      deallog << "OK" << std::endl;
      tria_in.clear();
      tria_out.clear();
      surface_mesh.clear();
    }
}

int
main()
{
  initlog();
  test<2, 2>();
  test<3, 3>();
}

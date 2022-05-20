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

// Convert a deal.II triangulation to a CGAL surface mesh. In the 2D case,
// thw whole triangulation is a 2D surface mesh. In 3D, the surface mesh
// describes the boundary of the deal.II Triangulation.

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
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

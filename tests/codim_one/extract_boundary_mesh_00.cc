// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2018 by the deal.II authors
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


/*
 Code for testing the function
 GridGenerator::extract_boundary_mesh (...).
 We test that the order of cells and the orientation
 of the vertices is consistent between the two meshes.

 This test checks the whole thing for a 2d and 3d hypercube.
*/


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

using namespace std;


template <int dim, int spacedim>
void
save_mesh(const Triangulation<dim, spacedim> &tria)
{
  GridOut grid_out;
  grid_out.write_ucd(tria, deallog.get_file_stream());
}


int
main()
{
  initlog();

  {
    // Extract the whole boundary of a hyper-cube
    const int dim = 2;

    deallog << "Testing hyper_cube in dim: " << dim << "..." << endl;
    map<Triangulation<dim - 1, dim>::cell_iterator,
        Triangulation<dim, dim>::face_iterator>
      surface_to_volume_mapping;

    Triangulation<dim> volume_mesh;
    GridGenerator::hyper_cube(volume_mesh);
    volume_mesh.begin_active()->face(0)->set_boundary_id(1);
    volume_mesh.refine_global(1);

    save_mesh(volume_mesh);

    Triangulation<dim - 1, dim> boundary_mesh;

    surface_to_volume_mapping =
      GridGenerator::extract_boundary_mesh(volume_mesh, boundary_mesh);

    save_mesh(boundary_mesh);
  }

  {
    // Extract the whole boundary of a hyper-cube
    const int dim = 3;

    deallog << "Testing hyper_cube in dim: " << dim << "..." << endl;
    map<Triangulation<dim - 1, dim>::cell_iterator,
        Triangulation<dim, dim>::face_iterator>
      surface_to_volume_mapping;

    Triangulation<dim> volume_mesh;
    GridGenerator::hyper_cube(volume_mesh);
    volume_mesh.begin_active()->face(0)->set_boundary_id(1);
    volume_mesh.refine_global(1);

    save_mesh(volume_mesh);

    Triangulation<dim - 1, dim> boundary_mesh;

    surface_to_volume_mapping =
      GridGenerator::extract_boundary_mesh(volume_mesh, boundary_mesh);

    save_mesh(boundary_mesh);
  }


  {
    // Extract a piece of the boundary of a hyper-cube

    const int dim = 3;
    deallog << "Testing hyper_cube in dim: " << dim << "..." << endl;

    map<Triangulation<dim - 1, dim>::cell_iterator,
        Triangulation<dim, dim>::face_iterator>
      surface_to_volume_mapping;

    Triangulation<dim> volume_mesh;
    GridGenerator::hyper_cube(volume_mesh);
    volume_mesh.begin_active()->face(0)->set_boundary_id(1);
    volume_mesh.refine_global(1);

    save_mesh(volume_mesh);

    Triangulation<dim - 1, dim> boundary_mesh;
    set<types::boundary_id>     boundary_ids;
    boundary_ids.insert(0);

    surface_to_volume_mapping =
      GridGenerator::extract_boundary_mesh(volume_mesh,
                                           boundary_mesh,
                                           boundary_ids);

    save_mesh(boundary_mesh);
  }



  return 0;
}

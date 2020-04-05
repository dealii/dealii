// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2019 by the deal.II authors
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

 Like _02, but set volume boundary indicator to 1 and make sure that
 extract_boundary_mesh also copies this.
*/

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

using namespace std;


template <int s_dim, int spacedim>
void
test_vertices_orientation(
  const Triangulation<s_dim, spacedim> &boundary_mesh,
  map<typename Triangulation<s_dim, spacedim>::cell_iterator,
      typename Triangulation<s_dim + 1, spacedim>::face_iterator>
    &surface_to_volume_mapping)
{
  typename Triangulation<s_dim, spacedim>::active_cell_iterator
    cell = boundary_mesh.begin_active(),
    endc = boundary_mesh.end();
  typename Triangulation<s_dim + 1, spacedim>::face_iterator face;

  for (; cell != endc; ++cell)
    {
      face = surface_to_volume_mapping[cell];
      Assert(face->at_boundary(), ExcInternalError());

      deallog << "Surface cell: " << cell << " with vertices:" << std::endl;
      for (const unsigned int k : GeometryInfo<s_dim>::vertex_indices())
        {
          deallog << "  " << cell->vertex(k) << std::endl;
          Assert(std::fabs(cell->vertex(k).distance(Point<spacedim>()) - 1) <
                   1e-12,
                 ExcInternalError());
        }

      deallog << "Volume face: " << face << " with vertices:" << std::endl;
      for (const unsigned int k : GeometryInfo<s_dim>::vertex_indices())
        {
          deallog << "  " << face->vertex(k) << std::endl;
          Assert(std::fabs(face->vertex(k).distance(Point<spacedim>()) - 1) <
                   1e-12,
                 ExcInternalError());
        }


      Point<spacedim> diff(face->center());
      diff -= cell->center();
      AssertThrow(diff.square() < 1.e-15 * face->center().square(),
                  ExcInternalError());
    }
}

template <int dim, int spacedim>
void
save_mesh(const Triangulation<dim, spacedim> &tria)
{
  GridOut grid_out;
  grid_out.write_gnuplot(tria, deallog.get_file_stream());
}


int
main()
{
  initlog();


  {
    // Extract the boundary of a hyper-sphere

    const int dim = 3;
    deallog << "Testing hyper_ball in dim: " << dim << "..." << endl;

    map<Triangulation<dim - 1, dim>::cell_iterator,
        Triangulation<dim, dim>::face_iterator>
                                 surface_to_volume_mapping;
    const SphericalManifold<dim> boundary_description;
    Triangulation<dim>           volume_mesh;
    GridGenerator::hyper_ball(volume_mesh);
    for (Triangulation<dim>::active_cell_iterator cell =
           volume_mesh.begin_active();
         cell != volume_mesh.end();
         ++cell)
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        if (cell->at_boundary(f))
          cell->face(f)->set_all_boundary_ids(1);
    GridTools::copy_boundary_to_manifold_id(volume_mesh);
    volume_mesh.set_manifold(1, boundary_description);
    volume_mesh.refine_global(1);

    const SphericalManifold<dim - 1, dim> surface_description;
    Triangulation<dim - 1, dim>           boundary_mesh;
    boundary_mesh.set_manifold(1, surface_description);

    surface_to_volume_mapping =
      GridGenerator::extract_boundary_mesh(volume_mesh, boundary_mesh);
    deallog << volume_mesh.n_active_cells() << std::endl;
    deallog << boundary_mesh.n_active_cells() << std::endl;
    save_mesh(boundary_mesh);


    test_vertices_orientation(boundary_mesh, surface_to_volume_mapping);
    save_mesh(boundary_mesh);
  }


  return 0;
}

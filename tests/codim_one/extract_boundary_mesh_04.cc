// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


/*
 Code for testing the function
 GridGenerator::extract_boundary_mesh (...).
 We test that the order of cells and the orientation
 of the vertices is consistent between the two meshes.

 Like _03, but set volume boundary indicator to 1 only on the upper
 half of the sphere .
*/

#include "../tests.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

using namespace std;


template <int s_dim, int spacedim>
void test_vertices_orientation(const Triangulation<s_dim,spacedim> &boundary_mesh,
                               map< typename Triangulation<s_dim,spacedim>::cell_iterator,
                               typename Triangulation<s_dim+1,spacedim>::face_iterator >
                               &surface_to_volume_mapping)
{
  typename Triangulation<s_dim,spacedim>::active_cell_iterator
  cell = boundary_mesh.begin_active(),
  endc = boundary_mesh.end();
  typename Triangulation<s_dim+1,spacedim>::face_iterator face;

  for (; cell!=endc; ++cell)
    {

      face = surface_to_volume_mapping [cell];
      Assert (face->at_boundary(), ExcInternalError());

      deallog << "Surface cell: " << cell << " with vertices:" << std::endl;
      for (unsigned int k=0; k<GeometryInfo<s_dim>::vertices_per_cell; ++k)
        {
          deallog << "  " << cell->vertex(k) << std::endl;
          Assert (std::fabs(cell->vertex(k).distance (Point<spacedim>()) - 1) < 1e-12,
                  ExcInternalError());
        }

      deallog << "Volume face: " << face << " with vertices:" << std::endl;
      for (unsigned int k=0; k<GeometryInfo<s_dim>::vertices_per_cell; ++k)
        {
          deallog << "  " << face->vertex(k) << std::endl;
          Assert (std::fabs(face->vertex(k).distance (Point<spacedim>()) - 1) < 1e-12,
                  ExcInternalError());
        }


      for (unsigned int k=0; k<GeometryInfo<s_dim>::vertices_per_cell; ++k)
        {
          Point<spacedim> diff(face->vertex(k));
          diff -= cell->vertex(k);
          Assert (diff.square() == 0, ExcInternalError());
        }
    }
}

template <int dim, int spacedim>
void save_mesh(const Triangulation<dim,spacedim> &tria)
{
  GridOut grid_out;
  grid_out.write_gnuplot (tria, deallog.get_file_stream());
}


int main ()
{

  ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);


  {
    // Extract the boundary of a hyper-sphere

    const int dim = 3;
    deallog << "Testing hyper_cube in dim: " << dim << "..."<< endl;

    map< Triangulation<dim-1,dim>::cell_iterator,
         Triangulation<dim,dim>::face_iterator>
         surface_to_volume_mapping;
    const HyperBallBoundary<dim> boundary_description;
    Triangulation<dim> volume_mesh;
    GridGenerator::hyper_ball(volume_mesh);
    for (Triangulation<dim>::active_cell_iterator
         cell = volume_mesh.begin_active();
         cell != volume_mesh.end(); ++cell)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->at_boundary(f) && (cell->center()[dim-1]>0))
          cell->face(f)->set_all_boundary_indicators (1);

    volume_mesh.set_boundary (1, boundary_description);
    volume_mesh.set_boundary (0, boundary_description);
    volume_mesh.refine_global (1);

    const HyperBallBoundary<dim-1,dim> surface_description;
    Triangulation<dim-1,dim> boundary_mesh;
    boundary_mesh.set_boundary (1, surface_description);
    boundary_mesh.set_boundary (0, surface_description);

    set<types::boundary_id> boundary_ids;
    boundary_ids.insert(1);

    surface_to_volume_mapping
      = GridGenerator::extract_boundary_mesh (volume_mesh, boundary_mesh,
                                          boundary_ids);

    deallog << volume_mesh.n_active_cells () << std::endl;
    deallog << boundary_mesh.n_active_cells () << std::endl;
    deallog << "Surface mesh:" << std::endl;
    save_mesh(boundary_mesh);

    test_vertices_orientation(boundary_mesh, surface_to_volume_mapping);
  }


  return 0;
}

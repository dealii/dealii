// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2015 by the deal.II authors
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

 Compared to the _00 test, we here test things with the surface of a
 3d sphere. There is no surface object attached, so the surface is
 topologically a sphere but not geometrically. The _02 test checks the
 same with a boundary object attached.
*/

#include "../tests.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

using namespace std;


template <int s_dim, int spacedim>
bool test_vertices_orientation(const Triangulation<s_dim,spacedim> &boundary_mesh,
                               map< typename Triangulation<s_dim,spacedim>::cell_iterator,
                               typename Triangulation<s_dim+1,spacedim>::face_iterator >
                               &surface_to_volume_mapping,
                               const int verbosity = 1)
{
  typename Triangulation<s_dim,spacedim>::active_cell_iterator
  cell = boundary_mesh.begin_active(),
  endc = boundary_mesh.end();
  typename Triangulation<s_dim+1,spacedim>::face_iterator face;

  bool success = true;

  if (verbosity>1)
    {
      deallog << "The last column should be 0" <<endl;
      deallog << "Vol faces" << "\t\t" << "Surf cell" <<
              "\t\t" << "Distance" <<endl<<endl;
    }

  for (; cell!=endc; ++cell)
    {

      face = surface_to_volume_mapping [cell];

      for (unsigned int k=0; k<GeometryInfo<s_dim>::vertices_per_cell; ++k)
        {
          Point<spacedim> diff(face->vertex(k));
          diff -= cell->vertex(k);
          if (verbosity>1)
            {
              deallog << face->vertex(k) << "\t\t";
              deallog << cell->vertex(k) << "\t\t\t" << diff.square() << endl;
            }
          if (diff.square()>0)
            {
              success = false;
              break;
            }
        }
      if (verbosity>1) deallog << endl;
    }
  return success;
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


  {
    // Extract the boundary of a hyper-sphere

    const int dim = 3;
    deallog << "Testing hyper_ball in dim: " << dim << "..."<< endl;

    map< Triangulation<dim-1,dim>::cell_iterator,
         Triangulation<dim,dim>::face_iterator>
         surface_to_volume_mapping;
    Triangulation<dim> volume_mesh;
    GridGenerator::hyper_ball(volume_mesh);
    volume_mesh.reset_manifold (0);

    volume_mesh.refine_global (1);

    Triangulation<dim-1,dim> boundary_mesh;

    surface_to_volume_mapping
      = GridGenerator::extract_boundary_mesh (volume_mesh, boundary_mesh);

    if (!test_vertices_orientation(boundary_mesh, surface_to_volume_mapping,2))
      Assert (false, ExcInternalError());
    save_mesh(boundary_mesh);
  }


  return 0;
}

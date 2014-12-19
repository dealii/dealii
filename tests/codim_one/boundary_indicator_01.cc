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



// for surfaces, we need some sort of mapping also for interior cells
// and faces, but at the time of writing this, we can only set
// boundary_indicators of faces and edges that are truly at the
// boundary of the domain. to work around this, we use the material_id
// field that one can set on cells, and copy it also to the adjacent
// faces. initially, however, we also copied the material_id to the
// boundary_indicator of adjacent faces that truly were at the
// boundary of the domain, and which might have had something
// purposefully set already
//
// this test verifies that this is now fixed.

#include "../tests.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

using namespace std;



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

  // Extract the boundary of 3/4 of a sphere
  {
    const int dim = 3;
    deallog << "Testing hyper_cube in dim: " << dim << "..."<< endl;

    const HyperBallBoundary<dim> boundary_description;
    Triangulation<dim> volume_mesh;
    GridGenerator::hyper_ball(volume_mesh);
    volume_mesh.set_boundary (0, boundary_description);

    // exclude one of the 6 faces
    // from the surface mesh
    // extraction
    for (Triangulation<dim>::active_cell_iterator
         cell = volume_mesh.begin_active();
         cell != volume_mesh.end(); ++cell)
      {
        bool done = false;
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->at_boundary(f))
            {
              cell->face(f)->set_boundary_indicator(1);
              done = true;
              break;
            }
        if (done)
          break;
      }

    const HyperBallBoundary<dim-1,dim> surface_description;
    Triangulation<dim-1,dim> boundary_mesh;
    boundary_mesh.set_boundary (0, surface_description);

    // now extract a mesh of the 5
    // surface faces
    std::set<types::boundary_id> boundary_indicators;
    boundary_indicators.insert (0);
    GridGenerator::extract_boundary_mesh (volume_mesh, boundary_mesh,
                                      boundary_indicators);
    deallog << volume_mesh.n_active_cells () << std::endl;
    deallog << boundary_mesh.n_active_cells () << std::endl;

    // at this point, all cells and
    // edges of the surface mesh
    // should have boundary indicator
    // 0. set those at the boundary
    // of the mesh to 1 to force
    // straight line refinement, then
    // refine
    for (Triangulation<dim-1,dim>::active_cell_iterator
         cell = boundary_mesh.begin_active();
         cell != boundary_mesh.end(); ++cell)
      for (unsigned int f=0; f<GeometryInfo<dim-1>::faces_per_cell; ++f)
        if (cell->at_boundary(f))
          cell->face(f)->set_boundary_indicator(1);

    boundary_mesh.refine_global (2);

    save_mesh(boundary_mesh);
  }


  return 0;
}

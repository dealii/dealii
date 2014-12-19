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



// an extract of the _01 test that shows the essence of what's going
// wrong. compared to the _03 testcase, upon refinement of the middle
// cell (0.0), cell 0.1 forgets who some of its neighbors are. this is
// clearly not good
//
// the underlying reason was that when we update neighborship
// information in the triangulation class, we have to take into
// account the direction_flag of cells

#include "../tests.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  const unsigned int spacedim = 3;
  const unsigned int dim = spacedim-1;

  Triangulation<dim,spacedim> boundary_mesh;

  /*****************************************************************/
  // Create Surface Mesh:  Boundary of hypercube without one face
  {
    Triangulation<spacedim> volume_mesh;
    GridGenerator::hyper_cube(volume_mesh);
    Triangulation<spacedim>::active_cell_iterator
    cell = volume_mesh.begin_active();

    cell->face(0)->set_all_boundary_indicators (1);
    std::set<types::boundary_id> boundary_ids;
    boundary_ids.insert(0);
    GridGenerator::extract_boundary_mesh (volume_mesh, boundary_mesh, boundary_ids);
  }
  boundary_mesh.begin_active()->set_refine_flag ();
  boundary_mesh.execute_coarsening_and_refinement ();

  Triangulation<dim,spacedim>::active_cell_iterator
  cell = boundary_mesh.begin_active();
  for (; cell!=boundary_mesh.end(); ++cell)
    {
      deallog << "Cell = " << cell << std::endl;
      deallog << "  direction_flag = " << cell->direction_flag() << std::endl;

      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
        {
          deallog << "  face = " << face
                  << "  (neighbor = " << cell->neighbor(face) << ")"
                  << std::endl;

          if (cell->face(face)->has_children())
            for (unsigned int c=0; c<cell->face(face)->n_children(); ++c)
              {
                deallog << "    subface = " << c << std::endl;
                deallog << "              "
                        << cell->neighbor_child_on_subface(face, c)
                        << std::endl;
              }
        }
    }

  return 0;
}

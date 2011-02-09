//----------------------------  hanging_nodes_02.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  hanging_nodes_02.cc  ---------------------------


// an extract of the _01 test that shows the essence of what's going
// wrong. compared to the _03 testcase, upon refinement of the middle
// cell (0.0), cell 0.1 forgets who some of its neighbors are. this is
// clearly not good
//
// the underlying reason was that when we update neighborship
// information in the triangulation class, we have to take into
// account the direction_flag of cells

#include "../tests.h"

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_out.h>
#include <grid/grid_tools.h>


int main ()
{
  std::ofstream logfile("hanging_nodes_02/output");
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
    std::set<unsigned char> boundary_ids;
    boundary_ids.insert(0);
    GridTools::extract_boundary_mesh (volume_mesh, boundary_mesh, boundary_ids);
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

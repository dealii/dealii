//----------------------------  hanging_nodes_01.cc  ---------------------------
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
//----------------------------  hanging_nodes_01.cc  ---------------------------


// we crash when building hanging nodes. however, as elucidated by the
// _02 testcase, the reason is actually somewhere entirely different:
// upon refinement of one cell, we forget who the neighbors of another
// cell are. from there, no recovery is possible

#include "../tests.h"

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_out.h>
#include <grid/grid_tools.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <lac/constraint_matrix.h>


int main ()
{
  std::ofstream logfile("hanging_nodes_01/output");
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

				   /*****************************************************************/
  FE_Q<dim,spacedim> fe(1);
  DoFHandler<dim,spacedim>  dh (boundary_mesh);
  ConstraintMatrix     hanging_node_constraints;

  dh.distribute_dofs(fe);
  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dh,hanging_node_constraints);
  hanging_node_constraints.close ();

  hanging_node_constraints.print (deallog.get_file_stream());

  return 0;
}

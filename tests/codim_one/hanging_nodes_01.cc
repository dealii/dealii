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



// we crash when building hanging nodes. however, as elucidated by the
// _02 testcase, the reason is actually somewhere entirely different:
// upon refinement of one cell, we forget who the neighbors of another
// cell are. from there, no recovery is possible

#include "../tests.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/constraint_matrix.h>


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

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// we crash when building hanging nodes. however, as elucidated by the
// _02 testcase, the reason is actually somewhere entirely different:
// upon refinement of one cell, we forget who the neighbors of another
// cell are. from there, no recovery is possible

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>

#include "../tests.h"


int
main()
{
  initlog();

  const unsigned int spacedim = 3;
  const unsigned int dim      = spacedim - 1;

  Triangulation<dim, spacedim> boundary_mesh;

  /*****************************************************************/
  // Create Surface Mesh:  Boundary of hypercube without one face
  {
    Triangulation<spacedim> volume_mesh;
    GridGenerator::hyper_cube(volume_mesh);
    Triangulation<spacedim>::active_cell_iterator cell =
      volume_mesh.begin_active();

    cell->face(0)->set_all_boundary_ids(1);
    const std::set<types::boundary_id> boundary_ids = {0};
    GridGenerator::extract_boundary_mesh(volume_mesh,
                                         boundary_mesh,
                                         boundary_ids);
  }
  boundary_mesh.begin_active()->set_refine_flag();
  boundary_mesh.execute_coarsening_and_refinement();

  /*****************************************************************/
  FE_Q<dim, spacedim>       fe(1);
  DoFHandler<dim, spacedim> dh(boundary_mesh);
  AffineConstraints<double> hanging_node_constraints;

  dh.distribute_dofs(fe);
  hanging_node_constraints.clear();
  DoFTools::make_hanging_node_constraints(dh, hanging_node_constraints);
  hanging_node_constraints.close();

  hanging_node_constraints.print(deallog.get_file_stream());

  return 0;
}

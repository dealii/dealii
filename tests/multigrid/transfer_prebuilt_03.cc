// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Output transfer matrices on locally refined meshes with hanging node
// constraints

#include <deal.II/base/function_lib.h>
#include <deal.II/base/mg_level_object.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>

#include <algorithm>

#include "../tests.h"



template <int dim>
void
check_simple(const FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;

  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();
  tr.begin(2)->set_refine_flag();
  tr.execute_coarsening_and_refinement();

  DoFHandler<dim> mgdof(tr);
  mgdof.distribute_dofs(fe);
  mgdof.distribute_mg_dofs();

  MGConstrainedDoFs mg_constrained_dofs;
  mg_constrained_dofs.initialize(mgdof);
  mg_constrained_dofs.make_zero_boundary_constraints(mgdof, {0});

  MGTransferPrebuilt<Vector<double>> transfer(mg_constrained_dofs);
  transfer.build(mgdof);

  transfer.print_matrices(deallog.get_file_stream());
  transfer.print_indices(deallog.get_file_stream());
}


int
main()
{
  initlog(__FILE__);

  check_simple(FE_DGP<2>(0));
  check_simple(FE_DGP<2>(1));
  check_simple(FE_DGQ<2>(1));
  check_simple(FE_Q<2>(1));
  check_simple(FE_Q<2>(2));
  check_simple(FESystem<2>(FE_DGQ<2>(1), 2));

  check_simple(FE_RaviartThomas<2>(1));
  check_simple(FESystem<2>(FE_RaviartThomas<2>(1), 1, FE_DGQ<2>(1), 1));

  check_simple(FE_DGQ<3>(1));
  check_simple(FE_Q<3>(2));
}

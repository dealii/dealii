// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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


// Output transfer matrices on locally refined meshes with hanging node constraints

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>

#include <fstream>
#include <iomanip>
#include <iomanip>
#include <algorithm>

using namespace std;

template <int dim>
void check_simple(const FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;

  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();
  tr.begin(2)->set_refine_flag();
  tr.execute_coarsening_and_refinement();

  MGDoFHandler<dim> mgdof(tr);
  mgdof.distribute_dofs(fe);

  ConstraintMatrix     hanging_node_constraints;
  DoFTools::make_hanging_node_constraints (mgdof, hanging_node_constraints);
  hanging_node_constraints.close ();

  typename FunctionMap<dim>::type dirichlet_boundary;
  ZeroFunction<dim> homogeneous_dirichlet_bc (1);
  dirichlet_boundary[0] = &homogeneous_dirichlet_bc;

  MGConstrainedDoFs mg_constrained_dofs;
  mg_constrained_dofs.initialize(mgdof, dirichlet_boundary);

  MGTransferPrebuilt<Vector<double> > transfer(hanging_node_constraints, mg_constrained_dofs);
  transfer.build_matrices(mgdof);

  transfer.print_matrices(deallog.get_file_stream());
  transfer.print_indices(deallog.get_file_stream());
}


int main()
{
  initlog(__FILE__);

  check_simple (FE_DGP<2>(0));
  check_simple (FE_DGP<2>(1));
  check_simple (FE_DGQ<2>(1));
  check_simple (FE_Q<2>(1));
  check_simple (FE_Q<2>(2));
  check_simple (FESystem<2>(FE_DGQ<2>(1), 2));

  check_simple (FE_RaviartThomas<2>(1));
  check_simple (FESystem<2>(FE_RaviartThomas<2>(1), 1, FE_DGQ<2>(1), 1));

  check_simple (FE_DGQ<3>(1));
  check_simple (FE_Q<3>(2));
}

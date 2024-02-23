// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// have a 2x1 coarse mesh (or 2x1x1) and verify DoF indices in the hp-
// case with an FECollection that contains multiple copies of the same
// FE_Q(2) element. in the sequential case, the hp-code will unify DoF
// indices on boundaries between locally owned subdomains; in early
// versions of the parallel hp-support, we don't do that, but the
// final version now does
//
// this test gives a different perspective on this issue. in the _01
// test, we output the DoF indices on all cells as viewed by different
// processors. here, we also compute hanging node constraints which
// however should all have been resolved to identities during DoF
// index unification -- so there should be no constraints at all here


#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>

#include <numeric>

#include "../tests.h"

#include "../test_grids.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> triangulation(
    MPI_COMM_WORLD, Triangulation<dim>::limit_level_difference_at_vertices);
  TestGrids::hyper_line(triangulation, 2);
  Assert(triangulation.n_active_cells() == 2, ExcInternalError());

  hp::FECollection<dim> fe(FE_Q<dim>(2), FE_Q<dim>(2));

  DoFHandler<dim> dof_handler(triangulation);
  if (dof_handler.begin_active()->is_locally_owned())
    dof_handler.begin_active()->set_active_fe_index(0);
  if ((std::next(dof_handler.begin_active()))->is_locally_owned())
    (std::next(dof_handler.begin_active()))->set_active_fe_index(1);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  deallog << "n_dofs: " << dof_handler.n_dofs() << std::endl;
  deallog << "n_constraints: " << constraints.n_constraints() << std::endl;
  constraints.print(deallog.get_file_stream());
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}

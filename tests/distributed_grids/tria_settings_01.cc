// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test the setting
// parallel::distributed::Triangulation<dim>::mesh_reconstruction_after_repartitioning
// by looking which order the cells get enumerated. With the flag, the cells
// stay in z-order. Without it, the order in which the cells are created
// matters.

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

#include "coarse_grid_common.h"


template <int dim>
void
testit(parallel::distributed::Triangulation<dim> &tr)
{
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);

  (std::next(tr.begin_active()))->set_refine_flag();
  tr.execute_coarsening_and_refinement();

  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();

  typename parallel::distributed::Triangulation<dim>::active_cell_iterator it =
    tr.begin_active();
  for (; it != tr.end(); ++it)
    {
      deallog << it->center() << ", ";
    }
  deallog << std::endl;
}


template <int dim>
void
test(std::ostream & /*out*/)
{
  {
    parallel::distributed::Triangulation<dim> tr(
      MPI_COMM_WORLD,
      dealii::Triangulation<dim>::none,
      parallel::distributed::Triangulation<dim>::default_setting);
    testit(tr);
  }

  {
    parallel::distributed::Triangulation<dim> tr(
      MPI_COMM_WORLD,
      dealii::Triangulation<dim>::none,
      parallel::distributed::Triangulation<
        dim>::mesh_reconstruction_after_repartitioning);
    testit(tr);
  }
}


int
main(int argc, char *argv[])
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  deallog.push("2d");
  test<2>(deallog.get_file_stream());
  deallog.pop();
}

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test that DofRenumbering::hierarchical() also works for a
// parallel::shared::Triangulation


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/trilinos_vector.h>

#include "../tests.h"



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();

  const int dim = 2;

  dealii::FE_DGQ<dim> fe(0);

  // create a p::s::Triangulation. specifically choose z-order, for
  // which the algorithm previously worked by accident (because the
  // partitioning of cells was contiguous in z-order)

  dealii::parallel::shared::Triangulation<dim> tria_world(
    MPI_COMM_WORLD,
    dealii::Triangulation<dim>::none,
    false,
    parallel::shared::Triangulation<dim>::partition_zorder);
  dealii::GridGenerator::hyper_cube(tria_world, -1, 1);
  tria_world.refine_global(2);

  dealii::DoFHandler<dim> dh_world(tria_world);
  dh_world.distribute_dofs(fe);
  dealii::DoFRenumbering::hierarchical(dh_world);

  dealii::parallel::shared::Triangulation<dim> tria_self(MPI_COMM_SELF);
  dealii::GridGenerator::hyper_cube(tria_self, -1, 1);
  tria_self.refine_global(2);

  dealii::DoFHandler<dim> dh_self(tria_self);
  dh_self.distribute_dofs(fe);
  dealii::DoFRenumbering::hierarchical(dh_self);

  typename dealii::DoFHandler<dim>::active_cell_iterator
    cell_world = dh_world.begin_active(),
    endc_world = dh_world.end(), cell_self = dh_self.begin_active(),
    endc_self = dh_self.end();
  for (; cell_world != endc_world && cell_self != endc_self;
       ++cell_world, ++cell_self)
    if (cell_world->is_locally_owned())
      Assert(cell_world->dof_index(0) == cell_self->dof_index(0),
             ExcInternalError());

  deallog << "OK" << std::endl;
}

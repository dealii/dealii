// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// document a bug in MGLevelGlobalTransfer in 3d:

/*
An error occurred in line <341> of file
<../source/multigrid/mg_level_global_transfer.cc> in function void
dealii::MGLevelGlobalTransfer<dealii::LinearAlgebra::distributed::Vector<Number,
dealii::MemorySpace::Host> >::fill_and_communicate_copy_indices(const
dealii::DoFHandler<dh_dim, spacedim>&) [with int dim = 3; int spacedim = 3;
Number = double] The violated condition was:
    ::dealii::deal_II_exceptions::internals::compare_for_equality(this->copy_indices_level_mine.back().n_rows(),
0) Additional information: Dimension 2 not equal to 0.
*/

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include "../tests.h"

template <int dim>
void
check()
{
  parallel::distributed::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    typename parallel::distributed::Triangulation<dim>::Settings(
      parallel::distributed::Triangulation<
        dim>::mesh_reconstruction_after_repartitioning |
      parallel::distributed::Triangulation<
        dim>::construct_multigrid_hierarchy));

  GridGenerator::hyper_cube(tr);

  tr.refine_global(2);

  // coarsen a few cells:
  for (const auto &cell : tr.active_cell_iterators())
    {
      if (cell->is_artificial())
        continue;
      const std::string id = cell->parent()->id().to_string();
      if (id == "0_1:0" || id == "0_1:1" || id == "0_1:2" || id == "0_1:3")
        cell->set_coarsen_flag();
    }

  tr.execute_coarsening_and_refinement();

  FE_Q<dim>                 fe(1);
  DoFHandler<dim>           mgdof(tr);
  MGTransferMF<dim, double> transfer_mf;

  mgdof.distribute_dofs(fe);
  mgdof.distribute_mg_dofs();

  deallog << "n_dofs = " << mgdof.n_dofs() << std::endl;

  transfer_mf.build(mgdof);

  tr.coarsen_global();
  deallog << "OK" << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    all;

  check<2>();
  check<3>();
}

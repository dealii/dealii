// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Distribute DoFs as we do in hp_unify_dof_indices_02, but then
// renumber DoFs. In this test, we use the identity renumbering, i.e.,
// this does not actually change any DoF indices, but it still runs
// through the renumbering machinery to verify that things work
// correctly.


#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>

#include <numeric>

#include "../tests.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> triangulation(
    MPI_COMM_WORLD, Triangulation<dim>::limit_level_difference_at_vertices);

  std::vector<unsigned int> reps(dim, 1U);
  reps[0] = 2;
  reps[1] = 2;
  Point<dim> top_right;
  for (unsigned int d = 0; d < dim; ++d)
    top_right[d] = (d == 0 ? 2 : 1);
  GridGenerator::subdivided_hyper_rectangle(triangulation,
                                            reps,
                                            Point<dim>(),
                                            top_right);
  Assert(triangulation.n_global_active_cells() == 4, ExcInternalError());
  Assert(triangulation.n_active_cells() == 4, ExcInternalError());

  hp::FECollection<dim> fe;
  fe.push_back(FE_Q<dim>(2));
  fe.push_back(FE_Q<dim>(2));
  fe.push_back(FE_Q<dim>(2));

  DoFHandler<dim> dof_handler(triangulation);
  for (auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          if (cell->id().to_string() == "0_0:")
            cell->set_active_fe_index(0);
          if (cell->id().to_string() == "1_0:")
            cell->set_active_fe_index(1);
          if (cell->id().to_string() == "2_0:")
            cell->set_active_fe_index(2);
          if (cell->id().to_string() == "3_0:")
            cell->set_active_fe_index(0);
        }
    }
  dof_handler.distribute_dofs(fe);

  // now do the identity renumbering
  const IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
  std::vector<types::global_dof_index> new_indices(
    locally_owned_dofs.n_elements());
  for (unsigned int i = 0; i < new_indices.size(); ++i)
    new_indices[i] = locally_owned_dofs.nth_index_in_set(i);
  dof_handler.renumber_dofs(new_indices);

  deallog << "Processor: " << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
          << std::endl;
  for (auto &cell : dof_handler.active_cell_iterators())
    {
      deallog << "  Cell: " << cell << std::endl;

      std::vector<types::global_dof_index> dof_indices(
        cell->get_fe().dofs_per_cell);
      cell->get_dof_indices(dof_indices);
      deallog << "    ";
      for (auto i : dof_indices)
        deallog << i << ' ';
      deallog << std::endl;
    }
  deallog << "  n_locally_owned_dofs: " << dof_handler.n_locally_owned_dofs()
          << std::endl;
  deallog << "  n_global_dofs: " << dof_handler.n_dofs() << std::endl;
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

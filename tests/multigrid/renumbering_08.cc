// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check renumbering the degrees of freedom on the multigrid levels in
// parallel for FE_Q, otherwise similar to renumbering_05

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <algorithm>

#include "../tests.h"



template <int dim>
void
print_dof_numbers(const DoFHandler<dim> &dof)
{
  std::vector<types::global_dof_index> dof_indices(dof.get_fe().dofs_per_cell);
  deallog << "DoF numbers on active cells" << std::endl;
  for (auto &cell : dof.active_cell_iterators())
    if (!cell->is_artificial())
      {
        cell->get_dof_indices(dof_indices);
        deallog << "cell " << cell->id() << ": ";
        for (const types::global_dof_index i : dof_indices)
          deallog << i << ' ';
        deallog << std::endl;
      }
  for (unsigned int l = 0; l < dof.get_triangulation().n_global_levels(); ++l)
    {
      deallog << "DoF numbers on level " << l << std::endl;
      for (auto &cell : dof.cell_iterators_on_level(l))
        if (cell->level_subdomain_id() != numbers::artificial_subdomain_id)
          {
            cell->get_mg_dof_indices(dof_indices);
            deallog << "cell " << cell->id() << ": ";
            for (const types::global_dof_index i : dof_indices)
              deallog << i << ' ';
            deallog << std::endl;
          }
    }
}

template <int dim>
void
check()
{
  FE_Q<dim> fe(2);

  parallel::distributed::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  for (unsigned int cycle = 0; cycle < 5; ++cycle)
    {
      tr.clear();
      const unsigned int n_refine  = cycle / 3;
      const unsigned int remainder = cycle % 3;
      Point<dim>         p1;
      for (unsigned int d = 0; d < dim; ++d)
        p1[d] = -1;
      Point<dim> p2;
      for (unsigned int d = 0; d < remainder; ++d)
        p2[d] = 2.8;
      for (unsigned int d = remainder; d < dim; ++d)
        p2[d] = 1;
      std::vector<unsigned int> subdivisions(dim, 1);
      for (unsigned int d = 0; d < remainder; ++d)
        subdivisions[d] = 2;
      GridGenerator::subdivided_hyper_rectangle(tr, subdivisions, p1, p2);

      DoFHandler<dim> mgdof(tr);
      mgdof.distribute_dofs(fe);
      mgdof.distribute_mg_dofs();

      print_dof_numbers(mgdof);

      // compute a renumbering on the level degrees of freedom by simply
      // flipping the local numbers
      for (unsigned int l = 0; l < tr.n_global_levels(); ++l)
        {
          std::vector<types::global_dof_index> new_indices;
          if (mgdof.locally_owned_mg_dofs(l).n_elements() > 0)
            {
              const types::global_dof_index first =
                mgdof.locally_owned_mg_dofs(l).nth_index_in_set(0);
              const types::global_dof_index last =
                first + mgdof.locally_owned_mg_dofs(l).n_elements();
              for (unsigned int i = 0; i < last - first; ++i)
                new_indices.push_back(last - 1 - i);
            }
          mgdof.renumber_dofs(l, new_indices);
        }

      print_dof_numbers(mgdof);
    }
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv);
  MPILogInitAll                    log;

  check<2>();
  check<3>();
}

// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// Testcase from https://github.com/dealii/dealii/issues/7174 that shows
// broken DoF renumbering when using GMG

// Note: This test currently only works in serial and asserts when run with
// more than one processor.

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <algorithm>
#include <memory>

#include "../tests.h"

using namespace std;


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
        for (types::global_dof_index i : dof_indices)
          deallog << i << " ";
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
            for (types::global_dof_index i : dof_indices)
              deallog << i << " ";
            deallog << std::endl;
          }
    }
}

template <int dim>
void
check()
{
  dealii::parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    decltype(tria)::limit_level_difference_at_vertices,
    decltype(tria)::construct_multigrid_hierarchy);
  dealii::GridGenerator::hyper_cube(tria);
  tria.refine_global(3);
  /*
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
   tr.begin_active(tr.n_global_levels() - 1)->set_refine_flag();
 tr.execute_coarsening_and_refinement();
  */

  auto fe_scalar = std::make_unique<dealii::FE_Q<dim>>(1);
  auto fe        = std::make_unique<dealii::FESystem<dim>>(*fe_scalar, dim);

  dealii::DoFHandler<dim> dofhandler(tria);
  dofhandler.distribute_dofs(*fe);
  dofhandler.distribute_mg_dofs();

  print_dof_numbers(dofhandler);

  dealii::DoFRenumbering::component_wise(dofhandler);
  deallog << "Finished fine lvl renumbering" << std::endl;

  for (unsigned int lvl = 0; lvl < tria.n_global_levels(); lvl++)
    {
      dealii::DoFRenumbering::component_wise(dofhandler, lvl);
      deallog << "Finished renumbering on lvl " << lvl << std::endl;
    }

  print_dof_numbers(dofhandler);
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv);
  MPILogInitAll                    log;
  deallog.depth_console(5); ///////!

  check<2>();
  //  check<3>();
}

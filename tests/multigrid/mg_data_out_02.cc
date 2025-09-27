// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test parallel DataOut::set_cell_selection for multilevel cells

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"


template <int dim>
void
print(DataOut<dim> &data_out, const Triangulation<dim> &tria)
{
  auto p = data_out.get_cell_selection();

  for (auto cell = p.first(tria); cell.state() == IteratorState::valid;
       cell      = p.second(tria, cell))
    {
      deallog << " cell lvl=" << cell->level()
              << " active=" << cell->is_active()
              << " owner=" << static_cast<int>(cell->level_subdomain_id())
              << " id=" << cell->id().to_string() << std::endl;
    }
}

template <int dim>
void
do_test()
{
  parallel::distributed::Triangulation<dim> triangulation(
    MPI_COMM_WORLD,
    dealii::Triangulation<dim>::none,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);
  triangulation.begin_active()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();

  deallog << "dim= " << dim << " cells=" << triangulation.n_cells()
          << std::endl;

  DataOut<dim> data_out;
  data_out.attach_triangulation(triangulation);

  deallog << "* default:" << std::endl;
  print(data_out, triangulation);

  deallog << "* all cells:" << std::endl;
  data_out.set_cell_selection(
    [](const typename Triangulation<dim>::cell_iterator &cell) {
      return true;
    });
  print(data_out, triangulation);

  deallog << "* LocallyOwnedLevelCell:" << std::endl;
  data_out.set_cell_selection(IteratorFilters::LocallyOwnedLevelCell());
  print(data_out, triangulation);

  for (unsigned int level = 0; level < triangulation.n_global_levels(); ++level)
    {
      DataOut<dim> data_out;
      data_out.attach_triangulation(triangulation);

      deallog << "* LevelEqualTo " << level << std::endl;
      data_out.set_cell_selection(IteratorFilters::LevelEqualTo(level));
      print(data_out, triangulation);

      deallog << "* owned on level " << level << std::endl;
      data_out.set_cell_selection(
        [level](const typename Triangulation<dim>::cell_iterator &cell) {
          return (cell->level() == static_cast<int>(level) &&
                  cell->is_locally_owned_on_level());
        });
      print(data_out, triangulation);
    }
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    log;

  do_test<2>();
  do_test<3>();
  return 0;
}

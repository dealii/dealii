// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// validate algorithms that will flag cells for p-adaptivity:
// - hp::Refinement::p_adaptivity_fixed_number


#include <deal.II/base/geometry_info.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/refinement.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"



template <int dim>
void
validate(const DoFHandler<dim> &dh)
{
  deallog << " (cellid,feidx):";
  for (const auto &cell : dh.active_cell_iterators())
    if (!cell->is_artificial() && cell->is_locally_owned())
      {
        const std::string  cellid = cell->id().to_string();
        const unsigned int coarse_cellid =
          std::stoul(cellid.substr(0, cellid.find("_")));

        deallog << " (" << coarse_cellid << ',' << cell->future_fe_index()
                << ')';
      }
  deallog << std::endl;
}



template <int dim>
void
setup(Triangulation<dim> &tria, const DoFHandler<dim> &dh)
{
  // Initialize triangulation.
  GridGenerator::subdivided_hyper_cube(tria, 4);
  Assert(tria.n_cells(0) == tria.n_global_active_cells(), ExcInternalError());

  // Set all active FE indices to 1.
  // Flag first half of cells for refinement, and the other half for coarsening.
  for (const auto &cell : dh.active_cell_iterators())
    if (!cell->is_artificial() && cell->is_locally_owned())
      {
        cell->set_active_fe_index(1);

        const std::string  cellid = cell->id().to_string();
        const unsigned int coarse_cellid =
          std::stoul(cellid.substr(0, cellid.find('_')));

        if (coarse_cellid < 0.5 * tria.n_global_active_cells())
          cell->set_refine_flag();
        else
          cell->set_coarsen_flag();
      }
}



template <int dim>
void
test()
{
  hp::FECollection<dim> fes;
  for (unsigned int d = 1; d <= 3; ++d)
    fes.push_back(FE_Q<dim>(d));

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  DoFHandler<dim>                           dh(tria);
  setup(tria, dh);
  dh.distribute_dofs(fes);

  deallog << "starting situation" << std::endl;
  validate(dh);


  // We flag the first half of all cells to be refined and the last half of all
  // cells to be coarsened for p-adapativity. Ultimately, the first quarter of
  // all cells will be flagged for p-refinement, and the last quarter for p-
  // coarsening.

  Vector<double> indicators(tria.n_active_cells());
  for (const auto &cell : tria.active_cell_iterators())
    if (!cell->is_artificial() && cell->is_locally_owned())
      {
        const std::string  cellid = cell->id().to_string();
        const unsigned int coarse_cellid =
          std::stoul(cellid.substr(0, cellid.find('_')));

        if (coarse_cellid < .25 * tria.n_global_active_cells())
          indicators[cell->active_cell_index()] = 2.;
        else if (coarse_cellid < .75 * tria.n_global_active_cells())
          indicators[cell->active_cell_index()] = 1.;
        else
          indicators[cell->active_cell_index()] = 0.;
      }

  hp::Refinement::p_adaptivity_fixed_number(dh, indicators);

  deallog << "p-adaptivity fixed number" << std::endl;
  validate(dh);


  deallog << "OK" << std::endl;
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

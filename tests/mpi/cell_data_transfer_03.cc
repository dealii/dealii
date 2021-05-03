// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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



// p::d::CellDataTransfer test for an arbitrary value type
// and variable data sizes


#include <deal.II/distributed/cell_data_transfer.templates.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>

#include <numeric>

#include "../tests.h"


template <int dim, int spacedim>
std::vector<int>
get_data_of_first_child(
  const typename dealii::Triangulation<dim, spacedim>::cell_iterator &,
  const std::vector<std::vector<int>> &children_values)
{
  return children_values[0];
}



template <int dim, int spacedim>
void
test()
{
  // ------ setup ------
  parallel::distributed::Triangulation<dim, spacedim> tria(MPI_COMM_WORLD);
  GridGenerator::subdivided_hyper_cube(tria, 2);
  tria.refine_global(1);
  deallog << "cells before: " << tria.n_global_active_cells() << std::endl;

  // ----- prepare -----
  // set refinement/coarsening flags
  for (auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        if (cell->id().to_string() == "0_1:0")
          cell->set_refine_flag();
        else if (cell->parent()->id().to_string() ==
                 ((dim == 2) ? "3_0:" : "7_0:"))
          cell->set_coarsen_flag();
      }

  // ----- gather -----
  // store increasing amount of integers on all cells
  std::vector<std::vector<int>> cell_data(tria.n_active_cells());
  int                           i = 0;
  for (auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        std::vector<int> this_cell_data(i++);
        std::iota(this_cell_data.begin(), this_cell_data.end(), 0);

        deallog << "cellid=" << cell->id()
                << " vector_size=" << this_cell_data.size() << " vector_sum="
                << std::accumulate(this_cell_data.begin(),
                                   this_cell_data.end(),
                                   0);
        if (cell->coarsen_flag_set())
          deallog << " coarsening";
        else if (cell->refine_flag_set())
          deallog << " refining";
        deallog << std::endl;

        cell_data[cell->active_cell_index()] = std::move(this_cell_data);
      }

  // ----- transfer -----
  parallel::distributed::
    CellDataTransfer<dim, spacedim, std::vector<std::vector<int>>>
    cell_data_transfer(
      tria,
      /*transfer_variable_size_data=*/true,
      /*refinement_strategy=*/
      &dealii::AdaptationStrategies::Refinement::
        preserve<dim, spacedim, std::vector<int>>,
      /*coarsening_strategy=*/&get_data_of_first_child<dim, spacedim>);

  cell_data_transfer.prepare_for_coarsening_and_refinement(cell_data);
  tria.execute_coarsening_and_refinement();
  deallog << "cells after: " << tria.n_global_active_cells() << std::endl;

  cell_data.resize(tria.n_active_cells());
  cell_data_transfer.unpack(cell_data);

  // ------ verify ------
  // check if all children adopted the correct data
  for (auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        const std::vector<int> &this_cell_data =
          cell_data[(cell->active_cell_index())];

        deallog << "cellid=" << cell->id()
                << " vector_size=" << this_cell_data.size() << " vector_sum="
                << std::accumulate(this_cell_data.begin(),
                                   this_cell_data.end(),
                                   0)
                << std::endl;
      }

  // make sure no processor is hanging
  MPI_Barrier(MPI_COMM_WORLD);

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog.push("2d");
  test<2, 2>();
  deallog.pop();
  deallog.push("3d");
  test<3, 3>();
  deallog.pop();
}

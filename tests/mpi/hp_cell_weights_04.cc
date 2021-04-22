// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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



// Cell Weights Test
// -----------------
// Create a 4x4(x4) grid, on which all cells are associated with a Q1
// element besides the very first one, which has a Q5 element.
// We choose a cell weighting algorithm based on the number of degrees
// of freedom and check if load is balanced as expected after
// repartitioning the triangulation. The expected accumulated weight on
// each processor should correlate to the sum of all degrees of
// freedom on all cells of the corresponding subdomain.
// We employ a large proportionality factor on our weighting function
// to neglect the standard weight of '1000' per cell.
//
// This test works on a parallel::shared::Triangulation with METIS
// as a partitioner. Cell weighting with ZOLTAN was not available
// during the time this test was written.
//
// We consider aritifical cells in this case.


#include <deal.II/distributed/cell_weights.h>
#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"


template <int dim>
void
test()
{
  parallel::shared::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    true,
    parallel::shared::Triangulation<dim>::Settings::partition_metis);

  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  // Apply ndof cell weights.
  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FE_Q<dim>(1));
  fe_collection.push_back(FE_Q<dim>(5));

  DoFHandler<dim> dh(tria);

  // default: active_fe_index = 0
  for (auto &cell : dh.active_cell_iterators())
    if (cell->is_locally_owned())
      if (cell->id().to_string() == "0_2:00")
        cell->set_active_fe_index(1);

  dh.distribute_dofs(fe_collection);

  deallog << "Number of cells before repartitioning: "
          << tria.n_locally_owned_active_cells() << std::endl;
  {
    unsigned int dof_counter = 0;
    for (auto &cell : dh.active_cell_iterators())
      if (cell->is_locally_owned())
        dof_counter += cell->get_fe().dofs_per_cell;
    deallog << "  Cumulative dofs per cell: " << dof_counter << std::endl;
  }


  const parallel::CellWeights<dim> cell_weights(
    dh, parallel::CellWeights<dim>::ndofs_weighting({100000, 1}));

  // we didn't mark any cells, but we want to repartition our domain
  tria.execute_coarsening_and_refinement();


  deallog << "Number of cells after repartitioning: "
          << tria.n_locally_owned_active_cells() << std::endl;
  {
    unsigned int dof_counter = 0;
    for (auto &cell : dh.active_cell_iterators())
      if (cell->is_locally_owned())
        dof_counter += cell->get_fe().dofs_per_cell;
    deallog << "  Cumulative dofs per cell: " << dof_counter << std::endl;
  }

#ifdef DEBUG
  parallel::shared::Triangulation<dim> other_tria(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    false,
    parallel::shared::Triangulation<dim>::Settings::partition_metis);
  GridGenerator::hyper_cube(other_tria);
  other_tria.refine_global(3);

  dh.reinit(other_tria);
  dh.distribute_dofs(fe_collection);

  try
    {
      tria.execute_coarsening_and_refinement();
    }
  catch (const ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }
#else
  deallog
    << "ExcMessage(\"Triangulation associated with the DoFHandler has changed!\")"
    << std::endl;
#endif

  // make sure no processor is hanging
  MPI_Barrier(MPI_COMM_WORLD);

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deal_II_exceptions::disable_abort_on_exception();

  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();
}

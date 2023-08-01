// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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


// Similar as copy_serial_tria_01 but permute the cells on the coarse grid
// s.t. the cells in the auto-created TriangulationDescription::cell_infos
// does not match the order of coarse cells and sorting is required in
// Triangulation::create_triangulation() (see also PR #10092).

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include "../grid/tests.h"

using namespace dealii;

template <int dim>
void
test(int n_refinements, MPI_Comm comm)
{
  // create serial triangulation
  Triangulation<dim> basetria;
  GridGenerator::hyper_L(basetria);
  basetria.refine_global(n_refinements);

  GridTools::partition_triangulation_zorder(
    Utilities::MPI::n_mpi_processes(comm), basetria);

  // create instance of pft
  parallel::fullydistributed::Triangulation<dim> tria_pft(comm);

  // extract relevant information form serial triangulation
  auto construction_data =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      basetria, comm);

  for (unsigned int i = 0; i < construction_data.coarse_cells.size() / 2; ++i)
    {
      std::swap(construction_data.coarse_cells[i],
                construction_data
                  .coarse_cells[construction_data.coarse_cells.size() - 1 - i]);
      std::swap(construction_data.coarse_cell_index_to_coarse_cell_id[i],
                construction_data.coarse_cell_index_to_coarse_cell_id
                  [construction_data.coarse_cells.size() - 1 - i]);
    }

  // actually create triangulation
  tria_pft.create_triangulation(construction_data);

  // test triangulation
  FE_Q<dim>       fe(2);
  DoFHandler<dim> dof_handler(tria_pft);
  dof_handler.distribute_dofs(fe);

  // print statistics
  print_statistics(tria_pft);
  print_statistics(dof_handler);
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  const int      n_refinements = 3;
  const MPI_Comm comm          = MPI_COMM_WORLD;

  {
    deallog.push("2d");
    test<2>(n_refinements, comm);
    deallog.pop();
  }
  {
    deallog.push("3d");
    test<3>(n_refinements, comm);
    deallog.pop();
  }
}

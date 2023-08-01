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


// Create a serial triangulation with multigrid levels and copy it with .
// copy_triangulation.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../grid/tests.h"

using namespace dealii;

template <int dim>
void
test(int n_refinements, const int n_subdivisions, MPI_Comm comm)
{
  // create serial triangulation
  Triangulation<dim> basetria(
    Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::subdivided_hyper_cube(basetria, n_subdivisions);
  basetria.refine_global(n_refinements);

  // create instance of pft
  parallel::fullydistributed::Triangulation<dim> tria_pft(comm);

  // register dynamic construction data
  tria_pft.set_partitioner(
    [](dealii::Triangulation<dim> &tria, const unsigned int n_partitions) {
      GridTools::partition_triangulation_zorder(n_partitions, tria);
    },
    TriangulationDescription::Settings::construct_multigrid_hierarchy);

  // actually create triangulation
  tria_pft.copy_triangulation(basetria);

  // test triangulation
  FE_Q<dim>       fe(2);
  DoFHandler<dim> dof_handler(tria_pft);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  // print statistics
  print_statistics(tria_pft, true);
  print_statistics(dof_handler, true);
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  const MPI_Comm comm = MPI_COMM_WORLD;


  {
    deallog.push("1d");
    const int n_refinements  = 6;
    const int n_subdivisions = 8;
    test<1>(n_refinements, n_subdivisions, comm);
    deallog.pop();
  }
  {
    deallog.push("2d");
    const int n_refinements  = 3;
    const int n_subdivisions = 8;
    test<2>(n_refinements, n_subdivisions, comm);
    deallog.pop();
  }

  {
    deallog.push("3d");
    const int n_refinements  = 3;
    const int n_subdivisions = 4;
    test<3>(n_refinements, n_subdivisions, comm);
    deallog.pop();
  }
}

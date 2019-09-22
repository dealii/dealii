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


// create a tria mesh and copy it

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/fully_distributed_tria_util.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

using namespace dealii;

template <int dim>
void
test(const int /*n_refinements*/, const int /*n_subdivisions*/, MPI_Comm comm)
{
  parallel::distributed::Triangulation<dim> tria_pdt(
    comm,
    dealii::Triangulation<dim>::none,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::subdivided_hyper_cube(tria_pdt, 2);
  tria_pdt.refine_global(2);

  GridOut grid_out;
  grid_out.write_mesh_per_processor_as_vtu(tria_pdt, "tria_pdt", true, true);

  // extract relevant information form pdt
  auto construction_data = parallel::fullydistributed::Utilities::
    create_construction_data_from_triangulation(tria_pdt, comm, true);

  // create instance of pft
  parallel::fullydistributed::Triangulation<dim> tria_pft(comm);
  tria_pft.create_triangulation(construction_data);
  grid_out.write_mesh_per_processor_as_vtu(tria_pft, "tria_pft_0", true, true);

  // perform global coarsening (1x)
  for (auto cell : tria_pft.active_cell_iterators())
    cell->set_coarsen_flag();
  tria_pft.execute_coarsening_and_refinement();
  grid_out.write_mesh_per_processor_as_vtu(tria_pft, "tria_pft_1", true, true);
  for (auto cell : tria_pft.active_cell_iterators())
    cell->set_coarsen_flag();
  tria_pft.execute_coarsening_and_refinement();
  grid_out.write_mesh_per_processor_as_vtu(tria_pft, "tria_pft_2", true, true);

  // perform global refinement (4x)
  tria_pft.refine_global(1);
  grid_out.write_mesh_per_processor_as_vtu(tria_pft, "tria_pft_3", true, true);
  tria_pft.refine_global(1);
  grid_out.write_mesh_per_processor_as_vtu(tria_pft, "tria_pft_4", true, true);
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  const int      dim            = 2;
  const int      n_refinements  = 1;
  const int      n_subdivisions = 1;
  const MPI_Comm comm           = MPI_COMM_WORLD;

  if (dim == 1)
    test<1>(n_refinements, n_subdivisions, comm);
  else if (dim == 2)
    test<2>(n_refinements, n_subdivisions, comm);
  else if (dim == 3)
    test<3>(n_refinements, n_subdivisions, comm);
}

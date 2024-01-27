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
  // create instance of pft
  parallel::fullydistributed::Triangulation<dim> tria_pft(comm);

  // extract relevant information form serial triangulation
  parallel::fullydistributed::ConstructionData<dim, dim> construction_data;
  construction_data.cell_infos.resize(1);

  // create vertices
  for (unsigned int j = 0; j < 3; j++)
    for (unsigned int i = 0; i < 3; i++)
      construction_data.coarse_cell_vertices.emplace_back(i, j);

  // create cells
  for (unsigned int j = 0; j < 2; j++)
    for (unsigned int i = 0; i < 2; i++)
      {
        unsigned int gid = i + j * 2;

        // mapping
        construction_data.coarse_cell_index_to_coarse_cell_id.push_back(gid);

        // cell data
        CellData<dim> cell_data;
        cell_data.vertices[0] = i + 0 + (j + 0) * 3;
        cell_data.vertices[1] = i + 1 + (j + 0) * 3;
        cell_data.vertices[2] = i + 0 + (j + 1) * 3;
        cell_data.vertices[3] = i + 1 + (j + 1) * 3;
        construction_data.coarse_cells.push_back(cell_data);

        // cell info
        parallel::fullydistributed::CellData<dim> cell_info;
        cell_info.subdomain_id       = gid;
        cell_info.level_subdomain_id = gid;
        cell_info.id[0]              = gid;
        cell_info.id[1]              = dim;
        cell_info.id[2]              = 0;
        cell_info.id[3]              = 0;

        construction_data.cell_infos[0].push_back(cell_info);
      }

  construction_data.comm = comm;
  construction_data.settings =
    parallel::fullydistributed::Settings::construct_multigrid_hierarchy;

  // actually create triangulation
  tria_pft.create_triangulation(construction_data);

  GridOut grid_out;
  grid_out.write_mesh_per_processor_as_vtu(tria_pft, "tria_pft_0", true, true);
  tria_pft.execute_coarsening_and_refinement();
  grid_out.write_mesh_per_processor_as_vtu(tria_pft, "tria_pft_1", true, true);

  // perform global refinement (4x)
  // tria_pft.refine_global(1);
  // grid_out.write_mesh_per_processor_as_vtu(tria_pft, "tria_pft_1", true,
  // true); tria_pft.refine_global(1);
  // grid_out.write_mesh_per_processor_as_vtu(tria_pft, "tria_pft_2", true,
  // true); tria_pft.refine_global(1);
  // grid_out.write_mesh_per_processor_as_vtu(tria_pft, "tria_pft_3", true,
  // true);

  // perform global coarsening (1x)
  // for (auto cell : tria_pft.active_cell_iterators())
  //  cell->set_coarsen_flag();
  // tria_pft.execute_coarsening_and_refinement();
  // grid_out.write_mesh_per_processor_as_vtu(tria_pft, "tria_pft_4", true,
  // true); for (auto cell : tria_pft.active_cell_iterators())
  //  cell->set_coarsen_flag();
  // tria_pft.execute_coarsening_and_refinement();
  // grid_out.write_mesh_per_processor_as_vtu(tria_pft, "tria_pft_5", true,
  // true); for (auto cell : tria_pft.active_cell_iterators())
  //  cell->set_coarsen_flag();
  // tria_pft.execute_coarsening_and_refinement();
  // grid_out.write_mesh_per_processor_as_vtu(tria_pft, "tria_pft_6", true,
  // true);
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

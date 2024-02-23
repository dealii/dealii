// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Create a parallel::fullydistributed::Triangulation manually.

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


template <int dim>
void
test(MPI_Comm comm)
{
  // subdivided hypercube settings
  std::array<unsigned int, dim> n_subdivision = {{5, 6}};
  std::array<double, dim>       sizes         = {{1.0, 1.0}};

  // parallel info
  const auto n_ranks = Utilities::MPI::n_mpi_processes(comm);
  const auto my_rank = Utilities::MPI::this_mpi_process(comm);

  // create description
  TriangulationDescription::Description<dim, dim> construction_data;
  construction_data.comm = comm;
  construction_data.cell_infos.resize(1);

  const auto is_local =
    Utilities::create_evenly_distributed_partitioning(my_rank,
                                                      n_ranks,
                                                      n_subdivision[1]);

  if (!is_local.is_empty())
    {
      // 1) collect vertices
      for (unsigned int j = std::max<int>(0, is_local.nth_index_in_set(0) - 1);
           j <= std::min<unsigned int>(
                  is_local.nth_index_in_set(is_local.n_elements() - 1) + 3,
                  n_subdivision[1] + 1);
           ++j)
        for (unsigned int i = 0; i <= n_subdivision[0]; ++i)
          construction_data.coarse_cell_vertices.emplace_back(
            sizes[0] / n_subdivision[0] * i, sizes[1] / n_subdivision[1] * j);

      // 2) setup cells
      for (unsigned int j = std::max<int>(0, is_local.nth_index_in_set(0) - 1),
                        j_local = 0;
           j < std::min<unsigned int>(
                 is_local.nth_index_in_set(is_local.n_elements() - 1) + 2,
                 n_subdivision[1]);
           ++j, ++j_local)
        for (unsigned int i = 0; i < n_subdivision[0]; ++i)
          {
            const unsigned int c_global = i + j * n_subdivision[0];

            // 2a) coarse cell info
            CellData<dim> cell_data;

            const unsigned int stride = n_subdivision[0] + 1;
            cell_data.vertices[0]     = (i + 0) + (j_local + 0) * stride;
            cell_data.vertices[1]     = (i + 1) + (j_local + 0) * stride;
            cell_data.vertices[2]     = (i + 0) + (j_local + 1) * stride;
            cell_data.vertices[3]     = (i + 1) + (j_local + 1) * stride;

            construction_data.coarse_cells.emplace_back(cell_data);

            // 2b) local to global cell id map
            construction_data.coarse_cell_index_to_coarse_cell_id.emplace_back(
              c_global);

            // 2c) level cell info
            TriangulationDescription::CellData<dim> cell_data_;

            cell_data_.id = CellId(c_global, {}).template to_binary<dim>();
            cell_data_.subdomain_id =
              is_local.is_element(j) ?
                my_rank :
                (j < is_local.nth_index_in_set(0) ? (my_rank - 1) :
                                                    (my_rank + 1));
            cell_data_.level_subdomain_id = cell_data_.subdomain_id;

            if (i == 0)
              cell_data_.boundary_ids.emplace_back(0, 0);
            else if (i + 1 == n_subdivision[0])
              cell_data_.boundary_ids.emplace_back(1, 1);
            else if (j == 0)
              cell_data_.boundary_ids.emplace_back(2, 2);
            else if (j + 1 == n_subdivision[1])
              cell_data_.boundary_ids.emplace_back(3, 3);

            construction_data.cell_infos[0].emplace_back(cell_data_);
          }
    }

  // actually create triangulation
  parallel::fullydistributed::Triangulation<dim> tria_pft(comm);
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

  const MPI_Comm comm = MPI_COMM_WORLD;

  {
    deallog.push("2d");
    test<2>(comm);
    deallog.pop();
  }
}

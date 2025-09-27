// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/distributed/tria.h>

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <map>
#include <vector>

#include "../tests.h"

// test compute_local_to_global_vertex_index_map()

void
test()
{
  const MPI_Comm mpi_comm = MPI_COMM_WORLD;

  parallel::distributed::Triangulation<2> tria(mpi_comm);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);
  unsigned int rank = Utilities::MPI::this_mpi_process(mpi_comm);
  if (rank == 0)
    {
      typename Triangulation<2>::active_cell_iterator cell =
                                                        tria.begin_active(),
                                                      end_cell = tria.end();
      for (; cell != end_cell; ++cell)
        if (cell->is_locally_owned())
          cell->set_refine_flag();
    }

  tria.execute_coarsening_and_refinement();


  std::map<unsigned int, types::global_vertex_index> local_to_global_id =
    GridTools::compute_local_to_global_vertex_index_map(tria);

  std::vector<types::global_vertex_index>         vertex_global_index;
  typename Triangulation<2>::active_cell_iterator cell = tria.begin_active(),
                                                  endc = tria.end();
  for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
        {
          for (unsigned int i = 0; i < GeometryInfo<2>::lines_per_cell; ++i)
            {
              vertex_global_index.push_back(
                local_to_global_id[cell->line(i)->vertex_index(0)]);
              vertex_global_index.push_back(
                local_to_global_id[cell->line(i)->vertex_index(1)]);
              if (cell->line(i)->has_children())
                vertex_global_index.push_back(
                  local_to_global_id[cell->line(i)->child(0)->vertex_index(1)]);
            }
        }
    }
  std::sort(vertex_global_index.begin(), vertex_global_index.end());
  std::vector<types::global_vertex_index>::iterator it;
  it = std::unique(vertex_global_index.begin(), vertex_global_index.end());
  vertex_global_index.resize(it - vertex_global_index.begin());

  if (rank == 0)
    {
      std::vector<types::global_vertex_index> reference;
      for (unsigned int i = 0; i < 31; ++i)
        reference.push_back(i);
      for (unsigned int i = 0; i < vertex_global_index.size(); ++i)
        AssertThrow(reference[i] == vertex_global_index[i],
                    ExcMessage("Wrong global index"));
    }
  if (rank == 1)
    {
      std::vector<types::global_vertex_index> reference;
      reference.push_back(14);
      reference.push_back(18);
      reference.push_back(19);
      reference.push_back(20);
      reference.push_back(22);
      reference.push_back(23);
      reference.push_back(24);
      for (unsigned int i = 27; i < 55; ++i)
        reference.push_back(i);
      for (unsigned int i = 0; i < vertex_global_index.size(); ++i)
        AssertThrow(reference[i] == vertex_global_index[i],
                    ExcMessage("Wrong global index"));
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test();

  deallog << "OK" << std::endl;

  return 0;
}

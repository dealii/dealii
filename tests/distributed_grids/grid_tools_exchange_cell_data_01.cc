// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test GridTools::exchange_cell_data_to_ghosts

#include <deal.II/base/logstream.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <algorithm>

#include "../tests.h"

template <int dim>
void
test()
{
  const MPI_Comm mpi_communicator = MPI_COMM_WORLD;
  deallog << "dim = " << dim << std::endl;

  parallel::distributed::Triangulation<dim> tria(mpi_communicator);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  std::set<std::string> input, output;


  using DT = double;
  std::map<CellId, DT> map;
  DT                   counter = 0.0;

  std::map<unsigned int, std::set<dealii::types::subdomain_id>>
    vertices_with_ghost_neighbors =
      GridTools::compute_vertices_with_ghost_neighbors(tria);

  for (const auto &cell :
       tria.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
    {
      for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
        {
          const std::map<unsigned int,
                         std::set<dealii::types::subdomain_id>>::const_iterator
            neighbor_subdomains_of_vertex =
              vertices_with_ghost_neighbors.find(cell->vertex_index(v));

          if (neighbor_subdomains_of_vertex !=
              vertices_with_ghost_neighbors.end())
            {
              map[cell->id()] = ++counter;
              break;
            }
        }
    }

  using cell_iterator =
    typename parallel::distributed::Triangulation<dim>::active_cell_iterator;
  GridTools::
    exchange_cell_data_to_ghosts<DT, parallel::distributed::Triangulation<dim>>(
      tria,
      [&](const cell_iterator &cell) {
        DT                 value = map[cell->id()];
        std::ostringstream oss;
        oss << "pack " << cell->id() << ' ' << value;
        input.insert(oss.str());

        return value;
      },
      [&](const cell_iterator &cell, const DT &data) {
        std::ostringstream oss;
        oss << "unpack " << cell->id() << ' ' << data << " from "
            << cell->subdomain_id();

        output.insert(oss.str());
      });

  // sort the output because it will come in in random order
  for (auto &it : input)
    deallog << it << std::endl;
  for (auto &it : output)
    deallog << it << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test<2>();
  test<3>();

  return 0;
}

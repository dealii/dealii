// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// test GridTools::exchange_cell_data_to_ghosts

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <algorithm>

template <int dim>
void
test()
{
  const MPI_Comm& mpi_communicator = MPI_COMM_WORLD;
  deallog << "dim = " << dim << std::endl;

  parallel::distributed::Triangulation<dim> tria(mpi_communicator);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  std::set<std::string> output;

  typedef
    typename parallel::distributed::Triangulation<dim>::active_cell_iterator
                 cell_iterator;
  typedef double DT;
  DT             counter = 0.0;
  GridTools::
    exchange_cell_data_to_ghosts<DT, parallel::distributed::Triangulation<dim>>(
      tria,
      [&](const cell_iterator& cell) {
        DT value = ++counter;

        deallog << "pack " << cell->id() << " " << value << std::endl;
        return value;
      },
      [&](const cell_iterator& cell, const DT& data) {
        std::ostringstream oss;
        oss << "unpack " << cell->id() << " " << data << " from "
            << cell->subdomain_id();

        output.insert(oss.str());
      });

  // sort the output because it will come in in random order
  for(auto& it : output)
    deallog << it << std::endl;
}

int
main(int argc, char* argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test<2>();
  test<3>();

  return 0;
}

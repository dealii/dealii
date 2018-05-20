// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2017 by the deal.II authors
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

// save a tria mesh and load it

#include "../tests.h"
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <sstream>

template <int dim>
void
test()
{
  Triangulation<dim> tr1;

  GridGenerator::hyper_cube(tr1);
  tr1.refine_global(2);

  deallog << " n_active_cells: " << tr1.n_active_cells() << "\n" << std::endl;

  std::ostringstream oss;
  {
    boost::archive::text_oarchive oa(oss, boost::archive::no_header);
    tr1.save(oa, 0);
  }

  parallel::shared::Triangulation<dim> tr2(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    false,
    parallel::shared::Triangulation<dim>::partition_metis);
  {
    std::istringstream            iss(oss.str());
    boost::archive::text_iarchive ia(iss, boost::archive::no_header);
    tr2.load(ia, 0);
  }

  deallog << " n_active_cells: " << tr2.n_active_cells() << "\n"
          << " locally_owned_subdomain(): " << tr2.locally_owned_subdomain()
          << "\n"
          << " n_locally_owned_active_cells: "
          << tr2.n_locally_owned_active_cells() << "\n"
          << " n_global_active_cells: " << tr2.n_global_active_cells() << "\n"
          << std::endl;
}

int
main(int argc, char* argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();
}

// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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



// check that p::s::Tria.locally_owned_subdomain() works in 1d.

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  {
    parallel::shared::Triangulation<1> tria(
      MPI_COMM_WORLD,
      ::Triangulation<1>::none,
      false,
      parallel::shared::Triangulation<1>::partition_metis);

    deallog << tria.locally_owned_subdomain() << std::endl;
  }
}

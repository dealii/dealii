// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

// Test for the internal::send_and_receive function in the case of
// collective communication

#include "../tests.h"
#include "../../include/deal.II/base/mpi.templates.h"

#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>

#include <mpi.h>
#include <vector>



int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();

  auto n_procs = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  auto my_proc = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  // Creating the map array of points to be sent
  std::map<unsigned int, std::vector< Point<2> > > m;
  for (unsigned int i=0; i<n_procs; ++i)
    if (i != my_proc)
      m[i].push_back(Point<2>(my_proc, -2.0*my_proc));

  auto received_pts = dealii::Utilities::MPI::internal::send_and_receive< Point<2> >(MPI_COMM_WORLD, m);

  bool test_passed = true;
  for (auto const &pt: received_pts)
    if ( std::abs(pt.first - pt.second[0][0]) > 1e-12 ||
         std::abs(2.0*pt.first + pt.second[0][1]) > 1e-12 )
      {
        test_passed = false;
        deallog << "Error with point " << pt.second[0] << " received from rank " << pt.first << std::endl;
      }
  if (test_passed)
    deallog << "Test: ok" << std::endl;
  else
    deallog << "Test: FAILED" << std::endl;
}

// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2023 by the deal.II authors
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


// check that handling of ghost elements in parallel distributed vectors works
// appropriately when assigning from ghosted to non-ghosted vectors

#include <deal.II/base/cuda.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/read_write_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "numproc=" << numproc << std::endl;


  // processor 0 and 1 own 2 indices each, higher processors nothing, all are
  // ghosting global elements 1 and 3
  IndexSet local_owned(std::min(numproc * 2, 4U));
  if (myid < 2)
    local_owned.add_range(myid * 2, myid * 2 + 2);
  IndexSet local_relevant(local_owned.size());
  local_relevant = local_owned;
  local_relevant.add_range(1, 2);
  if (numproc > 1)
    local_relevant.add_range(3, 4);

  // run this twice, once where the vectors have called update_ghost_values
  // and once where they have not
  for (unsigned int run = 0; run < 2; ++run)
    {
      LinearAlgebra::distributed::Vector<double, MemorySpace::Default> v(
        local_owned, local_relevant, MPI_COMM_WORLD);

      // set local values
      LinearAlgebra::ReadWriteVector<double> rw_vector(local_owned);
      if (myid < 2)
        {
          rw_vector(myid * 2)     = myid * 2.0;
          rw_vector(myid * 2 + 1) = myid * 2.0 + 1.0;
        }

      v.import_elements(rw_vector, VectorOperation::insert);

      LinearAlgebra::distributed::Vector<double, MemorySpace::Default> w(v),
        u(v);
      u = 0;

      v *= 2.0;
      v.add(1.0);

      if (run == 1)
        {
          v.update_ghost_values();
          w.update_ghost_values();
          u.update_ghost_values();
        }

      rw_vector.import_elements(v, VectorOperation::insert);
      if (myid < 2)
        {
          Assert(rw_vector(myid * 2) == myid * 4.0 + 1, ExcInternalError());
          Assert(rw_vector(myid * 2 + 1) == myid * 4.0 + 3.0,
                 ExcInternalError());
        }

      // copy vector content to non-ghosted vectors, manually created.
      LinearAlgebra::distributed::Vector<double, MemorySpace::Default> v_dist(
        local_owned, MPI_COMM_WORLD),
        w_dist(v_dist), u_dist(v_dist);

      v_dist = v;
      w_dist = w;
      u_dist = u;

      u_dist.add(1.0, v_dist, -1.0, w);

      // copy back to a ghosted vector and update ghost values there
      u = u_dist;
      u.update_ghost_values();

      rw_vector.import_elements(u_dist, VectorOperation::insert);
      if (myid < 2)
        {
          Assert(rw_vector(myid * 2) == myid * 2.0 + 1, ExcInternalError());
          Assert(rw_vector(myid * 2 + 1) == myid * 2.0 + 2.0,
                 ExcInternalError());
        }

      rw_vector.import_elements(u, VectorOperation::insert);
      if (myid < 2)
        {
          Assert(rw_vector(myid * 2) == myid * 2.0 + 1, ExcInternalError());
          Assert(rw_vector(myid * 2 + 1) == myid * 2.0 + 2.0,
                 ExcInternalError());
        }

      IndexSet u_ghost_set(local_owned.size());
      u_ghost_set.add_index(1);
      u_ghost_set.add_index(3);
      LinearAlgebra::ReadWriteVector<double> u_ghost_vector(u_ghost_set);
      u_ghost_vector.import_elements(u, VectorOperation::insert);
      Assert(u_ghost_vector(1) == 2., ExcInternalError());
      if (numproc > 1)
        {
          if (run == 1)
            {
              IndexSet v_ghost_set(local_owned.size());
              v_ghost_set.add_index(3);
              LinearAlgebra::ReadWriteVector<double> v_ghost_vector(
                v_ghost_set);
              v_ghost_vector.import_elements(v, VectorOperation::insert);
              Assert(v_ghost_vector(3) == 7., ExcInternalError());
            }
          Assert(u_ghost_vector(3) == 4., ExcInternalError());
        }

      // check l2 norm
      const double l2_norm = u.l2_norm();
      if (myid == 0 && run == 1)
        deallog << "L2 norm: " << l2_norm << std::endl;
    }

  if (myid == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();
      deallog << std::setprecision(4);

      test();
    }
  else
    test();
}

// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2014 by the deal.II authors
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


// check that handling of ghost elements in parallel distributed vectors works
// appropriately when assigning from ghosted to non-ghosted vectors

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/parallel_vector.h>
#include <fstream>
#include <iostream>
#include <vector>


void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0) deallog << "numproc=" << numproc << std::endl;


  // processor 0 and 1 own 2 indices each, higher processors nothing, all are
  // ghosting global elements 1 and 3
  IndexSet local_owned(std::min(numproc*2, 4U));
  if (myid < 2)
    local_owned.add_range(myid*2,myid*2+2);
  IndexSet local_relevant(local_owned.size());
  local_relevant = local_owned;
  local_relevant.add_range(1,2);
  if (numproc > 1)
    local_relevant.add_range(3,4);

  // run this twice, once where the vectors have called update_ghost_values
  // and once where they have not
  for (unsigned int run = 0; run < 2; ++run)
    {
      parallel::distributed::Vector<double> v(local_owned, local_relevant, MPI_COMM_WORLD);

      // set local values
      if (myid < 2)
        {
          v(myid*2)=myid*2.0;
          v(myid*2+1)=myid*2.0+1.0;
        }

      v.compress(VectorOperation::insert);

      parallel::distributed::Vector<double> w(v), u(v);
      u = 0;

      v*=2.0;
      v.add(1.0);

      if (run == 1)
        {
          v.update_ghost_values();
          w.update_ghost_values();
          u.update_ghost_values();
        }

      if (myid < 2)
        {
          Assert(v(myid*2) == myid*4.0+1, ExcInternalError());
          Assert(v(myid*2+1) == myid*4.0+3.0, ExcInternalError());
        }

      // copy vector content to non-ghosted vectors, manually created.
      parallel::distributed::Vector<double> v_dist(local_owned, MPI_COMM_WORLD),
        w_dist(v_dist), u_dist(v_dist);

      v_dist = v;
      w_dist = w;
      u_dist = u;

      u_dist.add(1.0, v_dist, -1.0, w);

      // copy back to a ghosted vector and update ghost values there
      u = u_dist;
      u.update_ghost_values();

      if (myid < 2)
        {
          Assert(u_dist(myid*2) == myid*2.0+1, ExcInternalError());
          Assert(u_dist(myid*2+1) == myid*2.0+2.0, ExcInternalError());
          Assert(u(myid*2) == myid*2.0+1, ExcInternalError());
          Assert(u(myid*2+1) == myid*2.0+2.0, ExcInternalError());
        }

      Assert(u(1) == 2., ExcInternalError());
      if (numproc > 1)
        {
          if (run == 1)
            Assert(v(3) == 7., ExcInternalError());
          Assert(u(3) == 4., ExcInternalError());
        }

      // check l2 norm
      const double l2_norm = u.l2_norm();
      if (myid == 0 && run == 1)
        deallog << "L2 norm: " << l2_norm << std::endl;
    }

  if (myid == 0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();

}

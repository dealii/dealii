// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2013 by the deal.II authors
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


// check n_ghost_entries() and is_ghost_entry()

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

  const unsigned int set = 200;
  AssertIndexRange (numproc, set-2);
  const unsigned int local_size = set - myid;
  unsigned int global_size = 0;
  unsigned int my_start = 0;
  for (unsigned int i=0; i<numproc; ++i)
    {
      global_size += set - i;
      if (i<myid)
        my_start += set - i;
    }
  // each processor owns some indices and all
  // are ghosting elements from three
  // processors (the second). some entries
  // are right around the border between two
  // processors
  IndexSet local_owned(global_size);
  local_owned.add_range(my_start, my_start + local_size);
  IndexSet local_relevant(global_size);
  local_relevant = local_owned;
  unsigned int ghost_indices [10] = {1, 2, 13, set-2, set-1, set, set+1, 2*set,
                                     2*set+1, 2*set+3
                                    };
  local_relevant.add_indices (&ghost_indices[0], &ghost_indices[0]+10);

  parallel::distributed::Vector<double> v(local_owned, local_relevant, MPI_COMM_WORLD);

  // check number of ghosts everywhere (counted
  // the above)
  if (myid == 0)
    {
      AssertDimension (v.n_ghost_entries(), 5);
    }
  else if (myid == 1)
    {
      AssertDimension (v.n_ghost_entries(), 8);
    }
  else if (myid == 2)
    {
      AssertDimension (v.n_ghost_entries(), 7);
    }
  else
    {
      AssertDimension (v.n_ghost_entries(), 10);
    }

  // count that 13 is ghost only on non-owning
  // processors
  if (myid == 0)
    {
      Assert (v.is_ghost_entry (13) == false, ExcInternalError());
    }
  else
    {
      Assert (v.is_ghost_entry (13) == true, ExcInternalError());
    }

  // count that 27 is ghost nowhere
  Assert (v.is_ghost_entry (27) == false, ExcInternalError());
  if (myid == 0)
    {
      Assert (v.in_local_range (27) == true, ExcInternalError());
    }
  else
    {
      Assert (v.in_local_range (27) == false, ExcInternalError());
    }

  // element with number set is ghost
  if (myid == 1)
    {
      Assert (v.is_ghost_entry (set) == false, ExcInternalError());
    }
  else
    {
      Assert (v.is_ghost_entry (set) == true, ExcInternalError());
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

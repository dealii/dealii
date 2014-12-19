// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2013 by the deal.II authors
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


// let processors write ghost_targets, import_indices, import_targets to file

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/partitioner.h>
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

  Utilities::MPI::Partitioner v(local_owned, local_relevant, MPI_COMM_WORLD);

  // write the info on ghost processors and import indices to file
  {
    std::ofstream file((std::string("dat.") + Utilities::int_to_string(myid)).c_str());
    file << "**** proc " << myid << std::endl;
    file << "ghost targets: ";
    for (unsigned int i=0; i<v.ghost_targets().size(); ++i)
      file << "[" << v.ghost_targets()[i].first << "/"
           << v.ghost_targets()[i].second << "] ";
    file << std::endl;
    file << "import targets: ";
    for (unsigned int i=0; i<v.import_targets().size(); ++i)
      file << "[" << v.import_targets()[i].first << "/"
           << v.import_targets()[i].second << "] ";
    file << std::endl;
    file << "import indices:" << std::endl;
    for (unsigned int i=0; i<v.import_indices().size(); ++i)
      file << "[" << v.import_indices()[i].first << "/"
           << v.import_indices()[i].second << ")" << std::endl;
    file << "****" << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (myid==0)
    {
      for (unsigned int i=0; i<numproc; ++i)
        {
          cat_file((std::string("dat.") + Utilities::int_to_string(i)).c_str());
        }

    }

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

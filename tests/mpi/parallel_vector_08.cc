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


// check parallel_vector::copy_from to update ghost values. Same vector layout
// as in parallel_vector_07.cc

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

  // v has ghosts, w has none. set some entries
  // on w, copy into v and check if they are
  // there
  parallel::distributed::Vector<double> v(local_owned, local_relevant, MPI_COMM_WORLD);
  parallel::distributed::Vector<double> w(local_owned, local_owned, MPI_COMM_WORLD);

  // set a few of the local elements
  for (unsigned i=0; i<local_size; ++i)
    w.local_element(i) = 2.0 * (i + my_start);

  v.copy_from(w);
  v.update_ghost_values();

  // check local values for correctness
  for (unsigned int i=0; i<local_size; ++i)
    Assert (v.local_element(i) == 2.0 * (i + my_start), ExcInternalError());

  // check local values with two different
  // access operators
  for (unsigned int i=0; i<local_size; ++i)
    Assert (v.local_element(i) == v(local_owned.nth_index_in_set (i)), ExcInternalError());
  for (unsigned int i=0; i<local_size; ++i)
    Assert (v.local_element(i) == v(i+my_start), ExcInternalError());

  // check non-local entries on all processors
  for (unsigned int i=0; i<10; ++i)
    Assert (v(ghost_indices[i])== 2. * ghost_indices[i], ExcInternalError());

  // compare direct access local_element with access ()
  for (unsigned int i=0; i<10; ++i)
    if (ghost_indices[i] < my_start)
      Assert (v(ghost_indices[i])==v.local_element(local_size+i), ExcInternalError());

  if (myid == 0)
    for (unsigned int i=5; i<10; ++i)
      Assert (v(ghost_indices[i])==v.local_element(local_size+i-5), ExcInternalError());


  // now the same again, but import ghosts
  // through the call to copy_from
  v.reinit (local_owned, local_relevant, MPI_COMM_WORLD);
  v.copy_from(w, true);

  // check local values for correctness
  for (unsigned int i=0; i<local_size; ++i)
    Assert (v.local_element(i) == 2.0 * (i + my_start), ExcInternalError());

  // check local values with two different
  // access operators
  for (unsigned int i=0; i<local_size; ++i)
    Assert (v.local_element(i) == v(local_owned.nth_index_in_set (i)), ExcInternalError());
  for (unsigned int i=0; i<local_size; ++i)
    Assert (v.local_element(i) == v(i+my_start), ExcInternalError());

  // check non-local entries on all processors
  for (unsigned int i=0; i<10; ++i)
    Assert (v(ghost_indices[i])== 2. * ghost_indices[i], ExcInternalError());

  // compare direct access [] with access ()
  for (unsigned int i=0; i<10; ++i)
    if (ghost_indices[i] < my_start)
      Assert (v(ghost_indices[i])==v.local_element(local_size+i), ExcInternalError());

  if (myid == 0)
    for (unsigned int i=5; i<10; ++i)
      Assert (v(ghost_indices[i])==v.local_element(local_size+i-5), ExcInternalError());


  // now do not call import_ghosts and check
  // whether ghosts really are zero
  v.reinit (local_owned, local_relevant, MPI_COMM_WORLD);
  v.copy_from(w, false);

  // check non-local entries on all processors
  for (unsigned int i=0; i<10; ++i)
    if (local_owned.is_element (ghost_indices[i]) == false)
      Assert (v(ghost_indices[i]) == 0., ExcInternalError());

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

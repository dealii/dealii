// ---------------------------------------------------------------------
//
// Copyright (C) 2014 by the deal.II authors
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



// check correct behavior of sadd of Trilinos vectors
// if they have different Epetra maps 

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/trilinos_vector.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

void test ()
{
  const int n_proc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const int my_id = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  //All processes should own 10 entries
  const int entries_per_process = 10;

  IndexSet locally_owned(entries_per_process*n_proc);
  const int begin_index = my_id*entries_per_process;
  const int end_index = (my_id+1)*entries_per_process;
  locally_owned.add_range(begin_index, end_index);

  IndexSet locally_relevant(entries_per_process*n_proc);
  const int local_begin = std::max(0, begin_index-entries_per_process/2);
  const int local_end = entries_per_process*n_proc;
  locally_relevant.add_range (local_begin, local_end);
 
  TrilinosWrappers::MPI::Vector ghosted, distributed;
  distributed.reinit(locally_owned, MPI_COMM_WORLD);
  ghosted.reinit (locally_owned, locally_relevant, MPI_COMM_WORLD);

  // set the 'distributed' vector to all ones and store its results in
  // 'ghosted'
  distributed=1.;
  ghosted=distributed;

  // then multiply 'distributed' by two and add 'ghosted' to it again
  distributed*=2;
  distributed.add(ghosted, true);

  // assign the result, which should contain all 3s to 'ghosted'
  ghosted = distributed;

  if (my_id==0)
    {
      deallog << "Distributed:" << std::endl;
      for (unsigned int i=begin_index; i<end_index; ++i)
	deallog << i << ": " << distributed(i) << std::endl; 
      
      deallog << "Ghosted:" << std::endl;
      for (unsigned int i=local_begin; i<local_end; ++i)
	deallog << i << ": " << ghosted(i) << std::endl;
    }

  // verify correct value
  for (unsigned int i=begin_index; i<end_index; ++i)
    Assert(distributed(i)==3, ExcInternalError());
      
  for (unsigned int i=local_begin; i<local_end; ++i)
    Assert(ghosted(i)==3, ExcInternalError());
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  
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

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



// check correct behaviour of sadd of Trilinos vectors
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
  TrilinosWrappers::MPI::Vector ghosted, distributed;
  //All processes should own 10 entries
  const int entries_per_process = 10;
  const int n_proc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  IndexSet locally_owned(entries_per_process*n_proc);
  IndexSet locally_relevant(entries_per_process*n_proc);
  const int my_id = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const int begin_index = my_id*entries_per_process;
  const int end_index = (my_id+1)*entries_per_process;
  const int local_begin = std::max(0, begin_index-entries_per_process/2);
  const int local_end = entries_per_process*n_proc;

  locally_owned.add_range(begin_index, end_index);
  locally_relevant.add_range
    (local_begin, local_end);
 
  distributed.reinit(locally_owned, MPI_COMM_WORLD);
  ghosted.reinit (locally_owned, locally_relevant, MPI_COMM_WORLD);

  distributed=1.;
  ghosted=distributed;
  distributed.sadd (2., ghosted);
  ghosted = distributed;
  deallog << "sadd (s, v)" <<std::endl;
  if (my_id==0)
    for (unsigned int i=local_begin; i<local_end; ++i)
      deallog << i << ": " << ghosted(i) << std::endl;

  distributed=1.;
  ghosted=distributed;
  distributed.sadd (2., 3., ghosted);
  ghosted = distributed;
  deallog << "sadd (s, a, v)" <<std::endl;
  if (my_id==0)
    for (unsigned int i=local_begin; i<local_end; ++i)
      deallog << i << ": " << ghosted(i) << std::endl;

  distributed=1.;
  ghosted=distributed;
  distributed.sadd (2., 3., ghosted, 4., ghosted);
  ghosted = distributed;
  deallog << "sadd (s, a, v, b, w)" <<std::endl;
  if (my_id==0)
    for (unsigned int i=local_begin; i<local_end; ++i)
      deallog << i << ": " << ghosted(i) << std::endl;

  distributed=1.;
  ghosted=distributed;
  distributed.sadd (2., 3., ghosted, 4., ghosted, 5., ghosted);
  ghosted = distributed;
  deallog << "sadd (s, a, v, b, w, c, x)" <<std::endl;
  if (my_id==0)
    for (unsigned int i=local_begin; i<local_end; ++i)
      deallog << i << ": " << ghosted(i) << std::endl;
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

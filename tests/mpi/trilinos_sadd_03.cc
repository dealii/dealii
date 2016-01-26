// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2015 by the deal.II authors
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
// if there are processes that doesn't own anything

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
  TrilinosWrappers::MPI::Vector ghosted, ghosted2, distributed, distributed2;
  const unsigned int my_id = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int n_proc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  //All processes except for the last one own elements
  const int entries_per_process = (my_id==n_proc-1)?0:10;
  IndexSet locally_owned(entries_per_process*(n_proc-1));
  IndexSet locally_relevant(entries_per_process*(n_proc-1));
  const int begin_index = my_id*entries_per_process;
  const unsigned int end_index = (my_id+1)*entries_per_process;
  const unsigned int local_begin = std::max(0, begin_index-entries_per_process/2);
  const unsigned int local_end = entries_per_process*(n_proc-1);

  if (my_id != n_proc-1)
    {
      locally_owned.add_range(begin_index, end_index);
      locally_relevant.add_range(local_begin, local_end);
    }

  distributed.reinit(locally_owned, MPI_COMM_WORLD);
  distributed2.reinit(locally_owned, MPI_COMM_WORLD);
  ghosted.reinit (locally_owned, locally_relevant, MPI_COMM_WORLD);
  ghosted2.reinit (locally_owned, locally_relevant, MPI_COMM_WORLD);

  //sadd(s,v) with ghosted vector
  distributed=1.;
  ghosted=distributed;
  ghosted *= 1.6;
  distributed.sadd (2., ghosted);
  ghosted = distributed;
  deallog << "sadd (s, v) ghosted" <<std::endl;
  if (my_id==0)
    for (unsigned int i=local_begin; i<local_end; ++i)
      deallog << i << ": " << ghosted(i) << std::endl;

  //sadd(s,v) with distributed vector
  distributed=1.;
  distributed.sadd (1.5, distributed);
  ghosted = distributed;
  deallog << "sadd (s, v) distributed" <<std::endl;
  if (my_id==0)
    for (unsigned int i=local_begin; i<local_end; ++i)
      deallog << i << ": " << ghosted(i) << std::endl;

  //sadd(s,a,v) with ghosted vector
  distributed=1.;
  ghosted=distributed;
  ghosted *= 1.6;
  distributed.sadd (2., 3., ghosted);
  ghosted = distributed;
  deallog << "sadd (s, a, v) ghosted" <<std::endl;
  if (my_id==0)
    for (unsigned int i=local_begin; i<local_end; ++i)
      deallog << i << ": " << ghosted(i) << std::endl;

  //sadd(s,a,v) with distributed vector
  distributed=1.;
  distributed.sadd (2., 3., distributed);
  ghosted = distributed;
  deallog << "sadd (s, a, v) distributed" <<std::endl;
  if (my_id==0)
    for (unsigned int i=local_begin; i<local_end; ++i)
      deallog << i << ": " << ghosted(i) << std::endl;

  //sadd(s,a,v,b,w) with ghosted vector
  distributed=1.;
  ghosted=distributed;
  ghosted *= 1.6;
  ghosted2 = distributed;
  ghosted2 *= 0.354;
  distributed.sadd (2., 3., ghosted, 4., ghosted2);
  ghosted = distributed;
  deallog << "sadd (s, a, v, b, w) ghosted" <<std::endl;
  if (my_id==0)
    for (unsigned int i=local_begin; i<local_end; ++i)
      deallog << i << ": " << ghosted(i) << std::endl;

  //sadd(s,a,v,b,w) with distributed vector
  distributed=1.;
  distributed2=1.6;
  distributed.sadd (2., 3., distributed, 4., distributed2);
  ghosted = distributed;
  deallog << "sadd (s, a, v, b, w) distributed" <<std::endl;
  if (my_id==0)
    for (unsigned int i=local_begin; i<local_end; ++i)
      deallog << i << ": " << ghosted(i) << std::endl;

  //sadd(s,a,v,b,w,c,x) with ghosted vector
  distributed=1.;
  ghosted=distributed;
  distributed2 = 7.22;
  ghosted2=distributed2;
  distributed.sadd (2., 3., ghosted, 4., ghosted2, 5., distributed2);
  ghosted = distributed;
  deallog << "sadd (s, a, v, b, w, c, x) ghosted" <<std::endl;
  if (my_id==0)
    for (unsigned int i=local_begin; i<local_end; ++i)
      deallog << i << ": " << ghosted(i) << std::endl;

  //sadd(s,a,v,b,w,c,x) with distributed vector
  distributed=1.;
  distributed2 = 7.22;
  distributed.sadd (2., 3., distributed, 4., distributed, 5., distributed2);
  ghosted = distributed;
  deallog << "sadd (s, a, v, b, w, c, x) distributed" <<std::endl;
  if (my_id==0)
    for (unsigned int i=local_begin; i<local_end; ++i)
      deallog << i << ": " << ghosted(i) << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, testing_max_num_threads());

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();

}

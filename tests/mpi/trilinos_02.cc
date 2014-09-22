// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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



// check Trilinos has_ghost_elements()

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/trilinos_vector.h>
#include <fstream>
#include <iostream>
#include <vector>

void check(TrilinosWrappers::MPI::Vector &v, bool ghost)
{
  Assert(v.has_ghost_elements()==ghost, ExcMessage("wrong ghost elements"));
}

void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0) deallog << "numproc=" << numproc << std::endl;

  // each processor owns 2 indices and all
  // are ghosting element 1 (the second)
  IndexSet local_active(numproc*2);
  local_active.add_range(myid*2,myid*2+2);
  IndexSet local_relevant(numproc*2);
  local_relevant = local_active;
  local_relevant.add_range(1,2);

  TrilinosWrappers::MPI::Vector v(local_active, MPI_COMM_WORLD);
  check(v, false);

  TrilinosWrappers::MPI::Vector v2(local_relevant, MPI_COMM_WORLD);
  check(v2, true);

  TrilinosWrappers::MPI::Vector v3=v;
  check(v3, false);

  v3.reinit(v2);
  check(v3, true);

  v3.reinit(local_active);
  check(v3, false);

  v3.reinit(local_relevant);
  check(v3, true);

  TrilinosWrappers::MPI::Vector v4=v2;
  check(v4, true);

  v4 = v;
  check(v4, true); //this only copies contents!

  v4 = v2;
  check(v4, true);

  TrilinosWrappers::MPI::Vector v5(v2);
  check(v5, true);

  if (myid==0)
    deallog << "OK" << std::endl;
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

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



// creation and size of LA::MPI::Vector

#include "../tests.h"
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/base/index_set.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

#include "gla.h"

template <class LA>
void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0)
    deallog << "numproc=" << numproc << std::endl;

  // each processor owns 2 indices and all
  // are ghosting Element 1 (the second)

  IndexSet local_active(numproc*2);
  local_active.add_range(myid*2,myid*2+2);
  IndexSet local_relevant(numproc*2);
  local_relevant.add_range(1,2);


  IndexSet something(100);
  something.add_range(myid,myid+1);
  if (myid==numproc-1)
    something.add_range(numproc,100);


  {
    typename LA::MPI::Vector v1;

    v1.reinit(something, MPI_COMM_WORLD);
    Assert(!v1.has_ghost_elements(), ExcInternalError());
    Assert(v1.size()==100, ExcInternalError());

    typename LA::MPI::Vector v2(local_active, local_relevant, MPI_COMM_WORLD);
    Assert(v2.has_ghost_elements(), ExcInternalError());
    Assert(v2.size()==numproc*2, ExcInternalError());

    v2.reinit(local_active,MPI_COMM_WORLD);
    Assert(!v2.has_ghost_elements(), ExcInternalError());
    Assert(v2.size()==numproc*2, ExcInternalError());

    v2.reinit(local_active, local_relevant, MPI_COMM_WORLD);
    Assert(v2.has_ghost_elements(), ExcInternalError());
    Assert(v2.size()==numproc*2, ExcInternalError());

  }

  // done
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

      {
        deallog.push("PETSc");
        test<LA_PETSc>();
        deallog.pop();
        deallog.push("Trilinos");
        test<LA_Trilinos>();
        deallog.pop();
      }

    }
  else
    {
      deallog.push("PETSc");
      test<LA_PETSc>();
      deallog.pop();
      deallog.push("Trilinos");
      test<LA_Trilinos>();
      deallog.pop();
    }

  if (myid==9999)
    test<LA_Dummy>();


}

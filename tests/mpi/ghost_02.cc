// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2015 by the deal.II authors
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



// Test correct copying  of ghosted vectors in PETScWrappers::mpi::vectors

#include "../tests.h"
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/base/index_set.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>


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

  PETScWrappers::MPI::Vector vb(local_active, MPI_COMM_WORLD);
  PETScWrappers::MPI::Vector v(local_active, local_relevant, MPI_COMM_WORLD);
  PETScWrappers::MPI::Vector v2(local_active, local_relevant, MPI_COMM_WORLD);

  vb = 1.5;
  v2 = vb;
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "ghost: " << get_real_assert_zero_imag(v2(1)) << std::endl;
  Assert(get_real_assert_zero_imag(v2(1)) == 1.5, ExcInternalError());
  Assert(get_real_assert_zero_imag(v2(myid*2)) == 1.5, ExcInternalError());
  Assert(get_real_assert_zero_imag(v2(myid*2+1)) == 1.5, ExcInternalError());

  // set local values
  vb(myid*2)=myid*2.0;
  vb(myid*2+1)=myid*2.0+1.0;

  vb.compress(VectorOperation::insert);
  vb*=2.0;
  v=vb;


  // check local values
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      deallog << myid*2 << ":" << get_real_assert_zero_imag(v(myid*2)) << std::endl;
      deallog << myid*2+1 << ":" << get_real_assert_zero_imag(v(myid*2+1)) << std::endl;
    }

  Assert(get_real_assert_zero_imag(v(myid*2)) == myid*4.0, ExcInternalError());
  Assert(get_real_assert_zero_imag(v(myid*2+1)) == myid*4.0+2.0, ExcInternalError());


  // check ghost values
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "ghost: " << get_real_assert_zero_imag(v(1)) << std::endl;
  Assert(get_real_assert_zero_imag(v(1)) == 2.0, ExcInternalError());

  //assignment from ghosted to ghosted
  v2 = v;
  Assert(get_real_assert_zero_imag(v2(1)) == 2.0, ExcInternalError());
  Assert(get_real_assert_zero_imag(v2(myid*2)) == myid*4.0, ExcInternalError());
  Assert(get_real_assert_zero_imag(v2(myid*2+1)) == myid*4.0+2.0, ExcInternalError());
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "ghost: " << get_real_assert_zero_imag(v2(1)) << std::endl;

  // done
  if (myid==0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
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

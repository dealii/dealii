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


// build a vector whose elements exceed the size of unsigned int in case of 64
// bit indices. To avoid excessive memory consumption, let the vector start at
// a number close to the maximum of unsigned int but extend past the last
// index

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <fstream>
#include <iostream>
#include <vector>


void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0) deallog << "numproc=" << numproc << std::endl;

  types::global_dof_index min_index = 0xffffffffU - 39;
  types::global_dof_index local_size = 42;
  IndexSet local_owned(min_index+numproc*local_size);
  local_owned.add_range(min_index+myid*local_size, min_index+(myid+1)*local_size);

  // all processors ghost some entries around invalid_unsigned_int and on the
  // border between two processors
  IndexSet local_relevant(local_owned.size());
  local_relevant = local_owned;
  local_relevant.add_range(min_index+38,min_index+40);
  local_relevant.add_range(min_index+41,min_index+43);

  LinearAlgebra::distributed::Vector<double> v(local_owned, local_relevant, MPI_COMM_WORLD);

  deallog << "Local range of proc 0: " << v.local_range().first << " "
          << v.local_range().second << std::endl;

  // set local values
  for (types::global_dof_index i=min_index+myid*local_size;
       i<min_index+(myid+1)*local_size; ++i)
    v(i) = (double)i;

  deallog << "vector norm: " << v.l2_norm() << std::endl;

  // check ghost values
  deallog << "v(42) = " << v(min_index+42) << std::endl;
  v.update_ghost_values();
  deallog << "v(42) = " << v(min_index+42) << std::endl;
  Assert(v(min_index+41) == min_index+41, ExcInternalError());
  Assert(v(min_index+39) == min_index+39, ExcInternalError());
  Assert(v(min_index+38) == min_index+38, ExcInternalError());

  v.zero_out_ghosts();
  v(min_index+38) = min_index;
  v(min_index+39) = min_index*2;
  v(min_index+41) = min_index+7;
  v(min_index+42) = -static_cast<double>(min_index);
  v.compress(VectorOperation::add);
  v.update_ghost_values();
  deallog << "v(38) = " << v(min_index+38) << std::endl;
  deallog << "v(39) = " << v(min_index+39) << std::endl;
  deallog << "v(41) = " << v(min_index+41) << std::endl;
  deallog << "v(42) = " << v(min_index+42) << std::endl;

  deallog << "OK" << std::endl;
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
      deallog << std::setprecision(12);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();

}

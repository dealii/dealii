// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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



// document bug in the new ConstraintMatrix::distribute() from r29593
// if one vector owns all DoFs, then distribute() hangs.

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/constraint_matrix.h>

#include <fstream>
#include <sstream>



void test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  const unsigned int n_processes = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  // create a vector that consists of elements indexed from 0 to 100
  // everything belongs to proc 0
  IndexSet idx(100);
  if (myid==0)
    idx.add_range(0,100);

  TrilinosWrappers::MPI::Vector vec (idx, MPI_COMM_WORLD);

  for (unsigned int i=vec.local_range().first; i<vec.local_range().second; ++i)
    vec(i) = i;
  vec.compress(VectorOperation::insert);

  ConstraintMatrix cm;
  cm.close ();

  cm.distribute (vec);

  if (myid == 0)
    deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();

}

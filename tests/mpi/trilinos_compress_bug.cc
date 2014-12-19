// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2014 by the deal.II authors
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



// document a bug in Trilinos regarding compress()
// this fails in 10.6.4, but works in 10.4.2!

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/constraint_matrix.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/trilinos_vector.h>


void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0) deallog << "numproc=" << numproc << std::endl;

  IndexSet locally_owned(21);
  if (myid==0)
    locally_owned.add_range(0,9);
  else
    locally_owned.add_range(9,21);
//
  TrilinosWrappers::MPI::Vector test1(locally_owned);
  if (myid==0)
    {
      types::global_dof_index idx[]= {0,1,2,3,4,5,6,7,8};
      double val[]= {0,1,2,3,4,5,6,7,8};
      test1.add(9,idx,val);
    }
  else
    {
      {
        types::global_dof_index idx[]= {1,9,3,10,5,11,12,13,14};
        double val[]= {1,9,3,10,5,11,12,13,14};
        test1.add(9,idx,val);
      }
    }

  test1.compress(VectorOperation::add);

  //TrilinosWrappers::MPI::Vector test(test1.vector_partitioner()); // works
  //TrilinosWrappers::MPI::Vector test(locally_owned); // works
  TrilinosWrappers::MPI::Vector test(test1); // fails

  test = 0;

  if (myid==0)
    test(locally_owned.nth_index_in_set(5))=7;

  if (myid==0)
    deallog << "before compress: " << test(locally_owned.nth_index_in_set(5)) << std::endl;

  test.compress(VectorOperation::insert);

  if (myid==0)
    deallog << "after compress: " << test(locally_owned.nth_index_in_set(5)) << std::endl;

  // Trilinos produces a 0 instead of a 7 here. Why?
  if (myid==0)
    {
      Assert(test(locally_owned.nth_index_in_set(5)) == 7, ExcInternalError());
    }

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
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

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



// Test if Add() and Insert work correctly in Trilinos in parallel
// it turns out Trilinos always adds

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/parallel_vector.h>
#include <fstream>
#include <iostream>
#include <vector>


void test ()
{
  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const unsigned int my_id = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  const unsigned int n_rows = 3;
  const unsigned int n_cols = 4;

  IndexSet row_partitioning (n_rows);
  IndexSet col_partitioning (n_cols);

  if (n_procs == 1)
    {
      row_partitioning.add_range(0, n_rows);
      col_partitioning.add_range(0, n_cols);
    }
  else if (n_procs == 2)
    {
      // row_partitioning should be { [0, 2), [2, n_rows) }
      // col_partitioning should be { [0, 2), [2, n_cols) }
      // col_relevant_set should be { [0, 3), [1, n_cols) }
      if (my_id == 0)
        {
          row_partitioning.add_range(0, 2);
          col_partitioning.add_range(0, 2);
        }
      else if (my_id == 1)
        {
          row_partitioning.add_range(2, n_rows);
          col_partitioning.add_range(2, n_cols);
        }
    }
  else
    Assert (false, ExcNotImplemented());

  TrilinosWrappers::SparsityPattern sp (row_partitioning,
                                        col_partitioning, MPI_COMM_WORLD);
  if (my_id == 0)
    {
      sp.add (0, 0);
      sp.add (0, 2);
    }
  if ((n_procs == 1) || (my_id == 1))
    sp.add(2,3);
  sp.compress();

  TrilinosWrappers::SparseMatrix A;
  A.reinit (sp);
  if (my_id==0)
    {
      A.set(0, 0, 0.1);
      A.set(0, 2, 0.2);
    }
  if ((n_procs == 1) || (my_id == 1))
    {
      A.set(0, 0, 0.1);
      A.set(0, 2, 0.2);
      A.set(2,3, 0.3);
    }

  A.compress(VectorOperation::insert);

  if (my_id == 0)
    {
      //A.add (0, 0, 1.33);
      //A.add (0, 2, 1.0);
    }
  if ((n_procs == 1) || (my_id == 1))
    A.add(0,0, 1.67);
  //A.begin()->value() += 1.67; // crashes

  A.compress(VectorOperation::add);

  if (my_id==0)
    {
      deallog << "A(0,0)=" << A(0,0) << std::endl;
    }


  if (my_id == 0) deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, numbers::invalid_unsigned_int);

  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
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

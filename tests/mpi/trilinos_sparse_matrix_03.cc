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



// test TrilinosWrappers::SparseMatrix::iterator semantics

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
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

  // now access elements by iterator. we know what should be in the
  // matrix, just make sure
  if (my_id == 0)
    {
      for (TrilinosWrappers::SparseMatrix::iterator p=A.begin(0);
           p != A.end(0);
           ++p)
        {
          AssertThrow (p->row() == 0, ExcInternalError());
          if (p->column() == 0)
            AssertThrow (p->value() == 0.1, ExcInternalError())
            else if (p->column() == 2)
              AssertThrow (p->value() == 0.2, ExcInternalError())
              else
                // well, we didn't write here, so the only thing that
                // should be in there is a zero
                AssertThrow (p->value() == 0.0, ExcInternalError());
        }
    }
  else
    {
      for (TrilinosWrappers::SparseMatrix::iterator p=A.begin(2);
           p != A.end(2);
           ++p)
        {
          AssertThrow (p->row() == 2, ExcInternalError());
          if (p->column() == 3)
            AssertThrow (p->value() == 0.3, ExcInternalError())
            else
              // well, we didn't write here, so the only thing that
              // should be in there is a zero
              AssertThrow (p->value() == 0.0, ExcInternalError());
        }
    }

  // now let each of the processors write something into the memory of
  // the other. the values we can locally access are of course the
  // same ones as before because we haven't communicated the values
  // yet:
  if (my_id == 0)
    {
      A.set (2, 3, 42.);
      for (TrilinosWrappers::SparseMatrix::iterator p=A.begin(0);
           p != A.end(0);
           ++p)
        {
          AssertThrow (p->row() == 0, ExcInternalError());
          if (p->column() == 0)
            AssertThrow (p->value() == 0.1, ExcInternalError())
            else if (p->column() == 2)
              AssertThrow (p->value() == 0.2, ExcInternalError())
              else
                // well, we didn't write here, so the only thing that
                // should be in there is a zero
                AssertThrow (p->value() == 0.0, ExcInternalError());
        }
    }
  else
    {
      A.set (0, 0, 108.);
      for (TrilinosWrappers::SparseMatrix::iterator p=A.begin(2);
           p != A.end(2);
           ++p)
        {
          AssertThrow (p->row() == 2, ExcInternalError());
          if (p->column() == 3)
            AssertThrow (p->value() == 0.3, ExcInternalError())
            else
              // well, we didn't write here, so the only thing that
              // should be in there is a zero
              AssertThrow (p->value() == 0.0, ExcInternalError());
        }
    }

  // then call compress() and ensure that we get the correct values
  A.compress (VectorOperation::insert);
  if (my_id == 0)
    {
      for (TrilinosWrappers::SparseMatrix::iterator p=A.begin(0);
           p != A.end(0);
           ++p)
        {
          AssertThrow (p->row() == 0, ExcInternalError());
          if (p->column() == 0)
            AssertThrow (p->value() == 108, ExcInternalError())
            else if (p->column() == 2)
              AssertThrow (p->value() == 0.2, ExcInternalError())
              else
                // well, we didn't write here, so the only thing that
                // should be in there is a zero
                AssertThrow (p->value() == 0.0, ExcInternalError());
        }
    }
  else
    {
      for (TrilinosWrappers::SparseMatrix::iterator p=A.begin(2);
           p != A.end(2);
           ++p)
        {
          AssertThrow (p->row() == 2, ExcInternalError());
          if (p->column() == 3)
            AssertThrow (p->value() == 42, ExcInternalError())
            else
              // well, we didn't write here, so the only thing that
              // should be in there is a zero
              AssertThrow (p->value() == 0.0, ExcInternalError());
        }
    }

  if (my_id == 0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, testing_max_num_threads());

  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
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

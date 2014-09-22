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



// TrilinosWrappers::SparseMatrix::print got column indices wrong

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <fstream>
#include <iostream>
#include <vector>


void test ()
{
  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const unsigned int my_id = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  const unsigned int n_rows = 2;
  const unsigned int n_cols = 2;

  IndexSet row_partitioning (n_rows);
  IndexSet col_partitioning (n_cols);

  if (n_procs == 1)
    {
      row_partitioning.add_range(0, n_rows);
      col_partitioning.add_range(0, n_cols);
    }
  else if (n_procs == 2)
    {
      if (my_id == 0)
        {
          row_partitioning.add_range(0, 1);
          col_partitioning.add_range(0, 1);
        }
      else if (my_id == 1)
        {
          row_partitioning.add_range(1, n_rows);
          col_partitioning.add_range(1, n_cols);
        }
    }
  else
    Assert (false, ExcNotImplemented());

  /* A is

     0     1
     0     1
  */
  const unsigned int n_entries = 2;
  const unsigned int line [n_entries] = {0, 1};
  const unsigned int local_index [n_entries] = {1, 1};
  const double local_value [n_entries] = {1.0, 1.0};

  TrilinosWrappers::SparsityPattern sp (row_partitioning, col_partitioning, MPI_COMM_WORLD);
  for (unsigned int i = 0; i < n_entries; ++i)
    if (row_partitioning.is_element(line[i]))
      sp.add(line[i], local_index[i]);
  sp.compress();

  TrilinosWrappers::SparseMatrix A;
  A.clear ();
  A.reinit (sp);
  for (unsigned int i = 0; i<n_entries; ++i)
    if (row_partitioning.is_element(line[i]))
      A.add(line[i], local_index[i], local_value[i]);
  A.compress(VectorOperation::add);

  if (my_id == 0)
    {
      Assert(A.el(0, 0) == 0, ExcMessage("Wrong element in A!"));
      Assert(A.el(0, 1) == 1, ExcMessage("Wrong element in A!"));
    }
  if ((n_procs == 1) || (my_id == 1))
    {
      Assert(A.el(1, 0) == 0, ExcMessage("Wrong element in A!"));
      Assert(A.el(1, 1) == 1, ExcMessage("Wrong element in A!"));
    }

  TrilinosWrappers::SparseMatrix AtA;
  A.Tmmult (AtA, A);

  /* AtA should be

     0     0
     0     2
  */

  // checking AtA row partitioning
  if (n_procs == 2)
    {
      if (my_id == 0)
        {
          Assert(AtA.local_range().first == 0,
                 ExcMessage("AtA Local Range is not as expected."));
          Assert(AtA.local_range().second == 1,
                 ExcMessage("AtA Local Range is not as expected."));
        }
      if (my_id == 1)
        {
          Assert(AtA.local_range().first == 1,
                 ExcMessage("AtA Local Range is not as expected."));
          Assert(AtA.local_range().second == 2,
                 ExcMessage("AtA Local Range is not as expected."));
        }
    }

  // checking AtA elements. note that
  // the el() function either returns
  // the correct value, or zero in
  // case the element is stored on
  // other processors. consequently,
  // in the following we only test
  // that the values that *must* be
  // zero in AtA are indeed zero, but
  // not that the others have the
  // correct value
  Assert(AtA.el(0, 0) == 0, ExcMessage("Wrong element in AtA!"));
  Assert(AtA.el(0, 1) == 0, ExcMessage("Wrong element in AtA!"));
  Assert(AtA.el(1, 0) == 0, ExcMessage("Wrong element in AtA!"));

  // now also check the one nonzero
  // element
  if ((n_procs == 1) || (my_id == 1))
    Assert(AtA.el(1, 1) == 2, ExcMessage("Wrong element in AtA!"));

  deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  try
    {
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
  catch (const char *p)
    {
      std::cerr << "Uncaught exception: " << p << std::endl;
      std::exit (1);
    }
}

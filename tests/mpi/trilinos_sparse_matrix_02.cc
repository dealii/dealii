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



// This demonstrates a bug where copying a Trilinos matrix and then modifying
// the source will edit both matrices.

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

  const unsigned int n_rows = 4;

  IndexSet locally_owned (n_rows);

  if (n_procs == 1)
    {
      locally_owned.add_range(0, n_rows);
    }
  else if (n_procs == 2)
    {
      // should be { [0, 2), [2, n_rows) }
      if (my_id == 0)
        locally_owned.add_range(0, 2);
      else if (my_id == 1)
        locally_owned.add_range(2, n_rows);
    }
  else
    Assert (false, ExcNotImplemented());

  TrilinosWrappers::SparsityPattern sp (locally_owned, locally_owned, MPI_COMM_WORLD);
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
  TrilinosWrappers::SparseMatrix B;
  B.reinit (A);

  A.add(0, 0, 0.1);
  A.add(0, 2, 0.2);
  if ((n_procs == 1) || (my_id == 1))
    A.add(2,3, 0.3);

  double l1a = (n_procs==1) ? 0.3 : 0.4;
  double l1b = n_procs*1.2;

  A.compress(VectorOperation::add);
  deallog << "1: " << A.l1_norm() << " " << B.l1_norm()
          << " (should be " << l1a << " 0.0)" << std::endl;

  deallog << "set B=A..." << std::endl;

  B.copy_from(A);

  deallog << "2: " << A.l1_norm() << " " << B.l1_norm()
          << " (should be " << l1a << " " << l1a << ")" << std::endl;

  if (my_id==0)
    {
      deallog << "A(0,0)=" << A(0,0) << std::endl;
      deallog << "B(0,0)=" << B(0,0) << std::endl;
    }

  deallog << "reassemble A..." << std::endl;

  A = 0;
  A.add(0, 0, -1.2);
  A.compress(VectorOperation::add);
  deallog << "3: " << A.l1_norm() << " " << B.l1_norm()
          << " (should be " << l1b << " " << l1a << ")" << std::endl;

  if (my_id==0)
    {
      deallog << "A(0,0)=" << A(0,0) << std::endl;
      deallog << "B(0,0)=" << B(0,0) << std::endl;
    }

  if (my_id == 0) deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
  MPILogInitAll log;
  test();
}

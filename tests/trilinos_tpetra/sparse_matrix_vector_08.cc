// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check SparseMatrix::vmult, vmult_add with deal.II vector

#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test(Vector<double> &v, Vector<double> &w)
{
  LinearAlgebra::TpetraWrappers::SparseMatrix<double, MemorySpace::Default> m(
    w.size(), v.size(), v.size());
  for (unsigned int i = 0; i < m.m(); ++i)
    for (unsigned int j = 0; j < m.n(); ++j)
      m.set(i, j, i + 2 * j);

  for (unsigned int i = 0; i < v.size(); ++i)
    v(i) = i;

  m.compress(VectorOperation::insert);

  // w:=Mv
  m.vmult(w, v);

  // make sure we get the expected result
  for (unsigned int i = 0; i < m.m(); ++i)
    {
      double result = 0;
      for (unsigned int j = 0; j < m.n(); ++j)
        result += (i + 2 * j) * j;
      AssertThrow(w(i) == result, ExcInternalError());
    }

  m.vmult_add(w, v);
  // make sure we get the expected result
  for (unsigned int i = 0; i < m.m(); ++i)
    {
      double result = 0;
      for (unsigned int j = 0; j < m.n(); ++j)
        result += (i + 2 * j) * j;
      AssertThrow(w(i) == result + result, ExcInternalError());
    }

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());


  try
    {
      {
        Vector<double> v(100);
        Vector<double> w(95);
        test(v, w);
      }
    }
  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}

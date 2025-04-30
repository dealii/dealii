// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// check SparseMatrix::vmult, vmult_add with distributed deal.II vector (but
// without using distributed things)

#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test(LinearAlgebra::distributed::Vector<double> &v,
     LinearAlgebra::distributed::Vector<double> &w)
{
  LinearAlgebra::TpetraWrappers::SparseMatrix<double> m(w.size(),
                                                        v.size(),
                                                        v.size());
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
        LinearAlgebra::distributed::Vector<double> v(100);
        LinearAlgebra::distributed::Vector<double> w(95);
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

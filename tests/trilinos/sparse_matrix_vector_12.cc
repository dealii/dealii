// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2019 by the deal.II authors
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



// check SparseMatrix::vmult, vmult_add with
// LinearAlgebra::EpetraWrappers::Vector

#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test(LinearAlgebra::EpetraWrappers::Vector &v,
     LinearAlgebra::EpetraWrappers::Vector &w)
{
  LinearAlgebra::ReadWriteVector<double> read_write_v(v.size());
  LinearAlgebra::ReadWriteVector<double> read_write_w(w.size());

  TrilinosWrappers::SparseMatrix m(w.size(), v.size(), v.size());
  for (unsigned int i = 0; i < m.m(); ++i)
    for (unsigned int j = 0; j < m.n(); ++j)
      m.set(i, j, i + 2 * j);

  for (unsigned int i = 0; i < v.size(); ++i)
    read_write_v(i) = i;

  m.compress(VectorOperation::insert);
  v.import(read_write_v, VectorOperation::insert);

  // w:=Mv
  m.vmult(w, v);

  // make sure we get the expected result
  read_write_w.import(w, VectorOperation::insert);
  for (unsigned int i = 0; i < m.m(); ++i)
    {
      double result = 0;
      for (unsigned int j = 0; j < m.n(); ++j)
        result += (i + 2 * j) * j;
      AssertThrow(read_write_w(i) == result, ExcInternalError());
    }

  m.vmult_add(w, v);

  // make sure we get the expected result
  read_write_w.import(w, VectorOperation::insert);
  for (unsigned int i = 0; i < m.m(); ++i)
    {
      double result = 0;
      for (unsigned int j = 0; j < m.n(); ++j)
        result += (i + 2 * j) * j;
      AssertThrow(read_write_w(i) == result + result, ExcInternalError());
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
        LinearAlgebra::EpetraWrappers::Vector v(complete_index_set(100),
                                                MPI_COMM_SELF);
        LinearAlgebra::EpetraWrappers::Vector w(complete_index_set(95),
                                                MPI_COMM_SELF);
        test(v, w);
      }
    }
  catch (std::exception &exc)
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

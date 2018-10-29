// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2017 by the deal.II authors
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



// check SparseMatrix::matrix_norm_square

#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test(TrilinosWrappers::MPI::Vector &v)
{
  TrilinosWrappers::SparseMatrix m(v.size(), v.size(), v.size());
  for (unsigned int i = 0; i < m.m(); ++i)
    for (unsigned int j = 0; j < m.m(); ++j)
      m.set(i, j, i + 2 * j);

  for (unsigned int i = 0; i < v.size(); ++i)
    v(i) = i;

  m.compress(VectorOperation::insert);
  v.compress(VectorOperation::insert);

  // <w,Mv>
  const TrilinosScalar s = m.matrix_norm_square(v);

  // make sure we get the expected result
  for (unsigned int i = 0; i < v.size(); ++i)
    AssertThrow(v(i) == i, ExcInternalError());

  TrilinosScalar result = 0;
  for (unsigned int i = 0; i < m.m(); ++i)
    for (unsigned int j = 0; j < m.m(); ++j)
      result += (i + 2 * j) * j * i;

  AssertThrow(s == result, ExcInternalError());

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
        TrilinosWrappers::MPI::Vector v;
        v.reinit(complete_index_set(30), MPI_COMM_WORLD);
        test(v);
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

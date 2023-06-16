// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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

// Check copy_to() and copy_from() between tensors and full matrices.

#include <deal.II/base/tensor.h>

#include <deal.II/lac/full_matrix.h>

#include "../tests.h"

int
main()
{
  initlog();

  const unsigned int matrix_dimension = 4;

  Tensor<2, matrix_dimension> myTensor;

  FullMatrix<double> myMatrix(matrix_dimension, matrix_dimension);

  for (unsigned int i = 0; i < matrix_dimension; i++)
    for (unsigned int j = 0; j < matrix_dimension; j++)
      {
        myMatrix(i, j) = i + j;
      }

  myMatrix.copy_to(myTensor);

  deallog << matrix_dimension << std::endl;
  deallog.get_file_stream() << myTensor << std::endl;

  myTensor *= 2;
  myMatrix.copy_from(myTensor);
  for (unsigned int i = 0; i < matrix_dimension; i++)
    for (unsigned int j = 0; j < matrix_dimension; j++)
      {
        Assert(myMatrix(i, j) == 2. * (i + j), ExcInternalError());
      }
}

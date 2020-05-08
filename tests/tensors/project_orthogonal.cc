// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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

// Test the function projection_onto_orthogonal_tensors

#include <deal.II/base/tensor.h>

#include <deal.II/lac/full_matrix.h>

#include "../tests.h"

template <int dim, typename number>
void
print_type(const Tensor<2, dim, number> &)
{
  if (std::is_same<number, float>::value)
    deallog << " Tensor<2, " << dim << ", float>" << std::endl;
  else if (std::is_same<number, double>::value)
    deallog << " Tensor<2, " << dim << ", double>" << std::endl;
  else
    Assert(false, ExcNotImplemented());
}

template <int dim, typename number>
void
test(const Tensor<2, dim, number> &A)
{
  print_type(A);
  deallog << "original tensor A = " << A << std::endl;
  deallog << "det A = " << determinant(A) << std::endl;
  deallog << "A A^T = " << (A * transpose(A)) << std::endl;
  const Tensor<2, dim, number> Q = project_onto_orthogonal_tensors(A);
  deallog << "projected to Q = " << Q << std::endl;
  deallog << "det Q = " << determinant(Q) << std::endl;
  deallog << "Q Q^T = " << (Q * transpose(Q)) << std::endl;
  deallog << std::endl;
}

int
main()
{
  initlog();
  // not orthogonal
  test(Tensor<2, 3, double>{{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}});
  // already orthogonal
  test(Tensor<2, 3, double>{{{1, 0, 0}, {0, 0, 1}, {0, 1, 0}}});
  // not orthogonal, but det = 1
  test(Tensor<2, 3, double>{{{1, 2, 0}, {0, 1, 0}, {0, 0, 1}}});
  // 2D not orthogonal
  test(Tensor<2, 2, double>{{{1, 2}, {3, 4}}});
  // 2D already orthogonal
  test(Tensor<2, 2, double>{{{0, 1}, {1, 0}}});
  // 1D not orthogonal
  test(Tensor<2, 1, double>{{{2.5}}});
  // 1D already orthogonal
  test(Tensor<2, 1, double>{{{1.}}});

  deallog << "OK" << std::endl;
}

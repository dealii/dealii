// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

// create and manipulate constexpr SymmetricTensor objects

#include <deal.II/base/symmetric_tensor.h>

#include "../tests.h"

DEAL_II_CONSTEXPR SymmetricTensor<2, 2>
                  get_tensor_2()
{
  SymmetricTensor<2, 2> A;
  A[0][0] = 1.;
  A[1][1] = 3.;
  A[0][1] = -5.;
  return A;
}


DEAL_II_CONSTEXPR SymmetricTensor<4, 2>
                  get_tensor_4()
{
  SymmetricTensor<4, 2> B;
  B[0][0][0][0] = 1.;
  B[1][1][1][1] = 2.5;
  B[0][1][0][1] = 0.2;
  return B;
}


int
main()
{
  initlog();

  DEAL_II_CONSTEXPR const auto A = get_tensor_2();
  {
    LogStream::Prefix              p("SymmetricTensor<2,2>");
    DEAL_II_CONSTEXPR const double calculation  = A[0][0] * A[1][0] - A[1][1];
    DEAL_II_CONSTEXPR const double invariants[] = {first_invariant(A),
                                                   second_invariant(A),
                                                   third_invariant(A)};
    DEAL_II_CONSTEXPR const double norm_square  = scalar_product(A, A);
    deallog << "calculation result = " << calculation << std::endl;
    deallog << "invariants = " << invariants[0] << ", " << invariants[1] << ", "
            << invariants[2] << std::endl;
    deallog << "norm square = " << norm_square << std::endl;
  }

  DEAL_II_CONSTEXPR const auto B = get_tensor_4();
  {
    LogStream::Prefix              p("SymmetricTensor<4,2>");
    DEAL_II_CONSTEXPR const double calculation =
      B[0][0][0][0] * B[1][0][0][1] - B[1][1][0][0];
    DEAL_II_CONSTEXPR const auto B_times_A = B * A;
    deallog << "calculation result = " << calculation << std::endl;
    deallog << "B times A (0,0) = " << B_times_A[0][0] << std::endl;
  }
  deallog << "OK" << std::endl;
  return 0;
}
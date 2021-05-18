//-----------------------------------------------------------
//
//    Copyright (C) 2020 - 2021 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//-----------------------------------------------------------

// Test SUNDIALS' clones on our own vectors.

#include <deal.II/lac/vector.h>

#include <deal.II/sundials/n_vector.h>

#include "../tests.h"

using namespace SUNDIALS::internal;

int
main()
{
  initlog();

  const unsigned int N = 4;
  using VectorType     = Vector<double>;

  VectorType y_dealii(N);
  y_dealii = 1.0;

  // The commented version works, the other one segfaults.
  // auto y_sundials = make_nvector_view(y_dealii);
  N_Vector y_sundials     = make_nvector_view(y_dealii);
  N_Vector scale_sundials = N_VClone(y_sundials);

  // scale_sundials = 2*y_sundials
  N_VScale(2.0, y_sundials, scale_sundials);

  const auto *scale_dealii = unwrap_nvector<VectorType>(scale_sundials);

  deallog << "Deal.II Norm: " << y_dealii.l2_norm() << std::endl;
  deallog << "Twice deal.II Norm: " << scale_dealii->l2_norm() << std::endl;

  return 0;
}
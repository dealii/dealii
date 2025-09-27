// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test the invariants of tensors using the Cayley-Hamilton theorem,
// see http://en.wikipedia.org/wiki/Invariants_of_tensors

#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"


template <int dim>
void
check()
{
  for (unsigned int round = 0; round < 10; ++round)
    {
      SymmetricTensor<2, dim> S;
      for (unsigned int i = 0; i < S.n_independent_components; ++i)
        S[S.unrolled_to_component_indices(i)] = Testing::rand() % 10;

      deallog << "S = " << S << std::endl;
      deallog << "first invariant  = " << first_invariant(S) << std::endl;
      deallog << "second invariant = " << second_invariant(S) << std::endl;
      deallog << "third invariant  = " << third_invariant(S) << std::endl;

      Tensor<2, dim> S_cubed;
      for (unsigned int d = 0; d < dim; ++d)
        for (unsigned int e = 0; e < dim; ++e)
          for (unsigned int f = 0; f < dim; ++f)
            for (unsigned int g = 0; g < dim; ++g)
              S_cubed[d][e] += S[d][f] * S[f][g] * S[g][e];

      Tensor<2, dim> S_squared;
      for (unsigned int d = 0; d < dim; ++d)
        for (unsigned int e = 0; e < dim; ++e)
          for (unsigned int f = 0; f < dim; ++f)
            S_squared[d][e] += S[d][f] * S[f][e];

      Tensor<2, dim> R = S_cubed - first_invariant(S) * S_squared +
                         second_invariant(S) * S -
                         third_invariant(S) * unit_symmetric_tensor<dim>();
      deallog << R << std::endl;

      AssertThrow(R.norm() < 1e-10, ExcInternalError());
    }
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  check<3>();

  deallog << "OK" << std::endl;
}

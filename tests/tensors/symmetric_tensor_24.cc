// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check SymmetricTensor<2,dim>::component_to_unrolled_index and the
// other way round

#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"


template <int dim>
void
check()
{
  using S = SymmetricTensor<2, dim>;
  for (unsigned int i = 0; i < S::n_independent_components; ++i)
    {
      deallog << i << "  --  " << S::unrolled_to_component_indices(i)
              << std::endl;
      AssertThrow(S::component_to_unrolled_index(
                    S::unrolled_to_component_indices(i)) == i,
                  ExcInternalError());
    }
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  check<1>();
  check<2>();
  check<3>();
  check<4>();
  check<5>();
  check<6>();

  deallog << "OK" << std::endl;
}

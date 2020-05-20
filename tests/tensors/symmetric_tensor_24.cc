// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
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


// check SymmetricTensor<2,dim>::component_to_unrolled_index and the
// other way round

#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"


template <int dim>
void
check()
{
  typedef SymmetricTensor<2, dim> S;
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

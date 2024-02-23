// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test SymmetricTensor::norm() for complex-valued tensors

#include <deal.II/base/symmetric_tensor.h>

#include "../tests.h"


template <int rank, int dim>
void
check()
{
  // build a regular tensor
  SymmetricTensor<rank, dim> t;

  // build one in which all numbers are the same but purely imaginary
  SymmetricTensor<rank, dim, std::complex<double>> ti;

  // build one in which all numbers have both real and imaginary components
  SymmetricTensor<rank, dim, std::complex<double>> tc;

  for (unsigned int i = 0; i < t.n_independent_components; ++i)
    {
      t.access_raw_entry(i)  = 1.0 * (i + 1);
      ti.access_raw_entry(i) = std::complex<double>(0, 1.0 * (i + 1));
      tc.access_raw_entry(i) =
        std::complex<double>(1.0 * (i + 1), 1.0 * (i + 1));
    }

  deallog << t.norm() << ' ' << ti.norm() << ' ' << tc.norm() << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  deallog << "check rank 2 tensors" << std::endl;
  check<2, 1>();
  check<2, 2>();
  check<2, 3>();

  deallog << "check rank 4 tensors" << std::endl;
  check<4, 1>();
  check<4, 2>();
  check<4, 3>();

  deallog << "OK" << std::endl;
}

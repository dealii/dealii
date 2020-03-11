// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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

// Check SymmetricTensor::norm() for complex-valued rank-2 tensors by
// computing it via multiplication

#include <deal.II/base/symmetric_tensor.h>

#include "../tests.h"


template <int rank, int dim>
void
check()
{
  // build a regular tensor
  SymmetricTensor<rank, dim> t, t_conj;

  // build one in which all numbers are the same but purely imaginary
  SymmetricTensor<rank, dim, std::complex<double>> ti, ti_conj;

  // build one in which all numbers have both real and imaginary components
  SymmetricTensor<rank, dim, std::complex<double>> tc, tc_conj;

  for (unsigned int i = 0; i < t.n_independent_components; ++i)
    {
      t.access_raw_entry(i)  = 1.0 * (i + 1);
      ti.access_raw_entry(i) = std::complex<double>(0, 1.0 * (i + 1));
      tc.access_raw_entry(i) =
        std::complex<double>(1.0 * (i + 1), 1.0 * (i + 1));

      t_conj.access_raw_entry(i)  = 1.0 * (i + 1);
      ti_conj.access_raw_entry(i) = std::complex<double>(0, -1.0 * (i + 1));
      tc_conj.access_raw_entry(i) =
        std::complex<double>(1.0 * (i + 1), -1.0 * (i + 1));
    }

  deallog << t.norm() << " vs " << std::sqrt(t_conj * t) << std::endl
          << ti.norm() << " vs " << std::sqrt(ti_conj * ti) << std::endl
          << tc.norm() << " vs " << std::sqrt(tc_conj * tc) << std::endl
          << std::endl;
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

  deallog << "OK" << std::endl;
}

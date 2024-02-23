// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check Tensor<0,dim>

#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"

template <typename U, typename V>
void
compare(const U &u, const V &v)
{
  AssertThrow(static_cast<double>(u) == static_cast<double>(v),
              ExcInternalError());
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);

  using T = Tensor<0, 1>;
  T t1(13.), t2(42);

  compare(T(), 0.);
  compare(T(13.), 13.);
  compare(T(t1), 13.);
  compare(static_cast<double>(t1), 13.);
  compare(static_cast<double &>(t1), 13.);
  compare((T() = t1), 13.);
  {
    T t;
    t = 13.;
    compare(t, 13.);
  }
  compare((t1 == t1), true);
  compare((t1 == t2), false);
  compare((t1 != t2), true);
  compare((t1 != t1), false);
  compare((T() += t1), t1);
  compare((T() -= t1), -t1);
  compare(T(13.) *= 3., 39.);
  compare(T(39) /= 3., 13.);
  compare((t1 * t2), 13 * 42);
  compare((t1 + t2), 13 + 42);
  compare((t1 - t2), 13 - 42);
  compare(-t1, -13.);
  compare(T(-12).norm(), 12.);
  compare(T(-12).norm_square(), 12 * 12.);

  t1.clear();
  compare(t1, 0.);

  deallog << "OK" << std::endl;
}

// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2021 by the deal.II authors
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
  compare((T() = 13.), 13.);
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

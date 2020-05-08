// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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


#include <deal.II/base/numbers.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/point.h>

#include <boost/core/demangle.hpp>

#include <memory>

#include "../tests.h"

using namespace Patterns::Tools;

// Try conversion on non elementary types

template <class T>
void
test(T t)
{
  auto p = Convert<T>::to_pattern();
  deallog << "Pattern  : " << p->description() << std::endl;
  auto s = Convert<T>::to_string(t);
  deallog << "To String: " << s << std::endl;
  deallog << "To value : " << Convert<T>::to_string(Convert<T>::to_value(s))
          << std::endl;
}

int
main()
{
  initlog();

  Point<1> t0(1);
  Point<2> t1(2, 3);
  Point<3> t2(3, 4, 5);

  Tensor<1, 1> t3;
  Tensor<1, 2> t4;
  Tensor<1, 3> t5;

  Tensor<2, 1> t32;
  Tensor<2, 2> t42;
  Tensor<2, 3> t52;

  Tensor<3, 1> t33;
  Tensor<3, 2> t43;
  Tensor<3, 3> t53;

  std::complex<double> t6(1, 2);
  auto                 t7 = std::make_pair(1, 2.0);

  test(t0);
  test(t1);
  test(t2);
  test(t3);
  test(t4);
  test(t5);
  test(t6);
  test(t7);
  test(t32);
  test(t42);
  test(t52);
  test(t33);
  test(t43);
  test(t53);

  return 0;
}

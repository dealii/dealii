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

// Try conversion on complex map types

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

  Point<3>                          t0(1, 2, 3);
  std::complex<double>              t1(4, 5);
  std::vector<Point<3>>             t2(2, t0);
  std::vector<std::complex<double>> t3(2, t1);

  std::map<unsigned int, Point<3>>                          t10;
  std::map<unsigned int, std::complex<double>>              t11;
  std::map<unsigned int, std::vector<Point<3>>>             t12;
  std::map<unsigned int, std::vector<std::complex<double>>> t13;

  t10[0] = t0;
  t11[0] = t1;
  t12[0] = t2;
  t13[0] = t3;

  t10[2] = t0;
  t11[2] = t1;
  t12[2] = t2;
  t13[2] = t3;

  test(t10);
  test(t11);
  test(t12);
  test(t13);

  return 0;
}

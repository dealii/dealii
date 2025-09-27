// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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

  std::map<int, double>                t0;
  std::unordered_map<int, double>      t1;
  std::multimap<int, double>           t2;
  std::unordered_multimap<int, double> t3;

  auto p  = std::make_pair(5, 1.0);
  auto p2 = std::make_pair(5, 2.0);
  auto p3 = std::make_pair(1, 3.0);

  t0.insert(p);
  t1.insert(p);
  t2.insert(p);
  t3.insert(p);

  t0.insert(p3);
  t1.insert(p3);
  t2.insert(p3);
  t3.insert(p3);

  t2.insert(p2);
  t3.insert(p2);

  test(t0);
  test(t1);
  test(t2);
  test(t3);

  return 0;
}

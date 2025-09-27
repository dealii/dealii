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

// Try conversion on arbitrary container types

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

  std::vector<unsigned int>             t00;
  std::deque<unsigned int>              t01;
  std::list<unsigned int>               t02;
  std::set<unsigned int>                t03;
  std::multiset<unsigned int>           t04;
  std::unordered_set<unsigned int>      t05;
  std::unordered_multiset<unsigned int> t06;

  t00.insert(t00.end(), 0);
  t01.insert(t01.end(), 1);
  t02.insert(t02.end(), 2);
  t03.insert(t03.end(), 3);
  t04.insert(t04.end(), 4);
  t05.insert(t05.end(), 5);
  t06.insert(t06.end(), 6);


  t00.insert(t00.end(), 1);
  t01.insert(t01.end(), 2);
  t02.insert(t02.end(), 3);
  t03.insert(t03.end(), 4);
  t04.insert(t04.end(), 5);
  t05.insert(t05.end(), 6);
  t06.insert(t06.end(), 7);

  test(t00);
  test(t01);
  test(t02);
  test(t03);
  test(t04);
  test(t05);
  test(t06);

  return 0;
}

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

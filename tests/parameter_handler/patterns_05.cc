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


#include <deal.II/base/patterns.h>
#include <deal.II/base/point.h>

#include <memory>

#include "../tests.h"

using namespace Patterns::Tools;

// Try conversion on elementary types

template <class T>
void
test(T t, std::string s)
{
  auto p = Convert<T>::to_pattern();
  deallog << "Pattern  : " << p->description() << std::endl;
  deallog << "String   : " << s << std::endl;
  deallog << "To string: " << Convert<T>::to_string(t) << std::endl;
  deallog << "To value : " << Convert<T>::to_string(Convert<T>::to_value(s))
          << std::endl;
}

int
main()
{
  initlog();


  std::vector<int>          l0(2, -1);
  std::vector<unsigned int> l1(3, 2);
  std::vector<double>       l2(4, 3.14);
  std::vector<std::string>  l3(2, "bar");
  Point<3>                  l4(3, 2, 1);
  std::complex<double>      l5(2.0, 1.0);

  std::vector<std::vector<int>>          vl0(2, std::vector<int>(3, 1));
  std::vector<std::vector<unsigned int>> vl1(2,
                                             std::vector<unsigned int>(3, 2));
  std::vector<std::vector<double>>       vl2(2, std::vector<double>(3, 3.14));
  std::vector<std::vector<std::string>>  vl3(2,
                                            std::vector<std::string>(3, "foo"));
  std::vector<Point<3>>                  vl4(2, Point<3>(4, 3, 2));
  std::vector<std::complex<double>>      vl5(2, std::complex<double>(1.0, 3.0));

  deallog << "List of int         : " << Convert<decltype(l0)>::to_string(l0)
          << std::endl;
  deallog << "List of unsigned int: " << Convert<decltype(l1)>::to_string(l1)
          << std::endl;
  deallog << "List of double      : " << Convert<decltype(l2)>::to_string(l2)
          << std::endl;
  deallog << "List of string      : " << Convert<decltype(l3)>::to_string(l3)
          << std::endl;
  deallog << "Point<3>            : " << Convert<decltype(l4)>::to_string(l4)
          << std::endl;
  deallog << "std::complex        : " << Convert<decltype(l5)>::to_string(l5)
          << std::endl;

  deallog << "List of lists of int         : "
          << Convert<decltype(vl0)>::to_string(vl0) << std::endl;
  deallog << "List of lists of unsigned int: "
          << Convert<decltype(vl1)>::to_string(vl1) << std::endl;
  deallog << "List of lists of double      : "
          << Convert<decltype(vl2)>::to_string(vl2) << std::endl;
  deallog << "List of lists of string      : "
          << Convert<decltype(vl3)>::to_string(vl3) << std::endl;
  deallog << "List of Point<3>             : "
          << Convert<decltype(vl4)>::to_string(vl4) << std::endl;
  deallog << "List of std::complex         : "
          << Convert<decltype(vl5)>::to_string(vl5) << std::endl;

  deallog << "=============================" << std::endl;

  l0 = Convert<decltype(l0)>::to_value("1,2,3");
  l1 = Convert<decltype(l1)>::to_value("3,4,5");
  l2 = Convert<decltype(l2)>::to_value("5,6,7");
  l3 = Convert<decltype(l3)>::to_value("8,9,8.5");
  l4 = Convert<decltype(l4)>::to_value("8,9,8.5");
  l5 = Convert<decltype(l5)>::to_value("2.0,9.0");

  vl0 = Convert<decltype(vl0)>::to_value("1,2,3  ; 1,2,3  ");
  vl1 = Convert<decltype(vl1)>::to_value("3,4,5  ; 3,4,5  ");
  vl2 = Convert<decltype(vl2)>::to_value("5,6,7  ; 5,6,7  ");
  vl3 = Convert<decltype(vl3)>::to_value("8,9,8.5; 8,9,8.5");
  vl4 = Convert<decltype(vl4)>::to_value("8,9,8.5; 8,9,8.5");
  vl5 = Convert<decltype(vl5)>::to_value("1.0,2.0; 6.0,7.0");

  deallog << "List of int         : " << Convert<decltype(l0)>::to_string(l0)
          << std::endl;
  deallog << "List of unsigned int: " << Convert<decltype(l1)>::to_string(l1)
          << std::endl;
  deallog << "List of double      : " << Convert<decltype(l2)>::to_string(l2)
          << std::endl;
  deallog << "List of string      : " << Convert<decltype(l3)>::to_string(l3)
          << std::endl;
  deallog << "Point<3>            : " << Convert<decltype(l4)>::to_string(l4)
          << std::endl;
  deallog << "std::complex        : " << Convert<decltype(l5)>::to_string(l5)
          << std::endl;

  deallog << "List of lists of int         : "
          << Convert<decltype(vl0)>::to_string(vl0) << std::endl;
  deallog << "List of lists of unsigned int: "
          << Convert<decltype(vl1)>::to_string(vl1) << std::endl;
  deallog << "List of lists of double      : "
          << Convert<decltype(vl2)>::to_string(vl2) << std::endl;
  deallog << "List of lists of string      : "
          << Convert<decltype(vl3)>::to_string(vl3) << std::endl;
  deallog << "List of Point<3>             : "
          << Convert<decltype(vl4)>::to_string(vl4) << std::endl;
  deallog << "List of std::complex         : "
          << Convert<decltype(vl5)>::to_string(vl5) << std::endl;

  return 0;
}

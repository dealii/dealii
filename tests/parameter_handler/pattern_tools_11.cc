// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check Convert<...>::to_string() for small numbers, and make sure
// we don't break anything for integral types

#include <deal.II/base/parameter_handler.h>

#include <memory>
#include <tuple>

#include "../tests.h"

int
main()
{
  initlog();

  double        a  = 1e-12;
  int           i  = -1;
  unsigned int  j  = 3;
  unsigned char c  = 3;
  signed char   c1 = -3;
  char          c2 = -3;
  bool          b  = false;

  auto t = std::make_tuple(a, i, j, c, c1, c2, b);

  auto s  = Patterns::Tools::Convert<decltype(t)>::to_string(t);
  auto tt = Patterns::Tools::Convert<decltype(t)>::to_value(s);
  auto r  = Patterns::Tools::Convert<decltype(t)>::to_string(tt);

  deallog << "String: " << s << ", roundtrip: " << r << std::endl;
}

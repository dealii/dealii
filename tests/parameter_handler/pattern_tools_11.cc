// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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

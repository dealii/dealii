// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test expression only constructor with vector based function

#include <deal.II/base/function_parser.h>
#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>

#include <map>

#include "../tests.h"


double
eval(const std::string &exp, unsigned int comp)
{
  FunctionParser<2> fp(exp);

  Point<2> p;

  return fp.value(p, comp);
}


int
main()
{
  initlog();

  for (unsigned int i = 0; i < 2; ++i)
    {
      double random = eval("rand(); rand()", i);
      if (0.0 <= random && random <= 1.0)
        deallog << "OK" << std::endl;
    }
}

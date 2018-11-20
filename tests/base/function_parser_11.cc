// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

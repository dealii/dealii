// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test that ParameterAcceptor::initialize() works as expected


#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/point.h>

#include "../tests.h"

// Test subsectioning

class Test : public ParameterAcceptor
{
public:
  Test()
  {
    add_parameter("A point", a_point);
  };

private:
  Point<3> a_point;
};

void
test_ext(const std::string &ext)
{
  Test a;

  deallog << "Generate and read input." << ext << std::endl;
  try
    {
      ParameterAcceptor::initialize("input." + ext);
    }
  catch (...)
    {
      // The above call must have created a file named input.ext
      cat_file(("input." + ext).c_str());
    }
}

int
main()
{
  initlog();
  test_ext("prm");
  test_ext("xml");
  test_ext("json");
}

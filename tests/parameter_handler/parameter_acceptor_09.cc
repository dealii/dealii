//-----------------------------------------------------------
//
//    Copyright (C) 2020 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//-----------------------------------------------------------

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

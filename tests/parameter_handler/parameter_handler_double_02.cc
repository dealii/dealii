// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2017 by the deal.II authors
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



// test ParameterHandler::Double::create()

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"

void
test(const std::string &desc)
{
  deallog << desc << " -> ";
  std::unique_ptr<dealii::Patterns::Double> c = Patterns::Double::create(desc);
  if (!c)
    {
      deallog << "NULL" << std::endl;
      return;
    }
  deallog << c->description() << std::endl;
}


int
main()
{
  initlog();

  // invalid:
  test("invalid");
  test("[Double invalid]");
  test("[Double");

  test(Patterns::Double().description());      // no limit
  test(Patterns::Double(-2.13).description()); // lower limit
  test(Patterns::Double(Patterns::Double::min_double_value, 42.0)
         .description());                          // upper limit
  test(Patterns::Double(0.2, 42.0).description()); // both limits
  test(Patterns::Double(1.0, -1.0).description()); // no limits

  return 0;
}

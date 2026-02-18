// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2003 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



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

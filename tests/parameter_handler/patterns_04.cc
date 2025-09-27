// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// pattern_factory tests

#include <deal.II/base/parameter_handler.h>

#include <memory>

#include "../tests.h"


void
test(const Patterns::PatternBase &p)
{
  // Print description
  auto desc = p.description();
  deallog << p.description() << std::endl;

  // Check that we can clone
  auto p2 = p.clone();
  deallog << p2->description() << std::endl;

  auto p3 = Patterns::pattern_factory(desc);
  deallog << p3->description() << std::endl;

  AssertThrow(desc == p2->description(),
              ExcInternalError("Cloned pattern does not have the same "
                               "description as the original pattern!"));

  AssertThrow(desc == p3->description(),
              ExcInternalError("Pattern created using factory "
                               "does not have the same "
                               "description as the original pattern!"));
}


using namespace Patterns;

int
main()
{
  initlog();
  test(Integer());
  test(Integer(-1));
  test(Integer(-1, 50));

  test(Double());
  test(Double(0.0, 50.0));

  test(Bool());

  test(List(Integer()));
  test(List(Double()));
  test(List(Bool()));

  test(Selection("alpha|beta"));
  test(List(Double(), 0, 3, ";"));

  test(Map(Integer(), Double()));
  test(Map(Integer(0, 10), Double(), 0, 10, ";", "="));

  deallog << "OK" << std::endl;
}

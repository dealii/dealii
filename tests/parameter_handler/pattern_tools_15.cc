// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check Patterns::Tools::Convert for custom enum types
// Test the conversion of a standard enum to string and back, and of an enum
// class to string and back

#include <deal.II/fe/fe_values.h>

#include <memory>

#include "../tests.h"

enum TestFlag
{
  red   = 1,
  green = 2,
  blue  = 4
};

enum class TestFlagClass
{
  yellow = 1,
  orange = 2,
  white  = 4
};

template <typename T>
using C = Patterns::Tools::Convert<T>;

int
main()
{
  initlog();
  {
    TestFlag flags = red;
    deallog << C<TestFlag>::to_string(flags) << std::endl;
    flags = C<TestFlag>::to_value("green|blue");
    deallog << C<TestFlag>::to_string(flags) << std::endl;
  }
  {
    TestFlagClass flags = TestFlagClass::yellow;
    deallog << C<TestFlagClass>::to_string(flags) << std::endl;
    flags = C<TestFlagClass>::to_value("orange|white");
    deallog << C<TestFlagClass>::to_string(flags) << std::endl;
  }
}

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


// check comparators in numbers namespace

#include <deal.II/base/numbers.h>

#include "../tests.h"


template <typename Number>
void
check(const Number &a, const Number &b)
{
  AssertThrow((numbers::values_are_equal(a, b)) == (a == b),
              ExcMessage("Comparator is incorrect."));
  AssertThrow((numbers::values_are_not_equal(a, b)) == (a != b),
              ExcMessage("Comparator is incorrect."));
  AssertThrow((numbers::value_is_zero(a)) == (a == 0.0),
              ExcMessage("Comparator is incorrect."));
  AssertThrow((numbers::value_is_zero(b)) == (b == 0.0),
              ExcMessage("Comparator is incorrect."));
  AssertThrow((numbers::value_is_less_than(a, b)) == (a < b),
              ExcMessage("Comparator is incorrect."));
  AssertThrow((numbers::value_is_less_than_or_equal_to(a, b)) == (a <= b),
              ExcMessage("Comparator is incorrect."));
  AssertThrow((numbers::value_is_greater_than(a, b)) == (a > b),
              ExcMessage("Comparator is incorrect."));
  AssertThrow((numbers::value_is_greater_than_or_equal_to(a, b)) == (a >= b),
              ExcMessage("Comparator is incorrect."));
}



int
main()
{
  initlog();

  check(0.0, 2.5);
  check(2.5, 1.5);
  check(1.5, 1.5);

  deallog << "OK" << std::endl;

  return 0;
}

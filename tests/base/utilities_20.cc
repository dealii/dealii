// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test functions in namespace Utilities

#include <deal.II/base/utilities.h>

#include "../tests.h"


void
test()
{
  double number        = 1.23456789e5;
  double number_1digit = Utilities::truncate_to_n_digits(number, 1);
  double number_2digit = Utilities::truncate_to_n_digits(number, 2);
  double number_3digit = Utilities::truncate_to_n_digits(number, 3);

  deallog << std::endl;
  deallog << "number = " << std::scientific << number << std::endl;
  deallog << "number_1digit = " << std::scientific << number_1digit
          << std::endl;
  deallog << "number_2digit = " << std::scientific << number_2digit
          << std::endl;
  deallog << "number_3digit = " << std::scientific << number_3digit
          << std::endl;

  number        = 0.999999999;
  number_1digit = Utilities::truncate_to_n_digits(number, 1);
  number_2digit = Utilities::truncate_to_n_digits(number, 2);
  number_3digit = Utilities::truncate_to_n_digits(number, 3);

  deallog << std::endl;
  deallog << "number = " << std::scientific << number << std::endl;
  deallog << "number_1digit = " << std::scientific << number_1digit
          << std::endl;
  deallog << "number_2digit = " << std::scientific << number_2digit
          << std::endl;
  deallog << "number_3digit = " << std::scientific << number_3digit
          << std::endl;

  number        = 0.0;
  number_1digit = Utilities::truncate_to_n_digits(number, 1);
  number_2digit = Utilities::truncate_to_n_digits(number, 2);
  number_3digit = Utilities::truncate_to_n_digits(number, 3);

  deallog << std::endl;
  deallog << "number = " << std::scientific << number << std::endl;
  deallog << "number_1digit = " << std::scientific << number_1digit
          << std::endl;
  deallog << "number_2digit = " << std::scientific << number_2digit
          << std::endl;
  deallog << "number_3digit = " << std::scientific << number_3digit
          << std::endl;

  float number_float        = -9.87654321e-6;
  float number_1digit_float = Utilities::truncate_to_n_digits(number_float, 1);
  float number_2digit_float = Utilities::truncate_to_n_digits(number_float, 2);
  float number_3digit_float = Utilities::truncate_to_n_digits(number_float, 3);

  deallog << std::endl;
  deallog << "number = " << std::scientific << number_float << std::endl;
  deallog << "number_1digit = " << std::scientific << number_1digit_float
          << std::endl;
  deallog << "number_2digit = " << std::scientific << number_2digit_float
          << std::endl;
  deallog << "number_3digit = " << std::scientific << number_3digit_float
          << std::endl;
}



int
main()
{
  initlog();

  test();
}

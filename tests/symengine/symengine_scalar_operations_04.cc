// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// Check that the std::complex overload for the division and multiplication
// operator works.

#include <deal.II/differentiation/sd.h>

#include "../tests.h"

namespace SD = Differentiation::SD;

template <typename Number1, typename Number2>
void
test()
{
  using SD_number_t = SD::Expression;

  const Number1 a = 1.0;
  const Number2 b = 2.0;

  const SD_number_t symb_a = SD::make_symbol("a");
  const SD_number_t symb_b = SD::make_symbol("b");

  {
    const SD_number_t symb_a_division_b = symb_a / b;
    deallog << "symbolic a/b: " << symb_a_division_b << std::endl;

    const SD_number_t a_division_symb_b = a / symb_b;
    ;
    deallog << "a/symbolic b: " << a_division_symb_b << std::endl;

    const SD_number_t substitute_symb_a_division_b =
      SD::substitute(symb_a_division_b, symb_a, a);
    deallog << "a/b: " << substitute_symb_a_division_b << std::endl;

    const SD_number_t a_division_substitute_symb_b =
      SD::substitute(a_division_symb_b, symb_b, b);
    deallog << "a/b: " << a_division_substitute_symb_b << std::endl;
  }

  {
    const SD_number_t symb_a_product_b = symb_a * b;
    deallog << "symbolic a*b: " << symb_a_product_b << std::endl;

    const SD_number_t a_product_symb_b = a * symb_b;
    ;
    deallog << "a*symbolic b: " << a_product_symb_b << std::endl;

    const SD_number_t substitute_symb_a_product_b =
      SD::substitute(symb_a_product_b, symb_a, a);
    deallog << "a*b: " << substitute_symb_a_product_b << std::endl;

    const SD_number_t a_product_substitute_symb_b =
      SD::substitute(a_product_symb_b, symb_b, b);
    deallog << "a*b: " << a_product_substitute_symb_b << std::endl;
  }
}

int
main()
{
  initlog();
  test<std::complex<float>, std::complex<double>>();
  test<std::complex<double>, std::complex<float>>();
  test<std::complex<double>, std::complex<double>>();
  test<std::complex<float>, std::complex<float>>();
  deallog << "OK" << std::endl;
}

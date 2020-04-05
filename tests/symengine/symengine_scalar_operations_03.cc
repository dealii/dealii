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


// Check that symbol substitution works for scalars

#include <deal.II/differentiation/sd.h>

#include "../tests.h"

namespace SD = Differentiation::SD;


int
main()
{
  initlog();

  using SD_number_t = SD::Expression;

  const double a = 1.0;
  const double b = 2.0;

  const SD_number_t symb_a = SD::make_symbol("a");
  const SD_number_t symb_b = SD::make_symbol("b");

  deallog.push("Substitution: Individual");
  {
    const SD_number_t symb_a_plus_b = symb_a + symb_b;
    deallog << "Symbolic a+b: " << symb_a_plus_b << std::endl;

    const SD_number_t symb_a_plus_b_1 =
      SD::substitute(symb_a_plus_b, symb_a, a);
    deallog << "Symbolic a+b (a=1): " << symb_a_plus_b_1 << std::endl;

    const SD_number_t symb_a_plus_b_subs =
      SD::substitute(symb_a_plus_b_1, symb_b, b);
    deallog << "Symbolic a+b (a=1, b=2): " << symb_a_plus_b_subs << std::endl;
  }
  deallog.pop();

  deallog.push("Substitution: Vector");
  {
    const SD_number_t symb_a_plus_b = symb_a + symb_b;

    const std::vector<std::pair<SD_number_t, double>> symbol_values{
      std::make_pair(symb_a, a), std::make_pair(symb_b, b)};

    const SD_number_t symb_a_plus_b_subs =
      SD::substitute(symb_a_plus_b, symbol_values);
    deallog << "Symbolic a+b (a=1, b=2): " << symb_a_plus_b_subs << std::endl;
  }
  deallog.pop();

  deallog.push("Substitution: Map");
  {
    const SD_number_t symb_a_plus_b = symb_a + symb_b;

    const SD::types::substitution_map substitution_map =
      SD::make_substitution_map(std::make_pair(symb_a, a),
                                std::make_pair(symb_b, b));

    const SD_number_t symb_a_plus_b_subs =
      SD::substitute(symb_a_plus_b, substitution_map);
    deallog << "Symbolic a+b (a=1, b=2): " << symb_a_plus_b_subs << std::endl;
  }
  deallog.pop();

  deallog.push("Substitution with evaluation");
  {
    const SD_number_t symb_a_plus_b = symb_a + symb_b;

    const SD::types::substitution_map substitution_map =
      SD::make_substitution_map(std::make_pair(symb_a, a),
                                std::make_pair(symb_b, b));

    const double symb_a_plus_b_subs =
      SD::substitute_and_evaluate<double>(symb_a_plus_b, substitution_map);
    Assert(symb_a_plus_b_subs == (a + b), ExcInternalError());
  }
  deallog.pop();

  deallog << "OK" << std::endl;
}

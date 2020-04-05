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


// Check that symbols and symbolic functions can be created using the
// utility functions

#include <deal.II/differentiation/sd.h>

#include "../tests.h"

namespace SD = Differentiation::SD;


int
main()
{
  initlog();

  using SD_number_t = SD::Expression;

  const SD::types::substitution_map sub_map{
    {SD_number_t("c"), SD_number_t(1.0)},
    {SD_number_t("b"), SD_number_t(2)},
    {SD_number_t("a"), SD_number_t(3.0f)}};

  const SD::types::symbol_vector symbols{SD_number_t("c"),
                                         SD_number_t("b"),
                                         SD_number_t("a")};

  const SD_number_t symb        = SD::make_symbol("d");
  const SD_number_t symb_func_1 = SD::make_symbolic_function("f", symbols);
  const SD_number_t symb_func_2 = SD::make_symbolic_function("g", sub_map);

  deallog << "Symbol: " << symb << std::endl;
  deallog << "Symbolic function (vector of arguments): " << symb_func_1
          << std::endl;
  deallog << "Symbolic function (map of arguments): " << symb_func_2
          << std::endl;

  deallog << "OK" << std::endl;
}

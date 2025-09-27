// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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

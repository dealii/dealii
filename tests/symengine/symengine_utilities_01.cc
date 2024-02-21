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


// Check that substitution map decomposition and printing work as expected

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

  const SD::types::symbol_vector symbols =
    SD::Utilities::extract_symbols(sub_map);
  const std::vector<double> values =
    SD::Utilities::extract_values<double>(sub_map);
  Assert(values.size() == symbols.size(),
         ExcDimensionMismatch(values.size(), symbols.size()));

  // Print the map itself (this should be dictionary ordered)
  deallog << "Print substitution map" << std::endl;
  SD::Utilities::print_substitution_map(deallog, sub_map);

  // Print the extracted symbol-value pairs (the ordering should match that of
  // the map)
  deallog << "Print extracted symbol-value pairs" << std::endl;
  for (unsigned int i = 0; i < symbols.size(); ++i)
    deallog << symbols[i] << " = " << values[i] << std::endl;

  deallog << "OK" << std::endl;
}

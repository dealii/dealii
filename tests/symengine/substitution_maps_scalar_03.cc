// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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


// Check that functions to create substitution maps from scalar variables work
// correctly.

#include <deal.II/differentiation/sd.h>

#include <symengine/basic.h>
#include <symengine/symengine_rcp.h>

#include <complex>
#include <iostream>
#include <string>

#include "../tests.h"

namespace SD = Differentiation::SD;
namespace SE = ::SymEngine;

int
main()
{
  initlog();


  const SD::types::substitution_map substitution_map_1 =
    SD::make_substitution_map(SD::Expression("x1"), SD::Expression(1));
  const SD::types::substitution_map substitution_map_2 =
    SD::make_substitution_map(SD::Expression("x3"), 3.0);
  const SD::types::substitution_map substitution_map_3 =
    SD::make_substitution_map(SD::types::symbol_vector{SD::Expression("x6"),
                                                       SD::Expression("x7")},
                              std::vector<int>{6, 7});
  const SD::types::substitution_map substitution_map_4{
    {SD::Expression("x8"), SD::Expression(8.0)},
    {SD::Expression("x9"), SD::Expression(9.0f)},
    {SD::Expression("x10"), SD::Expression(10u)}};
  const SD::types::substitution_map substitution_map_5 =
    SD::make_substitution_map(std::make_pair(SD::Expression("x11"), 11));
  const SD::types::substitution_map substitution_map_6 =
    SD::make_substitution_map(std::vector<std::pair<SD::Expression, int>>{
      std::make_pair(SD::Expression("x12"), 12),
      std::make_pair(SD::Expression("x13"), 13)});
  const SD::types::substitution_map substitution_map_7 =
    SD::make_substitution_map(std::make_pair(SD::Expression("x14"), 14.0f),
                              std::make_pair(SD::Expression("x15"), 15),
                              std::make_pair(SD::Expression("x16"), 16.0),
                              SD::types::substitution_map{
                                {SD::Expression("x17"), SD::Expression(17.0)},
                                {SD::Expression("x18"),
                                 SD::Expression(std::complex<double>(18.0))}});
  const SD::types::substitution_map substitution_map =
    SD::merge_substitution_maps(substitution_map_1,
                                substitution_map_2,
                                substitution_map_3,
                                substitution_map_4,
                                substitution_map_5,
                                substitution_map_6,
                                substitution_map_7);

  SD::Utilities::print_substitution_map(deallog, substitution_map);

  deallog << "OK" << std::endl;
}

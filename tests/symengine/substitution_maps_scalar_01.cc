// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that functions to add scalar variables to substitution maps work
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


  SD::types::substitution_map substitution_map;
  SD::add_to_substitution_map(substitution_map,
                              SD::Expression("x1"),
                              SD::Expression(1));
  SD::add_to_substitution_map(substitution_map, SD::Expression("x3"), 3.0);
  SD::add_to_substitution_map(substitution_map,
                              SD::types::symbol_vector{SD::Expression("x6"),
                                                       SD::Expression("x7")},
                              std::vector<int>{6, 7});
  SD::add_to_substitution_map(
    substitution_map,
    SD::types::substitution_map{{SD::Expression("x8"), SD::Expression(8.0)},
                                {SD::Expression("x9"), SD::Expression(9.0f)},
                                {SD::Expression("x10"), SD::Expression(10u)}});
  SD::add_to_substitution_map(substitution_map,
                              std::make_pair(SD::Expression("x11"), 11));
  SD::add_to_substitution_map(substitution_map,
                              std::vector<std::pair<SD::Expression, int>>{
                                std::make_pair(SD::Expression("x12"), 12),
                                std::make_pair(SD::Expression("x13"), 13)});
  SD::add_to_substitution_map(
    substitution_map,
    std::make_pair(SD::Expression("x14"), 14.0f),
    std::make_pair(SD::Expression("x15"), 15),
    std::make_pair(SD::Expression("x16"), 16.0),
    SD::types::substitution_map{{SD::Expression("x17"), SD::Expression(17.0)},
                                {SD::Expression("x18"),
                                 SD::Expression(std::complex<double>(18.0))}});

  SD::Utilities::print_substitution_map(deallog, substitution_map);

  if constexpr (running_in_debug_mode())
    {
      // Check that exceptions are raised when duplicate symbols
      // are found in a substitution map
      deal_II_exceptions::disable_abort_on_exception();
      try
        {
          SD::add_to_substitution_map(substitution_map,
                                      std::make_pair(SD::Expression("x14"),
                                                     14));

          deallog << "Duplicate symbol in map did not raise an error."
                  << std::endl;
        }
      catch (const ExcMessage &)
        {}
    }

  deallog << "OK" << std::endl;
}

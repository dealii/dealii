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

#ifdef DEBUG
  // Check that exceptions are raised when duplicate symbols
  // are found in a substitution map
  deal_II_exceptions::disable_abort_on_exception();
  try
    {
      SD::add_to_substitution_map(substitution_map,
                                  std::make_pair(SD::Expression("x14"), 14));

      deallog << "Duplicate symbol in map did not raise an error." << std::endl;
    }
  catch (const ExcMessage &)
    {}
#endif

  deallog << "OK" << std::endl;
}

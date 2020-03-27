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


// Check that functions to merge substitution maps work correctly.

#include <deal.II/differentiation/sd.h>

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

  deallog.push("2 maps");
  {
    SD::types::substitution_map substitution_map;
    SD::add_to_substitution_map(substitution_map,
                                SD::Expression("x1"),
                                SD::Expression(1));

    SD::types::substitution_map substitution_map_1;
    SD::add_to_substitution_map(substitution_map_1,
                                SD::Expression("x2"),
                                SD::Expression(2));

    SD::merge_substitution_maps(substitution_map, substitution_map_1);
    SD::Utilities::print_substitution_map(deallog, substitution_map);
    deallog << std::endl;
  }
  deallog.pop();

  deallog.push("Multiple maps");
  {
    SD::types::substitution_map substitution_map;
    SD::add_to_substitution_map(substitution_map,
                                SD::Expression("x1"),
                                SD::Expression(1));

    SD::types::substitution_map substitution_map_1;
    SD::add_to_substitution_map(substitution_map_1,
                                SD::Expression("x2"),
                                SD::Expression(2));

    SD::types::substitution_map substitution_map_2;
    SD::add_to_substitution_map(substitution_map_2,
                                SD::Expression("x3"),
                                SD::Expression(3));

    SD::types::substitution_map substitution_map_3;
    SD::add_to_substitution_map(substitution_map_3,
                                SD::Expression("x4"),
                                SD::Expression(4));

    SD::types::substitution_map substitution_map_4;
    SD::add_to_substitution_map(substitution_map_4,
                                SD::Expression("x5"),
                                SD::Expression(5));

    SD::merge_substitution_maps(substitution_map,
                                substitution_map_1,
                                substitution_map_2,
                                substitution_map_3,
                                substitution_map_4);
    SD::Utilities::print_substitution_map(deallog, substitution_map);
    deallog << std::endl;
  }
  deallog.pop();



  deallog.push("Const maps");
  {
    SD::types::substitution_map substitution_map_1;
    SD::add_to_substitution_map(substitution_map_1,
                                SD::Expression("x1"),
                                SD::Expression(1));

    SD::types::substitution_map substitution_map_2;
    SD::add_to_substitution_map(substitution_map_2,
                                SD::Expression("x2"),
                                SD::Expression(2));

    const SD::types::substitution_map substitution_map_1c = substitution_map_1;
    const SD::types::substitution_map substitution_map_2c = substitution_map_2;

    const SD::types::substitution_map substitution_map_out =
      SD::merge_substitution_maps(substitution_map_1c, substitution_map_2c);
    SD::Utilities::print_substitution_map(deallog, substitution_map_out);
    deallog << std::endl;
  }

  // Check that no error is raised when a duplicate entries are found
  // in two substitution maps
  {
    SD::types::substitution_map substitution_map;
    SD::add_to_substitution_map(substitution_map,
                                SD::Expression("x1"),
                                SD::Expression(1));

    SD::types::substitution_map substitution_map_1;
    SD::add_to_substitution_map(substitution_map_1,
                                SD::Expression("x1"),
                                SD::Expression(1));

    SD::merge_substitution_maps(substitution_map, substitution_map_1);
  }

#ifdef DEBUG
  // Check that exceptions are raised when duplicate symbols
  // associated with unequal values are found in a substitution map
  deal_II_exceptions::disable_abort_on_exception();
  try
    {
      {
        SD::types::substitution_map substitution_map;
        SD::add_to_substitution_map(substitution_map,
                                    SD::Expression("x1"),
                                    SD::Expression(1));

        SD::types::substitution_map substitution_map_1;
        SD::add_to_substitution_map(substitution_map_1,
                                    SD::Expression("x1"),
                                    SD::Expression(2));

        SD::merge_substitution_maps(substitution_map, substitution_map_1);
      }

      deallog
        << "Duplicate symbol with non-equal value in map did not raise an error."
        << std::endl;
    }
  catch (const ExcMessage &)
    {}
#endif

  deallog << "OK" << std::endl;
}

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


// Check that functions to add tensor variables to symbol maps, and subsequently
// set their associated values, work correctly.

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
  const unsigned int dim = 2;

  Tensor<2, dim>          t;
  SymmetricTensor<2, dim> st;
  for (unsigned int i = 0; i < t.n_independent_components; ++i)
    t[t.unrolled_to_component_indices(i)] = i;
  for (unsigned int i = 0; i < st.n_independent_components; ++i)
    st[st.unrolled_to_component_indices(i)] = i;

  deallog << "Construct symbol map" << std::endl;
  SD::types::substitution_map symbol_map;
  SD::add_to_symbol_map(symbol_map, SD::make_tensor_of_symbols<2, dim>("t1"));
  SD::add_to_symbol_map(symbol_map,
                        SD::make_symmetric_tensor_of_symbols<2, dim>("st1"));
  SD::add_to_symbol_map(symbol_map,
                        SD::make_tensor_of_symbols<2, dim>("t2"),
                        SD::make_symmetric_tensor_of_symbols<2, dim>("st2"));

  SD::Utilities::print_substitution_map(deallog, symbol_map);


  deallog << "Set values in symbol map" << std::endl;
  SD::set_value_in_symbol_map(symbol_map,
                              SD::make_tensor_of_symbols<2, dim>("t1"),
                              t);
  SD::set_value_in_symbol_map(
    symbol_map, SD::make_symmetric_tensor_of_symbols<2, dim>("st1"), st);
  SD::set_value_in_symbol_map(
    symbol_map,
    std::make_pair(SD::make_tensor_of_symbols<2, dim>("t2"), t),
    std::make_pair(SD::make_symmetric_tensor_of_symbols<2, dim>("st2"), st));

  SD::Utilities::print_substitution_map(deallog, symbol_map);

  deallog << "OK" << std::endl;
}

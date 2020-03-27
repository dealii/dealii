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


// Check that functions to create substitution maps from tensor variables work
// correctly.

#include <deal.II/differentiation/sd.h>

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


  SD::types::substitution_map substitution_map;
  SD::add_to_substitution_map(substitution_map,
                              SD::Expression("x1"),
                              SD::Expression(1));

  SD::merge_substitution_maps(
    substitution_map,
    SD::make_substitution_map(SD::make_tensor_of_symbols<2, dim>("t1"), t),
    SD::make_substitution_map(
      SD::make_symmetric_tensor_of_symbols<2, dim>("st1"), st),
    SD::make_substitution_map(SD::make_tensor_of_symbols<2, dim>("t2"), t),
    SD::make_substitution_map(
      SD::make_symmetric_tensor_of_symbols<2, dim>("st2"), st));

  SD::Utilities::print_substitution_map(deallog, substitution_map);

  deallog << "OK" << std::endl;
}

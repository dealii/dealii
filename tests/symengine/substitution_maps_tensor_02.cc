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

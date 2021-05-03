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


// Check that the wrapper for symengine expressions can be use with some
// of the algorithms defined in the standard namespace

#include <deal.II/differentiation/sd.h>

#include <algorithm>

#include "../tests.h"

namespace SD = Differentiation::SD;

int
main()
{
  initlog();

  const SD::types::symbol_vector copy_from{SD::Expression(1),
                                           SD::Expression(true),
                                           SD::Expression(2.5),
                                           SD::Expression("a")};
  SD::types::symbol_vector       copy_to(copy_from.size());

  // Mimics Tensor operator=
  std::copy(copy_from.begin(), copy_from.end(), copy_to.begin());

  deallog << "Copy" << std::endl;
  for (unsigned int i = 0; i < copy_from.size(); ++i)
    deallog << i << ": " << copy_from[i] << " -> " << copy_to[i] << std::endl;


  SD::types::symbol_vector swap_to;
  std::swap(copy_to, swap_to);

  deallog << "Swap" << std::endl;
  for (unsigned int i = 0; i < swap_to.size(); ++i)
    deallog << i << ": " << swap_to[i] << std::endl;

  deallog << "OK" << std::endl;
}

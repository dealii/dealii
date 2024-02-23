// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// The test is used to check the restriction_is_additive flags. The
// face degrees of freedom of an RT_Bubbles element must be non-additive
// as they have continuity requirements. The interior DoFs however must
// be additive.

#include <deal.II/fe/fe_rt_bubbles.h>

#include <string>

#include "../tests.h"


template <int dim>
void
test(const unsigned int degree)
{
  FE_RT_Bubbles<dim> fe_rt_bubbles(degree);

  deallog << "Degree=" << degree
          << ", restriction is additive flags:" << std::endl;

  for (unsigned int i = 0; i < fe_rt_bubbles.dofs_per_cell; ++i)
    deallog << fe_rt_bubbles.restriction_is_additive(i) << ' ';

  deallog << std::endl;
}



int
main()
{
  initlog();

  deallog << "Dimension 2: " << std::endl;
  for (unsigned int i = 1; i < 4; ++i)
    test<2>(i);

  deallog << "Dimension 3: " << std::endl;
  for (unsigned int i = 1; i < 4; ++i)
    test<3>(i);

  return 0;
}

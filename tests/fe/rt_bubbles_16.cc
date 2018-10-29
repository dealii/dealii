// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// The test is used to check the restriction_is_additive flags. The
// face degrees of freedom of an RT_Bubbles element must be non-additive
// as they have continuity requrements. The interior DoFs however must
// be additive.

#include <deal.II/fe/fe_rt_bubbles.h>

#include <string>

#include "../tests.h"


std::ofstream logfile("output");

template <int dim>
void
test(const unsigned int degree)
{
  FE_RT_Bubbles<dim> fe_rt_bubbles(degree);

  deallog << "Degree=" << degree
          << ", restriction is additive flags:" << std::endl;

  for (unsigned int i = 0; i < fe_rt_bubbles.dofs_per_cell; ++i)
    deallog << fe_rt_bubbles.restriction_is_additive(i) << " ";

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

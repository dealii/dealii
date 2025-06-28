// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// it turns out that FE_Q::face_to_cell_index() had a bug for elements beyond
// Q2 when using the face orientation flag. this test is for the 2d case for the
// Q4 case

#include <deal.II/fe/fe_q.h>

#include <iostream>

#include "../tests.h"

template <int dim>
void
test()
{
  FE_Q<dim>          fe(4);
  const unsigned int dofs_per_face = fe.dofs_per_face;

  for (unsigned int face = 0; face < 4; ++face)
    {
      deallog << "Face=" << face << std::endl;

      for (types::geometric_orientation orientation = 0; orientation < 2;
           ++orientation)
        {
          deallog << "  orientation=" << (orientation == 0 ? "false" : "true")
                  << std::endl
                  << "    ";
          for (unsigned int i = 0; i < dofs_per_face; ++i)
            deallog << fe.face_to_cell_index(i, face, orientation) << " - ";
          deallog << std::endl;
        }
    }
}

int
main()
{
  initlog();

  test<2>();
}

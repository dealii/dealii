// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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


// it turns out that FE_Q::face_to_cell_index() had a bug for elements beyond
// Q2 when using the face flip flag. this test is for the 2d case for the Q4
// case

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

      for (int flip = 0; flip < 2; ++flip)
        {
          deallog << "  flip=" << (flip == 0 ? "false" : "true") << std::endl
                  << "    ";
          for (unsigned int i = 0; i < dofs_per_face; ++i)
            deallog << fe.face_to_cell_index(
                         i, face, true, (flip == 0 ? false : true), false)
                    << " - ";
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

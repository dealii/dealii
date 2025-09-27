// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that computation of hp-constraints works for RT elements correctly

char logname[] = "output";


#include <deal.II/fe/fe_raviart_thomas.h>

#include "../hp/hp_constraints_common.h"


template <int dim>
void
test()
{
  if (dim == 1)
    return;

  deallog << "Test for dim = " << dim << std::endl << std::endl;
  hp::FECollection<dim>     fe;
  std::vector<unsigned int> degrees;
  for (unsigned int i = 0; i < 7 - dim; ++i)
    {
      fe.push_back(FE_RaviartThomasNodal<dim>(i));
      degrees.push_back(i);
    }

  deallog << "No hanging nodes test" << std::endl;
  test_no_hanging_nodes(fe);
  deallog << std::endl << std::endl;

  deallog << "Hanging nodes test" << std::endl;
  test_with_hanging_nodes(fe);
  deallog << std::endl << std::endl;

  deallog << "Wrong face orientation test" << std::endl;
  test_with_wrong_face_orientation(fe);
  deallog << std::endl << std::endl;

  deallog << "2d deformed mesh test" << std::endl;
  test_with_2d_deformed_mesh(fe);
  deallog << std::endl << std::endl;

  deallog << "2d deformed refined mesh test" << std::endl;
  test_with_2d_deformed_refined_mesh(fe);
  deallog << std::endl << std::endl;

  deallog << "Interpolation test" << std::endl;
  test_interpolation(fe, degrees);

  deallog << std::endl << std::endl;
  deallog << std::endl << std::endl;
}

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
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
// face degrees of freedom of a BDM element must be non-additive as
// they have continuity requirements, however the interior DOFs must
// be additive, e.g. for order 1 elements all DOFs are non-additive,
// while for the order 2 element in 2d we have 12 non-additive face DOFs
// and 2 additive interior ones. The test should output a vector
// consisting of faces_per_cell * dofs_per_face zeros, followed by
// interior_dofs ones.

#include <deal.II/fe/fe_bdm.h>

#include <string>

#include "../tests.h"


template <int dim>
void
test(const unsigned int degree)
{
  FE_BDM<dim> fe_bdm(degree);

  deallog << "Degree=" << degree
          << ", restriction is additive flags:" << std::endl;

  for (unsigned int i = 0; i < fe_bdm.dofs_per_cell; ++i)
    deallog << fe_bdm.restriction_is_additive(i) << ' ';

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

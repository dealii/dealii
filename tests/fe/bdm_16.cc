// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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


// The test is used to check the restriction_is_additive flags. The
// face degrees of freedom of a BDM element must be non-additive as
// they have continuity requirements, however the interior DOFs must
// be additive, e.g. for order 1 elements all DOFs are non-additive,
// while for the order 2 element in 2d we have 12 non-additive face DOFs
// and 2 additive interior ones. The test should output a vector
// consisting of faces_per_cell * dofs_per_face zeros, followed by
// interior_dofs ones.

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/fe/fe_bdm.h>

#include <fstream>
#include <string>


std::ofstream logfile ("output");

template<int dim>
void
test (const unsigned int degree)
{
  FE_BDM<dim> fe_bdm(degree);

  deallog << "Degree=" << degree
          << ", restriction is additive flags:"
          << std::endl;

  for (unsigned int i=0; i<fe_bdm.dofs_per_cell; ++i)
    std::cout << fe_bdm.restriction_is_additive(i) << " ";

  deallog << std::endl;
}


int
main()
{
  initlog();

  deallog << "Dimension 2: " << std::endl;
  for (unsigned int i=1; i<4; ++i)
    test<2>(i);

  deallog << "Dimension 3: " << std::endl;
  for (unsigned int i=1; i<4; ++i)
    test<3>(i);

  return 0;
}




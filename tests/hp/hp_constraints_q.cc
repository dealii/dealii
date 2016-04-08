// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2014 by the deal.II authors
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



// check that computation of hp constraints works for Q elements correctly

char logname[] = "output";


#include "hp_constraints_common.h"


template <int dim>
void test ()
{
  deallog << "Test for dim = " << dim << std::endl << std::endl;
  hp::FECollection<dim> fe;
  std::vector<unsigned int> degrees;
  for (unsigned int i=1; i<4; ++i)
    {
      fe.push_back (FE_Q<dim>(QIterated<1>(QTrapez<1>(),i)));
      degrees.push_back (i);
    }

  deallog << "No hanging nodes test" << std::endl;
  test_no_hanging_nodes  (fe);
  deallog << std::endl << std::endl;

  deallog << "Hanging nodes test" << std::endl;
  test_with_hanging_nodes  (fe);
  deallog << std::endl << std::endl;

  deallog << "Wrong face orientation test" << std::endl;
  test_with_wrong_face_orientation  (fe);
  deallog << std::endl << std::endl;

  deallog << "2d deformed mesh test" << std::endl;
  test_with_2d_deformed_mesh  (fe);
  deallog << std::endl << std::endl;

  deallog << "2d deformed refined mesh test" << std::endl;
  test_with_2d_deformed_refined_mesh  (fe);
  deallog << std::endl << std::endl;

  deallog << "Interpolation test" << std::endl;
  test_interpolation  (fe, degrees);
  deallog << std::endl << std::endl;

  deallog << "Random hanging nodes" << std::endl;
  test_with_hanging_nodes_random (fe);

  deallog << std::endl << std::endl;
  deallog << std::endl << std::endl;
}

//----------------------------  hp_constraints_q_system_x.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  hp_constraints_q_system_x.cc  ---------------------------


// check that computation of hp constraints works for FESystem(FE_Q) elements correctly
// on a uniformly refined mesh for functions of degree q

// these tests check that we can deal with FESystem(FE_Q(p),FE_DGQ(q)) for
// different p,q

char logname[] = "hp_constraints_q_system_x/output";


#include "hp_constraints_common.h"


template <int dim>
void test ()
{
  deallog << "Test for dim = " << dim << std::endl << std::endl;
  hp::FECollection<dim> fe;
  for (unsigned int i=1; i<4; ++i)
    for (unsigned int j=0; j<4; ++j)
      fe.push_back (FESystem<dim>(FE_Q<dim>(i), 1,
				  FE_DGQ<dim>(j), 1));

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

  deallog << "Random hanging nodes" << std::endl;
  test_with_hanging_nodes_random (fe);

  deallog << std::endl << std::endl;
  deallog << std::endl << std::endl;
}

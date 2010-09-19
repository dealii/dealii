//----------------------------  project_q_system_03.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  project_q_system_03.cc  ---------------------------


// check that VectorTools::project works for FESystem(FE_Q) elements correctly
// on a uniformly refined mesh for functions of degree q

char logname[] = "project_q_system_03/output";


#include "project_common.h"


template <int dim>
void test ()
{
  for (unsigned int p=1; p<5-dim; ++p)
    test_with_wrong_face_orientation (FESystem<dim>(FE_Q<dim>(p), 1,
					 FE_DGQ<dim>(p+1), 1),
			   p);
}

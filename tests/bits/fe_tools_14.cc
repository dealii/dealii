//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_tools_14.cc  ---------------------------

#include "../tests.h"
#include "fe_tools_common.cc"

// comparison of the results of
//   FE::face_to_equivalent_cell_index (i)
// and
// FE::face_to_cell_index (i, 0)
// According to the documentation they should be the same!


std::string output_file_name = "fe_tools_14/output";


template <int dim>
void
check_this (const FiniteElement<dim> &fe1,
            const FiniteElement<dim> &fe2)
{
  bool ok = true;
  
  for (unsigned int i = 0; i < fe1.dofs_per_face; ++i)
    {
      if (fe1.face_to_equivalent_cell_index (i) !=
	  fe1.face_to_cell_index (i, 0))
	ok = false;

    }

  for (unsigned int i = 0; i < fe2.dofs_per_face; ++i)
    {
      if (fe2.face_to_equivalent_cell_index (i) !=
	  fe2.face_to_cell_index (i, 0))
	ok = false;
    }
  
  ok ? deallog << "OK" << std::endl :
    deallog << "Error!" << std::endl;
}


//    $Id$
//    Version: $Name: $ 
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
#include "dof_tools_common.h"

// comparison of the results of
//   FE::face_to_equivalent_cell_index (i)
// and
// FE::face_to_cell_index (i, 0)
// According to the documentation they should be the same!


std::string output_file_name = "fe_face_to_cell_index/output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();
  
  for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
    {
      deallog << fe.face_to_equivalent_cell_index (i)
	      << ' '
	      << fe.face_to_cell_index (i, 0)
	      << std::endl;
      Assert (fe.face_to_equivalent_cell_index (i) ==
	      fe.face_to_cell_index (i, 0),
	      ExcInternalError());
    }
  
  deallog << "OK" << std::endl;
}


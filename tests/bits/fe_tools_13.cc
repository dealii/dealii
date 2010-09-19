//----------------------------  fe_tools_13.cc  ---------------------------
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
//----------------------------  fe_tools_13.cc  ---------------------------

#include "../tests.h"
#include "fe_tools_common.h"

// check
//   FE::face_to_equivalent_cell_index


std::string output_file_name = "fe_tools_13/output";



template <int dim>
void
check_this (const FiniteElement<dim> &fe1,
            const FiniteElement<dim> &fe2)
{
  for (unsigned int i=0; i<fe1.dofs_per_face; ++i)
    deallog << fe1.face_to_equivalent_cell_index(i) << std::endl;

  for (unsigned int i=0; i<fe2.dofs_per_face; ++i)
    deallog << fe2.face_to_equivalent_cell_index(i) << std::endl;
}


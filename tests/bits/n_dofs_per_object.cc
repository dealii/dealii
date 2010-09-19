//----------------------------  n_dofs_per_object.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  n_dofs_per_object.cc  ---------------------------

#include "../tests.h"
#include "dof_tools_common.h"
#include <lac/vector.h>

// check
//   FiniteElement::n_dofs_per_object



std::string output_file_name = "n_dofs_per_object/output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();
  deallog << fe.dofs_per_vertex << ' '
	  << fe.dofs_per_line << ' '
	  << fe.dofs_per_quad << ' '
	  << fe.dofs_per_hex << std::endl;
  deallog << fe.template n_dofs_per_object<0>() << ' '
	  << fe.template n_dofs_per_object<1>() << ' '
	  << fe.template n_dofs_per_object<2>() << ' '
	  << fe.template n_dofs_per_object<3>() << std::endl;  
}

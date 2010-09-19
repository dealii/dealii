//----------------------------  hp_constraints_are_implemented.cc  ---------------------------
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
//----------------------------  hp_constraints_are_implemented.cc  ---------------------------

#include "../tests.h"
#include "dof_tools_common.h"
#include <lac/vector.h>

// check
//   FE::hp_constraints_are_implemented
// a bit like fe_tools_14, but works on a different set of elements



std::string output_file_name = "hp_constraints_are_implemented/output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  deallog << dof_handler.get_fe().get_name()
	  << ": "
	  << (dof_handler.get_fe().hp_constraints_are_implemented() ? "true" : "false")
	  << std::endl;
}

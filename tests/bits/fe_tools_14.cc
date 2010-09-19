//----------------------------  fe_tools_14.cc  ---------------------------
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
#include "fe_tools_common.h"

// check
//   FE::hp_constraints_are_implemented ()
// a bit like hp_constraints_are_implemented, but with a different
// set of elements


std::string output_file_name = "fe_tools_14/output";



template <int dim>
void
check_this (const FiniteElement<dim> &fe1,
            const FiniteElement<dim> &fe2)
{
  deallog << (fe1.hp_constraints_are_implemented () ? "true" : "false")
	  << std::endl;
  deallog << (fe2.hp_constraints_are_implemented () ? "true" : "false")
	  << std::endl;
}


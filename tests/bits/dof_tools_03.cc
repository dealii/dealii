//----------------------------  dof_tools_3.cc  ---------------------------
//    $Id$
//    Version: 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_tools_3.cc  ---------------------------

#include "../tests.h"
#include "dof_tools_common.cc"
#include <dofs/dof_constraints.h>

// check
//   DoFTools::
//   make_hanging_node_constraints (const DoFHandler<dim> &,
//	                            ConstraintMatrix      &);


std::string output_file_name = "dof_tools_03.output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
                                   // don't run this test if hanging
                                   // nodes are not implemented
  if (dof_handler.get_fe().constraints_are_implemented() == false)
    return;
  
  ConstraintMatrix cm;
  DoFTools::make_hanging_node_constraints (dof_handler, cm);
  cm.close ();

  deallog << cm.n_constraints () << std::endl;
  deallog << cm.max_constraint_indirections () << std::endl;

  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    deallog << (cm.is_constrained(i) ? '0' : '1');
  deallog << std::endl;
  
  deallog << cm.n_constraints () << std::endl;
  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    deallog << (cm.is_identity_constrained(i) ? '0' : '1');
  deallog << std::endl;  
}

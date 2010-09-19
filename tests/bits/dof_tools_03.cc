//----------------------------  dof_tools_3.cc  ---------------------------
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
//----------------------------  dof_tools_3.cc  ---------------------------

#include "../tests.h"
#include "dof_tools_common.h"
#include <lac/constraint_matrix.h>
#include <lac/vector.h>

// check
//   DoFTools::
//   make_hanging_node_constraints (const DoFHandler<dim> &,
//	                            ConstraintMatrix      &);


std::string output_file_name = "dof_tools_03/output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
				   // use a higher output accuracy for this
				   // test. the reason is that many of the
				   // constraints are negative powers of 2,
				   // which have exact representations with 3
				   // or 4 digits of accuracy, but not with
				   // the usual 2 digits (for example, 0.375,
				   // which sometimes rounds to 0.38 and
				   // sometimes to 0.37, depending on how
				   // intermediate errors have accumulated)
  deallog.get_file_stream().precision (5);
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

  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    if (cm.is_constrained(i))
      {
	deallog << "not identity constrained: " << i << std::endl;
	Vector<double> v (dof_handler.n_dofs());
	v(i) = 1;
	cm.condense (v);
	for (unsigned int j=0; j<dof_handler.n_dofs(); ++j)
	  if (v(j) != 0)
	    deallog << "  line " << j << ": " << v(j) << std::endl;
      }
}

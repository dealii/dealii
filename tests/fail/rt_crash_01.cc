//----------------------------  rt_crash_01.cc  ---------------------------
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
//----------------------------  rt_crash_01.cc  ---------------------------

#include "../tests.h"
#include "../bits/dof_tools_common.cc"
#include <lac/constraint_matrix.h>
#include <base/function.h>
#include <base/quadrature_lib.h>
#include <numerics/vectors.h>
#include <lac/vector.h>
#include <dofs/dof_tools.h>
#include <numerics/data_out.h>

// check an abort in FEPolyTensor when used with RaviartThomas
// elements, when computing some sort of edge directions. This is the
// case left over from dof_tools_19

std::string output_file_name = "rt_crash_01/output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
                                   // there's presently a crash in the
                                   // Raviart-Thomas element. only
                                   // check for these elements,
                                   // everything else is handled in
                                   // dof_tools_19
  if (dof_handler.get_fe().get_name().find ("RaviartThomas") ==
      std::string::npos)
    return;
  
  ConstantFunction<dim> test_func (1, dof_handler.get_fe().n_components ()); 

                                   // don't run this test if hanging
                                   // nodes are not implemented
  if (dof_handler.get_fe().constraints_are_implemented() == false)
    return;

  ConstraintMatrix cm;
  DoFTools::make_hanging_node_constraints (dof_handler, cm);
  cm.close ();

  deallog << cm.n_constraints () << std::endl;
  deallog << cm.max_constraint_indirections () << std::endl;

  // L_2 project constant function onto field
  QGauss<dim> quadrature (6);
  Vector<double> solution (dof_handler.n_dofs ());

  VectorTools::project (dof_handler, cm,
			quadrature, test_func,
			solution);
  cm.distribute (solution);

				   // all values of the projected solution
				   // should be close around 1 (the value of
				   // the original function)
  for (unsigned int i=0; i<solution.size(); ++i)
    Assert (std::fabs (solution(i)-1) < 1e-6,
	    ExcInternalError());
  
  // Evaluate error
  Vector<double> cellwise_errors (dof_handler.get_tria ().n_active_cells());
  VectorTools::integrate_difference (dof_handler, solution, test_func,
                                     cellwise_errors, quadrature,
                                     VectorTools::L2_norm);
  const double p_l2_error = cellwise_errors.l2_norm();

  deallog << "L2_Error : " << p_l2_error << std::endl;
  
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "u");
    data_out.build_patches ();
    data_out.write_gnuplot (deallog.get_file_stream());
  }
  
  Assert (p_l2_error < 1e-12, ExcInternalError());  
}

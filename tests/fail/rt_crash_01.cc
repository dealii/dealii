// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include "../tests.h"
#include "../bits/dof_tools_common.h"
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>

// check an abort in FEPolyTensor when used with RaviartThomas
// elements, when computing some sort of edge directions. This is the
// case left over from dof_tools_19

std::string output_file_name = "output";


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

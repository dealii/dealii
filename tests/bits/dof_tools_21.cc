//----------------------------  grid_tools.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006, 2008, 2011, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_tools_21.cc  ---------------------------

#include "../tests.h"
#include "dof_tools_periodic.h"
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/numerics/vectors.h>
#include <deal.II/lac/vector.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <sstream>

// A simple test for
//   DoFTools::
//   make_periodicity_constraints (const FaceIterator       &,
//                                 const FaceIterator       &,
//                                 dealii::ConstraintMatrix &,
//                                 const std::vector<bool>  &)
//
// We project an already periodic function onto the FE space of
// periodic functions and store the resulting L2 difference
// between the projected and unprojected function.
// This should reveal any errors in the constraint_matrix:

std::string output_file_name("dof_tools_21/output");


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  Functions::CosineFunction<dim> test_func (dof_handler.get_fe().n_components ());

  ConstraintMatrix cm;

  // Apply periodic boundary conditions only in the one direction where
  // we can match the (locally refined) faces:
  DoFTools::make_periodicity_constraints (dof_handler.begin(0)->face(0),
                                          dof_handler.begin(0)->face(1),
                                          cm);
  cm.close ();

  deallog << cm.n_constraints () << std::endl;
  deallog << cm.max_constraint_indirections () << std::endl;

  QGauss<dim> quadrature (6);

  Vector<double> unconstrained (dof_handler.n_dofs ());
  Vector<double> constrained (dof_handler.n_dofs ());

  // Ensure that we can interpolate:
  if (dof_handler.get_fe().get_unit_support_points().size() == 0)
    return;

  VectorTools::interpolate (dof_handler,
                            test_func,
                            unconstrained);

  constrained = unconstrained;
  cm.distribute (constrained);

  constrained -= unconstrained;

  const double p_l2_error = constrained.l2_norm();

  Assert (p_l2_error < 1e-11, ExcInternalError());

  deallog << "L2_Error : " << p_l2_error << std::endl;
}

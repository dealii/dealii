// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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
#include "fe_tools_common.h"
#include <deal.II/base/quadrature_lib.h>

// check
//   FETools::compute_projection_from_quadrature_points_matrix
// we put this into the fe_tools_common framework for simplicity, but
// in fact we ignore the second FE it passes to the check_this() function
// and we can only test as well for scalar elements, since this is all
// the function presently supports.
//
// the test makes sure that if the FE has support points and we use
// these support points as quadrature points, that the resulting
// matrix is the unit matrix


std::string output_file_name = "output";


template <int dim>
void
check_this (const FiniteElement<dim> &fe,
            const FiniteElement<dim> &/*fe2*/)
{
  // only check if both elements have
  // support points. otherwise,
  // interpolation doesn't really
  // work
  if (fe.n_components() != 1)
    return;

  // ignore this check if this fe has already
  // been treated
  static std::set<std::string> already_checked;
  if (already_checked.find(fe.get_name()) != already_checked.end())
    return;
  already_checked.insert (fe.get_name());

  // only test elements with support
  // points
  if (fe.has_support_points() == false)
    return;

  // test with different quadrature formulas
  Quadrature<dim> q_rhs(fe.get_unit_support_points(),
                        std::vector<double> (fe.dofs_per_cell,
                                             1./fe.dofs_per_cell));

  FullMatrix<double> X (fe.dofs_per_cell,
                        q_rhs.size());

  Assert (X.m() == X.n(), ExcInternalError());

  FETools::compute_projection_from_quadrature_points_matrix (fe,
                                                             q_rhs, q_rhs,
                                                             X);

  for (unsigned int i=0; i<X.m(); ++i)
    X(i,i) -= 1;

  Assert (X.frobenius_norm() < 1e-10, ExcInternalError());
}


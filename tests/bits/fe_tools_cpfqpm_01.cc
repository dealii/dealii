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
// this test simply computes the matrix and outputs some of its
// characteristics


std::string output_file_name = "output";


template <int dim>
void
check_this (const FiniteElement<dim> &fe,
            const FiniteElement<dim> &/*fe2*/)
{
  deallog << std::setprecision (9);

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


  // test with different quadrature formulas
  QGauss<dim> q_lhs(fe.degree+1);
  QGauss<dim> q_rhs(fe.degree+1>2 ? fe.degree+1-2 : 1);

  FullMatrix<double> X (fe.dofs_per_cell,
                        q_rhs.size());

  FETools::compute_projection_from_quadrature_points_matrix (fe,
                                                             q_lhs, q_rhs,
                                                             X);

  output_matrix (X);
}


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
// this test makes sure that projecting onto a finite element space
// sufficiently fine to hold the quadrature point data, then interpolating
// back to the quadrature points is an identity operation


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


  // test with the same quadrature formulas
  // of a degree that is high enough to
  // exactly capture the data
  QGauss<dim> q_lhs(fe.degree+1);
  QGauss<dim> q_rhs(fe.degree+1);

  // this test can only succeed if there are
  // at least as many degrees of freedom in
  // the finite element as there are
  // quadrature points
  if (fe.dofs_per_cell < q_rhs.size())
    return;

  deallog << "dofs_per_cell=" << fe.dofs_per_cell
          << ", n_q_points=" << q_rhs.size()
          << std::endl;

  FullMatrix<double> X (fe.dofs_per_cell,
                        q_rhs.size());

  FETools::compute_projection_from_quadrature_points_matrix (fe,
                                                             q_lhs, q_rhs,
                                                             X);

  // then compute the matrix that
  // interpolates back to the quadrature
  // points
  FullMatrix<double> I_q (q_rhs.size(), fe.dofs_per_cell);
  FETools::compute_interpolation_to_quadrature_points_matrix (fe, q_rhs,
                                                              I_q);

  FullMatrix<double> product (q_rhs.size(),
                              q_rhs.size());
  I_q.mmult (product, X);

  // the product should be the identity
  // matrix now. make sure that this is
  // indeed the case
  for (unsigned int i=0; i<product.m(); ++i)
    product(i,i) -= 1;

  output_matrix (product);
  Assert (product.frobenius_norm() < 1e-10, ExcInternalError());
}


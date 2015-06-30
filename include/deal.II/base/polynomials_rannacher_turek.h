// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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


#ifndef dealii__polynomials_rannacher_turek_h
#define dealii__polynomials_rannacher_turek_h

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <vector>

DEAL_II_NAMESPACE_OPEN


/**
 * Basis for polynomial space on the unit square used for lowest order
 * Rannacher Turek element.
 *
 * The i-th basis function is the dual basis element corresponding to
 * the dof which evaluates the function's mean value across the i-th
 * face. The numbering can be found in GeometryInfo.
 *
 * @ingroup Polynomials
 * @author Patrick Esser
 * @date 2015
 **/
template <int dim>
class PolynomialsRannacherTurek
{
public:
  /**
   * Dimension we are working in.
   */
  static const unsigned int dimension = dim;

  /**
   * Constructor, checking that the basis is implemented in this
   * dimension.
   */
  PolynomialsRannacherTurek();

  /** Value of basis function @p i at @p p.
  */
  double compute_value(const unsigned int i,
                       const Point<dim> &p) const;
  /** Gradient of basis function @p i at @p p.
  */
  Tensor<1, dim> compute_grad(const unsigned int i,
                              const Point<dim> &p) const;
  /** Gradient of gradient of basis function @p i at @p p.
  */
  Tensor<2, dim> compute_grad_grad(const unsigned int i,
                                   const Point<dim> &p) const;

  /**
   * Compute values and derivatives of all basis functions at @p
   * unit_point.
   *
   * Size of the vectors must be either equal to the number of
   * polynomials or zero. A size of zero means that we are not
   * computing the vector entries.
   */
  void compute(const Point<dim> &unit_point,
               std::vector<double> &values,
               std::vector<Tensor<1, dim> > &grads,
               std::vector<Tensor<2, dim> > &grad_grads) const;
};


DEAL_II_NAMESPACE_CLOSE

#endif

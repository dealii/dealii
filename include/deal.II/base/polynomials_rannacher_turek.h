// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#ifndef dealii_polynomials_rannacher_turek_h
#define dealii_polynomials_rannacher_turek_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>
#include <deal.II/base/scalar_polynomials_base.h>
#include <deal.II/base/tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN


/**
 * Basis for polynomial space on the unit square used for lowest order
 * Rannacher Turek element.
 *
 * The i-th basis function is the dual basis element corresponding to the dof
 * which evaluates the function's mean value across the i-th face. The
 * numbering can be found in GeometryInfo.
 *
 * @ingroup Polynomials
 * @author Patrick Esser
 * @date 2015
 */
template <int dim>
class PolynomialsRannacherTurek : public ScalarPolynomialsBase<dim>
{
public:
  /**
   * Dimension we are working in.
   */
  static const unsigned int dimension = dim;

  /**
   * Constructor, checking that the basis is implemented in this dimension.
   */
  PolynomialsRannacherTurek();

  /**
   * Value of basis function @p i at @p p.
   */
  double
  compute_value(const unsigned int i, const Point<dim> &p) const override;

  /**
   * <tt>order</tt>-th of basis function @p i at @p p.
   *
   * Consider using evaluate() instead.
   */
  template <int order>
  Tensor<order, dim>
  compute_derivative(const unsigned int i, const Point<dim> &p) const;

  /**
   * @copydoc ScalarPolynomialsBase::compute_1st_derivative()
   */
  virtual Tensor<1, dim>
  compute_1st_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_2nd_derivative()
   */
  virtual Tensor<2, dim>
  compute_2nd_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_3rd_derivative()
   */
  virtual Tensor<3, dim>
  compute_3rd_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_4th_derivative()
   */
  virtual Tensor<4, dim>
  compute_4th_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * Gradient of basis function @p i at @p p.
   */
  Tensor<1, dim>
  compute_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * Gradient of gradient of basis function @p i at @p p.
   */
  Tensor<2, dim>
  compute_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * Compute values and derivatives of all basis functions at @p unit_point.
   *
   * Size of the vectors must be either equal to the number of polynomials or
   * zero. A size of zero means that we are not computing the vector entries.
   */
  void
  evaluate(const Point<dim> &           unit_point,
           std::vector<double> &        values,
           std::vector<Tensor<1, dim>> &grads,
           std::vector<Tensor<2, dim>> &grad_grads,
           std::vector<Tensor<3, dim>> &third_derivatives,
           std::vector<Tensor<4, dim>> &fourth_derivatives) const override;

  /**
   * Return the name of the space, which is <tt>RannacherTurek</tt>.
   */
  std::string
  name() const override;

  /**
   * @copydoc ScalarPolynomialsBase::clone()
   */
  virtual std::unique_ptr<ScalarPolynomialsBase<dim>>
  clone() const override;
};


namespace internal
{
  namespace PolynomialsRannacherTurekImplementation
  {
    template <int order, int dim>
    inline Tensor<order, dim>
    compute_derivative(const unsigned int, const Point<dim> &)
    {
      Assert(dim == 2, ExcNotImplemented());
      return Tensor<order, dim>();
    }


    template <int order>
    inline Tensor<order, 2>
    compute_derivative(const unsigned int i, const Point<2> &p)
    {
      const unsigned int dim = 2;

      Tensor<order, dim> derivative;
      switch (order)
        {
          case 1:
            {
              Tensor<1, dim> &grad =
                *reinterpret_cast<Tensor<1, dim> *>(&derivative);
              if (i == 0)
                {
                  grad[0] = -2.5 + 3 * p(0);
                  grad[1] = 1.5 - 3 * p(1);
                }
              else if (i == 1)
                {
                  grad[0] = -0.5 + 3.0 * p(0);
                  grad[1] = 1.5 - 3.0 * p(1);
                }
              else if (i == 2)
                {
                  grad[0] = 1.5 - 3.0 * p(0);
                  grad[1] = -2.5 + 3.0 * p(1);
                }
              else if (i == 3)
                {
                  grad[0] = 1.5 - 3.0 * p(0);
                  grad[1] = -0.5 + 3.0 * p(1);
                }
              else
                {
                  Assert(false, ExcNotImplemented());
                }
              return derivative;
            }
          case 2:
            {
              Tensor<2, dim> &grad_grad =
                *reinterpret_cast<Tensor<2, dim> *>(&derivative);
              if (i == 0)
                {
                  grad_grad[0][0] = 3;
                  grad_grad[0][1] = 0;
                  grad_grad[1][0] = 0;
                  grad_grad[1][1] = -3;
                }
              else if (i == 1)
                {
                  grad_grad[0][0] = 3;
                  grad_grad[0][1] = 0;
                  grad_grad[1][0] = 0;
                  grad_grad[1][1] = -3;
                }
              else if (i == 2)
                {
                  grad_grad[0][0] = -3;
                  grad_grad[0][1] = 0;
                  grad_grad[1][0] = 0;
                  grad_grad[1][1] = 3;
                }
              else if (i == 3)
                {
                  grad_grad[0][0] = -3;
                  grad_grad[0][1] = 0;
                  grad_grad[1][0] = 0;
                  grad_grad[1][1] = 3;
                }
              return derivative;
            }
          default:
            {
              // higher derivatives are all zero
              return Tensor<order, dim>();
            }
        }
    }
  } // namespace PolynomialsRannacherTurekImplementation
} // namespace internal



// template functions
template <int dim>
template <int order>
Tensor<order, dim>
PolynomialsRannacherTurek<dim>::compute_derivative(const unsigned int i,
                                                   const Point<dim> & p) const
{
  return internal::PolynomialsRannacherTurekImplementation::compute_derivative<
    order>(i, p);
}



template <int dim>
inline Tensor<1, dim>
PolynomialsRannacherTurek<dim>::compute_1st_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  return compute_derivative<1>(i, p);
}



template <int dim>
inline Tensor<2, dim>
PolynomialsRannacherTurek<dim>::compute_2nd_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  return compute_derivative<2>(i, p);
}



template <int dim>
inline Tensor<3, dim>
PolynomialsRannacherTurek<dim>::compute_3rd_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  return compute_derivative<3>(i, p);
}



template <int dim>
inline Tensor<4, dim>
PolynomialsRannacherTurek<dim>::compute_4th_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  return compute_derivative<4>(i, p);
}



template <int dim>
inline std::string
PolynomialsRannacherTurek<dim>::name() const
{
  return "RannacherTurek";
}


DEAL_II_NAMESPACE_CLOSE

#endif

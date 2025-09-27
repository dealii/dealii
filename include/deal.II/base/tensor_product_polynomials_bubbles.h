// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_tensor_product_polynomials_bubbles_h
#define dealii_tensor_product_polynomials_bubbles_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/utilities.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup Polynomials
 * @{
 */

/**
 * A class that represents a space of tensor product polynomials, augmented
 * by $dim$ (non-normalized) bubble functions of form
 * $\varphi_j(\mathbf x)
 * = 2^{\text{degree}-1}\left(x_j-frac 12\right)^{\text{degree}-1}
 * \left[\prod_{i=0}^{dim-1}(x_i(1-x_i))\right]$
 * for $j=0,\ldots,dim-1$. If `degree` is one, then the first factor
 * disappears and one receives the usual bubble function centered
 * at the mid-point of the cell.
 *
 * This class inherits
 * most of its functionality from TensorProductPolynomials. The bubble
 * enrichments are added for the last index.
 */
template <int dim>
class TensorProductPolynomialsBubbles : public ScalarPolynomialsBase<dim>
{
public:
  /**
   * Access to the dimension of this object, for checking and automatic
   * setting of dimension in other classes.
   */
  static constexpr unsigned int dimension = dim;

  /**
   * Constructor. <tt>pols</tt> is a vector of objects that should be derived
   * or otherwise convertible to one-dimensional polynomial objects. It will
   * be copied element by element into a private variable.
   */
  template <class Pol>
  TensorProductPolynomialsBubbles(const std::vector<Pol> &pols);

  /**
   * Print the list of <tt>tensor_polys</tt> indices to <tt>out</tt>.
   */
  void
  output_indices(std::ostream &out) const;

  /**
   * Set the ordering of the polynomials. Requires
   * <tt>renumber.size()==tensor_polys.n()</tt>.  Stores a copy of
   * <tt>renumber</tt>.
   */
  void
  set_numbering(const std::vector<unsigned int> &renumber);

  /**
   * Give read access to the renumber vector.
   */
  const std::vector<unsigned int> &
  get_numbering() const;

  /**
   * Give read access to the inverse renumber vector.
   */
  const std::vector<unsigned int> &
  get_numbering_inverse() const;

  /**
   * Compute the value and the first and second derivatives of each tensor
   * product polynomial at <tt>unit_point</tt>.
   *
   * The size of the vectors must either be equal 0 or equal n(). In the first
   * case, the function will not compute these values.
   *
   * If you need values or derivatives of all tensor product polynomials then
   * use this function, rather than using any of the compute_value(),
   * compute_grad() or compute_grad_grad() functions, see below, in a loop
   * over all tensor product polynomials.
   */
  void
  evaluate(const Point<dim>            &unit_point,
           std::vector<double>         &values,
           std::vector<Tensor<1, dim>> &grads,
           std::vector<Tensor<2, dim>> &grad_grads,
           std::vector<Tensor<3, dim>> &third_derivatives,
           std::vector<Tensor<4, dim>> &fourth_derivatives) const override;

  /**
   * Compute the value of the <tt>i</tt>th tensor product polynomial at
   * <tt>unit_point</tt>. Here <tt>i</tt> is given in tensor product
   * numbering.
   *
   * Note, that using this function within a loop over all tensor product
   * polynomials is not efficient, because then each point value of the
   * underlying (one-dimensional) polynomials is (unnecessarily) computed
   * several times.  Instead use the evaluate() function with
   * <tt>values.size()==</tt>n() to get the point values of all tensor
   * polynomials all at once and in a much more efficient way.
   */
  double
  compute_value(const unsigned int i, const Point<dim> &p) const override;

  /**
   * Compute the order @p order derivative of the <tt>i</tt>th tensor product
   * polynomial at <tt>unit_point</tt>. Here <tt>i</tt> is given in tensor
   * product numbering.
   *
   * Note, that using this function within a loop over all tensor product
   * polynomials is not efficient, because then each derivative value of the
   * underlying (one-dimensional) polynomials is (unnecessarily) computed
   * several times.  Instead use the evaluate() function, see above, with the
   * size of the appropriate parameter set to n() to get the point value of
   * all tensor polynomials all at once and in a much more efficient way.
   */
  template <int order>
  Tensor<order, dim>
  compute_derivative(const unsigned int i, const Point<dim> &p) const;

  /**
   * @copydoc ScalarPolynomialsBase::compute_1st_derivative()
   */
  virtual Tensor<1, dim>
  compute_1st_derivative(const unsigned int i,
                         const Point<dim>  &p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_2nd_derivative()
   */
  virtual Tensor<2, dim>
  compute_2nd_derivative(const unsigned int i,
                         const Point<dim>  &p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_3rd_derivative()
   */
  virtual Tensor<3, dim>
  compute_3rd_derivative(const unsigned int i,
                         const Point<dim>  &p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_4th_derivative()
   */
  virtual Tensor<4, dim>
  compute_4th_derivative(const unsigned int i,
                         const Point<dim>  &p) const override;

  /**
   * Compute the grad of the <tt>i</tt>th tensor product polynomial at
   * <tt>unit_point</tt>. Here <tt>i</tt> is given in tensor product
   * numbering.
   *
   * Note, that using this function within a loop over all tensor product
   * polynomials is not efficient, because then each derivative value of the
   * underlying (one-dimensional) polynomials is (unnecessarily) computed
   * several times.  Instead use the evaluate() function, see above, with
   * <tt>grads.size()==</tt>n() to get the point value of all tensor
   * polynomials all at once and in a much more efficient way.
   */
  Tensor<1, dim>
  compute_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * Compute the second derivative (grad_grad) of the <tt>i</tt>th tensor
   * product polynomial at <tt>unit_point</tt>. Here <tt>i</tt> is given in
   * tensor product numbering.
   *
   * Note, that using this function within a loop over all tensor product
   * polynomials is not efficient, because then each derivative value of the
   * underlying (one-dimensional) polynomials is (unnecessarily) computed
   * several times.  Instead use the evaluate() function, see above, with
   * <tt>grad_grads.size()==</tt>n() to get the point value of all tensor
   * polynomials all at once and in a much more efficient way.
   */
  Tensor<2, dim>
  compute_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * Return the number of tensor product polynomials plus the bubble
   * enrichments. For <i>n</i> 1d polynomials this is <i>n<sup>dim</sup>+1</i>
   * if the maximum degree of the polynomials is one and
   * <i>n<sup>dim</sup>+dim</i> otherwise.
   */
  unsigned int
  n() const;

  /**
   * Return the name of the space, which is
   * <tt>TensorProductPolynomialsBubbles</tt>.
   */
  std::string
  name() const override;

  /**
   * @copydoc ScalarPolynomialsBase::clone()
   */
  virtual std::unique_ptr<ScalarPolynomialsBase<dim>>
  clone() const override;

private:
  /**
   * The TensorProductPolynomials object
   */
  TensorProductPolynomials<dim> tensor_polys;

  /**
   * Index map for reordering the polynomials.
   */
  std::vector<unsigned int> index_map;

  /**
   * Index map for reordering the polynomials.
   */
  std::vector<unsigned int> index_map_inverse;
};

/** @} */


/* ---------------- template and inline functions ---------- */

#ifndef DOXYGEN

template <int dim>
template <class Pol>
inline TensorProductPolynomialsBubbles<dim>::TensorProductPolynomialsBubbles(
  const std::vector<Pol> &pols)
  : ScalarPolynomialsBase<dim>(1,
                               Utilities::fixed_power<dim>(pols.size()) + dim)
  , tensor_polys(pols)
  , index_map(tensor_polys.n() +
              ((tensor_polys.polynomials.size() <= 2) ? 1 : dim))
  , index_map_inverse(tensor_polys.n() +
                      ((tensor_polys.polynomials.size() <= 2) ? 1 : dim))
{
  const unsigned int q_degree  = tensor_polys.polynomials.size() - 1;
  const unsigned int n_bubbles = ((q_degree <= 1) ? 1 : dim);
  // append index for renumbering
  for (unsigned int i = 0; i < tensor_polys.n() + n_bubbles; ++i)
    {
      index_map[i]         = i;
      index_map_inverse[i] = i;
    }
}


template <int dim>
inline unsigned int
TensorProductPolynomialsBubbles<dim>::n() const
{
  return tensor_polys.n() + dim;
}


template <>
inline unsigned int
TensorProductPolynomialsBubbles<0>::n() const
{
  return numbers::invalid_unsigned_int;
}


template <int dim>
inline const std::vector<unsigned int> &
TensorProductPolynomialsBubbles<dim>::get_numbering() const
{
  return index_map;
}


template <int dim>
inline const std::vector<unsigned int> &
TensorProductPolynomialsBubbles<dim>::get_numbering_inverse() const
{
  return index_map_inverse;
}


template <int dim>
inline std::string
TensorProductPolynomialsBubbles<dim>::name() const
{
  return "TensorProductPolynomialsBubbles";
}


template <int dim>
template <int order>
Tensor<order, dim>
TensorProductPolynomialsBubbles<dim>::compute_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  const unsigned int q_degree      = tensor_polys.polynomials.size() - 1;
  const unsigned int max_q_indices = tensor_polys.n();
  Assert(i < max_q_indices + /* n_bubbles= */ ((q_degree <= 1) ? 1 : dim),
         ExcInternalError());

  // treat the regular basis functions
  if (i < max_q_indices)
    return tensor_polys.template compute_derivative<order>(i, p);

  [[maybe_unused]] const unsigned int comp = i - tensor_polys.n();

  if constexpr (order == 1)
    {
      Tensor<1, dim> derivative;
      for (unsigned int d = 0; d < dim; ++d)
        {
          derivative[d] = 1.;
          // compute grad(4*\prod_{i=1}^d (x_i(1-x_i)))(p)
          for (unsigned j = 0; j < dim; ++j)
            derivative[d] *=
              (d == j ? 4 * (1 - 2 * p[j]) : 4 * p[j] * (1 - p[j]));
          // and multiply with (2*x_i-1)^{r-1}
          for (unsigned int i = 0; i < q_degree - 1; ++i)
            derivative[d] *= 2 * p[comp] - 1;
        }

      if (q_degree >= 2)
        {
          // add \prod_{i=1}^d 4*(x_i(1-x_i))(p)
          double value = 1.;
          for (unsigned int j = 0; j < dim; ++j)
            value *= 4 * p[j] * (1 - p[j]);
          // and multiply with grad(2*x_i-1)^{r-1}
          double tmp = value * 2 * (q_degree - 1);
          for (unsigned int i = 0; i < q_degree - 2; ++i)
            tmp *= 2 * p[comp] - 1;
          derivative[comp] += tmp;
        }

      return derivative;
    }
  else if constexpr (order == 2)
    {
      Tensor<2, dim> derivative;

      double v[dim + 1][3];
      {
        for (unsigned int c = 0; c < dim; ++c)
          {
            v[c][0] = 4 * p[c] * (1 - p[c]);
            v[c][1] = 4 * (1 - 2 * p[c]);
            v[c][2] = -8;
          }

        double tmp = 1.;
        for (unsigned int i = 0; i < q_degree - 1; ++i)
          tmp *= 2 * p[comp] - 1;
        v[dim][0] = tmp;

        if (q_degree >= 2)
          {
            double tmp = 2 * (q_degree - 1);
            for (unsigned int i = 0; i < q_degree - 2; ++i)
              tmp *= 2 * p[comp] - 1;
            v[dim][1] = tmp;
          }
        else
          v[dim][1] = 0.;

        if (q_degree >= 3)
          {
            double tmp = 4 * (q_degree - 2) * (q_degree - 1);
            for (unsigned int i = 0; i < q_degree - 3; ++i)
              tmp *= 2 * p[comp] - 1;
            v[dim][2] = tmp;
          }
        else
          v[dim][2] = 0.;
      }

      // calculate (\partial_j \partial_k \psi) * monomial
      Tensor<2, dim> grad_grad_1;
      for (unsigned int d1 = 0; d1 < dim; ++d1)
        for (unsigned int d2 = 0; d2 < dim; ++d2)
          {
            grad_grad_1[d1][d2] = v[dim][0];
            for (unsigned int x = 0; x < dim; ++x)
              {
                unsigned int derivative = 0;
                if (d1 == x || d2 == x)
                  {
                    if (d1 == d2)
                      derivative = 2;
                    else
                      derivative = 1;
                  }
                grad_grad_1[d1][d2] *= v[x][derivative];
              }
          }

      // calculate (\partial_j  \psi) *(\partial_k monomial)
      // and (\partial_k  \psi) *(\partial_j monomial)
      Tensor<2, dim> grad_grad_2;
      Tensor<2, dim> grad_grad_3;
      for (unsigned int d = 0; d < dim; ++d)
        {
          grad_grad_2[d][comp] = v[dim][1];
          grad_grad_3[comp][d] = v[dim][1];
          for (unsigned int x = 0; x < dim; ++x)
            {
              grad_grad_2[d][comp] *= v[x][d == x];
              grad_grad_3[comp][d] *= v[x][d == x];
            }
        }

      // calculate \psi *(\partial j \partial_k monomial) and sum
      double psi_value = 1.;
      for (unsigned int x = 0; x < dim; ++x)
        psi_value *= v[x][0];

      for (unsigned int d1 = 0; d1 < dim; ++d1)
        for (unsigned int d2 = 0; d2 < dim; ++d2)
          derivative[d1][d2] =
            grad_grad_1[d1][d2] + grad_grad_2[d1][d2] + grad_grad_3[d1][d2];
      derivative[comp][comp] += psi_value * v[dim][2];

      return derivative;
    }
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
      return {};
    }
}



template <int dim>
inline Tensor<1, dim>
TensorProductPolynomialsBubbles<dim>::compute_1st_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  return compute_derivative<1>(i, p);
}



template <int dim>
inline Tensor<2, dim>
TensorProductPolynomialsBubbles<dim>::compute_2nd_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  return compute_derivative<2>(i, p);
}



template <int dim>
inline Tensor<3, dim>
TensorProductPolynomialsBubbles<dim>::compute_3rd_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  return compute_derivative<3>(i, p);
}



template <int dim>
inline Tensor<4, dim>
TensorProductPolynomialsBubbles<dim>::compute_4th_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  return compute_derivative<4>(i, p);
}

#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif

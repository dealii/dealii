// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2000 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_tensor_product_polynomials_h
#define dealii_tensor_product_polynomials_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/ndarray.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/scalar_polynomials_base.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations for friends
// TODO: We may be able to modify these classes so they aren't
// required to be friends
template <int dim>
class TensorProductPolynomialsBubbles;
template <int dim>
class TensorProductPolynomialsConst;

/**
 * @addtogroup Polynomials
 * @{
 */

/**
 * Tensor product of given polynomials.
 *
 * Given a vector of <i>n</i> one-dimensional polynomials <i>P<sub>1</sub></i>
 * to <i>P<sub>n</sub></i>, this class generates <i>n<sup>dim</sup></i>
 * polynomials of the form <i>Q<sub>ijk</sub>(x,y,z) =
 * P<sub>i</sub>(x)P<sub>j</sub>(y)P<sub>k</sub>(z)</i>. If the base
 * polynomials are mutually orthogonal on the interval [-1,1] or [0,1], then
 * the tensor product polynomials are orthogonal on [-1,1]<sup>dim</sup> or
 * [0,1]<sup>dim</sup>, respectively.
 *
 * Indexing is as follows: the order of dim-dimensional polynomials is
 * x-coordinates running fastest, then y-coordinate, etc. The first few
 * polynomials are thus <i>P<sub>1</sub>(x)P<sub>1</sub>(y),
 * P<sub>2</sub>(x)P<sub>1</sub>(y), P<sub>3</sub>(x)P<sub>1</sub>(y), ...,
 * P<sub>1</sub>(x)P<sub>2</sub>(y), P<sub>2</sub>(x)P<sub>2</sub>(y),
 * P<sub>3</sub>(x)P<sub>2</sub>(y), ...</i> and likewise in 3d.
 *
 * The output_indices() function prints the ordering of the dim-dimensional
 * polynomials, i.e. for each polynomial in the polynomial space it gives the
 * indices i,j,k of the one-dimensional polynomials in x,y and z direction.
 * The ordering of the dim-dimensional polynomials can be changed by using the
 * set_numbering() function.
 *
 * @tparam PolynomialType A class that satisfies the required interface for computing
 *   tensor products. Typical choices for this template argument are
 *   Polynomials::Polynomial and Polynomials::PiecewisePolynomial.
 */
template <int dim, typename PolynomialType = Polynomials::Polynomial<double>>
class TensorProductPolynomials : public ScalarPolynomialsBase<dim>
{
public:
  /**
   * Access to the dimension of this object, for checking and automatic
   * setting of dimension in other classes.
   */
  static constexpr unsigned int dimension = dim;

  /**
   * Constructor. <tt>pols</tt> is a vector of objects that should be derived
   * or otherwise convertible to one-dimensional polynomial objects of type
   * `PolynomialType` (template argument of class). It will be copied element
   * by element into a protected member variable.
   */
  template <class Pol>
  TensorProductPolynomials(const std::vector<Pol> &pols);

  /**
   * Print the list of the indices to <tt>out</tt>.
   */
  void
  output_indices(std::ostream &out) const;

  /**
   * Set the ordering of the polynomials. Requires
   * <tt>renumber.size()==n()</tt>.  Stores a copy of <tt>renumber</tt>.
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
   * Compute the <tt>order</tt>th derivative of the <tt>i</tt>th tensor
   * product polynomial at <tt>unit_point</tt>. Here <tt>i</tt> is given in
   * tensor product numbering.
   *
   * Note, that using this function within a loop over all tensor product
   * polynomials is not efficient, because then each derivative value of the
   * underlying (one-dimensional) polynomials is (unnecessarily) computed
   * several times.  Instead use the evaluate() function, see above, with the
   * size of the appropriate parameter set to n() to get the point value of
   * all tensor polynomials all at once and in a much more efficient way.
   *
   * @tparam order The derivative order.
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
   * Return the name of the space, which is <tt>TensorProductPolynomials</tt>.
   */
  std::string
  name() const override;

  /**
   * @copydoc ScalarPolynomialsBase::clone()
   */
  virtual std::unique_ptr<ScalarPolynomialsBase<dim>>
  clone() const override;

  /**
   * Return an estimate (in bytes) for the memory consumption of this object.
   */
  virtual std::size_t
  memory_consumption() const override;

  /**
   * Return a copy of the underlying one-dimensional polynomials given to the
   * constructor of this class.
   */
  std::vector<PolynomialType>
  get_underlying_polynomials() const;

protected:
  /**
   * Copy of the vector <tt>pols</tt> of polynomials given to the constructor.
   */
  std::vector<PolynomialType> polynomials;

  /**
   * Index map for reordering the polynomials.
   */
  std::vector<unsigned int> index_map;

  /**
   * Index map for reordering the polynomials.
   */
  std::vector<unsigned int> index_map_inverse;

  /**
   * Each tensor product polynomial <i>i</i> is a product of one-dimensional
   * polynomials in each space direction. Compute the indices of these
   * one-dimensional polynomials for each space direction, given the index
   * <i>i</i>.
   */
  void
  compute_index(const unsigned int             i,
                std::array<unsigned int, dim> &indices) const;

  /**
   * TensorProductPolynomialsBubbles has a TensorProductPolynomials class
   * so we declare it as a friend class.
   */
  friend class TensorProductPolynomialsBubbles<dim>;

  /**
   * TensorProductPolynomialsConst has a TensorProductPolynomials class
   * so we declare it as a friend class.
   */
  friend class TensorProductPolynomialsConst<dim>;
};



/**
 * Anisotropic tensor product of given polynomials.
 *
 * Given one-dimensional polynomials $P^x_1(x), P^x_2(x), \ldots$ in
 * $x$-direction, $P^y_1(y), P^y_2(y), \ldots$ in $y$-direction, and
 * so on, this class generates polynomials of the form $Q_{ijk}(x,y,z)
 * = P^x_i(x)P^y_j(y)P^z_k(z)$. (With obvious generalization if @p dim
 * is in fact only 2. If @p dim is in fact only 1, then the result is
 * simply the same set of one-dimensional polynomials passed to the
 * constructor.)
 *
 * If the elements of each set of base polynomials are mutually
 * orthogonal on the interval $[-1,1]$ or $[0,1]$, then the tensor
 * product polynomials are orthogonal on $[-1,1]^d$ or $[0,1]^d$,
 * respectively.
 *
 * The resulting @p dim-dimensional tensor product polynomials are
 * ordered as follows: We iterate over the $x$ coordinates running
 * fastest, then the $y$ coordinate, etc. For example, for @p dim==2,
 * the first few polynomials are thus
 * $P^x_1(x)P^y_1(y)$,
 * $P^x_2(x)P^y_1(y)$, $P^x_3(x)P^y_1(y)$, ...,
 * $P^x_1(x)P^y_2(y)$, $P^x_2(x)P^y_2(y)$,
 * $P^x_3(x)P^y_2(y)$, etc.
 */
template <int dim>
class AnisotropicPolynomials : public ScalarPolynomialsBase<dim>
{
public:
  /**
   * Constructor. @p base_polynomials is a table of one-dimensional
   * polynomials. The number of rows in this table (the first index
   * when indexing into @p base_polynomials) needs to be equal to the
   * space dimension, with the elements of each row (i.e., the second
   * index) giving the polynomials that shall be used in this
   * particular coordinate direction.
   *
   * Since we want to build <i>anisotropic</i> polynomials, the @p dim
   * sets of polynomials passed in as arguments may of course be
   * different, and may also vary in number.
   *
   * The number of tensor product polynomials is <tt>Nx*Ny*Nz</tt>, or with
   * terms dropped if the number of space dimensions is less than 3.
   */
  AnisotropicPolynomials(
    const std::vector<std::vector<Polynomials::Polynomial<double>>>
      &base_polynomials);

  /**
   * Set the ordering of the polynomials. Requires
   * <tt>renumber.size()==n()</tt>.  Stores a copy of <tt>renumber</tt>.
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
   * The size of the vectors must either be equal <tt>0</tt> or equal
   * <tt>this->n()</tt>.  In the first case, the function will not compute
   * these values.
   *
   * If you need values or derivatives of all tensor product polynomials then
   * use this function, rather than using any of the <tt>compute_value</tt>,
   * <tt>compute_grad</tt> or <tt>compute_grad_grad</tt> functions, see below,
   * in a loop over all tensor product polynomials.
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
   * several times.  Instead use the <tt>compute</tt> function, see above,
   * with <tt>values.size()==this->n()</tt> to get the point values of all
   * tensor polynomials all at once and in a much more efficient way.
   */
  double
  compute_value(const unsigned int i, const Point<dim> &p) const override;

  /**
   * Compute the <tt>order</tt>th derivative of the <tt>i</tt>th tensor
   * product polynomial at <tt>unit_point</tt>. Here <tt>i</tt> is given in
   * tensor product numbering.
   *
   * Note, that using this function within a loop over all tensor product
   * polynomials is not efficient, because then each derivative value of the
   * underlying (one-dimensional) polynomials is (unnecessarily) computed
   * several times.  Instead use the evaluate() function, see above, with the
   * size of the appropriate parameter set to n() to get the point value of
   * all tensor polynomials all at once and in a much more efficient way.
   *
   * @tparam order The derivative order.
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
   * several times.  Instead use the <tt>compute</tt> function, see above,
   * with <tt>grads.size()==this->n()</tt> to get the point value of all
   * tensor polynomials all at once and in a much more efficient way.
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
   * several times.  Instead use the <tt>compute</tt> function, see above,
   * with <tt>grad_grads.size()==this->n()</tt> to get the point value of
   * all tensor polynomials all at once and in a much more efficient way.
   */
  Tensor<2, dim>
  compute_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * Return the name of the space, which is <tt>AnisotropicPolynomials</tt>.
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
   * Copy of the vector <tt>pols</tt> of polynomials given to the constructor.
   */
  const std::vector<std::vector<Polynomials::Polynomial<double>>> polynomials;

  /**
   * Index map for reordering the polynomials.
   */
  std::vector<unsigned int> index_map;

  /**
   * Index map for reordering the polynomials.
   */
  std::vector<unsigned int> index_map_inverse;

  /**
   * Each tensor product polynomial $p_i$ is a product of one-dimensional
   * polynomials in each space direction. Compute the indices of these
   * one-dimensional polynomials for each space direction, given the index
   * <tt>i</tt>.
   */
  void
  compute_index(const unsigned int             i,
                std::array<unsigned int, dim> &indices) const;

  /**
   * Given the input to the constructor, compute <tt>n_pols</tt>.
   */
  static unsigned int
  get_n_tensor_pols(
    const std::vector<std::vector<Polynomials::Polynomial<double>>> &pols);
};

/** @} */

#ifndef DOXYGEN


/* ---------------- template and inline functions ---------- */


template <int dim, typename PolynomialType>
template <class Pol>
inline TensorProductPolynomials<dim, PolynomialType>::TensorProductPolynomials(
  const std::vector<Pol> &pols)
  : ScalarPolynomialsBase<dim>(1, Utilities::fixed_power<dim>(pols.size()))
  , polynomials(pols.begin(), pols.end())
  , index_map(this->n())
  , index_map_inverse(this->n())
{
  // per default set this index map to identity. This map can be changed by
  // the user through the set_numbering() function
  for (unsigned int i = 0; i < this->n(); ++i)
    {
      index_map[i]         = i;
      index_map_inverse[i] = i;
    }
}


template <int dim, typename PolynomialType>
inline const std::vector<unsigned int> &
TensorProductPolynomials<dim, PolynomialType>::get_numbering() const
{
  return index_map;
}


template <int dim, typename PolynomialType>
inline const std::vector<unsigned int> &
TensorProductPolynomials<dim, PolynomialType>::get_numbering_inverse() const
{
  return index_map_inverse;
}


template <int dim, typename PolynomialType>
inline std::string
TensorProductPolynomials<dim, PolynomialType>::name() const
{
  return "TensorProductPolynomials";
}


template <int dim, typename PolynomialType>
template <int order>
Tensor<order, dim>
TensorProductPolynomials<dim, PolynomialType>::compute_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  std::array<unsigned int, dim> indices;
  compute_index(i, indices);

  ndarray<double, dim, 5> v;
  {
    std::vector<double> tmp(5);
    for (unsigned int d = 0; d < dim; ++d)
      {
        polynomials[indices[d]].value(p[d], tmp);
        v[d][0] = tmp[0];
        v[d][1] = tmp[1];
        v[d][2] = tmp[2];
        v[d][3] = tmp[3];
        v[d][4] = tmp[4];
      }
  }

  if constexpr (order == 1)
    {
      Tensor<1, dim> derivative;
      for (unsigned int d = 0; d < dim; ++d)
        {
          derivative[d] = 1.;
          for (unsigned int x = 0; x < dim; ++x)
            {
              unsigned int x_order = 0;
              if (d == x)
                ++x_order;

              derivative[d] *= v[x][x_order];
            }
        }

      return derivative;
    }
  else if constexpr (order == 2)
    {
      Tensor<2, dim> derivative;
      for (unsigned int d1 = 0; d1 < dim; ++d1)
        for (unsigned int d2 = 0; d2 < dim; ++d2)
          {
            derivative[d1][d2] = 1.;
            for (unsigned int x = 0; x < dim; ++x)
              {
                unsigned int x_order = 0;
                if (d1 == x)
                  ++x_order;
                if (d2 == x)
                  ++x_order;

                derivative[d1][d2] *= v[x][x_order];
              }
          }

      return derivative;
    }
  else if constexpr (order == 3)
    {
      Tensor<3, dim> derivative;
      for (unsigned int d1 = 0; d1 < dim; ++d1)
        for (unsigned int d2 = 0; d2 < dim; ++d2)
          for (unsigned int d3 = 0; d3 < dim; ++d3)
            {
              derivative[d1][d2][d3] = 1.;
              for (unsigned int x = 0; x < dim; ++x)
                {
                  unsigned int x_order = 0;
                  if (d1 == x)
                    ++x_order;
                  if (d2 == x)
                    ++x_order;
                  if (d3 == x)
                    ++x_order;

                  derivative[d1][d2][d3] *= v[x][x_order];
                }
            }

      return derivative;
    }
  else if constexpr (order == 4)
    {
      Tensor<4, dim> derivative;
      for (unsigned int d1 = 0; d1 < dim; ++d1)
        for (unsigned int d2 = 0; d2 < dim; ++d2)
          for (unsigned int d3 = 0; d3 < dim; ++d3)
            for (unsigned int d4 = 0; d4 < dim; ++d4)
              {
                derivative[d1][d2][d3][d4] = 1.;
                for (unsigned int x = 0; x < dim; ++x)
                  {
                    unsigned int x_order = 0;
                    if (d1 == x)
                      ++x_order;
                    if (d2 == x)
                      ++x_order;
                    if (d3 == x)
                      ++x_order;
                    if (d4 == x)
                      ++x_order;

                    derivative[d1][d2][d3][d4] *= v[x][x_order];
                  }
              }

      return derivative;
    }
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
      return {};
    }
}



template <>
template <int order>
Tensor<order, 0>
TensorProductPolynomials<0, Polynomials::Polynomial<double>>::
  compute_derivative(const unsigned int, const Point<0> &) const
{
  AssertThrow(false, ExcNotImplemented());

  return {};
}



template <int dim, typename PolynomialType>
inline Tensor<1, dim>
TensorProductPolynomials<dim, PolynomialType>::compute_1st_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  return compute_derivative<1>(i, p);
}



template <int dim, typename PolynomialType>
inline Tensor<2, dim>
TensorProductPolynomials<dim, PolynomialType>::compute_2nd_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  return compute_derivative<2>(i, p);
}



template <int dim, typename PolynomialType>
inline Tensor<3, dim>
TensorProductPolynomials<dim, PolynomialType>::compute_3rd_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  return compute_derivative<3>(i, p);
}



template <int dim, typename PolynomialType>
inline Tensor<4, dim>
TensorProductPolynomials<dim, PolynomialType>::compute_4th_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  return compute_derivative<4>(i, p);
}



template <int dim>
template <int order>
Tensor<order, dim>
AnisotropicPolynomials<dim>::compute_derivative(const unsigned int i,
                                                const Point<dim>  &p) const
{
  std::array<unsigned int, dim> indices;
  compute_index(i, indices);

  std::vector<std::vector<double>> v(dim, std::vector<double>(order + 1));
  for (unsigned int d = 0; d < dim; ++d)
    polynomials[d][indices[d]].value(p[d], v[d]);

  if constexpr (order == 1)
    {
      Tensor<1, dim> derivative;
      for (unsigned int d = 0; d < dim; ++d)
        {
          derivative[d] = 1.;
          for (unsigned int x = 0; x < dim; ++x)
            {
              unsigned int x_order = 0;
              if (d == x)
                ++x_order;

              derivative[d] *= v[x][x_order];
            }
        }

      return derivative;
    }
  else if constexpr (order == 2)
    {
      Tensor<2, dim> derivative;
      for (unsigned int d1 = 0; d1 < dim; ++d1)
        for (unsigned int d2 = 0; d2 < dim; ++d2)
          {
            derivative[d1][d2] = 1.;
            for (unsigned int x = 0; x < dim; ++x)
              {
                unsigned int x_order = 0;
                if (d1 == x)
                  ++x_order;
                if (d2 == x)
                  ++x_order;

                derivative[d1][d2] *= v[x][x_order];
              }
          }

      return derivative;
    }
  else if constexpr (order == 3)
    {
      Tensor<3, dim> derivative;
      for (unsigned int d1 = 0; d1 < dim; ++d1)
        for (unsigned int d2 = 0; d2 < dim; ++d2)
          for (unsigned int d3 = 0; d3 < dim; ++d3)
            {
              derivative[d1][d2][d3] = 1.;
              for (unsigned int x = 0; x < dim; ++x)
                {
                  unsigned int x_order = 0;
                  if (d1 == x)
                    ++x_order;
                  if (d2 == x)
                    ++x_order;
                  if (d3 == x)
                    ++x_order;

                  derivative[d1][d2][d3] *= v[x][x_order];
                }
            }

      return derivative;
    }
  else if constexpr (order == 4)
    {
      Tensor<4, dim> derivative;
      for (unsigned int d1 = 0; d1 < dim; ++d1)
        for (unsigned int d2 = 0; d2 < dim; ++d2)
          for (unsigned int d3 = 0; d3 < dim; ++d3)
            for (unsigned int d4 = 0; d4 < dim; ++d4)
              {
                derivative[d1][d2][d3][d4] = 1.;
                for (unsigned int x = 0; x < dim; ++x)
                  {
                    unsigned int x_order = 0;
                    if (d1 == x)
                      ++x_order;
                    if (d2 == x)
                      ++x_order;
                    if (d3 == x)
                      ++x_order;
                    if (d4 == x)
                      ++x_order;

                    derivative[d1][d2][d3][d4] *= v[x][x_order];
                  }
              }

      return derivative;
    }
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
      return {};
    }
}



template <>
template <int order>
Tensor<order, 0>
AnisotropicPolynomials<0>::compute_derivative(const unsigned int,
                                              const Point<0> &) const
{
  AssertThrow(false, ExcNotImplemented());

  return {};
}



template <int dim>
inline Tensor<1, dim>
AnisotropicPolynomials<dim>::compute_1st_derivative(const unsigned int i,
                                                    const Point<dim>  &p) const
{
  return compute_derivative<1>(i, p);
}



template <int dim>
inline Tensor<2, dim>
AnisotropicPolynomials<dim>::compute_2nd_derivative(const unsigned int i,
                                                    const Point<dim>  &p) const
{
  return compute_derivative<2>(i, p);
}



template <int dim>
inline Tensor<3, dim>
AnisotropicPolynomials<dim>::compute_3rd_derivative(const unsigned int i,
                                                    const Point<dim>  &p) const
{
  return compute_derivative<3>(i, p);
}



template <int dim>
inline Tensor<4, dim>
AnisotropicPolynomials<dim>::compute_4th_derivative(const unsigned int i,
                                                    const Point<dim>  &p) const
{
  return compute_derivative<4>(i, p);
}



template <int dim>
inline std::string
AnisotropicPolynomials<dim>::name() const
{
  return "AnisotropicPolynomials";
}



#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif

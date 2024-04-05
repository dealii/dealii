// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_simplex_barycentric_polynomials_h
#define dealii_simplex_barycentric_polynomials_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/scalar_polynomials_base.h>
#include <deal.II/base/table.h>

#include <limits>

DEAL_II_NAMESPACE_OPEN

/**
 * Polynomial implemented in barycentric coordinates.
 *
 * Barycentric coordinates are a coordinate system defined on simplices that
 * are particularly easy to work with since they express coordinates in the
 * simplex as convex combinations of the vertices. For example, any point in a
 * triangle can be written as
 *
 * @f[
 *   (x, y) = c_0 (x_0, y_0) + c_1 (x_1, y_1) + c_2 (x_2, y_2).
 * @f]
 *
 * where each value $c_i$ is the relative weight of each vertex (so the
 * centroid is, in 2d, where each $c_i = 1/3$). Since we only consider convex
 * combinations we can rewrite this equation as
 *
 * @f[
 *   (x, y) = (1 - c_1 - c_2) (x_0, y_0) + c_1 (x_1, y_1) + c_2 (x_2, y_2).
 * @f]
 *
 * This results in three polynomials that are equivalent to $P^1$ in 2d. More
 * exactly, this class implements a polynomial space defined with the basis,
 * in 2d, of
 * @f{align*}{
 * t_0(x, y) &= 1 - x - y \\
 * t_1(x, y) &= x \\
 * t_2(x, y) &= y
 * @f}
 * and, in 3d,
 * @f{align*}{
 * t_0(x, y) &= 1 - x - y - z \\
 * t_1(x, y) &= x             \\
 * t_2(x, y) &= y             \\
 * t_2(x, y) &= z
 * @f}
 *
 * which is, in practice, a very convenient basis for defining simplex
 * polynomials: for example, the fourth basis function of a TRI6 element is
 *
 * @f[
 * 4 * t_1(x, y) * t_2(x, y).
 * @f]
 *
 * Barycentric polynomials in <code>dim</code>-dimensional space have
 * <code>dim + 1</code> variables in since <code>t_0</code> can be written in
 * terms of the other monomials.
 *
 * Monomials can be conveniently constructed with
 * BarycentricPolynomial::monomial().
 *
 * @ingroup Polynomials
 */
template <int dim, typename Number = double>
class BarycentricPolynomial
{
public:
  /**
   * Constructor for the zero polynomial.
   */
  BarycentricPolynomial();

  /**
   * Constructor for a monomial.
   */
  BarycentricPolynomial(const TableIndices<dim + 1> &powers,
                        const Number                 coefficient);

  /**
   * Return the specified monomial.
   */
  static BarycentricPolynomial<dim, Number>
  monomial(const unsigned int d);

  /**
   * Print the polynomial to the output stream with lowest-order terms first.
   * For example, the first P6 basis function is printed as
   * <code>-1 * t0^1 + 2 * t0^2</code>, where <code>t0</code> is the first
   * barycentric variable, <code>t1</code> is the second, etc.
   */
  void
  print(std::ostream &out) const;

  /**
   * Degree of each barycentric polynomial.
   */
  TableIndices<dim + 1>
  degrees() const;

  /**
   * Unary minus.
   */
  BarycentricPolynomial<dim, Number>
  operator-() const;

  /**
   * Add a scalar.
   */
  template <typename Number2>
  BarycentricPolynomial<dim, Number>
  operator+(const Number2 &a) const;

  /**
   * Subtract a scalar.
   */
  template <typename Number2>
  BarycentricPolynomial<dim, Number>
  operator-(const Number2 &a) const;

  /**
   * Multiply by a scalar.
   */
  template <typename Number2>
  BarycentricPolynomial<dim, Number>
  operator*(const Number2 &a) const;

  /**
   * Divide by a scalar.
   */
  template <typename Number2>
  BarycentricPolynomial<dim, Number>
  operator/(const Number2 &a) const;

  /**
   * Add another barycentric polynomial.
   */
  BarycentricPolynomial<dim, Number>
  operator+(const BarycentricPolynomial<dim, Number> &augend) const;

  /**
   * Subtract another barycentric polynomial.
   */
  BarycentricPolynomial<dim, Number>
  operator-(const BarycentricPolynomial<dim, Number> &augend) const;

  /**
   * Multiply by another barycentric polynomial.
   */
  BarycentricPolynomial<dim, Number>
  operator*(const BarycentricPolynomial<dim, Number> &multiplicand) const;

  /**
   * Differentiate in barycentric coordinates.
   */
  BarycentricPolynomial<dim, Number>
  barycentric_derivative(const unsigned int coordinate) const;

  /**
   * Differentiate in Cartesian coordinates.
   */
  BarycentricPolynomial<dim, Number>
  derivative(const unsigned int coordinate) const;

  /**
   * Evaluate the polynomial.
   */
  Number
  value(const Point<dim> &point) const;

  /**
   * Return an estimate, in bytes, of the memory usage of the object.
   */
  std::size_t
  memory_consumption() const;

protected:
  /**
   * Coefficients of the polynomial. The exponents are the integer indexes.
   */
  Table<dim + 1, Number> coefficients;

  /**
   * Utility function for barycentric polynomials: it is convenient to loop
   * over all the indices at once in a dimension-independent way, but we also
   * need to access the actual indices of the underlying Table object. This
   * utility function converts an integral index into the equivalent
   * TableIndices array (which are also the implicitly stored polynomial
   * exponents).
   */
  static TableIndices<dim + 1>
  index_to_indices(const std::size_t           &index,
                   const TableIndices<dim + 1> &extents);
};

/**
 * Scalar polynomial space based on barycentric polynomials.
 */
template <int dim>
class BarycentricPolynomials : public ScalarPolynomialsBase<dim>
{
public:
  /**
   * Alias for polynomial type.
   */
  using PolyType = BarycentricPolynomial<dim>;

  /**
   * Alias for polynomial gradient type.
   */
  using GradType = std::array<PolyType, dim>;

  /**
   * Alias for polynomial hessian type.
   */
  using HessianType = std::array<GradType, dim>;

  /**
   * Alias for polynomial third derivatives type.
   */
  using ThirdDerivativesType = std::array<HessianType, dim>;

  /**
   * Alias for polynomial fourth derivatives type.
   */
  using FourthDerivativesType = std::array<ThirdDerivativesType, dim>;

  /**
   * Make the dimension available to the outside.
   */
  static constexpr unsigned int dimension = dim;

  /**
   * Get the standard Lagrange basis for a specified degree.
   */
  static BarycentricPolynomials<dim>
  get_fe_p_basis(const unsigned int degree);

  /**
   * Constructor taking the polynomial @p degree as input.
   */
  BarycentricPolynomials(
    const std::vector<BarycentricPolynomial<dim>> &polynomials);

  /**
   * Access operator.
   */
  const BarycentricPolynomial<dim> &
  operator[](const std::size_t i) const;

  /**
   * @copydoc ScalarPolynomialsBase::evaluate()
   */
  void
  evaluate(const Point<dim>            &unit_point,
           std::vector<double>         &values,
           std::vector<Tensor<1, dim>> &grads,
           std::vector<Tensor<2, dim>> &grad_grads,
           std::vector<Tensor<3, dim>> &third_derivatives,
           std::vector<Tensor<4, dim>> &fourth_derivatives) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_value()
   */
  double
  compute_value(const unsigned int i, const Point<dim> &p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_1st_derivative()
   */
  Tensor<1, dim>
  compute_1st_derivative(const unsigned int i,
                         const Point<dim>  &p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_2nd_derivative()
   */
  Tensor<2, dim>
  compute_2nd_derivative(const unsigned int i,
                         const Point<dim>  &p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_3rd_derivative()
   */
  Tensor<3, dim>
  compute_3rd_derivative(const unsigned int i,
                         const Point<dim>  &p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_4th_derivative()
   */
  Tensor<4, dim>
  compute_4th_derivative(const unsigned int i,
                         const Point<dim>  &p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_grad()
   */
  Tensor<1, dim>
  compute_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_grad_grad()
   */
  Tensor<2, dim>
  compute_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::memory_consumption()
   */
  virtual std::size_t
  memory_consumption() const override;

  /**
   * @copydoc ScalarPolynomialsBase::name()
   */
  std::string
  name() const override;

  /**
   * @copydoc ScalarPolynomialsBase::clone()
   */
  virtual std::unique_ptr<ScalarPolynomialsBase<dim>>
  clone() const override;

protected:
  std::vector<PolyType>              polys;
  std::vector<GradType>              poly_grads;
  std::vector<HessianType>           poly_hessians;
  std::vector<ThirdDerivativesType>  poly_third_derivatives;
  std::vector<FourthDerivativesType> poly_fourth_derivatives;
};

// non-member template functions for algebra

/**
 * Multiply a BarycentricPolynomial by a constant.
 */
template <int dim, typename Number1, typename Number2>
BarycentricPolynomial<dim, Number1>
operator*(const Number2 &a, const BarycentricPolynomial<dim, Number1> &bp)
{
  return bp * Number1(a);
}

/**
 * Add a constant to a BarycentricPolynomial.
 */
template <int dim, typename Number1, typename Number2>
BarycentricPolynomial<dim, Number1>
operator+(const Number2 &a, const BarycentricPolynomial<dim, Number1> &bp)
{
  return bp + Number1(a);
}

/**
 * Subtract a BarycentricPolynomial from a constant.
 */
template <int dim, typename Number1, typename Number2>
BarycentricPolynomial<dim, Number1>
operator-(const Number2 &a, const BarycentricPolynomial<dim, Number1> &bp)
{
  return bp - Number1(a);
}

/**
 * Write a BarycentricPolynomial to the provided output stream.
 */
template <int dim, typename Number>
std::ostream &
operator<<(std::ostream &out, const BarycentricPolynomial<dim, Number> &bp)
{
  bp.print(out);
  return out;
}

// Template function definitions

// BarycentricPolynomial:
template <int dim, typename Number>
BarycentricPolynomial<dim, Number>::BarycentricPolynomial()
{
  TableIndices<dim + 1> extents;
  for (unsigned int d = 0; d < dim + 1; ++d)
    extents[d] = 1;
  coefficients.reinit(extents);

  coefficients(TableIndices<dim + 1>{}) = Number();
}



template <int dim, typename Number>
BarycentricPolynomial<dim, Number>::BarycentricPolynomial(
  const TableIndices<dim + 1> &powers,
  const Number                 coefficient)
{
  TableIndices<dim + 1> extents;
  for (unsigned int d = 0; d < dim + 1; ++d)
    extents[d] = powers[d] + 1;
  coefficients.reinit(extents);

  coefficients(powers) = coefficient;
}



template <int dim, typename Number>
BarycentricPolynomial<dim, Number>
BarycentricPolynomial<dim, Number>::monomial(const unsigned int d)
{
  AssertIndexRange(d, dim + 1);
  TableIndices<dim + 1> indices;
  indices[d] = 1;
  return BarycentricPolynomial<dim, Number>(indices, Number(1));
}



template <int dim, typename Number>
void
BarycentricPolynomial<dim, Number>::print(std::ostream &out) const
{
  const auto &coeffs     = this->coefficients;
  auto        first      = index_to_indices(0, coeffs.size());
  bool        print_plus = false;
  if (coeffs(first) != Number())
    {
      out << coeffs(first);
      print_plus = true;
    }
  for (std::size_t i = 1; i < coeffs.n_elements(); ++i)
    {
      const auto indices = index_to_indices(i, coeffs.size());
      if (coeffs(indices) == Number())
        continue;
      if (print_plus)
        out << " + ";
      out << coeffs(indices);
      for (unsigned int d = 0; d < dim + 1; ++d)
        {
          if (indices[d] != 0)
            out << " * t" << d << '^' << indices[d];
        }
      print_plus = true;
    }

  if (!print_plus)
    out << Number();
}



template <int dim, typename Number>
TableIndices<dim + 1>
BarycentricPolynomial<dim, Number>::degrees() const
{
  auto deg = coefficients.size();
  for (unsigned int d = 0; d < dim + 1; ++d)
    deg[d] -= 1;
  return deg;
}



template <int dim, typename Number>
BarycentricPolynomial<dim, Number>
BarycentricPolynomial<dim, Number>::operator-() const
{
  return *this * Number(-1);
}



template <int dim, typename Number>
template <typename Number2>
BarycentricPolynomial<dim, Number>
BarycentricPolynomial<dim, Number>::operator+(const Number2 &a) const
{
  BarycentricPolynomial<dim, Number> result(*this);
  result.coefficients(index_to_indices(0, result.coefficients.size())) += a;

  return result;
}



template <int dim, typename Number>
template <typename Number2>
BarycentricPolynomial<dim, Number>
BarycentricPolynomial<dim, Number>::operator-(const Number2 &a) const
{
  return *this + (-a);
}



template <int dim, typename Number>
template <typename Number2>
BarycentricPolynomial<dim, Number>
BarycentricPolynomial<dim, Number>::operator*(const Number2 &a) const
{
  if (a == Number2())
    {
      return BarycentricPolynomial<dim, Number>();
    }

  BarycentricPolynomial<dim, Number> result(*this);
  for (std::size_t i = 0; i < result.coefficients.n_elements(); ++i)
    {
      const auto index = index_to_indices(i, result.coefficients.size());
      result.coefficients(index) *= a;
    }

  return result;
}



template <int dim, typename Number>
template <typename Number2>
BarycentricPolynomial<dim, Number>
BarycentricPolynomial<dim, Number>::operator/(const Number2 &a) const
{
  Assert(a != Number2(), ExcDivideByZero());
  return *this * (Number(1) / Number(a));
}



template <int dim, typename Number>
BarycentricPolynomial<dim, Number>
BarycentricPolynomial<dim, Number>::operator+(
  const BarycentricPolynomial<dim, Number> &augend) const
{
  TableIndices<dim + 1> deg;
  for (unsigned int d = 0; d < dim + 1; ++d)
    {
      deg[d] = std::max(degrees()[d], augend.degrees()[d]);
    }

  BarycentricPolynomial<dim, Number> result(deg, Number());

  auto add_coefficients = [&](const Table<dim + 1, Number> &in) {
    for (std::size_t i = 0; i < in.n_elements(); ++i)
      {
        const auto index = index_to_indices(i, in.size());
        result.coefficients(index) += in(index);
      }
  };

  add_coefficients(this->coefficients);
  add_coefficients(augend.coefficients);
  return result;
}



template <int dim, typename Number>
BarycentricPolynomial<dim, Number>
BarycentricPolynomial<dim, Number>::operator-(
  const BarycentricPolynomial<dim, Number> &augend) const
{
  return *this + (-augend);
}



template <int dim, typename Number>
BarycentricPolynomial<dim, Number>
BarycentricPolynomial<dim, Number>::operator*(
  const BarycentricPolynomial<dim, Number> &multiplicand) const
{
  TableIndices<dim + 1> deg;
  for (unsigned int d = 0; d < dim + 1; ++d)
    {
      deg[d] = multiplicand.degrees()[d] + degrees()[d];
    }

  BarycentricPolynomial<dim, Number> result(deg, Number());

  const auto &coef_1   = this->coefficients;
  const auto &coef_2   = multiplicand.coefficients;
  auto       &coef_out = result.coefficients;

  for (std::size_t i1 = 0; i1 < coef_1.n_elements(); ++i1)
    {
      const auto index_1 = index_to_indices(i1, coef_1.size());
      for (std::size_t i2 = 0; i2 < coef_2.n_elements(); ++i2)
        {
          const auto index_2 = index_to_indices(i2, coef_2.size());

          TableIndices<dim + 1> index_out;
          for (unsigned int d = 0; d < dim + 1; ++d)
            index_out[d] = index_1[d] + index_2[d];
          coef_out(index_out) += coef_1(index_1) * coef_2(index_2);
        }
    }

  return result;
}



template <int dim, typename Number>
BarycentricPolynomial<dim, Number>
BarycentricPolynomial<dim, Number>::barycentric_derivative(
  const unsigned int coordinate) const
{
  AssertIndexRange(coordinate, dim + 1);

  if (degrees()[coordinate] == 0)
    return BarycentricPolynomial<dim, Number>();

  auto deg = degrees();
  deg[coordinate] -= 1;
  BarycentricPolynomial<dim, Number> result(deg,
                                            std::numeric_limits<Number>::max());
  const auto                        &coeffs_in  = coefficients;
  auto                              &coeffs_out = result.coefficients;
  for (std::size_t i = 0; i < coeffs_out.n_elements(); ++i)
    {
      const auto out_index   = index_to_indices(i, coeffs_out.size());
      auto       input_index = out_index;
      input_index[coordinate] += 1;

      coeffs_out(out_index) = coeffs_in(input_index) * input_index[coordinate];
    }

  return result;
}



template <int dim, typename Number>
BarycentricPolynomial<dim, Number>
BarycentricPolynomial<dim, Number>::derivative(
  const unsigned int coordinate) const
{
  AssertIndexRange(coordinate, dim);
  return -barycentric_derivative(0) + barycentric_derivative(coordinate + 1);
}



template <int dim, typename Number>
Number
BarycentricPolynomial<dim, Number>::value(const Point<dim> &point) const
{
  // TODO: this is probably not numerically stable for higher order.
  // We really need some version of Horner's method.
  Number result = {};

  // Begin by converting point (which is in Cartesian coordinates) to
  // barycentric coordinates:
  std::array<Number, dim + 1> b_point;
  b_point[0] = 1.0;
  for (unsigned int d = 0; d < dim; ++d)
    {
      b_point[0] -= point[d];
      b_point[d + 1] = point[d];
    }

  // Now evaluate the polynomial at the computed barycentric point:
  for (std::size_t i = 0; i < coefficients.n_elements(); ++i)
    {
      const auto indices = index_to_indices(i, coefficients.size());
      const auto coef    = coefficients(indices);
      if (coef == Number())
        continue;

      auto temp = Number(1);
      for (unsigned int d = 0; d < dim + 1; ++d)
        temp *= Utilities::pow(b_point[d], indices[d]);
      result += coef * temp;
    }

  return result;
}



template <int dim, typename Number>
std::size_t
BarycentricPolynomial<dim, Number>::memory_consumption() const
{
  return coefficients.memory_consumption();
}



template <int dim, typename Number>
TableIndices<dim + 1>
BarycentricPolynomial<dim, Number>::index_to_indices(
  const std::size_t           &index,
  const TableIndices<dim + 1> &extents)
{
  TableIndices<dim + 1> result;
  auto                  temp = index;

  for (unsigned int n = 0; n < dim + 1; ++n)
    {
      std::size_t slice_size = 1;
      for (unsigned int n2 = n + 1; n2 < dim + 1; ++n2)
        slice_size *= extents[n2];
      result[n] = temp / slice_size;
      temp %= slice_size;
    }
  return result;
}



template <int dim>
const BarycentricPolynomial<dim> &
BarycentricPolynomials<dim>::operator[](const std::size_t i) const
{
  AssertIndexRange(i, polys.size());
  return polys[i];
}

DEAL_II_NAMESPACE_CLOSE

#endif

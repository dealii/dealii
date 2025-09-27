// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_derivative_form_h
#define dealii_derivative_form_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/tensor.h>

#include <cstddef>


DEAL_II_NAMESPACE_OPEN

/**
 * This class represents the (tangential) derivatives of a function $ \mathbf F:
 * {\mathbb R}^{\text{dim}} \rightarrow {\mathbb R}^{\text{spacedim}}$. Such
 * functions are always used to map the reference dim-dimensional cell into
 * spacedim-dimensional space. For such objects, the first derivative of the
 * function is a linear map from ${\mathbb R}^{\text{dim}}$ to ${\mathbb
 * R}^{\text{spacedim}}$, i.e., it can be represented as a matrix in ${\mathbb
 * R}^{\text{spacedim}\times \text{dim}}$. This makes sense since one would
 * represent the first derivative, $\nabla \mathbf F(\mathbf x)$ with $\mathbf
 * x\in
 * {\mathbb R}^{\text{dim}}$, in such a way that the directional derivative in
 * direction $\mathbf d\in {\mathbb R}^{\text{dim}}$ so that
 * @f{align*}{
 *   \nabla \mathbf F(\mathbf x) \mathbf d
 *   = \lim_{\varepsilon\rightarrow 0}
 *     \frac{\mathbf F(\mathbf x + \varepsilon \mathbf d) - \mathbf F(\mathbf
 * x)}{\varepsilon},
 * @f}
 * i.e., one needs to be able to multiply the matrix $\nabla \mathbf F(\mathbf
 * x)$ by a vector in ${\mathbb R}^{\text{dim}}$, and the result is a difference
 * of function values, which are in ${\mathbb R}^{\text{spacedim}}$.
 * Consequently, the matrix must be of size $\text{spacedim}\times\text{dim}$.
 *
 * Similarly, the second derivative is a bilinear map from  ${\mathbb
 * R}^{\text{dim}} \times  {\mathbb R}^{\text{dim}}$ to ${\mathbb
 * R}^{\text{spacedim}}$, which one can think of a rank-3 object of size
 * $\text{spacedim}\times\text{dim}\times\text{dim}$.
 *
 * In deal.II we represent these derivatives using objects of type
 * DerivativeForm@<1,dim,spacedim,Number@>,
 * DerivativeForm@<2,dim,spacedim,Number@> and so on.
 */
template <int order, int dim, int spacedim, typename Number = double>
class DerivativeForm
{
public:
  /**
   * Constructor. Initialize all entries to zero.
   */
  DerivativeForm() = default;

  /**
   * Constructor from a tensor.
   */
  DerivativeForm(const Tensor<order + 1, dim, Number> &);

  /**
   * Constructor. This constructor initializes the data stored by this
   * object from a `spacedim x dim^order` array that is represented
   * by a Tensor with `spacedim` components each of which is a
   * Tensor of rank `order` and size `dim`.
   */
  DerivativeForm(const Tensor<1, spacedim, Tensor<order, dim, Number>> &);

  /**
   * Read-Write access operator.
   */
  Tensor<order, dim, Number> &
  operator[](const unsigned int i);

  /**
   * Read-only access operator.
   */
  const Tensor<order, dim, Number> &
  operator[](const unsigned int i) const;

  /**
   * Assignment operator.
   */
  DerivativeForm &
  operator=(const Tensor<order + 1, dim, Number> &);

  /**
   * Assignment operator.
   */
  DerivativeForm &
  operator=(const Tensor<order, spacedim, Tensor<1, dim, Number>> &);

  /**
   * Assignment operator.
   */
  DerivativeForm &
  operator=(const Tensor<1, dim, Number> &);

  /**
   * Number conversion operator.
   */
  template <typename OtherNumber>
  DerivativeForm &
  operator=(const DerivativeForm<order, dim, spacedim, OtherNumber> &df);

  /**
   * Converts a DerivativeForm <order, dim, dim, Number> to Tensor<order+1, dim,
   * Number>. In particular, if order == 1 and the derivative is the Jacobian of
   * $\mathbf F(\mathbf x)$, then Tensor[i] = $\nabla F_i(\mathbf x)$.
   */
  operator Tensor<order + 1, dim, Number>() const;

  /**
   * Converts a DerivativeForm<1, dim, 1, Number> to Tensor<1, dim, Number>.
   */
  operator Tensor<1, dim, Number>() const;

  /**
   * Return the transpose of a rectangular DerivativeForm,
   * viewed as a two dimensional matrix.
   */
  DerivativeForm<1, spacedim, dim, Number>
  transpose() const;

  /**
   * Compute the Frobenius norm of this form, i.e., the expression
   * $\sqrt{\sum_{ij} |DF_{ij}|^2} =
   * \sqrt{\sum_{ij} |\frac{\partial F_i}{\partial x_j}|^2}$.
   */
  typename numbers::NumberTraits<Number>::real_type
  norm() const;

  /**
   * Compute the volume element associated with the Jacobian of the
   * transformation $\mathbf F$. That is to say if $DF$ is square, it computes
   * $\det(DF)$, in case DF is not square returns $\sqrt{\det(DF^T \,DF)}$.
   */
  Number
  determinant() const;

  /**
   * Assuming that the current object stores the Jacobian of a mapping
   * $\mathbf F$, then the current function computes the <i>covariant</i> form
   * of the derivative, namely $(\nabla \mathbf F) {\mathbf G}^{-1}$, where
   * $\mathbf G = (\nabla \mathbf F)^{T}(\nabla \mathbf F)$. If $\nabla \mathbf
   * F$ is a square matrix (i.e., $\mathbf F:
   * {\mathbb R}^n \mapsto {\mathbb R}^n$), then this function
   * simplifies to computing $\nabla {\mathbf F}^{-T}$.
   */
  DerivativeForm<1, dim, spacedim, Number>
  covariant_form() const;


  /**
   * Compute the first fundamental form. This is the tensor of dot products of
   * the columns of the current object.
   *
   * @note This function is only defined for `order == 1`.
   */
  Tensor<2, dim, Number>
  first_fundamental_form() const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  static std::size_t
  memory_consumption();

  /**
   * Exception.
   */
  DeclException1(ExcInvalidTensorIndex,
                 int,
                 << "Invalid DerivativeForm index " << arg1);

private:
  /**
   * Auxiliary function that computes $A T^{T}$ where A represents the current
   * object.
   */
  DerivativeForm<1, dim, spacedim, Number>
  times_T_t(const Tensor<2, dim, Number> &T) const;


  /**
   * Array of tensors holding the subelements.
   */
  Tensor<order, dim, Number> tensor[spacedim];
};


/*--------------------------- Inline functions -----------------------------*/

#ifndef DOXYGEN

template <int order, int dim, int spacedim, typename Number>
inline DerivativeForm<order, dim, spacedim, Number>::DerivativeForm(
  const Tensor<order + 1, dim, Number> &T)
{
  Assert((dim == spacedim),
         ExcMessage("Only allowed for forms with dim==spacedim."));
  if (dim == spacedim)
    for (unsigned int j = 0; j < dim; ++j)
      (*this)[j] = T[j];
}



template <int order, int dim, int spacedim, typename Number>
inline DerivativeForm<order, dim, spacedim, Number>::DerivativeForm(
  const Tensor<1, spacedim, Tensor<order, dim, Number>> &T)
{
  for (unsigned int j = 0; j < spacedim; ++j)
    (*this)[j] = T[j];
}



template <int order, int dim, int spacedim, typename Number>
inline DerivativeForm<order, dim, spacedim, Number> &
DerivativeForm<order, dim, spacedim, Number>::operator=(
  const Tensor<order + 1, dim, Number> &ta)
{
  Assert((dim == spacedim), ExcMessage("Only allowed when dim==spacedim."));

  if (dim == spacedim)
    for (unsigned int j = 0; j < dim; ++j)
      (*this)[j] = ta[j];
  return *this;
}



template <int order, int dim, int spacedim, typename Number>
inline DerivativeForm<order, dim, spacedim, Number> &
DerivativeForm<order, dim, spacedim, Number>::operator=(
  const Tensor<order, spacedim, Tensor<1, dim, Number>> &T)
{
  for (unsigned int j = 0; j < spacedim; ++j)
    (*this)[j] = T[j];
  return *this;
}



template <int order, int dim, int spacedim, typename Number>
inline DerivativeForm<order, dim, spacedim, Number> &
DerivativeForm<order, dim, spacedim, Number>::operator=(
  const Tensor<1, dim, Number> &T)
{
  Assert((1 == spacedim) && (order == 1),
         ExcMessage("Only allowed for spacedim==1 and order==1."));

  (*this)[0] = T;

  return *this;
}



template <int order, int dim, int spacedim, typename Number>
template <typename OtherNumber>
inline DerivativeForm<order, dim, spacedim, Number> &
DerivativeForm<order, dim, spacedim, Number>::operator=(
  const DerivativeForm<order, dim, spacedim, OtherNumber> &df)
{
  for (unsigned int j = 0; j < spacedim; ++j)
    (*this)[j] = df[j];
  return *this;
}



template <int order, int dim, int spacedim, typename Number>
inline Tensor<order, dim, Number> &
DerivativeForm<order, dim, spacedim, Number>::operator[](const unsigned int i)
{
  AssertIndexRange(i, spacedim);

  return tensor[i];
}



template <int order, int dim, int spacedim, typename Number>
inline const Tensor<order, dim, Number> &
DerivativeForm<order, dim, spacedim, Number>::operator[](
  const unsigned int i) const
{
  AssertIndexRange(i, spacedim);

  return tensor[i];
}



template <int order, int dim, int spacedim, typename Number>
inline DerivativeForm<order, dim, spacedim, Number>::
operator Tensor<1, dim, Number>() const
{
  Assert((1 == spacedim) && (order == 1),
         ExcMessage("Only allowed for spacedim==1."));

  return (*this)[0];
}



template <int order, int dim, int spacedim, typename Number>
inline DerivativeForm<order, dim, spacedim, Number>::
operator Tensor<order + 1, dim, Number>() const
{
  Assert((dim == spacedim), ExcMessage("Only allowed when dim==spacedim."));

  Tensor<order + 1, dim, Number> t;

  if (dim == spacedim)
    for (unsigned int j = 0; j < dim; ++j)
      t[j] = (*this)[j];

  return t;
}



template <int order, int dim, int spacedim, typename Number>
inline DerivativeForm<1, spacedim, dim, Number>
DerivativeForm<order, dim, spacedim, Number>::transpose() const
{
  Assert(order == 1, ExcMessage("Only for rectangular DerivativeForm."));
  DerivativeForm<1, spacedim, dim, Number> tt;

  for (unsigned int i = 0; i < spacedim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      tt[j][i] = (*this)[i][j];

  return tt;
}



template <int order, int dim, int spacedim, typename Number>
inline DerivativeForm<1, dim, spacedim, Number>
DerivativeForm<order, dim, spacedim, Number>::times_T_t(
  const Tensor<2, dim, Number> &T) const
{
  Assert(order == 1, ExcMessage("Only for order == 1."));
  DerivativeForm<1, dim, spacedim, Number> dest;
  for (unsigned int i = 0; i < spacedim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      dest[i][j] = (*this)[i] * T[j];

  return dest;
}



template <int order, int dim, int spacedim, typename Number>
inline typename numbers::NumberTraits<Number>::real_type
DerivativeForm<order, dim, spacedim, Number>::norm() const
{
  typename numbers::NumberTraits<Number>::real_type sum_of_squares = 0;
  for (unsigned int i = 0; i < spacedim; ++i)
    sum_of_squares += tensor[i].norm_square();
  return std::sqrt(sum_of_squares);
}



template <int order, int dim, int spacedim, typename Number>
inline Number
DerivativeForm<order, dim, spacedim, Number>::determinant() const
{
  Assert(order == 1, ExcMessage("Only for order == 1."));
  if (dim == spacedim)
    {
      const Tensor<2, dim, Number> T =
        static_cast<Tensor<2, dim, Number>>(*this);
      return dealii::determinant(T);
    }
  else
    {
      Assert(spacedim > dim, ExcMessage("Only for spacedim>dim."));
      return (std::sqrt(dealii::determinant(first_fundamental_form())));
    }
}



template <int order, int dim, int spacedim, typename Number>
inline Tensor<2, dim, Number>
DerivativeForm<order, dim, spacedim, Number>::first_fundamental_form() const
{
  Assert(order == 1, ExcMessage("Only for order == 1."));
  const DerivativeForm<1, spacedim, dim, Number> DF_t = this->transpose();

  Tensor<2, dim, Number> G;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      G[i][j] = DF_t[i] * DF_t[j];

  return G;
}



template <int order, int dim, int spacedim, typename Number>
inline DerivativeForm<1, dim, spacedim, Number>
DerivativeForm<order, dim, spacedim, Number>::covariant_form() const
{
  if (dim == spacedim)
    {
      const Tensor<2, dim, Number> DF_t =
        dealii::transpose(invert(static_cast<Tensor<2, dim, Number>>(*this)));
      return DerivativeForm<1, dim, spacedim, Number>(DF_t);
    }
  else
    {
      return (this->times_T_t(invert(first_fundamental_form())));
    }
}


template <int order, int dim, int spacedim, typename Number>
inline std::size_t
DerivativeForm<order, dim, spacedim, Number>::memory_consumption()
{
  return sizeof(DerivativeForm<order, dim, spacedim, Number>);
}

#endif // DOXYGEN



/**
 * Output operator for DerivativeForm. Print the elements consecutively, with a
 * space in between, two spaces between rank 1 subtensors, three between rank 2
 * and so on.
 *
 * @relatesalso DerivativeForm
 */
template <int order, int dim, int spacedim, typename Number>
inline std::ostream &
operator<<(std::ostream                                       &out,
           const DerivativeForm<order, dim, spacedim, Number> &df)
{
  for (unsigned int i = 0; i < spacedim; ++i)
    {
      out << df[i];
      if (i != spacedim - 1)
        for (unsigned int j = 0; j < order + 1; ++j)
          out << ' ';
    }

  return out;
}



/**
 * One of the uses of DerivativeForm is to apply it as a linear transformation.
 * This function returns $\nabla \mathbf F(\mathbf x) \Delta \mathbf x$, which
 * approximates the change in $\mathbf F(\mathbf x)$ when $\mathbf x$ is changed
 * by the amount $\Delta \mathbf x$
 * @f[
 *   \nabla \mathbf F(\mathbf x) \; \Delta \mathbf x
 *   \approx
 *   \mathbf F(\mathbf x + \Delta \mathbf x) - \mathbf F(\mathbf x).
 * @f]
 * The transformation corresponds to
 * @f[
 *   [\text{result}]_{i_1,\dots,i_k} = i\sum_{j}
 *   \left[\nabla \mathbf F(\mathbf x)\right]_{i_1,\dots,i_k, j}
 *   \Delta x_j
 * @f]
 * in index notation and corresponds to
 * $[\Delta \mathbf x] [\nabla \mathbf F(\mathbf x)]^T$ in matrix notation.
 *
 * @relatesalso DerivativeForm
 */
template <int spacedim, int dim, typename Number1, typename Number2>
inline Tensor<1, spacedim, typename ProductType<Number1, Number2>::type>
apply_transformation(const DerivativeForm<1, dim, spacedim, Number1> &grad_F,
                     const Tensor<1, dim, Number2>                   &d_x)
{
  Tensor<1, spacedim, typename ProductType<Number1, Number2>::type> dest;
  for (unsigned int i = 0; i < spacedim; ++i)
    dest[i] = grad_F[i] * d_x;
  return dest;
}



/**
 * Similar to the previous apply_transformation().
 * Each row of the result corresponds to one of the rows of @p D_X transformed
 * by @p grad_F, equivalent to $\mathrm{D\_X} \, \mathrm{grad\_F}^T$ in matrix
 * notation.
 *
 * @relatesalso DerivativeForm
 */
// rank=2
template <int spacedim, int dim, typename Number1, typename Number2>
inline DerivativeForm<1,
                      spacedim,
                      dim,
                      typename ProductType<Number1, Number2>::type>
apply_transformation(const DerivativeForm<1, dim, spacedim, Number1> &grad_F,
                     const Tensor<2, dim, Number2>                   &D_X)
{
  DerivativeForm<1, spacedim, dim, typename ProductType<Number1, Number2>::type>
    dest;
  for (unsigned int i = 0; i < dim; ++i)
    dest[i] = apply_transformation(grad_F, D_X[i]);

  return dest;
}



/**
 * Similar to the previous apply_transformation(), specialized for the case `dim
 * == spacedim` where we can return a rank-2 tensor instead of the more general
 * `DerivativeForm`.
 * Each row of the result corresponds to one of the rows of @p D_X transformed
 * by @p grad_F, equivalent to $\mathrm{D\_X} \, \mathrm{grad\_F}^T$ in matrix
 * notation.
 *
 * @relatesalso DerivativeForm
 */
// rank=2
template <int dim, typename Number1, typename Number2>
inline Tensor<2, dim, typename ProductType<Number1, Number2>::type>
apply_transformation(const DerivativeForm<1, dim, dim, Number1> &grad_F,
                     const Tensor<2, dim, Number2>              &D_X)
{
  Tensor<2, dim, typename ProductType<Number1, Number2>::type> dest;
  for (unsigned int i = 0; i < dim; ++i)
    dest[i] = apply_transformation(grad_F, D_X[i]);

  return dest;
}



/**
 * Similar to the previous apply_transformation().
 * Each row of the result corresponds to one of the rows of @p D_X transformed
 * by @p grad_F.
 *
 * @relatesalso DerivativeForm
 */
template <int spacedim,
          int dim,
          int n_components,
          typename Number1,
          typename Number2>
inline Tensor<1,
              n_components,
              Tensor<1, spacedim, typename ProductType<Number1, Number2>::type>>
apply_transformation(
  const DerivativeForm<1, dim, spacedim, Number1>        &grad_F,
  const Tensor<1, n_components, Tensor<1, dim, Number2>> &D_X)
{
  Tensor<1,
         n_components,
         Tensor<1, spacedim, typename ProductType<Number1, Number2>::type>>
    dest;
  for (unsigned int i = 0; i < n_components; ++i)
    dest[i] = apply_transformation(grad_F, D_X[i]);

  return dest;
}



/**
 * Similar to the previous apply_transformation(). In matrix notation, it
 * computes $DF2 \, DF1^{T}$. Moreover, the result of this operation $\mathbf A$
 * can be interpreted as a metric tensor in
 * ${\mathbb R}^\text{spacedim}$ which corresponds to the Euclidean metric
 * tensor in
 * ${\mathbb R}^\text{dim}$. For every pair of vectors
 * $\mathbf u, \mathbf v \in {\mathbb R}^\text{spacedim}$, we have:
 * @f[
 *   \mathbf u \cdot \mathbf A \mathbf v =
 *   \text{DF2}^{-1}(\mathbf u) \cdot \text{DF1}^{-1}(\mathbf v)
 * @f]
 *
 * @relatesalso DerivativeForm
 */
template <int spacedim, int dim, typename Number1, typename Number2>
inline Tensor<2, spacedim, typename ProductType<Number1, Number2>::type>
apply_transformation(const DerivativeForm<1, dim, spacedim, Number1> &DF1,
                     const DerivativeForm<1, dim, spacedim, Number2> &DF2)
{
  Tensor<2, spacedim, typename ProductType<Number1, Number2>::type> dest;

  for (unsigned int i = 0; i < spacedim; ++i)
    dest[i] = apply_transformation(DF1, DF2[i]);

  return dest;
}



/**
 * Transpose of a rectangular DerivativeForm DF, mostly for compatibility
 * reasons.
 *
 * @relatesalso DerivativeForm
 */
template <int dim, int spacedim, typename Number>
inline DerivativeForm<1, spacedim, dim, Number>
transpose(const DerivativeForm<1, dim, spacedim, Number> &DF)
{
  return DF.transpose();
}



/**
 * Specialization of apply_transformation() for a diagonal DerivativeForm.
 *
 * @relatesalso DerivativeForm
 */
template <int spacedim, int dim, typename Number1, typename Number2>
inline Tensor<1, spacedim, typename ProductType<Number1, Number2>::type>
apply_diagonal_transformation(
  const DerivativeForm<1, dim, spacedim, Number1> &grad_F,
  const Tensor<1, dim, Number2>                   &d_x)
{
  Assert(dim == spacedim,
         ExcMessage("Only dim = spacedim allowed for diagonal transformation"));
  Tensor<1, spacedim, typename ProductType<Number1, Number2>::type> dest;
  for (unsigned int i = 0; i < spacedim; ++i)
    dest[i] = grad_F[i][i] * d_x[i];
  return dest;
}


/**
 * Similar to the previous apply_diagonal_transformation(), specialized for the
 * case `dim == spacedim` where we can return a rank-2 tensor instead of the
 * more general `DerivativeForm`.
 * Each row of the result corresponds to one of the rows of @p D_X transformed
 * by @p grad_F, equivalent to $\mathrm{D\_X} \, \mathrm{grad\_F}^T$ in matrix
 * notation.
 *
 * @relatesalso DerivativeForm
 */
template <int dim, typename Number1, typename Number2>
inline Tensor<2, dim, typename ProductType<Number1, Number2>::type>
apply_diagonal_transformation(
  const DerivativeForm<1, dim, dim, Number1> &grad_F,
  const Tensor<2, dim, Number2>              &D_X)
{
  Tensor<2, dim, typename ProductType<Number1, Number2>::type> dest;
  for (unsigned int i = 0; i < dim; ++i)
    dest[i] = apply_diagonal_transformation(grad_F, D_X[i]);

  return dest;
}



/**
 * Similar to the previous apply_diagonal_transformation().
 * Each row of the result corresponds to one of the rows of @p D_X transformed
 * by @p grad_F.
 *
 * @relatesalso DerivativeForm
 */
template <int spacedim,
          int dim,
          int n_components,
          typename Number1,
          typename Number2>
inline Tensor<1,
              n_components,
              Tensor<1, spacedim, typename ProductType<Number1, Number2>::type>>
apply_diagonal_transformation(
  const DerivativeForm<1, dim, spacedim, Number1>        &grad_F,
  const Tensor<1, n_components, Tensor<1, dim, Number2>> &D_X)
{
  Tensor<1,
         n_components,
         Tensor<1, spacedim, typename ProductType<Number1, Number2>::type>>
    dest;
  for (unsigned int i = 0; i < n_components; ++i)
    dest[i] = apply_diagonal_transformation(grad_F, D_X[i]);

  return dest;
}



/**
 * Similar to the previous apply_transformation().
 * Each row of the result corresponds to one of the rows of @p D_X transformed
 * by @p grad_F, equivalent to $\mathrm{D\_X} \, \mathrm{grad\_F}^T$ in matrix
 * notation.
 *
 * @relatesalso DerivativeForm
 */
// rank=2
template <int spacedim, int dim, typename Number1, typename Number2>
inline DerivativeForm<1,
                      spacedim,
                      dim,
                      typename ProductType<Number1, Number2>::type>
apply_diagonal_transformation(
  const DerivativeForm<1, dim, spacedim, Number1> &grad_F,
  const Tensor<2, dim, Number2>                   &D_X)
{
  DerivativeForm<1, spacedim, dim, typename ProductType<Number1, Number2>::type>
    dest;
  for (unsigned int i = 0; i < dim; ++i)
    dest[i] = apply_diagonal_transformation(grad_F, D_X[i]);

  return dest;
}


DEAL_II_NAMESPACE_CLOSE

#endif

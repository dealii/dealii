// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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

#ifndef dealii_elasticity_standard_tensors_h
#define dealii_elasticity_standard_tensors_h


#include <deal.II/base/config.h>

#include <deal.II/base/numbers.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

DEAL_II_NAMESPACE_OPEN

namespace Physics
{
  namespace Elasticity
  {
    /**
     * A collection of tensor definitions that mostly conform to notation used
     * in standard scientific literature, in particular the book of
     * Wriggers (2008). The citation for this reference, as well as other
     * notation used here, can be found in the description for the
     * Physics::Elasticity namespace.
     *
     * @note These hold specifically for the codimension 0 case with a
     * Cartesian basis, where the metric tensor is the identity tensor.
     *
     * @relatesalso Tensor
     * @relatesalso SymmetricTensor
     * @author Jean-Paul Pelteret, Andrew McBride, 2016
     */
    template <int dim>
    class StandardTensors
    {
    public:
      /**
       * @name Metric tensors
       */
      //@{

      /**
       * The second-order referential/spatial symmetric identity (metric) tensor
       * $\mathbf{I}$.
       *
       * This is defined such that, for any rank-2 tensor or symmetric tensor,
       * the following holds:
       * @f[
       *  \mathbf{I} \cdot \{ \bullet \} = \{ \bullet \} \cdot \mathbf{I} =
       * \{ \bullet \}
       *    \qquad \text{and} \qquad
       * \mathbf{I} : \{ \bullet \} = \textrm{trace} \{ \bullet \} \, .
       * @f]
       *
       * This definition aligns with the rank-2 symmetric tensor returned by
       * unit_symmetric_tensor(). If one is to interpret the tensor as a
       * matrix, then this simply corresponds to the identity matrix.
       */
      static DEAL_II_CONSTEXPR const SymmetricTensor<2, dim> I
#ifdef DEAL_II_HAVE_CXX14_CONSTEXPR
        = unit_symmetric_tensor<dim>()
#endif
        ;

      /**
       * The fourth-order referential/spatial unit symmetric tensor
       * $\mathcal{S}$.
       *
       * This is defined such that for a general rank-2 tensor $\{ \hat{\bullet}
       * \}$ the following holds:
       * @f[
       *   \mathcal{S} : \{ \hat{\bullet} \}
       *   \dealcoloneq \dfrac{1}{2}
       *   \left[ \{ \hat{\bullet} \} + \{ \hat{\bullet} \}^T \right] \, .
       * @f]
       *
       * As a corollary to this, for any second-order symmetric tensor $\{
       * \bullet \}$
       * @f[
       *  \mathcal{S} : \{ \bullet \}
       *    = \{ \bullet \} : \mathcal{S} = \{ \bullet \} \, .
       * @f]
       *
       * This definition aligns with the fourth-order symmetric tensor
       * $\mathcal{S}$ introduced in the Physics::Elasticity namespace
       * description and that which is returned by identity_tensor().
       *
       * @note If you apply this to a standard tensor then it doesn't behave like
       * the fourth-order identity tensor, but rather as a symmetrization
       * operator.
       */
      static DEAL_II_CONSTEXPR const SymmetricTensor<4, dim> S
#ifdef DEAL_II_HAVE_CXX14_CONSTEXPR
        = identity_tensor<dim>()
#endif
        ;

      /**
       * The fourth-order referential/spatial tensor $\mathbf{I} \otimes
       * \mathbf{I}$.
       *
       * This is defined such that, for any rank-2 tensor, the following holds:
       * @f[
       *  [\mathbf{I} \otimes \mathbf{I}] : \{ \bullet \} =
       *  \textrm{trace}\{ \bullet \} \mathbf{I} \, .
       * @f]
       */
      static DEAL_II_CONSTEXPR const SymmetricTensor<4, dim> IxI
#ifdef DEAL_II_HAVE_CXX14_CONSTEXPR
        = outer_product(unit_symmetric_tensor<dim>(),
                        unit_symmetric_tensor<dim>())
#endif
        ;

      //@}

      /**
       * @name Projection operators
       */
      //@{

      /**
       * The fourth-order spatial deviatoric tensor. Also known as the
       * deviatoric operator, this tensor projects a second-order symmetric
       * tensor onto a deviatoric space (for which the hydrostatic component is
       * removed).
       *
       * This is defined as
       * @f[
       *   \mathcal{P}
       *     \dealcoloneq \mathcal{S} - \frac{1}{\textrm{dim}} \mathbf{I}
       *     \otimes \mathbf{I}
       * @f]
       * where $\mathcal{S}$ is the fourth-order unit symmetric tensor and
       * $\mathbf{I}$ is the second-order identity tensor.
       *
       * For any second-order (spatial) symmetric tensor the following holds:
       * @f[
       *  \mathcal{P} : \{ \bullet \}
       *  \dealcoloneq \{ \bullet \} - \frac{1}{\textrm{dim}}
       *  \left[ \{ \bullet \} : \mathbf{I} \right]\mathbf{I}
       *  = \mathcal{P}^{T} : \{ \bullet \}
       *  = \texttt{dev\_P} \left( \{ \bullet \} \right)
       * @f]
       * and, therefore,
       * @f[
       * \texttt{dev\_P} \left( \{ \bullet \} \right) : \mathbf{I}
       *   = \textrm{trace}(\texttt{dev\_P} \left( \{ \bullet \} \right)) = 0 \,
       * .
       * @f]
       *
       * This definition aligns with the fourth-order symmetric tensor that
       * is returned by deviator_tensor().
       *
       * @dealiiWriggersA{47,3.129}
       * @dealiiHolzapfelA{232,6.105}
       */
      static DEAL_II_CONSTEXPR const SymmetricTensor<4, dim> dev_P
#ifdef DEAL_II_HAVE_CXX14_CONSTEXPR
        = deviator_tensor<dim>()
#endif
        ;

      /**
       * Return the fourth-order referential deviatoric tensor, as constructed
       * from
       * the deformation gradient tensor @p F.
       * Also known as the deviatoric operator, this tensor projects a
       * second-order symmetric tensor onto a deviatoric space (for which the
       * hydrostatic component is removed).
       *
       * This referential isochoric projection tensor is defined as
       * @f[
       *   \hat{\mathcal{P}}
       *     \dealcoloneq \frac{\partial \bar{\mathbf{C}}}{\partial \mathbf{C}}
       * @f]
       * with
       * @f[
       *  \bar{\mathbf{C}} \dealcoloneq J^{-2/\textrm{dim}} \mathbf{C}
       *    \qquad \text{,} \qquad
       *  \mathbf{C} = \mathbf{F}^{T}\cdot\mathbf{F}
       *    \qquad \text{and} \qquad
       *  J = \textrm{det}\mathbf{F}
       * @f]
       * such that, for any second-order (referential) symmetric tensor,
       * the following holds:
       * @f[
       *  \{ \bullet \} : \hat{\mathcal{P}}
       *    \dealcoloneq J^{-2/\textrm{dim}} \left[ \{ \bullet \} -
       * \frac{1}{\textrm{dim}}\left[\mathbf{C} : \{ \bullet \}\right]
       * \mathbf{C}^{-1} \right] = \texttt{Dev\_P} \left( \{ \bullet \} \right)
       * \, .
       * @f]
       * It can therefore be readily shown that
       * @f[
       *  \texttt{Dev\_P} \left( \{ \bullet \} \right) : \mathbf{C} = 0 \, .
       * @f]
       *
       * @note It may be observed that we have defined the tensor as the
       * transpose of that adopted by Wriggers (2008). We have done this so that
       * it may be strictly applied through the chain rule to achieve the
       * definition of the second Piola-Kirchhoff stress, i.e.
       * @f[
       *   \mathbf{S}
       *     = 2\frac{\partial \psi \left( \bar{\mathbf{C}} \right)}{\partial
       * \mathbf{C}} = 2\frac{\partial \psi \left( \bar{\mathbf{C}}
       * \right)}{\partial \bar{\mathbf{C}}} : \frac{\partial
       * \bar{\mathbf{C}}}{\partial \mathbf{C}} = \bar{\mathbf{S}} :
       * \hat{\mathcal{P}} \equiv \hat{\mathcal{P}}^{T} : \bar{\mathbf{S}} \, .
       * @f]
       *
       * @note Comparing the definition of this tensor in Holzapfel (2001) to that
       * adopted here, the inclusion of the extra factor $J^{-2/\textrm{dim}}$
       * does not, at the outset, seem to be a reasonable choice. However, in
       * the author's view it makes direct implementation of the expressions for
       * isochoric (referential) stress contributions and their linearization
       * simpler in practise.
       *
       * @dealiiWriggersA{46,3.125}
       * @dealiiHolzapfelA{229,6.83}
       */
      template <typename Number>
      static DEAL_II_CONSTEXPR SymmetricTensor<4, dim, Number>
                               Dev_P(const Tensor<2, dim, Number> &F);

      /**
       * Return the transpose of the fourth-order referential deviatoric tensor,
       * as constructed from the deformation gradient tensor @p F.
       * The result performs the following operation:
       * @f[
       *  \hat{\mathcal{P}}^{T} : \{ \bullet \}
       *    = J^{-2/\textrm{dim}} \left[ \{ \bullet \} - \frac{1}{\textrm{dim}}
       * \left[\mathbf{C}^{-1} : \{ \bullet \}\right] \mathbf{C} \right] =
       * \texttt{Dev\_P\_T} \{ \bullet \}
       * @f]
       */
      template <typename Number>
      static DEAL_II_CONSTEXPR SymmetricTensor<4, dim, Number>
                               Dev_P_T(const Tensor<2, dim, Number> &F);

      //@}

      /**
       * @name Scalar derivatives
       */
      //@{
      /**
       * Return the derivative of the volumetric Jacobian
       * $J = \text{det} \mathbf{F}$ with respect to the right Cauchy-Green
       * tensor, as constructed from the deformation gradient tensor @p F.
       * The computed result is
       * @f[
       *  \frac{\partial J}{\partial \mathbf{C}}
       *   = \frac{1}{2} J \mathbf{C}^{-1}
       * @f]
       * with
       * @f[
       *  \mathbf{C} = \mathbf{F}^{T}\cdot\mathbf{F} \, .
       * @f]
       *
       * @dealiiWriggersA{46,3.124}
       * @dealiiHolzapfelA{228,6.82}
       */
      template <typename Number>
      static DEAL_II_CONSTEXPR SymmetricTensor<2, dim, Number>
                               ddet_F_dC(const Tensor<2, dim, Number> &F);

      //@}

      /**
       * @name Tensor derivatives
       */
      //@{

      /**
       * Return the derivative of the inverse of the right Cauchy-Green
       * tensor with respect to the right Cauchy-Green tensor itself,
       * as constructed from the deformation gradient tensor @p F.
       * The result, accounting for symmetry, is defined in index notation as
       * @f[
       *  \left[ \frac{\partial \mathbf{C}^{-1}}{\partial \mathbf{C}}
       * \right]_{IJKL}
       *    \dealcoloneq -\frac{1}{2}[ C^{-1}_{IK}C^{-1}_{JL}
       *     + C^{-1}_{IL}C^{-1}_{JK}  ]
       * @f]
       *
       * @dealiiWriggersA{76,3.255}
       */
      template <typename Number>
      static DEAL_II_CONSTEXPR SymmetricTensor<4, dim, Number>
                               dC_inv_dC(const Tensor<2, dim, Number> &F);

      //@}
    };

  } // namespace Elasticity
} // namespace Physics



#ifndef DOXYGEN

// --------------------- inline functions and constants -------------------


template <int dim>
template <typename Number>
DEAL_II_CONSTEXPR inline SymmetricTensor<4, dim, Number>
Physics::Elasticity::StandardTensors<dim>::Dev_P(
  const Tensor<2, dim, Number> &F)
{
  const Number det_F = determinant(F);
  Assert(numbers::value_is_greater_than(det_F, 0.0),
         ExcMessage("Deformation gradient has a negative determinant."));
  const Tensor<2, dim, Number>          C_ns  = transpose(F) * F;
  const SymmetricTensor<2, dim, Number> C     = symmetrize(C_ns);
  const SymmetricTensor<2, dim, Number> C_inv = symmetrize(invert(C_ns));

  // See Wriggers p46 equ 3.125 (but transpose indices)
  SymmetricTensor<4, dim, Number> Dev_P =
    outer_product(C, C_inv);                   // Dev_P = C_x_C_inv
  Dev_P /= -dim;                               // Dev_P = -[1/dim]C_x_C_inv
  Dev_P += SymmetricTensor<4, dim, Number>(S); // Dev_P = S - [1/dim]C_x_C_inv
  Dev_P *=
    std::pow(det_F, -2.0 / dim); // Dev_P = J^{-2/dim} [S - [1/dim]C_x_C_inv]

  return Dev_P;
}



template <int dim>
template <typename Number>
DEAL_II_CONSTEXPR inline SymmetricTensor<4, dim, Number>
Physics::Elasticity::StandardTensors<dim>::Dev_P_T(
  const Tensor<2, dim, Number> &F)
{
  const Number det_F = determinant(F);
  Assert(numbers::value_is_greater_than(det_F, 0.0),
         ExcMessage("Deformation gradient has a negative determinant."));
  const Tensor<2, dim, Number>          C_ns  = transpose(F) * F;
  const SymmetricTensor<2, dim, Number> C     = symmetrize(C_ns);
  const SymmetricTensor<2, dim, Number> C_inv = symmetrize(invert(C_ns));

  // See Wriggers p46 equ 3.125 (not transposed)
  SymmetricTensor<4, dim, Number> Dev_P_T =
    outer_product(C_inv, C);                     // Dev_P = C_inv_x_C
  Dev_P_T /= -dim;                               // Dev_P = -[1/dim]C_inv_x_C
  Dev_P_T += SymmetricTensor<4, dim, Number>(S); // Dev_P = S - [1/dim]C_inv_x_C
  Dev_P_T *=
    std::pow(det_F, -2.0 / dim); // Dev_P = J^{-2/dim} [S - [1/dim]C_inv_x_C]

  return Dev_P_T;
}



template <int dim>
template <typename Number>
DEAL_II_CONSTEXPR SymmetricTensor<2, dim, Number>
                  Physics::Elasticity::StandardTensors<dim>::ddet_F_dC(
  const Tensor<2, dim, Number> &F)
{
  return Number(0.5) * determinant(F) * symmetrize(invert(transpose(F) * F));
}



template <int dim>
template <typename Number>
DEAL_II_CONSTEXPR inline SymmetricTensor<4, dim, Number>
Physics::Elasticity::StandardTensors<dim>::dC_inv_dC(
  const Tensor<2, dim, Number> &F)
{
  const SymmetricTensor<2, dim, Number> C_inv =
    symmetrize(invert(transpose(F) * F));

  SymmetricTensor<4, dim, Number> dC_inv_dC;
  for (unsigned int A = 0; A < dim; ++A)
    for (unsigned int B = A; B < dim; ++B)
      for (unsigned int C = 0; C < dim; ++C)
        for (unsigned int D = C; D < dim; ++D)
          dC_inv_dC[A][B][C][D] -=
            0.5 * (C_inv[A][C] * C_inv[B][D] + C_inv[A][D] * C_inv[B][C]);

  return dC_inv_dC;
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

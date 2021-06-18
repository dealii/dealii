// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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

#ifndef dealii_physics_invariants_h
#define dealii_physics_invariants_h

#include <deal.II/base/config.h>

// #include <deal.II/base/exceptions.h>
// #include <deal.II/base/numbers.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/physics/elasticity/standard_tensors.h>

#include <set>
#include <type_traits>

DEAL_II_NAMESPACE_OPEN


namespace Physics
{
  /**
   * @brief
   *
   * References:
   *
   * Isotropic
   * Chadwick
   * Holzapfel
   * Ogden (e.g. RWO-talk.pdf)
   *
   * Transverse iso:
   *  - [470, 472]
   *  - [205, 416]
   *  - [68, 105, 470]
   *
   * Coupled:
   * Pelteret (p124)
   *   - [496]
   *   - [573]
   *   - Steinmann CISM notes?
   *    - Khoi habilitation?
   *
   * Coupled trans iso
   * Pelteret
   *   - [68]
   *
   * Irreducible invariants...
   *
   * Reducable invariants that are used sufficiently frequently that they
   * warrant inclusion here.
   */
  namespace Invariants
  {
    enum InvariantList
    {
      // === ISOTROPIC MATERIALS ===
      /**
       * Invariant $I_{1} \left( \mathbf{C} \right) = \mathbf{C} : \mathbf{I}
       * = \text{tr} \left( \mathbf{C} \right)$ where $\mathbf{C}$ is a rank-2
       * symmetric field tensor.
       */
      I1,
      /**
       * Invariant $I_{2} \left( \mathbf{C} \right) = \frac{1}{2} \left[
       * \left[ \mathbf{C} : \mathbf{I} \right]^{2} + \mathbf{C}^{2} :
       * \mathbf{I}  \right] = \frac{1}{2} \left[ \left[ \text{tr} \left(
       * \mathbf{C} \right) \right]^{2} + \left[ \mathbf{C} : \mathbf{C}
       * \right]  \right] $ where $\mathbf{C}$ is a rank-2 symmetric field
       * tensor.
       */
      I2,
      /**
       * Invariant $I_{3} \left( \mathbf{C} \right) = \text{det} \left(
       * \mathbf{C} \right)$ where $\mathbf{C}$ is a rank-2 symmetric field
       * tensor.
       *
       * If one interprets $\mathbf{C} = \mathbf{F}^{T} \cdot \mathbf{F}$ as
       * The rank-2 symmetric field tensor $\mathbf{C}$, where $\mathbf{F}$ is
       * the deformation gradient tensor, then $I_{3} \equiv \text{det} \left(
       * \mathbf{F} \right)^{2} = J^{2}$ which is the square of the volumetric
       * Jacobian. When a material exhibits a point-wise incompressible response
       * then $I_{3} \left( \mathbf{C} \right) = 1$.
       */
      I3,
      /**
       * Pseudo-invariant $I_{3} \left( \mathbf{C} \right) = \sqrt{\text{det}
       * \left( \mathbf{C} \right)}$ where $\mathbf{C}$ is a rank-2 symmetric
       * field tensor.
       *
       * If one interprets $\mathbf{C} = \mathbf{F}^{T} \cdot \mathbf{F}$ as
       * The rank-2 symmetric field tensor $\mathbf{C}$, where $\mathbf{F}$ is
       * the deformation gradient tensor, then $I_{3a} \equiv \text{det} \left(
       * \mathbf{F} \right) = J$ which is the volumetric Jacobian. When a
       * material exhibits a point-wise incompressible response then $I_{3a}
       * \left( \mathbf{C} \right) = 1$.
       */
      pI3,
      // === TRANSVERSE ISOTROPIC MATERIALS (1 preferred direction) ===
      /**
       * Invariant $I_{4} \left( \mathbf{C}, \mathbf{G} \right)
       * = \mathbf{C} : \mathbf{G}$  where both $\mathbf{C}$ and $\mathbf{G}$
       * are rank-2 symmetric tensors ($\mathbf{C}$ being a field tensor).
       *
       * The tensor $\mathbf{G}$ is known as the structure tensor that is
       * associated with some preferred (referential) direction
       * $\mathrm{N}$, with $\mathbf{G} = \mathrm{N} \otimes \mathrm{N}$ in
       * the simplest case of transversely isotropic materials with a strongly
       * preferential orientation. If $\mathbf{C}$ represents the right
       * Cauchy-Green deformation tensor, then $I_{4} \left( \mathbf{C},
       * \mathbf{G} \right) = \lambda\left( \mathrm{N} \right)^{2}$; that is to
       * say that this invariant measures the square of the stretch in the
       * preferred direction.
       *
       * For the case of dispersed transversely isotropic materials, then in 3-d
       * the structure tensor can be composed from $\mathbf{G} = \kappa
       * \mathbf{I} + \left( 1 - 3\kappa \right) \mathrm{N} \otimes \mathrm{N}$
       * where $\kappa \in \left[0, \frac{1}{3} \right]$ is the dispersion
       * parameter and $\mathrm{N}$ is the mean direction of anisotropy.
       * When $\kappa = 0$ then there is no dispersion and the transversely
       * isotropic behaviour is recovered, while $\kappa = \frac{1}{3}$ is
       * the limit at which the material behaves isotropically. Intermediate
       * values allow control over the degree of anisotropy exhibited by the
       * materials. This is useful in modelling, for instance, biological
       * materials or composite fibrous media that do not have a perfect
       * orientation of all fibers, but rather a statistical distribution away
       * from the mean direction.
       *
       * TODO[JPP]: References: Ogden, Saxena, ...
       */
      I4,
      /**
       * Invariant $I_{5} \left( \mathbf{C}, \mathbf{G} \right)
       * = \mathbf{C}^{2} : \mathbf{G}$  where both $\mathbf{C}$ and
       * $\mathbf{G}$ are rank-2 symmetric tensors ($\mathbf{C}$ being a field
       * tensor, $\mathbf{G}$ being a structure tensor).
       *
       * If $\mathbf{G} = \mathrm{N} \otimes \mathrm{N}$ then this can also be
       * expressed as $I_{5} \left( \mathbf{C}, \mathrm{N} \right) = \left[
       * \mathbf{C} \cdot \mathrm{N} \right] \cdot \left[ \mathbf{C} \cdot
       * \mathrm{N} \right]$.
       */
      I5,
      /**
       * Pseudo-invariant $I_{5} \left( \mathbf{C}, \mathbf{G} \right)
       * = \left[ I_{3} \mathbf{C}^{-1} \right] : \mathbf{G}$ where both
       * $\mathbf{C}$ and $\mathbf{G}$ are rank-2 symmetric tensors
       * ($\mathbf{C}$ being a field tensor, $\mathbf{G}$ being a structure
       * tensor).
       *
       * If $\mathbf{G} = \mathrm{N} \otimes \mathrm{N}$ then this can also be
       * expressed as $I_{5} \left( \mathbf{C}, \mathrm{N} \right) = \left[
       * I_{3} \mathbf{C}^{-1} \cdot \mathrm{N} \right] \cdot \mathrm{N}$, and
       * if $\mathbf{C}$ represents The rank-2 symmetric field tensor
       * $\mathbf{C}$ then this can be physically interpreted as the square of
       * the ratio of the deformed to undeformed area element initially normal
       * to the direction $\mathrm{N}$.
       *
       * TODO[JPP]: Cite Ogden slides: RWO-talk.pdf
       */
      pI5,
      // === ORTHOTROPIC MATERIALS (2 preferred directions) ===
      /**
       * Invariant $I_{6} \left( \mathbf{C}, \mathbf{G}' \right)
       * = \mathbf{C} : \mathbf{G}'$  where both $\mathbf{C}$ and
       * $\mathbf{G}'$ are rank-2 symmetric tensors ($\mathbf{C}$ being a
       * field tensor, $\mathbf{G}'$ being a structure tensor).
       *
       * The structure of $I_{6}$ is the same as $I_{4}$, but the structure
       * tensor $\mathbf{G}'$ is associated with another preferred direction
       * to that of $\mathbf{G}$.
       */
      I6,
      /**
       * Invariant $I_{7} \left( \mathbf{C}, \mathbf{G}' \right)
       * = \mathbf{C}^{2} : \mathbf{G}'$  where both $\mathbf{C}$ and
       * $\mathbf{G}'$ are rank-2 symmetric tensors ($\mathbf{C}$ being a field
       * tensor, $\mathbf{G}'$ being a structure tensor).
       *
       * The structure of $I_{7}$ is the same as $I_{5}$, but the structure
       * tensor $\mathbf{G}'$ is associated with another preferred direction to
       * that of $\mathbf{G}$.
       */
      I7,
      /**
       * Pseudo-invariant $I_{5} \left( \mathbf{C}, \mathbf{G}' \right)
       * = \left[ I_{3} \mathbf{C}^{-1} \right] : \mathbf{G}'$ where both
       * $\mathbf{C}$ and $\mathbf{G}'$ are rank-2 symmetric tensors
       * ($\mathbf{C}$ being a field tensor, $\mathbf{G}'$ being a structure
       * tensor).
       *
       * The structure of pseudo-invariant $I_{7}$ is the same as
       * pseudo-invariant $I_{5}$, but the structure tensor
       * $\mathbf{G}'$ is associated with another preferred direction to
       * that of $\mathbf{G}$.
       */
      pI7,
      /**
       * Invariant $I_{8} \left( \mathbf{C}, \mathbf{G}, \mathbf{G}' \right)
       * = \mathbf{I} : \left[ \mathbf{G} \cdot \mathbf{C} \cdot \mathbf{G}'
       * \right] = \text{tr} \left( \mathbf{G} \cdot \mathbf{C} \cdot \mathbf{G}
       * \right)$  where both $\mathbf{C}$, $\mathbf{G}$ and $\mathbf{G}'$ are
       * rank-2 symmetric tensors ($\mathbf{C}$ being a field tensor, while
       * $\mathbf{G}$ and $\mathbf{G}'$ are both structure tensors).
       *
       * If the two structure tensors can be decomposed as $\mathbf{G} =
       * \mathrm{N} \otimes \mathrm{N}$ and $\mathbf{G}' = \mathrm{N}' \otimes
       * \mathrm{N}'$, then $I_{8} \left( \mathbf{C}, \mathbf{N}, \mathbf{N}'
       * \right) = \left[ \mathbf{N} \cdot \mathbf{C} \cdot \mathbf{N}' \right]
       * \left[ \mathbf{N} \cdot \mathbf{N}' \right]$.
       */
      I8,
      // === COUPLED ISOTROPIC MATERIALS ===
      /**
       * Invariant $I_{9} \left( \mathrm{H} \right) = \left[ \mathrm{H}
       * \otimes \mathrm{H} \right] : \mathbf{I} = \mathrm{H} \cdot
       * \mathrm{H}$ where $\mathrm{H}$ is a rank-1 field tensor.
       */
      I9,
      /**
       * Coupled invariant $I_{10} \left( \mathrm{H}, \mathbf{C} \right) =
       * \left[ \mathrm{H} \otimes \mathrm{H} \right] : \mathbf{C}$ where
       * $\mathrm{H}$ is a rank-1 field tensor and $\mathbf{C}$ is a rank-2
       * symmetric field tensor.
       */
      I10,
      /**
       * Coupled invariant $I_{11} \left( \mathrm{H}, \mathbf{C} \right) =
       * \left[ \mathrm{H} \otimes \mathrm{H} \right] : \mathbf{C}^{2}
       * = \left[ \mathbf{C} \cdot \mathrm{H} \right] \cdot \left[ \mathbf{C}
       * \cdot \mathrm{H} \right]$ where
       * $\mathrm{H}$ is a rank-1 field tensor and $\mathbf{C}$ is a rank-2
       * symmetric field tensor.
       */
      I11,
      /**
       * Coupled pseudo-invariant $I_{11a} \left( \mathrm{H}, \mathbf{C}
       * \right) = \left[ \mathrm{H} \otimes \mathrm{H} \right] :
       * \mathbf{C}^{-1}$ where $\mathrm{H}$ is a rank-1 field tensor and
       * $\mathbf{C}$ is a rank-2 symmetric field tensor.
       *
       * This is considered a pseudo-invariant as, through the Cayley-Hamilton
       * theorem, this invariant can be expressed in terms of the original
       * irreducible invariants: $I_{11a} = \frac{1}{I_{3}} \left[
       * I_{11} - I_{1} I_{10} + I_{2} I_{9} \right]$
       */
      pI11a,
      /**
       * Coupled pseudo-invariant $I_{11b} \left( \mathrm{H}, \mathbf{C}
       * \right) = I_{3a} I_{11a} = \left[ \mathrm{H} \otimes \mathrm{H}
       * \right] : \sqrt{\text{det} \left( \mathbf{C} \right)} \mathbf{C}$
       * where $\mathrm{H}$ is a rank-1 field tensor and $\mathbf{C}$ is a
       * rank-2 symmetric field tensor.
       *
       * In magneto-mechanics, if one interprets $\mathrm{H}$ as the magnetic
       * field vector and $\mathbf{C}$ as the right Cauchy-Green deformation
       * tensor, then $I_{7b}$ is the field-dependent component of the
       * magnetic energy stored in the free field.
       */
      pI11b,
      // === COUPLED TRANSVERSE ISOTROPIC MATERIALS (1 preferred direction) ===
      /**
       * Invariant $I_{12} \left( \mathrm{H}, \mathbf{G} \right) = \left[
       * \mathrm{H} \otimes \mathrm{H} \right] : \mathbf{G} = \mathrm{H}
       * \cdot \mathbf{G} \cdot \mathrm{H}$ where $\mathrm{H}$ is a rank-1
       * field tensor and $\mathbf{G}$ is a rank-2 symmetric structure
       * tensor.
       */
      I12,
      /**
       * Invariant $I_{13} \left( \mathrm{H}, \mathbf{C}, \mathbf{G} \right) =
       * \left[ \mathrm{H} \otimes \mathrm{H} \right] : \left[ \mathbf{C} \cdot
       * \mathbf{G} \cdot \mathbf{C} \right] = \left[ \mathbf{C} \cdot
       * \mathrm{H} \right] \cdot \mathbf{G} \cdot \left[ \mathbf{C} \cdot
       * \mathrm{H} \right]$ where $\mathrm{H}$ is a rank-1 field tensor,
       * $\mathbf{C}$ is a rank-2 symmetric field tensor and $\mathbf{G}$ is a
       * rank-2 structure tensor.
       */
      I13
    };


    namespace internal
    {
      /**
       * A function to help decide whether to add the contribution
       * that has a given invariant contribution as a scalar
       * coefficient. This relates to the efficiency of
       * computation while not losing track of the sensitivities
       * necessary to perform efficient SD/AD computations.
       */
      template <typename ScalarType, typename T = void>
      bool
      add_invariant_contribution(const ScalarType &value);


#ifndef DOXYGEN

      template <
        typename ScalarType,
        typename std::enable_if<std::is_arithmetic<ScalarType>::value>::type>
      bool
      add_invariant_contribution(const ScalarType &value)
      {
        // Add the (floating point) contribution only if its
        // non-zero valued.
        return value != ScalarType(0.0);
      }


      template <
        typename ScalarType,
        typename std::enable_if<!std::is_arithmetic<ScalarType>::value>::type>
      bool
      add_invariant_contribution(const ScalarType &value)
      {
        // Always add the contribution if the ScalarType a not a
        // floating point type. This way, if ScalarType is an AD
        // or SD type, then we ensure that we always track the
        // sensitivities and later compute the correct derivatives
        // of this contribution.
        return true;
      }


      template <
        typename ScalarType,
        typename std::enable_if<std::is_arithmetic<ScalarType>::value>::type>
      bool
      add_invariant_contribution(const VectorizedArray<ScalarType> &values)
      {
        // Add the (floating point) contribution only all of the
        // data is non-zero valued.
        return values != VectorizedArray<ScalarType>(0.0);
      }

#endif

      /**
       * Returns the value of the tensor product a symmetric tensor
       * with another symmetric tensor.
       *
       * @param[in] S A rank-2 symmetric tensor
       * @param[in] T A rank-2 symmetric tensor
       * @return The symmetric tensor $\mathbf{S} \cdot \mathbf{T}$
       */
      template <int dim, typename ScalarType>
      Tensor<2, dim, ScalarType>
      compute_ST_dot_ST(const SymmetricTensor<2, dim, ScalarType> &S,
                        const SymmetricTensor<2, dim, ScalarType> &T)
      {
        Tensor<2, dim, ScalarType> S_dot_T;
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int j = 0; j < dim; ++j)
            for (unsigned int k = 0; k < dim; ++k)
              S_dot_T[i][j] += S[i][k] * T[k][j];
        return S_dot_T;
      }

      /**
       * Returns the value of the tensor product a symmetric tensor
       * with itself.
       *
       * @param[in] C A rank-2 symmetric tensor
       * @return The symmetric tensor $\mathbf{C} \cdot \mathbf{C}$
       */
      template <int dim, typename ScalarType>
      SymmetricTensor<2, dim, ScalarType>
      compute_ST_squared(const SymmetricTensor<2, dim, ScalarType> &C)
      {
        SymmetricTensor<2, dim, ScalarType> C_squared;
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int j = i; j < dim; ++j)
            for (unsigned int k = 0; k < dim; ++k)
              C_squared[i][j] += C[i][k] * C[k][j];
        return C_squared;
      }

      /**
       * Computes the value of the derivative $\frac{d \mathbf{C}^{-1}}{d
       * \mathbf{C}}$, given the inverse of a rank-2 symmetric tensor
       * $\mathbf{C}$.
       *
       * @param C_inv The tensor that is the inverse of the rank-2 symmetric field tensor $\mathbf{C}$.
       * @return The rank-4 symmetric tensor that represents $\frac{d \mathbf{C}^{-1}}{d \mathbf{C}}$.
       */
      template <int dim, typename ScalarType>
      SymmetricTensor<4, dim, ScalarType>
      get_dC_inv_dC(const SymmetricTensor<2, dim, ScalarType> &C_inv)
      {
        SymmetricTensor<4, dim, ScalarType> dC_inv_dC;

        for (unsigned int A = 0; A < dim; ++A)
          for (unsigned int B = A; B < dim; ++B)
            for (unsigned int C = 0; C < dim; ++C)
              for (unsigned int D = C; D < dim; ++D)
                dC_inv_dC[A][B][C][D] -=
                  0.5 * (C_inv[A][C] * C_inv[B][D] + C_inv[A][D] * C_inv[B][C]);

        return dC_inv_dC;
      }
    } // namespace internal



    // /**
    //  * A simple container to hold the first derivatives of some function,
    //  * expressed in terms of a set of invariants. Each directional derivative
    //  * $i$ is
    //  * paired against a @p coefficient value.
    //  *
    //  * Consider the function $f\left( I_{1}, I_{2}, ..., I_{n} \right)$ that
    //  is
    //  * parameterised by a set of invariants. The first derivative of this
    //  * function with respect to some general tensor $\mathbf{A}$ would then
    //  be
    //  * @f[
    //  *   \frac{d f}{d \mathbf{A}}
    //  *     = \frac{d f}{d I_{1}} \frac{d I_{1}}{d \mathbf{A}}
    //  *     + \frac{d f}{d I_{2}} \frac{d I_{2}}{d \mathbf{A}}
    //  *     + ...
    //  *     + \frac{d f}{d I_{n}} \frac{d I_{n}}{d \mathbf{A}}
    //  *     = \sum_{\alpha}
    //  *         \frac{d f}{d I_{\alpha}} \frac{d I_{\alpha}}{d \mathbf{A}} .
    //  * @f]
    //  *
    //  * This data structure may be used to hold the coefficients $\frac{d f}{d
    //  * I_{\alpha}}$ (associated to the derivative in the direction of the
    //  * $I_{\alpha}$-th invariant).
    //  *
    //  * @tparam InvariantType An enumeration type for the invariant.
    //  * @tparam ScalarType The scalar type for the derivative coefficient.
    //  */
    // template <enum InvariantType, typename ScalarType = double>
    // struct FirstDerivative
    // {
    //   /**
    //    * @brief Constructor.
    //    */
    //   FirstDerivative(const InvariantType direction_i,
    //                   const ScalarType    coefficient)
    //     : direction_i(direction_i)
    //     , coefficient(coefficient)
    //   {}

    //   /**
    //    * The directional derivative that this structure holds the coefficient
    //    * of, i.e., the selected "$\alpha$".
    //    */
    //   const InvariantType direction_i;
    //   /**
    //    * The coefficient of the first derivative term, i.e., the value of
    //    * $\frac{d f}{d I_{\alpha}}$.
    //    */
    //   const ScalarType coefficient;
    // };



    // /**
    //  * A simple container to hold the second derivatives of some function,
    //  * expressed in terms of a set of invariants. Each directional derivative
    //  * $ij$ is
    //  * paired against a @p coefficient value.
    //  *
    //  * Consider the function $f\left( I_{1}, I_{2}, ..., I_{n} \right)$ that
    //  is
    //  * parameterised by a set of invariants. The first derivative of this
    //  * function with respect to some general tensor $\mathbf{A}$ would then
    //  be
    //  * @f[
    //  *   \frac{d f}{d \mathbf{A}}
    //  *     = \frac{d f}{d I_{1}} \frac{d I_{1}}{d \mathbf{A}}
    //  *     + \frac{d f}{d I_{2}} \frac{d I_{2}}{d \mathbf{A}}
    //  *     + ...
    //  *     + \frac{d f}{d I_{n}} \frac{d I_{n}}{d \mathbf{A}}
    //  *     = \sum_{\alpha}
    //  *         \frac{d f}{d I_{\alpha}} \frac{d I_{\alpha}}{d \mathbf{A}} .
    //  * @f]
    //  * By application of the chain rule, the second derivatives of the
    //  * function would then be
    //  * @f[
    //  *   \frac{d^{2} f}{d \mathbf{A} \otimes d \mathbf{A}}
    //  *     = \frac{d I_{1}}{d \mathbf{A}} \otimes \left[
    //  *         \frac{d^{2} f}{d I_{1}^{2}} \frac{d I_{1}}{d \mathbf{A}}
    //  *       + \frac{d^{2} f}{d I_{1} d I_{1}} \frac{d I_{2}}{d \mathbf{A}}
    //  *       + ...
    //  *       + \frac{d^{2} f}{d I_{1} d I_{n}} \frac{d I_{n}}{d \mathbf{A}}
    //  *     \right]
    //  *     + \frac{d f}{d I_{1}}
    //  *         \frac{d^{2} I_{1}}{d \mathbf{A} \otimes d \mathbf{A}}
    //  *     + ...
    //  *     = \sum_{\alpha} \left[
    //  *         \sum_{\beta} \frac{d^{2} f}{d I_{\beta} d I_{\alpha}}
    //  *           \frac{d I_{\alpha}}{d \mathbf{A}}
    //  *             \otimes \frac{d I_{\beta}}{d \mathbf{A}}
    //  *         + \frac{d f}{d I_{\alpha}}
    //  *             \frac{d^{2} I_{\alpha}}{d \mathbf{A} \otimes d \mathbf{A}}
    //  *       \right] .
    //  * @f]
    //  *
    //  * This data structure may be used to hold the coefficients $\frac{d^{2}
    //  * f}{d I_{\beta} d I_{\alpha}} = \frac{d}{dI_{\beta}} \left[ \frac{d
    //  f}{d
    //  * I_{\alpha}} \right]$ (associated to the derivative in the direction of
    //  * the $I_{\beta}$-th invariant of the first derivative that is already
    //  * taken in the direction of the $I_{\alpha}$-th invariant).
    //  *
    //  * @tparam InvariantType An enumeration type for the invariant.
    //  * @tparam ScalarType The scalar type for the derivative coefficient.
    //  */
    // template <enum InvariantType, typename ScalarType = double>
    // struct SecondDerivative
    // {
    //   SecondDerivative(const InvariantType direction_i,
    //                    const InvariantType direction_j,
    //                    const ScalarType    coefficient)
    //     : direction_i(direction_i)
    //     , direction_j(direction_j)
    //     , coefficient(coefficient)
    //   {}

    //   /**
    //    * The first direction for the directional derivative that this
    //    structure
    //    * holds the coefficient of, i.e., the selected "$\beta$".
    //    */
    //   const InvariantType direction_i;
    //   /**
    //    * The second direction for the directional derivative that this
    //    structure
    //    * holds the coefficient of, i.e., the selected "$\alpha$".
    //    */
    //   const InvariantType direction_j;
    //   /**
    //    * The coefficient of the second derivative term, i.e., the value of
    //    * $\frac{d^{2} f}{d I_{\beta} d I_{\alpha}}$.
    //    */
    //   const ScalarType coefficient;
    // };


    namespace Isotropic
    {
      /**
       * Get the set of invariants that are valid for isotropic media.
       */
      std::set<InvariantList>
      valid_invariants()
      {
        return {I1, I2, I3, pI3};
      }


      /**
       * Returns the selected invariant/pseudo-invariant value.
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       */
      template <int dim, typename ScalarType>
      ScalarType
      Ii(const enum InvariantList &                 i,
         const SymmetricTensor<2, dim, ScalarType> &C)
      {
        Assert(valid_invariants().count(i) != 0,
               ExcMessage(
                 "The selected invariant is not a valid isotropic invariant."));

        switch (i)
          {
            case I1:
              {
                // I1 = tr(C)
                return trace(C);
              }
              break;
            case I2:
              {
                // I2 = 0.5*(tr(C)*tr(C) - tr(C*C))
                const ScalarType tr_C     = trace(C);
                const ScalarType C_ddot_C = scalar_product(C, C);
                return 0.5 * (tr_C * tr_C - C_ddot_C);
              }
              break;
            case I3:
              {
                // I3 = det(C)
                return determinant(C);
              }
              break;
            case pI3:
              {
                // pI3 = sqrt(det(C))
                using std::sqrt;
                return sqrt(determinant(C));
              }
              break;
            default:
              break;
          }

        AssertThrow(false,
                    ExcMessage("This (pseudo-)invariant is not defined."));
        return ScalarType();
      }


      /**
       * Returns the selected invariant/pseudo-invariant first
       * derivative with respect to the rank-2 symmetric field tensor
       * $\mathbf{C}$.
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] C_inv The inverse of the rank-2 symmetric field tensor
       * $\mathbf{C}$. If it is known that the required invariant derivative
       * does not require this value (that is typically expensive to compute),
       * then any other value (such as the identity tensor) can be used instead.
       */
      template <int dim, typename ScalarType>
      SymmetricTensor<2, dim, ScalarType>
      dIi_dC(const enum InvariantList &                 i,
             const SymmetricTensor<2, dim, ScalarType> &C,
             const SymmetricTensor<2, dim, ScalarType> &C_inv)
      {
        Assert(valid_invariants().count(i) != 0,
               ExcMessage(
                 "The selected invariant is not a valid isotropic invariant."));

        switch (i)
          {
            case I1:
              {
                // I1 = tr(C)
                return static_cast<SymmetricTensor<2, dim, ScalarType>>(
                  Physics::Elasticity::StandardTensors<dim>::I);
              }
              break;
            case I2:
              {
                // I2 = 0.5*(tr(C)*tr(C) - tr(C*C))
                const ScalarType I1 = Ii(Invariants::I1, C);
                return I1 * Physics::Elasticity::StandardTensors<dim>::I - C;
              }
              break;
            case I3:
              {
                // I3 = det(C)
                const ScalarType I3 = Ii(Invariants::I3, C);
                return I3 * C_inv;
              }
              break;
            case pI3:
              {
                // pI3 = sqrt(det(C)) = sqrt(I3)
                const ScalarType pI3 = Ii(Invariants::pI3, C);
                return (0.5 * pI3) * C_inv;
              }
              break;
            default:
              break;
          }

        AssertThrow(false,
                    ExcMessage("This (pseudo-)invariant is not defined."));
        return SymmetricTensor<2, dim, ScalarType>();
      }


      /**
       * Returns the selected invariant/pseudo-invariant second
       * derivative with respect to the rank-2 symmetric field tensor
       * $\mathbf{C}$.
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] C_inv The inverse of the rank-2 symmetric field tensor
       * $\mathbf{C}$. If it is known that the required invariant derivative
       * does not require this value (that is typically expensive to compute),
       * then any other value (such as the identity tensor) can be used instead.
       */
      template <int dim, typename ScalarType>
      SymmetricTensor<4, dim, ScalarType>
      d2Ii_dC_dC(const enum InvariantList &                 i,
                 const SymmetricTensor<2, dim, ScalarType> &C,
                 const SymmetricTensor<2, dim, ScalarType> &C_inv)
      {
        Assert(valid_invariants().count(i) != 0,
               ExcMessage(
                 "The selected invariant is not a valid isotropic invariant."));

        switch (i)
          {
            case I1:
              {
                // I1 = tr(C)
                return SymmetricTensor<4, dim, ScalarType>();
              }
              break;
            case I2:
              {
                // I2 = 0.5*(tr(C)*tr(C) - tr(C*C))
                return static_cast<SymmetricTensor<4, dim, ScalarType>>(
                  Physics::Elasticity::StandardTensors<dim>::IxI -
                  Physics::Elasticity::StandardTensors<dim>::S);
              }
              break;
            case I3:
              {
                // I3 = det(C)
                const ScalarType I3 = Ii(Invariants::I3, C);
                return I3 * (outer_product(C_inv, C_inv) +
                             internal::get_dC_inv_dC(C_inv));
              }
              break;
            case pI3:
              {
                // pI3 = sqrt(det(C)) = sqrt(I3)
                const ScalarType pI3 = Ii(Invariants::pI3, C);
                return (0.5 * pI3) * (0.5 * outer_product(C_inv, C_inv) +
                                      internal::get_dC_inv_dC(C_inv));
              }
              break;
            default:
              break;
          }

        AssertThrow(false,
                    ExcMessage("This (pseudo-)invariant is not defined."));
        return SymmetricTensor<4, dim, ScalarType>();
      }
    } // namespace Isotropic



    namespace Transverse_Isotropic
    {
      /**
       * Get the set of invariants that are valid for transverse isotropic
       * media.
       */
      std::set<InvariantList>
      valid_invariants()
      {
        return {I1, I2, I3, pI3, I4, I5, pI5};
      }


      /**
       * Returns the selected invariant/pseudo-invariant value.
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] G The rank-2 symmetric structure tensor $\mathbf{G}$
       */
      template <int dim, typename ScalarType, typename ScalarType2>
      ScalarType
      Ii(const enum InvariantList &                  i,
         const SymmetricTensor<2, dim, ScalarType> & C,
         const SymmetricTensor<2, dim, ScalarType> & C_inv,
         const SymmetricTensor<2, dim, ScalarType2> &G)
      {
        Assert(
          valid_invariants().count(i) != 0,
          ExcMessage(
            "The selected invariant is not a valid transverse isotropic invariant."));

        switch (i)
          {
            case I1:
              {
                return Isotropic::Ii(i, C);
              }
              break;
            case I2:
              {
                return Isotropic::Ii(i, C);
              }
              break;
            case I3:
              {
                return Isotropic::Ii(i, C);
              }
              break;
            case pI3:
              {
                return Isotropic::Ii(i, C);
              }
              break;
            case I4:
              {
                // I4 = C : G
                return C * G;
              }
              break;
            case I5:
              {
                // I5 = C^{2} : G
                return Invariants::internal::compute_ST_squared(C) * G;
              }
              break;
            case pI5:
              {
                // pI5 = I3 C^{-1} : G
                const ScalarType I3 = Isotropic::Ii(Invariants::I3, C);
                return I3 * (C_inv * G);
              }
              break;
            default:
              break;
          }

        AssertThrow(false,
                    ExcMessage("This (pseudo-)invariant is not defined."));
        return ScalarType();
      }


      /**
       * Returns the selected invariant/pseudo-invariant first
       * derivative with respect to the rank-2 symmetric field tensor
       * $\mathbf{C}$.
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] C_inv The inverse of the rank-2 symmetric field tensor
       * $\mathbf{C}$. If it is known that the required invariant derivative
       * does not require this value (that is typically expensive to compute),
       * then any other value (such as the identity tensor) can be used instead.
       * @param[in] G The rank-2 symmetric structure tensor $\mathbf{G}$
       */
      template <int dim, typename ScalarType, typename ScalarType2>
      SymmetricTensor<2, dim, ScalarType>
      dIi_dC(const enum InvariantList &                  i,
             const SymmetricTensor<2, dim, ScalarType> & C,
             const SymmetricTensor<2, dim, ScalarType> & C_inv,
             const SymmetricTensor<2, dim, ScalarType2> &G)
      {
        Assert(
          valid_invariants().count(i) != 0,
          ExcMessage(
            "The selected invariant is not a valid transverse isotropic invariant."));

        switch (i)
          {
            case I1:
              {
                return Isotropic::dIi_dC(i, C, C_inv);
              }
              break;
            case I2:
              {
                return Isotropic::dIi_dC(i, C, C_inv);
              }
              break;
            case I3:
              {
                return Isotropic::dIi_dC(i, C, C_inv);
              }
              break;
            case pI3:
              {
                return Isotropic::dIi_dC(i, C, C_inv);
              }
              break;
            case I4:
              {
                // I4 = C : G
                return G;
              }
              break;
            case I5:
              {
                // I5 = C^{2} : G
                SymmetricTensor<2, dim, ScalarType> dI5_dC;

                for (unsigned int A = 0; A < dim; ++A)
                  for (unsigned int B = A; B < dim; ++B)
                    for (unsigned int D = 0; D < dim; ++D)
                      dI5_dC[A][B] += C[A][D] * G[D][B] + G[A][D] * C[D][B];

                return dI5_dC; // contract<1,0>(C,G) + contract<1,0>(G,C);
              }
              break;
            case pI5:
              {
                // pI5 = I3 C^{-1} : G
                const ScalarType I3 = Isotropic::Ii(Invariants::I3, C);
                return (I3 * (C_inv * G)) * C_inv +
                       (internal::get_dC_inv_dC(C_inv) * (I3 * G));
              }
              break;
            default:
              break;
          }

        AssertThrow(false,
                    ExcMessage("This (pseudo-)invariant is not defined."));
        return SymmetricTensor<2, dim, ScalarType>();
      }


      /**
       * Returns the selected invariant/pseudo-invariant second
       * derivative with respect to the rank-2 symmetric field tensor
       * $\mathbf{C}$.
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] C_inv The inverse of the rank-2 symmetric field tensor
       * $\mathbf{C}$. If it is known that the required invariant derivative
       * does not require this value (that is typically expensive to compute),
       * then any other value (such as the identity tensor) can be used instead.
       * @param[in] G The rank-2 symmetric structure tensor $\mathbf{G}$
       */
      template <int dim, typename ScalarType, typename ScalarType2>
      SymmetricTensor<4, dim, ScalarType>
      d2Ii_dC_dC(const enum InvariantList &                  i,
                 const SymmetricTensor<2, dim, ScalarType> & C,
                 const SymmetricTensor<2, dim, ScalarType> & C_inv,
                 const SymmetricTensor<2, dim, ScalarType2> &G)
      {
        Assert(
          valid_invariants().count(i) != 0,
          ExcMessage(
            "The selected invariant is not a valid transverse isotropic invariant."));

        switch (i)
          {
            case I1:
              {
                return Isotropic::d2Ii_dC_dC(i, C, C_inv);
              }
              break;
            case I2:
              {
                return Isotropic::d2Ii_dC_dC(i, C, C_inv);
              }
              break;
            case I3:
              {
                return Isotropic::d2Ii_dC_dC(i, C, C_inv);
              }
              break;
            case pI3:
              {
                return Isotropic::d2Ii_dC_dC(i, C, C_inv);
              }
              break;
            case I4:
              {
                // I4 = C : G
                return SymmetricTensor<4, dim, ScalarType>();
              }
              break;
            case I5:
              {
                // I5 = C^{2} : G
                const SymmetricTensor<4, dim> S =
                  Physics::Elasticity::StandardTensors<dim>::S;
                SymmetricTensor<4, dim, ScalarType> d2I5_dC_dC;

                for (unsigned int A = 0; A < dim; ++A)
                  for (unsigned int B = A; B < dim; ++B)
                    for (unsigned int C = 0; C < dim; ++C)
                      for (unsigned int D = C; D < dim; ++D)
                        for (unsigned int E = 0; E < dim; ++E)
                          d2I5_dC_dC[A][B][C][D] +=
                            G[A][E] * S[E][B][C][D] + S[A][E][C][D] * G[E][B];

                return d2I5_dC_dC; // contract<1,0>(G,S) + contract<1,0>(S,G);
              }
              break;
            case pI5:
              {
                // pI5 = I3 C^{-1} : G
                const ScalarType I3 = Isotropic::Ii(Invariants::I3, C);
                const SymmetricTensor<4, dim, ScalarType> dC_inv_dC =
                  internal::get_dC_inv_dC(C_inv);
                const ScalarType C_inv_G = C_inv * G;

                // See internal::get_dC_inv_dC()
                // This is d/dC[C_inv : G]
                SymmetricTensor<2, dim, ScalarType> d_C_inv_ddot_G_dC;
                for (unsigned int C = 0; C < dim; ++C)
                  for (unsigned int D = C; D < dim; ++D)
                    for (unsigned int A = 0; A < dim; ++A)
                      for (unsigned int B = 0; B < dim; ++B)
                        d_C_inv_ddot_G_dC[C][D] -=
                          (0.5 * G[A][B]) * (C_inv[A][C] * C_inv[B][D] +
                                             C_inv[A][D] * C_inv[B][C]);

                // This is d/dC[d_C_inv_dC : G]
                SymmetricTensor<4, dim, ScalarType> d_dC_inv_dC_ddot_G_dC;
                for (unsigned int A = 0; A < dim; ++A)
                  for (unsigned int B = A; B < dim; ++B)
                    for (unsigned int E = 0; E < dim; ++E)
                      for (unsigned int F = E; F < dim; ++F)
                        for (unsigned int C = 0; C < dim; ++C)
                          for (unsigned int D = 0; D < dim; ++D)
                            d_dC_inv_dC_ddot_G_dC[A][B][E][F] -=
                              (0.5 * G[C][D]) *
                              (dC_inv_dC[A][C][E][F] * C_inv[B][D] +
                               C_inv[A][C] * dC_inv_dC[B][D][E][F] +
                               dC_inv_dC[A][D][E][F] * C_inv[B][C] +
                               C_inv[A][D] * dC_inv_dC[B][C][E][F]);

                // d/dC of [(I3 * (C_inv * G)) * C_inv + (dC_inv_dC * (I3*G))]
                // TODO[JPP]: Check if the analytical result is just the zero
                // tensor. The result of the test invariants_02 seems to be just
                // that, for a number of tested inputs.
                return I3 *
                       (outer_product(C_inv_G * C_inv + dC_inv_dC * G, C_inv) +
                        outer_product(C_inv, d_C_inv_ddot_G_dC) +
                        C_inv_G * dC_inv_dC + d_dC_inv_dC_ddot_G_dC);
              }
              break;
            default:
              break;
          }

        AssertThrow(false,
                    ExcMessage("This (pseudo-)invariant is not defined."));
        return SymmetricTensor<4, dim, ScalarType>();
      }
    } // namespace Transverse_Isotropic



    namespace Orthotropic
    {
      /**
       * Get the set of invariants that are valid for orthotropic media.
       */
      std::set<InvariantList>
      valid_invariants()
      {
        return {I1, I2, I3, pI3, I4, I5, pI5, I6, I7, pI7, I8};
      }


      /**
       * Returns the selected invariant/pseudo-invariant value.
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] G1 The rank-2 symmetric structure tensor $\mathbf{G}$
       * @param[in] G2 The rank-2 symmetric structure tensor $\mathbf{G}'$
       */
      template <int dim, typename ScalarType, typename ScalarType2>
      ScalarType
      Ii(const enum InvariantList &                  i,
         const SymmetricTensor<2, dim, ScalarType> & C,
         const SymmetricTensor<2, dim, ScalarType> & C_inv,
         const SymmetricTensor<2, dim, ScalarType2> &G1,
         const SymmetricTensor<2, dim, ScalarType2> &G2)
      {
        Assert(
          valid_invariants().count(i) != 0,
          ExcMessage(
            "The selected invariant is not a valid orthotropic invariant."));

        switch (i)
          {
            case I1:
              {
                return Isotropic::Ii(i, C);
              }
              break;
            case I2:
              {
                return Isotropic::Ii(i, C);
              }
              break;
            case I3:
              {
                return Isotropic::Ii(i, C);
              }
              break;
            case pI3:
              {
                return Isotropic::Ii(i, C);
              }
              break;
            case I4:
              {
                return Transverse_Isotropic::Ii(i, C, C_inv, G1);
              }
              break;
            case I5:
              {
                return Transverse_Isotropic::Ii(i, C, C_inv, G1);
              }
              break;
            case pI5:
              {
                return Transverse_Isotropic::Ii(i, C, C_inv, G1);
              }
              break;
            case I6:
              {
                // I6 = C : G2 == I4(G2)
                return Transverse_Isotropic::Ii(I4, C, C_inv, G2);
              }
              break;
            case I7:
              {
                // I7 = C^{2} : G2 == I5(G2)
                return Transverse_Isotropic::Ii(I5, C, C_inv, G2);
              }
              break;
            case pI7:
              {
                // pI7 = I3 C^{-1} : G2 == pI5(G2)
                return Transverse_Isotropic::Ii(pI5, C, C_inv, G2);
              }
              break;
            case I8:
              {
                // I8 = I : [G . C . G']
                ScalarType I8 =
                  dealii::internal::NumberType<ScalarType>::value(0.0);

                for (unsigned int A = 0; A < dim; ++A)
                  for (unsigned int B = 0; B < dim; ++B)
                    for (unsigned int D = 0; D < dim; ++D)
                      I8 += G1[D][A] * C[A][B] * G2[B][D];

                return I8;
              }
              break;
            default:
              break;
          }

        AssertThrow(false,
                    ExcMessage("This (pseudo-)invariant is not defined."));
        return ScalarType();
      }


      /**
       * Returns the selected invariant/pseudo-invariant first
       * derivative with respect to the rank-2 symmetric field tensor
       * $\mathbf{C}$.
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] C_inv The inverse of the rank-2 symmetric field tensor
       * $\mathbf{C}$. If it is known that the required invariant derivative
       * does not require this value (that is typically expensive to compute),
       * then any other value (such as the identity tensor) can be used instead.
       * @param[in] G1 The rank-2 symmetric structure tensor $\mathbf{G}$
       * @param[in] G2 The rank-2 symmetric structure tensor $\mathbf{G}'$
       */
      template <int dim, typename ScalarType, typename ScalarType2>
      SymmetricTensor<2, dim, ScalarType>
      dIi_dC(const enum InvariantList &                  i,
             const SymmetricTensor<2, dim, ScalarType> & C,
             const SymmetricTensor<2, dim, ScalarType> & C_inv,
             const SymmetricTensor<2, dim, ScalarType2> &G1,
             const SymmetricTensor<2, dim, ScalarType2> &G2)
      {
        Assert(
          valid_invariants().count(i) != 0,
          ExcMessage(
            "The selected invariant is not a valid orthotropic invariant."));

        switch (i)
          {
            case I1:
              {
                return Isotropic::dIi_dC(i, C, C_inv);
              }
              break;
            case I2:
              {
                return Isotropic::dIi_dC(i, C, C_inv);
              }
              break;
            case I3:
              {
                return Isotropic::dIi_dC(i, C, C_inv);
              }
              break;
            case pI3:
              {
                return Isotropic::dIi_dC(i, C, C_inv);
              }
              break;
            case I4:
              {
                return Transverse_Isotropic::dIi_dC(i, C, C_inv, G1);
              }
              break;
            case I5:
              {
                return Transverse_Isotropic::dIi_dC(i, C, C_inv, G1);
              }
              break;
            case pI5:
              {
                return Transverse_Isotropic::dIi_dC(i, C, C_inv, G1);
              }
              break;
            case I6:
              {
                // I6 = C : G2 == I4(G2)
                return Transverse_Isotropic::dIi_dC(I4, C, C_inv, G2);
              }
              break;
            case I7:
              {
                // I7 = C^{2} : G2 == I5(G2)
                return Transverse_Isotropic::dIi_dC(I5, C, C_inv, G2);
              }
              break;
            case pI7:
              {
                // pI7 = I3 C^{-1} : G2 == pI5(G2)
                return Transverse_Isotropic::dIi_dC(pI5, C, C_inv, G2);
              }
              break;
            case I8:
              {
                // I8 = I : [G . C . G']
                SymmetricTensor<2, dim, ScalarType> dI8_dC;

                for (unsigned int A = 0; A < dim; ++A)
                  for (unsigned int B = A; B < dim; ++B)
                    for (unsigned int D = 0; D < dim; ++D)
                      dI8_dC[A][B] +=
                        0.5 * (G1[D][A] * G2[B][D] + G1[D][B] * G2[A][D]);

                return dI8_dC;
              }
              break;
            default:
              break;
          }

        AssertThrow(false,
                    ExcMessage("This (pseudo-)invariant is not defined."));
        return SymmetricTensor<2, dim, ScalarType>();
      }


      /**
       * Returns the selected invariant/pseudo-invariant second
       * derivative with respect to the rank-2 symmetric field tensor
       * $\mathbf{C}$.
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] C_inv The inverse of the rank-2 symmetric field tensor
       * $\mathbf{C}$. If it is known that the required invariant derivative
       * does not require this value (that is typically expensive to compute),
       * then any other value (such as the identity tensor) can be used instead.
       * @param[in] G1 The rank-2 symmetric structure tensor $\mathbf{G}$
       * @param[in] G2 The rank-2 symmetric structure tensor $\mathbf{G}'$
       */
      template <int dim, typename ScalarType, typename ScalarType2>
      SymmetricTensor<4, dim, ScalarType>
      d2Ii_dC_dC(const enum InvariantList &                  i,
                 const SymmetricTensor<2, dim, ScalarType> & C,
                 const SymmetricTensor<2, dim, ScalarType> & C_inv,
                 const SymmetricTensor<2, dim, ScalarType2> &G1,
                 const SymmetricTensor<2, dim, ScalarType2> &G2)
      {
        Assert(
          valid_invariants().count(i) != 0,
          ExcMessage(
            "The selected invariant is not a valid orthotropic invariant."));

        switch (i)
          {
            case I1:
              {
                return Isotropic::d2Ii_dC_dC(i, C, C_inv);
              }
              break;
            case I2:
              {
                return Isotropic::d2Ii_dC_dC(i, C, C_inv);
              }
              break;
            case I3:
              {
                return Isotropic::d2Ii_dC_dC(i, C, C_inv);
              }
              break;
            case pI3:
              {
                return Isotropic::d2Ii_dC_dC(i, C, C_inv);
              }
              break;
            case I4:
              {
                return Transverse_Isotropic::d2Ii_dC_dC(i, C, C_inv, G1);
              }
              break;
            case I5:
              {
                return Transverse_Isotropic::d2Ii_dC_dC(i, C, C_inv, G1);
              }
              break;
            case pI5:
              {
                return Transverse_Isotropic::d2Ii_dC_dC(i, C, C_inv, G1);
              }
              break;
            case I6:
              {
                // I6 = C : G2 == I4(G2)
                return Transverse_Isotropic::d2Ii_dC_dC(I4, C, C_inv, G2);
              }
              break;
            case I7:
              {
                // I7 = C^{2} : G2 == I5(G2)
                return Transverse_Isotropic::d2Ii_dC_dC(I5, C, C_inv, G2);
              }
              break;
            case pI7:
              {
                // pI7 = I3 C^{-1} : G2 == pI5(G2)
                return Transverse_Isotropic::d2Ii_dC_dC(pI5, C, C_inv, G2);
              }
              break;
            case I8:
              {
                // I8 = I : [G . C . G']
                return SymmetricTensor<4, dim, ScalarType>();
              }
              break;
            default:
              break;
          }

        AssertThrow(false,
                    ExcMessage("This (pseudo-)invariant is not defined."));
        return SymmetricTensor<4, dim, ScalarType>();
      }
    } // namespace Orthotropic



    namespace Coupled_Isotropic
    {
      /**
       * Get the set of invariants that are valid for field-coupled isotropic
       * media.
       */
      std::set<InvariantList>
      valid_invariants()
      {
        return {I1, I2, I3, pI3, I9, I10, I11, pI11a, pI11b};
      }


      /**
       * Returns the selected invariant/pseudo-invariant value.
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] H The rank-1 field tensor $\mathrm{H}$
       */
      template <int dim, typename ScalarType>
      ScalarType
      Ii(const enum InvariantList &                 i,
         const SymmetricTensor<2, dim, ScalarType> &C,
         const SymmetricTensor<2, dim, ScalarType> &C_inv,
         const Tensor<1, dim, ScalarType> &         H)
      {
        Assert(
          valid_invariants().count(i) != 0,
          ExcMessage(
            "The selected invariant is not a valid coupled isotropic invariant."));

        switch (i)
          {
            case I1:
              {
                return Isotropic::Ii(i, C);
              }
              break;
            case I2:
              {
                return Isotropic::Ii(i, C);
              }
              break;
            case I3:
              {
                return Isotropic::Ii(i, C);
              }
              break;
            case pI3:
              {
                return Isotropic::Ii(i, C);
              }
              break;
            case I9:
              {
                // I9 = [HxH]:I
                return H * H;
              }
              break;
            case I10:
              {
                // I10 = [HxH]:C
                return contract3(H, C, H);
              }
              break;
            case I11:
              {
                // I11 = [HxH]:C^2 = [CxH].[CxH]
                const Tensor<1, dim, ScalarType> Z = C * H;
                return Z * Z;
              }
              break;
            case pI11a:
              {
                // pI11a = [HxH]:C^{-1}
                return contract3(H, C_inv, H);
              }
              break;
            case pI11b:
              {
                // pI11b = [HxH]:[sqrt(det(C)) C^{-1}]
                const ScalarType pI3 = Isotropic::Ii(Invariants::pI3, C);
                return pI3 * contract3(H, C_inv, H);
              }
              break;
            default:
              break;
          }

        AssertThrow(false,
                    ExcMessage("This (pseudo-)invariant is not defined."));
        return ScalarType();
      }


      /**
       * Returns the selected invariant/pseudo-invariant first
       * derivative with respect to the rank-2 symmetric field tensor
       * $\mathbf{C}$.
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] C_inv The inverse of the rank-2 symmetric field tensor
       * $\mathbf{C}$. If it is known that the required invariant derivative
       * does not require this value (that is typically expensive to compute),
       * then any other value (such as the identity tensor) can be used instead.
       * @param[in] H The rank-1 field tensor $\mathrm{H}$
       */
      template <int dim, typename ScalarType>
      SymmetricTensor<2, dim, ScalarType>
      dIi_dC(const enum InvariantList &                 i,
             const SymmetricTensor<2, dim, ScalarType> &C,
             const SymmetricTensor<2, dim, ScalarType> &C_inv,
             const Tensor<1, dim, ScalarType> &         H)
      {
        Assert(
          valid_invariants().count(i) != 0,
          ExcMessage(
            "The selected invariant is not a valid coupled isotropic invariant."));

        switch (i)
          {
            case I1:
              {
                return Isotropic::dIi_dC(i, C, C_inv);
              }
              break;
            case I2:
              {
                return Isotropic::dIi_dC(i, C, C_inv);
              }
              break;
            case I3:
              {
                return Isotropic::dIi_dC(i, C, C_inv);
              }
              break;
            case pI3:
              {
                return Isotropic::dIi_dC(i, C, C_inv);
              }
              break;
            case I9:
              {
                // I9 = [HxH]:I
                return SymmetricTensor<2, dim, ScalarType>();
              }
              break;
            case I10:
              {
                // I10 = [HxH]:C
                return symmetrize(outer_product(H, H));
              }
              break;
            case I11:
              {
                // I11 = [HxH]:C^2 = [CxH].[CxH]
                const Tensor<1, dim, ScalarType> Z = C * H;
                return 2.0 * symmetrize(outer_product(Z, H));
              }
              break;
            case pI11a:
              {
                // pI11a = [HxH]:C^{-1}
                const Tensor<1, dim, ScalarType> Y = C_inv * H;
                return -symmetrize(outer_product(Y, Y));
              }
              break;
            case pI11b:
              {
                // pI11b = [HxH]:[sqrt(det(C)) C^{-1}]
                const ScalarType pI3 = Isotropic::Ii(Invariants::pI3, C);
                const Tensor<1, dim, ScalarType> Y = C_inv * H;
                return pI3 * ((0.5 * (H * Y)) * C_inv -
                              symmetrize(outer_product(Y, Y)));
              }
              break;
            default:
              break;
          }

        AssertThrow(false,
                    ExcMessage("This (pseudo-)invariant is not defined."));
        return SymmetricTensor<2, dim, ScalarType>();
      }


      /**
       * Returns the selected invariant/pseudo-invariant first
       * derivative with respect to the rank-1 field tensor
       * $\mathrm{H}$.
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] C_inv The inverse of the rank-2 symmetric field tensor
       * $\mathbf{C}$. If it is known that the required invariant derivative
       * does not require this value (that is typically expensive to compute),
       * then any other value (such as the identity tensor) can be used instead.
       * @param[in] H The rank-1 field tensor $\mathrm{H}$
       */
      template <int dim, typename ScalarType>
      Tensor<1, dim, ScalarType>
      dIi_dH(const enum InvariantList &                 i,
             const SymmetricTensor<2, dim, ScalarType> &C,
             const SymmetricTensor<2, dim, ScalarType> &C_inv,
             const Tensor<1, dim, ScalarType> &         H)
      {
        Assert(
          valid_invariants().count(i) != 0,
          ExcMessage(
            "The selected invariant is not a valid coupled isotropic invariant."));

        switch (i)
          {
            case I1:
              {
                return Tensor<1, dim, ScalarType>();
              }
              break;
            case I2:
              {
                return Tensor<1, dim, ScalarType>();
              }
              break;
            case I3:
              {
                return Tensor<1, dim, ScalarType>();
              }
              break;
            case pI3:
              {
                return Tensor<1, dim, ScalarType>();
              }
              break;
            case I9:
              {
                // I9 = [HxH]:I
                return 2.0 * H;
              }
              break;
            case I10:
              {
                // I10 = [HxH]:C
                return 2.0 * (C * H);
              }
              break;
            case I11:
              {
                // I11 = [HxH]:C^2 = [CxH].[CxH]
                return 2.0 * (C * (C * H));
              }
              break;
            case pI11a:
              {
                // pI11a = [HxH]:C^{-1}
                return 2.0 * (C_inv * H);
              }
              break;
            case pI11b:
              {
                // pI11b = [HxH]:[sqrt(det(C)) C^{-1}]
                const ScalarType pI3 = Isotropic::Ii(Invariants::pI3, C);
                return (2.0 * pI3) * (C_inv * H);
              }
              break;
            default:
              break;
          }

        AssertThrow(false,
                    ExcMessage("This (pseudo-)invariant is not defined."));
        return Tensor<1, dim, ScalarType>();
      }


      /**
       * Returns the selected invariant/pseudo-invariant second
       * derivative with respect to the rank-2 symmetric field tensor
       * $\mathbf{C}$.
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] C_inv The inverse of the rank-2 symmetric field tensor
       * $\mathbf{C}$. If it is known that the required invariant derivative
       * does not require this value (that is typically expensive to compute),
       * then any other value (such as the identity tensor) can be used instead.
       * @param[in] H The rank-1 field tensor $\mathrm{H}$
       */
      template <int dim, typename ScalarType>
      SymmetricTensor<4, dim, ScalarType>
      d2Ii_dC_dC(const enum InvariantList &                 i,
                 const SymmetricTensor<2, dim, ScalarType> &C,
                 const SymmetricTensor<2, dim, ScalarType> &C_inv,
                 const Tensor<1, dim, ScalarType> &         H)
      {
        Assert(
          valid_invariants().count(i) != 0,
          ExcMessage(
            "The selected invariant is not a valid coupled isotropic invariant."));

        switch (i)
          {
            case I1:
              {
                return Isotropic::d2Ii_dC_dC(i, C, C_inv);
              }
              break;
            case I2:
              {
                return Isotropic::d2Ii_dC_dC(i, C, C_inv);
              }
              break;
            case I3:
              {
                return Isotropic::d2Ii_dC_dC(i, C, C_inv);
              }
              break;
            case pI3:
              {
                return Isotropic::d2Ii_dC_dC(i, C, C_inv);
              }
              break;
            case I9:
              {
                // I9 = [HxH]:I
                return SymmetricTensor<4, dim, ScalarType>();
              }
              break;
            case I10:
              {
                // I10 = [HxH]:C
                return SymmetricTensor<4, dim, ScalarType>();
              }
              break;
            case I11:
              {
                // I11 = [HxH]:C^2 = [CxH].[CxH]
                const SymmetricTensor<2, dim, double> &I =
                  Physics::Elasticity::StandardTensors<dim>::I;

                SymmetricTensor<4, dim, ScalarType> d2I11_dC_dC;
                for (unsigned int A = 0; A < dim; ++A)
                  for (unsigned int B = A; B < dim; ++B)
                    for (unsigned int C = 0; C < dim; ++C)
                      for (unsigned int D = C; D < dim; ++D)
                        {
                          // Need to ensure symmetries of (A,B) and (C,D)
                          d2I11_dC_dC[A][B][C][D] =
                            0.5 *
                            (H[A] * H[D] * I[B][C] + H[A] * H[C] * I[B][D] +
                             H[B] * H[D] * I[A][C] + H[B] * H[C] * I[A][D]);
                        }

                return d2I11_dC_dC;
              }
              break;
            case pI11a:
              {
                // pI11a = [HxH]:C^{-1}
                const Tensor<1, dim, ScalarType> Y = C_inv * H;

                SymmetricTensor<4, dim, ScalarType> d2I11a_dC_dC;
                for (unsigned int A = 0; A < dim; ++A)
                  for (unsigned int B = A; B < dim; ++B)
                    for (unsigned int C = 0; C < dim; ++C)
                      for (unsigned int D = C; D < dim; ++D)
                        {
                          // Need to ensure symmetries of (A,B) and (C,D)
                          d2I11a_dC_dC[A][B][C][D] =
                            0.5 * (Y[A] * Y[C] * C_inv[B][D] +
                                   Y[B] * Y[C] * C_inv[A][D] +
                                   Y[A] * Y[D] * C_inv[B][C] +
                                   Y[B] * Y[D] * C_inv[A][C]);
                        }

                return d2I11a_dC_dC;
              }
              break;
            case pI11b:
              {
                // pI11b = [HxH]:[sqrt(det(C)) C^{-1}]
                const ScalarType pI3 = Isotropic::Ii(Invariants::pI3, C);
                const Tensor<1, dim, ScalarType> Y                 = C_inv * H;
                const ScalarType                 H_dot_C_inv_dot_H = H * Y;
                const SymmetricTensor<2, dim, ScalarType> symm_Y_x_Y =
                  symmetrize(outer_product(Y, Y));
                const SymmetricTensor<4, dim, ScalarType> dC_inv_dC =
                  internal::get_dC_inv_dC(C_inv);
                const SymmetricTensor<4, dim, ScalarType> d2I11a_dC_dC =
                  Coupled_Isotropic::d2Ii_dC_dC(pI11a, C, C_inv, H);

                return pI3 * d2I11a_dC_dC +
                       outer_product((0.5 * H_dot_C_inv_dot_H) * C_inv -
                                       symm_Y_x_Y,
                                     (0.5 * pI3) * C_inv) +
                       (outer_product((0.5 * pI3) * C_inv, -symm_Y_x_Y) +
                        ((0.5 * pI3) * H_dot_C_inv_dot_H) * dC_inv_dC);
              }
              break;
            default:
              break;
          }

        AssertThrow(false,
                    ExcMessage("This (pseudo-)invariant is not defined."));
        return SymmetricTensor<4, dim, ScalarType>();
      }


      /**
       * Returns the selected invariant/pseudo-invariant second
       * derivative with respect to the rank-1 field tensor
       * $\mathrm{H}$.
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] C_inv The inverse of the rank-2 symmetric field tensor
       * $\mathbf{C}$. If it is known that the required invariant derivative
       * does not require this value (that is typically expensive to compute),
       * then any other value (such as the identity tensor) can be used instead.
       * @param[in] H The rank-1 field tensor $\mathrm{H}$
       */
      template <int dim, typename ScalarType>
      SymmetricTensor<2, dim, ScalarType>
      d2Ii_dH_dH(const enum InvariantList &                 i,
                 const SymmetricTensor<2, dim, ScalarType> &C,
                 const SymmetricTensor<2, dim, ScalarType> &C_inv,
                 const Tensor<1, dim, ScalarType> &         H)
      {
        Assert(
          valid_invariants().count(i) != 0,
          ExcMessage(
            "The selected invariant is not a valid coupled isotropic invariant."));

        switch (i)
          {
            case I1:
              {
                return SymmetricTensor<2, dim, ScalarType>();
              }
              break;
            case I2:
              {
                return SymmetricTensor<2, dim, ScalarType>();
              }
              break;
            case I3:
              {
                return SymmetricTensor<2, dim, ScalarType>();
              }
              break;
            case pI3:
              {
                return SymmetricTensor<2, dim, ScalarType>();
              }
              break;
            case I9:
              {
                // I9 = [HxH]:I
                return 2.0 * unit_symmetric_tensor<dim, ScalarType>();
              }
              break;
            case I10:
              {
                // I10 = [HxH]:C
                return 2.0 * C;
              }
              break;
            case I11:
              {
                // I11 = [HxH]:C^2 = [CxH].[CxH]
                return 2.0 * Invariants::internal::compute_ST_squared(C);
              }
              break;
            case pI11a:
              {
                // pI11a = [HxH]:C^{-1}
                return 2.0 * C_inv;
              }
              break;
            case pI11b:
              {
                // pI11b = [HxH]:[sqrt(det(C)) C^{-1}]
                const ScalarType pI3 = Isotropic::Ii(Invariants::pI3, C);
                return (2.0 * pI3) * C_inv;
              }
              break;
            default:
              break;
          }

        AssertThrow(false,
                    ExcMessage("This (pseudo-)invariant is not defined."));
        return SymmetricTensor<2, dim, ScalarType>();
      }


      /**
       * Returns the selected invariant/pseudo-invariant second
       * derivative, with the derivatives first taken with respect to the rank-2
       * symmetric field tensor
       * $\mathbf{C}$ followed by the rank-1 field tensor $\mathrm{H}$.
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] C_inv The inverse of the rank-2 symmetric field tensor
       * $\mathbf{C}$. If it is known that the required invariant derivative
       * does not require this value (that is typically expensive to compute),
       * then any other value (such as the identity tensor) can be used instead.
       * @param[in] H The rank-1 field tensor $\mathrm{H}$
       */
      template <int dim, typename ScalarType>
      Tensor<3, dim, ScalarType>
      d2Ii_dC_dH(const enum InvariantList &                 i,
                 const SymmetricTensor<2, dim, ScalarType> &C,
                 const SymmetricTensor<2, dim, ScalarType> &C_inv,
                 const Tensor<1, dim, ScalarType> &         H)
      {
        Assert(
          valid_invariants().count(i) != 0,
          ExcMessage(
            "The selected invariant is not a valid coupled isotropic invariant."));

        switch (i)
          {
            case I1:
              {
                return Tensor<3, dim, ScalarType>();
              }
              break;
            case I2:
              {
                return Tensor<3, dim, ScalarType>();
              }
              break;
            case I3:
              {
                return Tensor<3, dim, ScalarType>();
              }
              break;
            case pI3:
              {
                return Tensor<3, dim, ScalarType>();
              }
              break;
            case I9:
              {
                // I9 = [HxH]:I
                return Tensor<3, dim, ScalarType>();
              }
              break;
            case I10:
              {
                // I10 = [HxH]:C
                const SymmetricTensor<2, dim, double> &I =
                  Physics::Elasticity::StandardTensors<dim>::I;

                Tensor<3, dim, ScalarType> d2I10_dC_dH;
                for (unsigned int A = 0; A < dim; ++A)
                  for (unsigned int B = 0; B < dim; ++B)
                    for (unsigned int C = 0; C < dim; ++C)
                      d2I10_dC_dH[A][B][C] += I[C][A] * H[B] + I[C][B] * H[A];

                return d2I10_dC_dH;
              }
              break;
            case I11:
              {
                // I11 = [HxH]:C^2 = [CxH].[CxH]
                const SymmetricTensor<2, dim, double> &I =
                  Physics::Elasticity::StandardTensors<dim>::I;
                const Tensor<1, dim, ScalarType> Z = C * H;

                Tensor<3, dim, ScalarType> d2I11_dC_dH;
                for (unsigned int A = 0; A < dim; ++A)
                  for (unsigned int B = 0; B < dim; ++B)
                    for (unsigned int D = 0; D < dim; ++D)
                      d2I11_dC_dH[A][B][D] += I[A][D] * Z[B] + I[B][D] * Z[A] +
                                              C[B][D] * H[A] + C[A][D] * H[B];

                return d2I11_dC_dH;
              }
              break;
            case pI11a:
              {
                // pI11a = [HxH]:C^{-1}
                const Tensor<1, dim, ScalarType> Y = C_inv * H;

                Tensor<3, dim, ScalarType> d2pI11a_dC_dH;
                for (unsigned int A = 0; A < dim; ++A)
                  for (unsigned int B = 0; B < dim; ++B)
                    for (unsigned int C = 0; C < dim; ++C)
                      d2pI11a_dC_dH[A][B][C] -=
                        C_inv[C][A] * Y[B] + C_inv[C][B] * Y[A];

                return d2pI11a_dC_dH;
              }
              break;
            case pI11b:
              {
                // pI11b = [HxH]:[sqrt(det(C)) C^{-1}]
                const ScalarType pI3 = Isotropic::Ii(Invariants::pI3, C);
                const Tensor<1, dim, ScalarType> Y = C_inv * H;
                const Tensor<3, dim, ScalarType> d2pI11a_dC_dH =
                  Coupled_Isotropic::d2Ii_dC_dH(pI11a, C, C_inv, H);

                Tensor<3, dim, ScalarType> d2pI11b_dC_dH;
                for (unsigned int A = 0; A < dim; ++A)
                  for (unsigned int B = 0; B < dim; ++B)
                    for (unsigned int C = 0; C < dim; ++C)
                      d2pI11b_dC_dH[A][B][C] =
                        pI3 * (C_inv[A][B] * Y[C] + d2pI11a_dC_dH[A][B][C]);

                return d2pI11b_dC_dH;
              }
              break;
            default:
              break;
          }

        AssertThrow(false,
                    ExcMessage("This (pseudo-)invariant is not defined."));
        return Tensor<3, dim, ScalarType>();
      }


      /**
       * Returns the selected invariant/pseudo-invariant second
       * derivative, with the derivatives first taken with respect to the rank-1
       * field tensor $\mathrm{H}$ followed by the rank-2 symmetric field tensor
       * $\mathbf{C}$ .
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] C_inv The inverse of the rank-2 symmetric field tensor
       * $\mathbf{C}$. If it is known that the required invariant derivative
       * does not require this value (that is typically expensive to compute),
       * then any other value (such as the identity tensor) can be used instead.
       * @param[in] H The rank-1 field tensor $\mathrm{H}$
       */
      template <int dim, typename ScalarType>
      Tensor<3, dim, ScalarType>
      d2Ii_dH_dC(const enum InvariantList &                 i,
                 const SymmetricTensor<2, dim, ScalarType> &C,
                 const SymmetricTensor<2, dim, ScalarType> &C_inv,
                 const Tensor<1, dim, ScalarType> &         H)
      {
        Assert(
          valid_invariants().count(i) != 0,
          ExcMessage(
            "The selected invariant is not a valid coupled isotropic invariant."));

        const Tensor<3, dim, ScalarType> d2Ii_dC_dH =
          Coupled_Isotropic::d2Ii_dC_dH(i, C, C_inv, H);

        Tensor<3, dim, ScalarType> d2Ii_dH_dC;
        for (unsigned int A = 0; A < dim; ++A)
          for (unsigned int B = 0; B < dim; ++B)
            for (unsigned int C = 0; C < dim; ++C)
              d2Ii_dH_dC[A][B][C] = d2Ii_dC_dH[C][B][A];

        return d2Ii_dH_dC;
      }
    } // namespace Coupled_Isotropic



    namespace Coupled_Transverse_Isotropic
    {
      /**
       * Get the set of invariants that are valid for field-coupled transverse
       * isotropic media.
       */
      std::set<InvariantList>
      valid_invariants()
      {
        return {
          I1, I2, I3, pI3, I4, I5, pI5, I9, I10, I11, pI11a, pI11b, I12, I13};
      }


      /**
       * Returns the selected invariant/pseudo-invariant value.
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] H The rank-1 field tensor $\mathrm{H}$
       * @param[in] G The rank-2 symmetric structure tensor $\mathbf{G}$
       */
      template <int dim, typename ScalarType, typename ScalarType2>
      ScalarType
      Ii(const enum InvariantList &                  i,
         const SymmetricTensor<2, dim, ScalarType> & C,
         const SymmetricTensor<2, dim, ScalarType> & C_inv,
         const Tensor<1, dim, ScalarType> &          H,
         const SymmetricTensor<2, dim, ScalarType2> &G)
      {
        Assert(
          valid_invariants().count(i) != 0,
          ExcMessage(
            "The selected invariant is not a valid coupled transverse isotropic invariant."));

        switch (i)
          {
            case I1:
              {
                return Isotropic::Ii(i, C);
              }
              break;
            case I2:
              {
                return Isotropic::Ii(i, C);
              }
              break;
            case I3:
              {
                return Isotropic::Ii(i, C);
              }
              break;
            case pI3:
              {
                return Isotropic::Ii(i, C);
              }
              break;
            case I4:
              {
                return Transverse_Isotropic::Ii(i, C, C_inv, G);
              }
              break;
            case I5:
              {
                return Transverse_Isotropic::Ii(i, C, C_inv, G);
              }
              break;
            case pI5:
              {
                return Transverse_Isotropic::Ii(i, C, C_inv, G);
              }
              break;
            case I9:
              {
                return Coupled_Isotropic::Ii(i, C, C_inv, H);
              }
              break;
            case I10:
              {
                return Coupled_Isotropic::Ii(i, C, C_inv, H);
              }
              break;
            case I11:
              {
                return Coupled_Isotropic::Ii(i, C, C_inv, H);
              }
              break;
            case pI11a:
              {
                return Coupled_Isotropic::Ii(i, C, C_inv, H);
              }
              break;
            case pI11b:
              {
                return Coupled_Isotropic::Ii(i, C, C_inv, H);
              }
              break;
            case I12:
              {
                // I12 = [HxH]:G
                return contract3(H, G, H);
              }
              break;
            case I13:
              {
                // I13 = [HxH]:[C.G.C]
                const Tensor<1, dim, ScalarType> Z = C * H;
                return contract3(Z, G, Z);
              }
              break;
            default:
              break;
          }

        AssertThrow(false,
                    ExcMessage("This (pseudo-)invariant is not defined."));
        return ScalarType();
      }


      /**
       * Returns the selected invariant/pseudo-invariant first
       * derivative with respect to the rank-2 symmetric field tensor
       * $\mathbf{C}$.
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] C_inv The inverse of the rank-2 symmetric field tensor
       * $\mathbf{C}$. If it is known that the required invariant derivative
       * does not require this value (that is typically expensive to compute),
       * then any other value (such as the identity tensor) can be used instead.
       * @param[in] H The rank-1 field tensor $\mathrm{H}$
       * @param[in] G The rank-2 symmetric structure tensor $\mathbf{G}$
       */
      template <int dim, typename ScalarType, typename ScalarType2>
      SymmetricTensor<2, dim, ScalarType>
      dIi_dC(const enum InvariantList &                  i,
             const SymmetricTensor<2, dim, ScalarType> & C,
             const SymmetricTensor<2, dim, ScalarType> & C_inv,
             const Tensor<1, dim, ScalarType> &          H,
             const SymmetricTensor<2, dim, ScalarType2> &G)
      {
        Assert(
          valid_invariants().count(i) != 0,
          ExcMessage(
            "The selected invariant is not a valid coupled transverse isotropic invariant."));

        switch (i)
          {
            case I1:
              {
                return Isotropic::dIi_dC(i, C, C_inv);
              }
              break;
            case I2:
              {
                return Isotropic::dIi_dC(i, C, C_inv);
              }
              break;
            case I3:
              {
                return Isotropic::dIi_dC(i, C, C_inv);
              }
              break;
            case pI3:
              {
                return Isotropic::dIi_dC(i, C, C_inv);
              }
              break;
            case I4:
              {
                return Transverse_Isotropic::dIi_dC(i, C, C_inv, G);
              }
              break;
            case I5:
              {
                return Transverse_Isotropic::dIi_dC(i, C, C_inv, G);
              }
              break;
            case pI5:
              {
                return Transverse_Isotropic::dIi_dC(i, C, C_inv, G);
              }
              break;
            case I9:
              {
                return Coupled_Isotropic::dIi_dC(i, C, C_inv, H);
              }
              break;
            case I10:
              {
                return Coupled_Isotropic::dIi_dC(i, C, C_inv, H);
              }
              break;
            case I11:
              {
                return Coupled_Isotropic::dIi_dC(i, C, C_inv, H);
              }
              break;
            case pI11a:
              {
                return Coupled_Isotropic::dIi_dC(i, C, C_inv, H);
              }
              break;
            case pI11b:
              {
                return Coupled_Isotropic::dIi_dC(i, C, C_inv, H);
              }
              break;
            case I12:
              {
                // I12 = [HxH]:G
                return SymmetricTensor<2, dim, ScalarType>();
              }
              break;
            case I13:
              {
                // I13 = [HxH]:[C.G.C]
                const Tensor<1, dim, ScalarType> Z  = C * H;
                const Tensor<1, dim, ScalarType> GZ = G * Z;
                return 2.0 * symmetrize(outer_product(H, GZ));
              }
              break;
            default:
              break;
          }

        AssertThrow(false,
                    ExcMessage("This (pseudo-)invariant is not defined."));
        return SymmetricTensor<2, dim, ScalarType>();
      }


      /**
       * Returns the selected invariant/pseudo-invariant first
       * derivative with respect to the rank-1 field tensor
       * $\mathrm{H}$.
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] C_inv The inverse of the rank-2 symmetric field tensor
       * $\mathbf{C}$. If it is known that the required invariant derivative
       * does not require this value (that is typically expensive to compute),
       * then any other value (such as the identity tensor) can be used instead.
       * @param[in] H The rank-1 field tensor $\mathrm{H}$
       * @param[in] G The rank-2 symmetric structure tensor $\mathbf{G}$
       */
      template <int dim, typename ScalarType, typename ScalarType2>
      Tensor<1, dim, ScalarType>
      dIi_dH(const enum InvariantList &                  i,
             const SymmetricTensor<2, dim, ScalarType> & C,
             const SymmetricTensor<2, dim, ScalarType> & C_inv,
             const Tensor<1, dim, ScalarType> &          H,
             const SymmetricTensor<2, dim, ScalarType2> &G)
      {
        Assert(
          valid_invariants().count(i) != 0,
          ExcMessage(
            "The selected invariant is not a valid coupled transverse isotropic invariant."));

        switch (i)
          {
            case I1:
              {
                return Tensor<1, dim, ScalarType>();
              }
              break;
            case I2:
              {
                return Tensor<1, dim, ScalarType>();
              }
              break;
            case I3:
              {
                return Tensor<1, dim, ScalarType>();
              }
              break;
            case pI3:
              {
                return Tensor<1, dim, ScalarType>();
              }
              break;
            case I4:
              {
                return Tensor<1, dim, ScalarType>();
              }
              break;
            case I5:
              {
                return Tensor<1, dim, ScalarType>();
              }
              break;
            case pI5:
              {
                return Tensor<1, dim, ScalarType>();
              }
              break;
            case I9:
              {
                return Coupled_Isotropic::dIi_dH(i, C, C_inv, H);
              }
              break;
            case I10:
              {
                return Coupled_Isotropic::dIi_dH(i, C, C_inv, H);
              }
              break;
            case I11:
              {
                return Coupled_Isotropic::dIi_dH(i, C, C_inv, H);
              }
              break;
            case pI11a:
              {
                return Coupled_Isotropic::dIi_dH(i, C, C_inv, H);
              }
              break;
            case pI11b:
              {
                return Coupled_Isotropic::dIi_dH(i, C, C_inv, H);
              }
              break;
            case I12:
              {
                // I12 = [HxH]:G
                return 2.0 * (G * H);
              }
              break;
            case I13:
              {
                // I13 = [HxH]:[C.G.C]
                const Tensor<1, dim, ScalarType> Z = C * H;
                return 2.0 * ((Z * G) * C);
              }
              break;
            default:
              break;
          }

        AssertThrow(false,
                    ExcMessage("This (pseudo-)invariant is not defined."));
        return Tensor<1, dim, ScalarType>();
      }


      /**
       * Returns the selected invariant/pseudo-invariant second
       * derivative with respect to the rank-2 symmetric field tensor
       * $\mathbf{C}$.
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] C_inv The inverse of the rank-2 symmetric field tensor
       * $\mathbf{C}$. If it is known that the required invariant derivative
       * does not require this value (that is typically expensive to compute),
       * then any other value (such as the identity tensor) can be used instead.
       * @param[in] H The rank-1 field tensor $\mathrm{H}$
       * @param[in] G The rank-2 symmetric structure tensor $\mathbf{G}$
       */
      template <int dim, typename ScalarType, typename ScalarType2>
      SymmetricTensor<4, dim, ScalarType>
      d2Ii_dC_dC(const enum InvariantList &                  i,
                 const SymmetricTensor<2, dim, ScalarType> & C,
                 const SymmetricTensor<2, dim, ScalarType> & C_inv,
                 const Tensor<1, dim, ScalarType> &          H,
                 const SymmetricTensor<2, dim, ScalarType2> &G)
      {
        Assert(
          valid_invariants().count(i) != 0,
          ExcMessage(
            "The selected invariant is not a valid coupled transverse isotropic invariant."));

        switch (i)
          {
            case I1:
              {
                return Isotropic::d2Ii_dC_dC(i, C, C_inv);
              }
              break;
            case I2:
              {
                return Isotropic::d2Ii_dC_dC(i, C, C_inv);
              }
              break;
            case I3:
              {
                return Isotropic::d2Ii_dC_dC(i, C, C_inv);
              }
              break;
            case pI3:
              {
                return Isotropic::d2Ii_dC_dC(i, C, C_inv);
              }
              break;
            case I4:
              {
                return Transverse_Isotropic::d2Ii_dC_dC(i, C, C_inv, G);
              }
              break;
            case I5:
              {
                return Transverse_Isotropic::d2Ii_dC_dC(i, C, C_inv, G);
              }
              break;
            case pI5:
              {
                return Transverse_Isotropic::d2Ii_dC_dC(i, C, C_inv, G);
              }
              break;
            case I9:
              {
                return Coupled_Isotropic::d2Ii_dC_dC(i, C, C_inv, H);
              }
              break;
            case I10:
              {
                return Coupled_Isotropic::d2Ii_dC_dC(i, C, C_inv, H);
              }
              break;
            case I11:
              {
                return Coupled_Isotropic::d2Ii_dC_dC(i, C, C_inv, H);
              }
              break;
            case pI11a:
              {
                return Coupled_Isotropic::d2Ii_dC_dC(i, C, C_inv, H);
              }
              break;
            case pI11b:
              {
                return Coupled_Isotropic::d2Ii_dC_dC(i, C, C_inv, H);
              }
              break;
            case I12:
              {
                // I12 = [HxH]:G
                return SymmetricTensor<4, dim, ScalarType>();
              }
              break;
            case I13:
              {
                // I13 = [HxH]:[C.G.C]
                SymmetricTensor<4, dim, ScalarType> d2I13_dC_dC;
                for (unsigned int A = 0; A < dim; ++A)
                  for (unsigned int B = A; B < dim; ++B)
                    for (unsigned int C = 0; C < dim; ++C)
                      for (unsigned int D = C; D < dim; ++D)
                        {
                          d2I13_dC_dC[A][B][C][D] =
                            0.5 *
                            (H[A] * G[B][C] * H[D] + H[A] * G[B][D] * H[C]);
                          d2I13_dC_dC[A][B][C][D] += G[A][B] * H[C] * H[D];
                        }

                return d2I13_dC_dC;
              }
              break;
            default:
              break;
          }

        AssertThrow(false,
                    ExcMessage("This (pseudo-)invariant is not defined."));
        return SymmetricTensor<4, dim, ScalarType>();
      }


      /**
       * Returns the selected invariant/pseudo-invariant second
       * derivative with respect to the rank-1 field tensor
       * $\mathrm{H}$.
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] C_inv The inverse of the rank-2 symmetric field tensor
       * $\mathbf{C}$. If it is known that the required invariant derivative
       * does not require this value (that is typically expensive to compute),
       * then any other value (such as the identity tensor) can be used instead.
       * @param[in] H The rank-1 field tensor $\mathrm{H}$
       * @param[in] G The rank-2 symmetric structure tensor $\mathbf{G}$
       */
      template <int dim, typename ScalarType, typename ScalarType2>
      SymmetricTensor<2, dim, ScalarType>
      d2Ii_dH_dH(const enum InvariantList &                  i,
                 const SymmetricTensor<2, dim, ScalarType> & C,
                 const SymmetricTensor<2, dim, ScalarType> & C_inv,
                 const Tensor<1, dim, ScalarType> &          H,
                 const SymmetricTensor<2, dim, ScalarType2> &G)
      {
        Assert(
          valid_invariants().count(i) != 0,
          ExcMessage(
            "The selected invariant is not a valid coupled transverse isotropic invariant."));

        switch (i)
          {
            case I1:
              {
                return SymmetricTensor<2, dim, ScalarType>();
              }
              break;
            case I2:
              {
                return SymmetricTensor<2, dim, ScalarType>();
              }
              break;
            case I3:
              {
                return SymmetricTensor<2, dim, ScalarType>();
              }
              break;
            case pI3:
              {
                return SymmetricTensor<2, dim, ScalarType>();
              }
              break;
            case I4:
              {
                return SymmetricTensor<2, dim, ScalarType>();
              }
              break;
            case I5:
              {
                return SymmetricTensor<2, dim, ScalarType>();
              }
              break;
            case pI5:
              {
                return SymmetricTensor<2, dim, ScalarType>();
              }
              break;
            case I9:
              {
                return Coupled_Isotropic::d2Ii_dH_dH(i, C, C_inv, H);
              }
              break;
            case I10:
              {
                return Coupled_Isotropic::d2Ii_dH_dH(i, C, C_inv, H);
              }
              break;
            case I11:
              {
                return Coupled_Isotropic::d2Ii_dH_dH(i, C, C_inv, H);
              }
              break;
            case pI11a:
              {
                return Coupled_Isotropic::d2Ii_dH_dH(i, C, C_inv, H);
              }
              break;
            case pI11b:
              {
                return Coupled_Isotropic::d2Ii_dH_dH(i, C, C_inv, H);
              }
              break;
            case I12:
              {
                // I12 = [HxH]:G
                return 2.0 * G;
              }
              break;
            case I13:
              {
                // I13 = [HxH]:[C.G.C]
                const Tensor<2, dim, ScalarType> G_dot_C =
                  internal::compute_ST_dot_ST(G, C);
                return 2.0 * symmetrize(C * G_dot_C);
              }
              break;
            default:
              break;
          }

        AssertThrow(false,
                    ExcMessage("This (pseudo-)invariant is not defined."));
        return SymmetricTensor<2, dim, ScalarType>();
      }


      /**
       * Returns the selected invariant/pseudo-invariant second
       * derivative, with the derivatives first taken with respect to the rank-2
       * symmetric field tensor
       * $\mathbf{C}$ followed by the rank-1 field tensor $\mathrm{H}$.
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] C_inv The inverse of the rank-2 symmetric field tensor
       * $\mathbf{C}$. If it is known that the required invariant derivative
       * does not require this value (that is typically expensive to compute),
       * then any other value (such as the identity tensor) can be used instead.
       * @param[in] H The rank-1 field tensor $\mathrm{H}$
       * @param[in] G The rank-2 symmetric structure tensor $\mathbf{G}$
       */
      template <int dim, typename ScalarType, typename ScalarType2>
      Tensor<3, dim, ScalarType>
      d2Ii_dC_dH(const enum InvariantList &                  i,
                 const SymmetricTensor<2, dim, ScalarType> & C,
                 const SymmetricTensor<2, dim, ScalarType> & C_inv,
                 const Tensor<1, dim, ScalarType> &          H,
                 const SymmetricTensor<2, dim, ScalarType2> &G)
      {
        Assert(
          valid_invariants().count(i) != 0,
          ExcMessage(
            "The selected invariant is not a valid coupled transverse isotropic invariant."));

        switch (i)
          {
            case I1:
              {
                return Tensor<3, dim, ScalarType>();
              }
              break;
            case I2:
              {
                return Tensor<3, dim, ScalarType>();
              }
              break;
            case I3:
              {
                return Tensor<3, dim, ScalarType>();
              }
              break;
            case pI3:
              {
                return Tensor<3, dim, ScalarType>();
              }
              break;
            case I4:
              {
                return Tensor<3, dim, ScalarType>();
              }
              break;
            case I5:
              {
                return Tensor<3, dim, ScalarType>();
              }
              break;
            case pI5:
              {
                return Tensor<3, dim, ScalarType>();
              }
              break;
            case I9:
              {
                return Coupled_Isotropic::d2Ii_dC_dH(i, C, C_inv, H);
              }
              break;
            case I10:
              {
                return Coupled_Isotropic::d2Ii_dC_dH(i, C, C_inv, H);
              }
              break;
            case I11:
              {
                return Coupled_Isotropic::d2Ii_dC_dH(i, C, C_inv, H);
              }
              break;
            case pI11a:
              {
                return Coupled_Isotropic::d2Ii_dC_dH(i, C, C_inv, H);
              }
              break;
            case pI11b:
              {
                return Coupled_Isotropic::d2Ii_dC_dH(i, C, C_inv, H);
              }
              break;
            case I12:
              {
                // I12 = [HxH]:G
                return Tensor<3, dim, ScalarType>();
              }
              break;
            case I13:
              {
                // I13 = [HxH]:[C.G.C]
                const Tensor<1, dim, ScalarType> Z  = C * H;
                const Tensor<1, dim, ScalarType> GZ = G * Z;
                const Tensor<2, dim, ScalarType> GC =
                  internal::compute_ST_dot_ST(G, C);
                const SymmetricTensor<2, dim, double> &I =
                  Physics::Elasticity::StandardTensors<dim>::I;

                Tensor<3, dim, ScalarType> d2I13_dC_dH;
                for (unsigned int A = 0; A < dim; ++A)
                  for (unsigned int B = 0; B < dim; ++B)
                    for (unsigned int D = 0; D < dim; ++D)
                      d2I13_dC_dH[A][B][D] = I[A][D] * GZ[B] + GZ[A] * I[B][D] +
                                             GC[A][D] * H[B] + H[A] * GC[B][D];

                // return 2.0 * symmetrize(outer_product(H,GZ));
                return d2I13_dC_dH;
              }
              break;
            default:
              break;
          }

        AssertThrow(false,
                    ExcMessage("This (pseudo-)invariant is not defined."));
        return Tensor<3, dim, ScalarType>();
      }


      /**
       * Returns the selected invariant/pseudo-invariant second
       * derivative, with the derivatives first taken with respect to the rank-1
       * field tensor $\mathrm{H}$ followed by the rank-2 symmetric field tensor
       * $\mathbf{C}$ .
       *
       * @param[in] i The invariant value to be computed
       * @param[in] C The rank-2 symmetric field tensor $\mathbf{C}$
       * @param[in] C_inv The inverse of the rank-2 symmetric field tensor
       * $\mathbf{C}$. If it is known that the required invariant derivative
       * does not require this value (that is typically expensive to compute),
       * then any other value (such as the identity tensor) can be used instead.
       * @param[in] H The rank-1 field tensor $\mathrm{H}$
       * @param[in] G The rank-2 symmetric structure tensor $\mathbf{G}$
       */
      template <int dim, typename ScalarType, typename ScalarType2>
      Tensor<3, dim, ScalarType>
      d2Ii_dH_dC(const enum InvariantList &                  i,
                 const SymmetricTensor<2, dim, ScalarType> & C,
                 const SymmetricTensor<2, dim, ScalarType> & C_inv,
                 const Tensor<1, dim, ScalarType> &          H,
                 const SymmetricTensor<2, dim, ScalarType2> &G)
      {
        Assert(
          valid_invariants().count(i) != 0,
          ExcMessage(
            "The selected invariant is not a valid coupled transverse isotropic invariant."));

        const Tensor<3, dim, ScalarType> d2Ii_dC_dH =
          Coupled_Transverse_Isotropic::d2Ii_dC_dH(i, C, C_inv, H, G);

        Tensor<3, dim, ScalarType> d2Ii_dH_dC;
        for (unsigned int A = 0; A < dim; ++A)
          for (unsigned int B = 0; B < dim; ++B)
            for (unsigned int C = 0; C < dim; ++C)
              d2Ii_dH_dC[A][B][C] = d2Ii_dC_dH[C][B][A];

        return d2Ii_dH_dC;
      }
    } // namespace Coupled_Transverse_Isotropic

  } // namespace Invariants
} // namespace Physics


DEAL_II_NAMESPACE_CLOSE


#endif // dealii_physics_invariants_h

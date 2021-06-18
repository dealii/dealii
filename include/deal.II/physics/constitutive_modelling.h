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

#ifndef dealii_physics_constitutive_modelling_h
#define dealii_physics_constitutive_modelling_h

#include <deal.II/base/config.h>

// #include <deal.II/base/exceptions.h>
// #include <deal.II/base/numbers.h>
#include <deal.II/base/std_cxx20/iota_view.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/physics/elasticity/standard_tensors.h>

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
   * Ogden
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
   */
  namespace ConstitutiveModelling
  {



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
    } // namespace internal



    namespace Isotropic
    {
      SymmetricTensor<2, dim, ScalarType>
      get_dPsi_dC(const Coupled_Material_Values<dim, ScalarType> &values) const
      {
        SymmetricTensor<2, dim, ScalarType> dPsi_dC;

        for (typename std::vector<int>::const_iterator it =
               Coupled_Material_Invariants::invariants.begin();
             it != Coupled_Material_Invariants::invariants.end();
             ++it)
          {
            const int &      i        = *it;
            const ScalarType dPsi_dIi = get_dPsi_dIi(i, values);

            // Since computing tensor derivatives is really expensive, we
            // only do so if the scalar coefficient is non-zero. But we
            // can only do this optimisation for floating point types!
            // We need to ensure that we track ALL sensitivities for
            // AD and SD types
            if (add_invariant_contribution(dPsi_dIi) == true)
              dPsi_dC += dPsi_dIi * Coupled_Material_Invariants::dIi_dC<dim>(
                                      i, values.H, values.C, values.C_inv);
          }

        return dPsi_dC;
      }

      SymmetricTensor<4, dim, ScalarType>
      get_d2Psi_dC_dC(
        const Coupled_Material_Values<dim, ScalarType> &values) const
      {
        SymmetricTensor<4, dim, ScalarType> d2Psi_dC_dC;

        for (typename std::vector<int>::const_iterator it1 =
               Coupled_Material_Invariants::invariants.begin();
             it1 != Coupled_Material_Invariants::invariants.end();
             ++it1)
          {
            const int &      i        = *it1;
            const ScalarType dPsi_dIi = get_dPsi_dIi(i, values);

            if (add_invariant_contribution(dPsi_dIi) == true)
              d2Psi_dC_dC +=
                dPsi_dIi * Coupled_Material_Invariants::d2Ii_dC_dC<dim>(
                             i, values.H, values.C, values.C_inv);

            for (typename std::vector<int>::const_iterator it2 =
                   Coupled_Material_Invariants::invariants.begin();
                 it2 != Coupled_Material_Invariants::invariants.end();
                 ++it2)
              {
                const int &      j = *it2;
                const ScalarType d2Psi_dIi_dIj =
                  get_d2Psi_dIi_dIj(i, j, values);

                if (add_invariant_contribution(d2Psi_dIi_dIj) == true)
                  d2Psi_dC_dC +=
                    d2Psi_dIi_dIj *
                    outer_product(Coupled_Material_Invariants::dIi_dC<dim>(
                                    i, values.H, values.C, values.C_inv),
                                  Coupled_Material_Invariants::dIi_dC<dim>(
                                    j, values.H, values.C, values.C_inv));
              }
          }

        return d2Psi_dC_dC;
      }

    } // namespace Isotropic


    //   // virtual SymmetricTensor<2, dim, ScalarType>
    //   // get_dPsi_dC(const Coupled_Material_Values<dim, ScalarType> &values)
    //   // const
    //   // {
    //   //   SymmetricTensor<2, dim, ScalarType> dPsi_dC;

    //   //   for (typename std::vector<int>::const_iterator it =
    //   //   Coupled_Material_Invariants::invariants.begin();
    //   //        it != Coupled_Material_Invariants::invariants.end();
    //   //        ++it)
    //   //     {
    //   //       const int &      i        = *it;
    //   //       const ScalarType dPsi_dIi = get_dPsi_dIi(i, values);

    //   //       // Since computing tensor derivatives is really expensive, we
    //   //       // only do so if the scalar coefficient is non-zero. But we
    //   //       // can only do this optimisation for floating point types!
    //   //       // We need to ensure that we track ALL sensitivities for
    //   //       // AD and SD types
    //   //       if (add_invariant_contribution(dPsi_dIi) == true)
    //   //         dPsi_dC += dPsi_dIi *
    //   //         Coupled_Material_Invariants::dIi_dC<dim>(i, values.H,
    //   values.C,
    //   //         values.C_inv);
    //   //     }

    //   //   return dPsi_dC;
    //   // }

    //   // /**
    //   //  * @copydoc Coupled_ME_Constitutive_Law_Base::get_d2Psi_dC_dC()
    //   //  *
    //   //  * @note In the current implementation, it is assumed that there are no cross-terms
    //   //  * involving invariants (e.g. Psi += I1*I4*I7) and that the energy
    //   is
    //   //  linear
    //   //  * in terms of the invariants.
    //   //  */
    //   // virtual SymmetricTensor<4, dim, ScalarType>
    //   // get_d2Psi_dC_dC(const Coupled_Material_Values<dim, ScalarType>
    //   &values)
    //   // const
    //   // {
    //   //   SymmetricTensor<4, dim, ScalarType> d2Psi_dC_dC;

    //   //   for (typename std::vector<int>::const_iterator it1 =
    //   //   Coupled_Material_Invariants::invariants.begin();
    //   //        it1 != Coupled_Material_Invariants::invariants.end();
    //   //        ++it1)
    //   //     {
    //   //       const int &      i        = *it1;
    //   //       const ScalarType dPsi_dIi = get_dPsi_dIi(i, values);

    //   //       if (add_invariant_contribution(dPsi_dIi) == true)
    //   //         d2Psi_dC_dC += dPsi_dIi *
    //   //         Coupled_Material_Invariants::d2Ii_dC_dC<dim>(i, values.H,
    //   //         values.C, values.C_inv);

    //   //       for (typename std::vector<int>::const_iterator it2 =
    //   //       Coupled_Material_Invariants::invariants.begin();
    //   //            it2 != Coupled_Material_Invariants::invariants.end();
    //   //            ++it2)
    //   //         {
    //   //           const int &      j             = *it2;
    //   //           const ScalarType d2Psi_dIi_dIj = get_d2Psi_dIi_dIj(i, j,
    //   //           values);

    //   //           if (add_invariant_contribution(d2Psi_dIi_dIj) == true)
    //   //             d2Psi_dC_dC += d2Psi_dIi_dIj *
    //   // outer_product(Coupled_Material_Invariants::dIi_dC<dim>(i,
    //   //             values.H, values.C, values.C_inv),
    //   // Coupled_Material_Invariants::dIi_dC<dim>(j,
    //   //                                                          values.H,
    //   //                                                          values.C,
    //   // values.C_inv));
    //   //         }
    //   //     }

    //   //   return d2Psi_dC_dC;
    //   // }

    //   // /**
    //   //  * @copydoc Coupled_ME_Constitutive_Law_Base::get_dPsi_dH()
    //   //  */
    //   // virtual Tensor<1, dim, ScalarType>
    //   // get_dPsi_dH(const Coupled_Material_Values<dim, ScalarType> &values)
    //   // const
    //   // {
    //   //   Tensor<1, dim, ScalarType> dPsi_dH;

    //   //   for (typename std::vector<int>::const_iterator it =
    //   //   Coupled_Material_Invariants::invariants.begin();
    //   //        it != Coupled_Material_Invariants::invariants.end();
    //   //        ++it)
    //   //     {
    //   //       const int &      i        = *it;
    //   //       const ScalarType dPsi_dIi = get_dPsi_dIi(i, values);

    //   //       if (add_invariant_contribution(dPsi_dIi) == true)
    //   //         dPsi_dH += dPsi_dIi *
    //   //         Coupled_Material_Invariants::dIi_dH<dim>(i, values.H,
    //   values.C,
    //   //         values.C_inv);
    //   //     }

    //   //   return dPsi_dH;
    //   // }

    //   // /**
    //   //  * @copydoc Coupled_ME_Constitutive_Law_Base::get_d2Psi_dH_dH()
    //   //  */
    //   // virtual SymmetricTensor<2, dim, ScalarType>
    //   // get_d2Psi_dH_dH(const Coupled_Material_Values<dim, ScalarType>
    //   &values)
    //   // const
    //   // {
    //   //   SymmetricTensor<2, dim, ScalarType> d2Psi_dH_dH;

    //   //   for (typename std::vector<int>::const_iterator it1 =
    //   //   Coupled_Material_Invariants::invariants.begin();
    //   //        it1 != Coupled_Material_Invariants::invariants.end();
    //   //        ++it1)
    //   //     {
    //   //       const int &      i        = *it1;
    //   //       const ScalarType dPsi_dIi = get_dPsi_dIi(i, values);

    //   //       if (add_invariant_contribution(dPsi_dIi) == true)
    //   //         d2Psi_dH_dH += dPsi_dIi *
    //   //         Coupled_Material_Invariants::d2Ii_dH_dH<dim>(i, values.H,
    //   //         values.C, values.C_inv);

    //   //       for (typename std::vector<int>::const_iterator it2 =
    //   //       Coupled_Material_Invariants::invariants.begin();
    //   //            it2 != Coupled_Material_Invariants::invariants.end();
    //   //            ++it2)
    //   //         {
    //   //           const int &      j             = *it2;
    //   //           const ScalarType d2Psi_dIi_dIj = get_d2Psi_dIi_dIj(i, j,
    //   //           values);

    //   //           if (add_invariant_contribution(d2Psi_dIi_dIj) == true)
    //   //             d2Psi_dH_dH += d2Psi_dIi_dIj *
    //   // symmetrize(outer_product(Coupled_Material_Invariants::dIi_dH<dim>(i,
    //   //             values.H, values.C, values.C_inv),
    //   // Coupled_Material_Invariants::dIi_dH<dim>(j,
    //   // values.H,
    //   // values.C,
    //   // values.C_inv)));
    //   //         }
    //   //     }

    //   //   return d2Psi_dH_dH;
    //   // }

    //   // /**
    //   //  * @copydoc Coupled_ME_Constitutive_Law_Base::get_d2Psi_dH_dC()
    //   //  */
    //   // virtual Tensor<3, dim, ScalarType>
    //   // get_d2Psi_dH_dC(const Coupled_Material_Values<dim, ScalarType>
    //   &values)
    //   // const
    //   // {
    //   //   Tensor<3, dim, ScalarType> d2Psi_dH_dC;

    //   //   for (typename std::vector<int>::const_iterator it1 =
    //   //   Coupled_Material_Invariants::invariants.begin();
    //   //        it1 != Coupled_Material_Invariants::invariants.end();
    //   //        ++it1)
    //   //     {
    //   //       const int &      i        = *it1;
    //   //       const ScalarType dPsi_dIi = get_dPsi_dIi(i, values);

    //   //       if (add_invariant_contribution(dPsi_dIi) == true)
    //   //         d2Psi_dH_dC += dPsi_dIi *
    //   //         Coupled_Material_Invariants::d2Ii_dH_dC<dim>(i, values.H,
    //   //         values.C, values.C_inv);

    //   //       for (typename std::vector<int>::const_iterator it2 =
    //   //       Coupled_Material_Invariants::invariants.begin();
    //   //            it2 != Coupled_Material_Invariants::invariants.end();
    //   //            ++it2)
    //   //         {
    //   //           const int &      j             = *it2;
    //   //           const ScalarType d2Psi_dIi_dIj = get_d2Psi_dIi_dIj(i, j,
    //   //           values);

    //   //           if (add_invariant_contribution(d2Psi_dIi_dIj) == true)
    //   //             d2Psi_dH_dC +=
    //   //               d2Psi_dIi_dIj *
    //   //               outer_product(static_cast<Tensor<2, dim,
    //   // ScalarType>>(Coupled_Material_Invariants::dIi_dC<dim>(i,
    //   //               values.H, values.C, values.C_inv)),
    //   // Coupled_Material_Invariants::dIi_dH<dim>(j,
    //   //                             values.H, values.C, values.C_inv));
    //   //         }
    //   //     }

    //   //   return d2Psi_dH_dC;
    //   // }

    // } // namespace Coupled_Isotropic



    // namespace Coupled_Transverse_Isotropic
    // {} // namespace Coupled_Transverse_Isotropic

  } // namespace ConstitutiveModelling
} // namespace Physics


DEAL_II_NAMESPACE_CLOSE


#endif // dealii_physics_constitutive_modelling_h

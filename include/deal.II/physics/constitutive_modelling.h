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



    namespace Transverse_Isotropic
    {} // namespace Transverse_Isotropic



    // namespace Coupled_Isotropic
    // {



    //   /**
    //    * Get an object that returns an iterator to all
    //    * invariants.
    //    *
    //    * @return std_cxx20::ranges::iota_view
    //    */
    //   std_cxx20::ranges::iota_view
    //   get_invariant_range()
    //   {
    //     return {I1, pI7b + 1};
    //   }



    //   /**
    //    * Returns the selected invariant/pseudo-invariant value
    //    *
    //    * @param[in] i The invariant value to be computed
    //    * @param[in] H The magnetic field vector
    //    * @param[in] C The right Cauchy-Green deformation tensor
    //    * @param[in] C_inv The inverse of the right Cauchy-Green
    //    * deformation tensor
    //    */
    //   template <int dim, typename ScalarType>
    //   ScalarType
    //   Ii(const enum Invariants &                    i,
    //      const Tensor<1, dim, ScalarType> &         H,
    //      const SymmetricTensor<2, dim, ScalarType> &C,
    //      const SymmetricTensor<2, dim, ScalarType> &C_inv)
    //   {
    //     switch (i)
    //       {
    //         case I1:
    //           {
    //             // I1 = tr(C)
    //             return Isotropic::Ii(Isotropic::I1, C);
    //           }
    //           break;
    //         case I2:
    //           {
    //             // I2 = 0.5*(tr(C)*tr(C) - tr(C*C))
    //             return Isotropic::Ii(Isotropic::I2, C);
    //           }
    //           break;
    //         case I3:
    //           {
    //             // I3 = det(C)
    //             return Isotropic::Ii(Isotropic::I3, C);
    //           }
    //           break;
    //         case pI3a:
    //           {
    //             // pI3a = sqrt(det(C))
    //             return Isotropic::Ii(Isotropic::pI3a, C);
    //           }
    //           break;
    //         case I4:
    //           {
    //             // I4 = [HxH]:I
    //             return H * H;
    //           }
    //           break;
    //         case I5:
    //           {
    //             // I5 = [HxH]:C
    //             //          return H*C*H; // Breaks CSE!
    //             return contract3(H, C, H);
    //           }
    //           break;
    //         case I6:
    //           {
    //             // I6 = [HxH]:C^2 = [CxH].[CxH]
    //             //          return H*(C*C)*H; // Breaks CSE!
    //             const Tensor<1, dim, ScalarType> Z = C * H;
    //             return Z * Z;
    //           }
    //           break;
    //         case pI7a:
    //           {
    //             // I7 = [HxH]:C^{-1}
    //             //          return H*C_inv*H; // Breaks CSE!
    //             return contract3(H, C_inv, H);
    //           }
    //           break;
    //         default:
    //           break;
    //       }

    //     AssertThrow(false,
    //                 ExcMessage("This (pseudo-)invariant is not defined."));
    //     return ScalarType();
    //   }

    //   /**
    //    * @name Invariant first derivatives
    //    */

    //   /**
    //    * Returns the selected invariant/pseudo-invariant first
    //    * derivative with respect to the
    //    * right Cauchy-Green deformation tensor
    //    *
    //    * @param[in] i The invariant value to be computed
    //    * @param[in] H The magnetic field vector
    //    * @param[in] C The right Cauchy-Green deformation tensor
    //    * @param[in] C_inv The inverse of the right Cauchy-Green
    //    * deformation tensor
    //    */
    //   template <int dim, typename ScalarType>
    //   SymmetricTensor<2, dim, ScalarType>
    //   dIi_dC(const enum Invariants &                    i,
    //          const Tensor<1, dim, ScalarType> &         H,
    //          const SymmetricTensor<2, dim, ScalarType> &C,
    //          const SymmetricTensor<2, dim, ScalarType> &C_inv)
    //   {
    //     switch (i)
    //       {
    //         case I1:
    //           {
    //             // I1 = tr(C)
    //             return static_cast<SymmetricTensor<2, dim, ScalarType>>(
    //               Physics::Elasticity::StandardTensors<dim>::I);
    //           }
    //           break;
    //         case I2:
    //           {
    //             // I2 = 0.5*(tr(C)*tr(C) - tr(C*C))
    //             const ScalarType I1 = Ii(Invariants::I1, H, C, C_inv);
    //             return I1 * Physics::Elasticity::StandardTensors<dim>::I - C;
    //           }
    //           break;
    //         case I3:
    //           {
    //             // I3 = det(C)
    //             const ScalarType I3 = Ii(Invariants::I3, H, C, C_inv);
    //             return I3 * C_inv;
    //           }
    //           break;
    //         case pI3a:
    //           {
    //             // pI3a = sqrt(det(C)) = sqrt(I3)
    //             const ScalarType pI3a = Ii(Invariants::pI3a, H, C, C_inv);
    //             return 0.5 * pI3a * C_inv;
    //           }
    //           break;
    //         case I4:
    //           {
    //             // I4= [HxH]:I
    //             return SymmetricTensor<2, dim, ScalarType>();
    //           }
    //           break;
    //         case I5:
    //           {
    //             // I5 = [HxH]:C
    //             return symmetrize(outer_product(H, H));
    //           }
    //           break;
    //         case I6:
    //           {
    //             // I6 = [HxH]:C^2 = [CxH].[CxH]
    //             const Tensor<1, dim, ScalarType> Z = C * H;
    //             return symmetrize(outer_product(Z, H) + outer_product(H, Z));
    //           }
    //           break;
    //         case pI7a:
    //           {
    //             // I7 = [HxH]:C^{-1}
    //             const Tensor<1, dim, ScalarType> Y = C_inv * H;
    //             return -symmetrize(outer_product(Y, Y));
    //           }
    //           break;
    //         default:
    //           break;
    //       }

    //     AssertThrow(false,
    //                 ExcMessage("This (pseudo-)invariant is not defined."));
    //     return SymmetricTensor<2, dim, ScalarType>();
    //   }

    //   /**
    //    * Returns the selected invariant/pseudo-invariant first
    //    * derivative with respect to the
    //    * magnetic field vector
    //    *
    //    * @param[in] i The invariant value to be computed
    //    * @param[in] H The magnetic field vector
    //    * @param[in] C The right Cauchy-Green deformation tensor
    //    * @param[in] C_inv The inverse of the right Cauchy-Green
    //    * deformation tensor
    //    */
    //   template <int dim, typename ScalarType>
    //   Tensor<1, dim, ScalarType>
    //   dIi_dH(const enum Invariants &                    i,
    //          const Tensor<1, dim, ScalarType> &         H,
    //          const SymmetricTensor<2, dim, ScalarType> &C,
    //          const SymmetricTensor<2, dim, ScalarType> &C_inv)
    //   {
    //     switch (i)
    //       {
    //         case I1:
    //           {
    //             // I1 = tr(C)
    //             return Tensor<1, dim, ScalarType>();
    //           }
    //           break;
    //         case I2:
    //           {
    //             // I2 = 0.5*(tr(C)*tr(C) - tr(C*C))
    //             return Tensor<1, dim, ScalarType>();
    //           }
    //           break;
    //         case I3:
    //           {
    //             // I3 = det(C) = J*J
    //             return Tensor<1, dim, ScalarType>();
    //           }
    //           break;
    //         case I4:
    //           {
    //             // I4= [HxH]:I
    //             return 2.0 * H;
    //           }
    //           break;
    //         case I5:
    //           {
    //             // I5 = [HxH]:C
    //             return 2.0 * (C * H);
    //           }
    //           break;
    //         case I6:
    //           {
    //             // I6 = [HxH]:C^2 = [CxH].[CxH]
    //             return 2.0 * (C * (C * H));
    //           }
    //           break;
    //         case pI7a:
    //           {
    //             // I7 = [HxH]:C^{-1}
    //             return 2.0 * (C_inv * H);
    //           }
    //           break;
    //         default:
    //           break;
    //       }

    //     AssertThrow(false,
    //                 ExcMessage("This (pseudo-)invariant is not defined."));
    //     return Tensor<1, dim, ScalarType>();
    //   }

    //   /**
    //    * @name Invariant second derivatives
    //    */

    //   /**
    //    * Returns the selected invariant/pseudo-invariant second
    //    * derivative with respect to the
    //    * right Cauchy-Green deformation tensor
    //    *
    //    * @param[in] i The invariant value to be computed
    //    * @param[in] H The magnetic field vector
    //    * @param[in] C The right Cauchy-Green deformation tensor
    //    * @param[in] C_inv The inverse of the right Cauchy-Green
    //    * deformation tensor
    //    */
    //   template <int dim, typename ScalarType>
    //   SymmetricTensor<4, dim, ScalarType>
    //   d2Ii_dC_dC(const enum Invariants &                    i,
    //              const Tensor<1, dim, ScalarType> &         H,
    //              const SymmetricTensor<2, dim, ScalarType> &C,
    //              const SymmetricTensor<2, dim, ScalarType> &C_inv)
    //   {
    //     switch (i)
    //       {
    //         case I1:
    //           {
    //             // I1 = tr(C)
    //             return SymmetricTensor<4, dim, ScalarType>();
    //           }
    //           break;
    //         case I2:
    //           {
    //             // I2 = 0.5*(tr(C)*tr(C) - tr(C*C))
    //             return static_cast<SymmetricTensor<4, dim, ScalarType>>(
    //               Physics::Elasticity::StandardTensors<dim>::IxI -
    //               Physics::Elasticity::StandardTensors<dim>::S);
    //           }
    //           break;
    //         case I3:
    //           {
    //             // I3 = det(C) = J*J
    //             return determinant(C) *
    //                    (outer_product(C_inv, C_inv) +
    //                     Coupled_Material_Values<dim, ScalarType>::dC_inv_dC(
    //                       C_inv));
    //           }
    //           break;
    //         case I4:
    //           {
    //             // I4= [HxH]:I
    //             return SymmetricTensor<4, dim, ScalarType>();
    //           }
    //           break;
    //         case I5:
    //           {
    //             // I5 = [HxH]:C
    //             return SymmetricTensor<4, dim, ScalarType>();
    //           }
    //           break;
    //         case I6:
    //           {
    //             // I6 = [HxH]:C^2 = [CxH].[CxH]
    //             const SymmetricTensor<2, dim, double> &I =
    //               Physics::Elasticity::StandardTensors<dim>::I;

    //             SymmetricTensor<4, dim, ScalarType> d2I6_dC_dC;
    //             for (unsigned int A = 0; A < dim; ++A)
    //               for (unsigned int B = A; B < dim; ++B)
    //                 for (unsigned int C = 0; C < dim; ++C)
    //                   for (unsigned int D = C; D < dim; ++D)
    //                     {
    //                       // Need to ensure symmetries of (A,B) and (C,D)
    //                       d2I6_dC_dC[A][B][C][D] =
    //                         0.5 *
    //                         (H[A] * H[D] * I[B][C] + H[A] * H[C] * I[B][D] +
    //                          H[B] * H[D] * I[A][C] + H[B] * H[C] * I[A][D]);
    //                     }

    //             return d2I6_dC_dC;
    //           }
    //           break;
    //         case pI7a:
    //           {
    //             // I7 = [HxH]:C^{-1}
    //             const Tensor<1, dim, ScalarType> Y = C_inv * H;

    //             SymmetricTensor<4, dim, ScalarType> d2I7_dC_dC;
    //             for (unsigned int A = 0; A < dim; ++A)
    //               for (unsigned int B = A; B < dim; ++B)
    //                 for (unsigned int C = 0; C < dim; ++C)
    //                   for (unsigned int D = C; D < dim; ++D)
    //                     {
    //                       // Need to ensure symmetries of (A,B) and (C,D)
    //                       d2I7_dC_dC[A][B][C][D] =
    //                         0.5 * (Y[A] * Y[C] * C_inv[B][D] +
    //                                Y[B] * Y[C] * C_inv[A][D] +
    //                                Y[A] * Y[D] * C_inv[B][C] +
    //                                Y[B] * Y[D] * C_inv[A][C]);
    //                     }

    //             return d2I7_dC_dC;
    //           }
    //           break;
    //         default:
    //           break;
    //       }

    //     AssertThrow(false,
    //                 ExcMessage("This (pseudo-)invariant is not defined."));
    //     return SymmetricTensor<4, dim, ScalarType>();
    //   }

    //   /**
    //    * Returns the selected invariant/pseudo-invariant second
    //    * derivative with respect to the
    //    * magnetic field vector
    //    *
    //    * @param[in] i The invariant value to be computed
    //    * @param[in] H The magnetic field vector
    //    * @param[in] C The right Cauchy-Green deformation tensor
    //    * @param[in] C_inv The inverse of the right Cauchy-Green
    //    * deformation tensor
    //    */
    //   template <int dim, typename ScalarType>
    //   SymmetricTensor<2, dim, ScalarType>
    //   d2Ii_dH_dH(const enum Invariants &i,
    //              const Tensor<1, dim, ScalarType> & /*H*/,
    //              const SymmetricTensor<2, dim, ScalarType> &C,
    //              const SymmetricTensor<2, dim, ScalarType> &C_inv)
    //   {
    //     switch (i)
    //       {
    //         case I1:
    //           {
    //             // I1 = tr(C)
    //             return SymmetricTensor<2, dim, ScalarType>();
    //           }
    //           break;
    //         case I2:
    //           {
    //             // I2 = 0.5*(tr(C)*tr(C) - tr(C*C))
    //             return SymmetricTensor<2, dim, ScalarType>();
    //           }
    //           break;
    //         case I3:
    //           {
    //             // I3 = det(C) = J*J
    //             return SymmetricTensor<2, dim, ScalarType>();
    //           }
    //           break;
    //         case I4:
    //           {
    //             // I4= [HxH]:I
    //             return static_cast<SymmetricTensor<2, dim, ScalarType>>(
    //               2.0 * unit_symmetric_tensor<dim>());
    //           }
    //           break;
    //         case I5:
    //           {
    //             // I5 = [HxH]:C
    //             return 2.0 * C;
    //           }
    //           break;
    //         case I6:
    //           {
    //             // I6 = [HxH]:C^2 = [CxH].[CxH]
    //             return 2.0 *
    //                    dealii::Invariants::internal::compute_ST_squared(C);
    //           }
    //           break;
    //         case pI7a:
    //           {
    //             // I7 = [HxH]:C^{-1}
    //             return 2.0 * C_inv;
    //           }
    //           break;
    //         default:
    //           break;
    //       }

    //     AssertThrow(false,
    //                 ExcMessage("This (pseudo-)invariant is not defined."));
    //     return SymmetricTensor<2, dim, ScalarType>();
    //   }

    //   /**
    //    * Returns the selected invariant/pseudo-invariant second
    //    * derivative with respect to the
    //    * right Cauchy-Green deformation tensor and the
    //    * magnetic field vector
    //    *
    //    * @param[in] i The invariant value to be computed
    //    * @param[in] H The magnetic field vector
    //    * @param[in] C The right Cauchy-Green deformation tensor
    //    * @param[in] C_inv The inverse of the right Cauchy-Green
    //    * deformation tensor
    //    */
    //   template <int dim, typename ScalarType>
    //   Tensor<3, dim, ScalarType>
    //   d2Ii_dH_dC(const enum Invariants &                    i,
    //              const Tensor<1, dim, ScalarType> &         H,
    //              const SymmetricTensor<2, dim, ScalarType> &C,
    //              const SymmetricTensor<2, dim, ScalarType> &C_inv)
    //   {
    //     switch (i)
    //       {
    //         case I1:
    //           {
    //             // I1 = tr(C)
    //             return Tensor<3, dim, ScalarType>();
    //           }
    //           break;
    //         case I2:
    //           {
    //             // I2 = 0.5*(tr(C)*tr(C) - tr(C*C))
    //             return Tensor<3, dim, ScalarType>();
    //           }
    //           break;
    //         case I3:
    //           {
    //             // I3 = det(C) = J*J
    //             return Tensor<3, dim, ScalarType>();
    //           }
    //           break;
    //         case I4:
    //           {
    //             // I4= [HxH]:I
    //             return Tensor<3, dim, ScalarType>();
    //           }
    //           break;
    //         case I5:
    //           {
    //             // I5 = [HxH]:C
    //             const SymmetricTensor<2, dim, double> &I =
    //               Physics::Elasticity::StandardTensors<dim>::I;

    //             Tensor<3, dim, ScalarType> d2I5_dH_dC;
    //             for (unsigned int A = 0; A < dim; ++A)
    //               for (unsigned int B = 0; B < dim; ++B)
    //                 for (unsigned int C = 0; C < dim; ++C)
    //                   d2I5_dH_dC[A][B][C] += I[C][A] * H[B] + I[C][B] * H[A];

    //             return d2I5_dH_dC;
    //           }
    //           break;
    //         case I6:
    //           {
    //             // I6 = [HxH]:C^2 = [CxH].[CxH]
    //             const SymmetricTensor<2, dim, double> &I =
    //               Physics::Elasticity::StandardTensors<dim>::I;
    //             const Tensor<1, dim, ScalarType>           Z = C * H;
    //             const SymmetricTensor<2, dim, ScalarType> &Ct =
    //               C; // Poor choice of indices used below...

    //             Tensor<3, dim, ScalarType> d2I6_dH_dC;
    //             for (unsigned int A = 0; A < dim; ++A)
    //               for (unsigned int B = 0; B < dim; ++B)
    //                 for (unsigned int C = 0; C < dim; ++C)
    //                   d2I6_dH_dC[A][B][C] += I[A][C] * Z[B] + I[B][C] * Z[A]
    //                   +
    //                                          Ct[B][C] * H[A] + Ct[A][C] *
    //                                          H[B];

    //             return d2I6_dH_dC;
    //           }
    //           break;
    //         case pI7a:
    //           {
    //             // I7 = [HxH]:C^{-1}
    //             const Tensor<1, dim, ScalarType> Y = C_inv * H;

    //             Tensor<3, dim, ScalarType> d2I7_dH_dC;
    //             for (unsigned int A = 0; A < dim; ++A)
    //               for (unsigned int B = 0; B < dim; ++B)
    //                 for (unsigned int C = 0; C < dim; ++C)
    //                   d2I7_dH_dC[A][B][C] -=
    //                     C_inv[C][A] * Y[B] + C_inv[C][B] * Y[A];

    //             return d2I7_dH_dC;
    //           }
    //           break;
    //         default:
    //           break;
    //       }

    //     AssertThrow(false,
    //                 ExcMessage("This (pseudo-)invariant is not defined."));
    //     return Tensor<3, dim, ScalarType>();
    //   }



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

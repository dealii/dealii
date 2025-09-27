// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_elasticity_kinematics_h
#define dealii_elasticity_kinematics_h


#include <deal.II/base/config.h>

#include <deal.II/base/numbers.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/physics/elasticity/standard_tensors.h>

DEAL_II_NAMESPACE_OPEN

namespace Physics
{
  namespace Elasticity
  {
    /**
     * A collection of tensor definitions for deformation and strain measures,
     * as well as a few special transformations, that conform to notation used
     * in standard scientific literature, in particular the books of
     * Holzapfel (2007) and Wriggers (2008). The citation for these references,
     * as well as other notation used here, can be found in the description for
     * the Physics::Elasticity namespace.
     *
     * @note These hold specifically for the codimension 0 case,
     * where the metric tensor is the identity tensor.
     */
    namespace Kinematics
    {
      /**
       * @name Deformation tensors
       */
      /** @{ */

      /**
       * Return the deformation gradient tensor,
       * as constructed from the material displacement gradient
       * tensor @p Grad_u.
       * The result is expressed as
       * @f[
       *  \mathbf{F}
       *    \dealcoloneq \nabla_{0} \boldsymbol{\varphi}
       *    \left( \mathbf{X} \right)
       *     =\mathbf{I} + \nabla_{0}\mathbf{u}
       * @f]
       * where $\mathbf{u} = \mathbf{u}\left(\mathbf{X}\right)$ is the
       * displacement at position $\mathbf{X}$ in the referential configuration.
       * The differential operator $\nabla_{0}$ is defined as
       * $\frac{\partial}{\partial \mathbf{X}}$.
       *
       * @dealiiWriggersA{23,3.14}
       * @dealiiHolzapfelA{71,2.39}
       */
      template <int dim, typename Number>
      DEAL_II_ALWAYS_INLINE Tensor<2, dim, Number>
                            F(const Tensor<2, dim, Number> &Grad_u);

      /**
       * Return the isochoric counterpart of the deformation gradient
       * tensor @p F .
       * The result is expressed as
       * @f[
       *  \mathbf{F}^{\text{iso}} \dealcoloneq J^{-1/\textrm{dim}} \mathbf{F}
       * @f]
       * where $J = \text{det}\left(\mathbf{F}\right)$.
       *
       * @dealiiWriggersA{29,3.28}
       * @dealiiHolzapfelA{228,6.79}
       */
      template <int dim, typename Number>
      DEAL_II_ALWAYS_INLINE Tensor<2, dim, Number>
                            F_iso(const Tensor<2, dim, Number> &F);

      /**
       * Return the volumetric counterpart of the deformation gradient
       * tensor @p F .
       * The result is expressed as
       * @f[
       *  \mathbf{F}^{\text{vol}} \dealcoloneq J^{1/\textrm{dim}} \mathbf{I}
       * @f]
       * where $J = \text{det}\left(\mathbf{F}\right)$.
       *
       * @dealiiWriggersA{29,3.28}
       * @dealiiHolzapfelA{228,6.79}
       */
      template <int dim, typename Number>
      DEAL_II_ALWAYS_INLINE SymmetricTensor<2, dim, Number>
                            F_vol(const Tensor<2, dim, Number> &F);

      /**
       * Return the symmetric right Cauchy-Green deformation tensor,
       * as constructed from the deformation gradient tensor @p F.
       * The result is expressed as
       * @f[
       *  \mathbf{C} \dealcoloneq \mathbf{F}^{T}\cdot\mathbf{F} \, .
       * @f]
       *
       * @dealiiWriggersA{23,3.15}
       * @dealiiHolzapfelA{78,2.65}
       */
      template <int dim, typename Number>
      DEAL_II_ALWAYS_INLINE SymmetricTensor<2, dim, Number>
                            C(const Tensor<2, dim, Number> &F);

      /**
       * Return the symmetric left Cauchy-Green deformation tensor,
       * as constructed from the deformation gradient tensor @p F.
       * The result is expressed as
       * @f[
       *  \mathbf{b} \dealcoloneq \mathbf{F}\cdot\mathbf{F}^{T} \, .
       * @f]
       *
       * @dealiiWriggersA{28,3.25}
       * @dealiiHolzapfelA{81,2.79}
       */
      template <int dim, typename Number>
      DEAL_II_ALWAYS_INLINE SymmetricTensor<2, dim, Number>
                            b(const Tensor<2, dim, Number> &F);

      /** @} */

      /**
       * @name Strain tensors
       */
      /** @{ */

      /**
       * Return the symmetric Green-Lagrange strain tensor,
       * as constructed from the deformation gradient tensor @p F.
       * The result is expressed as
       * @f[
       *  \mathbf{E} \dealcoloneq \frac{1}{2}
       *  \left[ \mathbf{F}^{T}\cdot\mathbf{F} - \mathbf{I} \right] \, .
       * @f]
       *
       * @dealiiWriggersA{23,3.15}
       * @dealiiHolzapfelA{79,6.29}
       */
      template <int dim, typename Number>
      DEAL_II_ALWAYS_INLINE SymmetricTensor<2, dim, Number>
                            E(const Tensor<2, dim, Number> &F);

      /**
       * Return the symmetric small strain tensor,
       * as constructed from the displacement gradient tensor @p Grad_u.
       * The result is expressed as
       * @f[
       *  \boldsymbol{\varepsilon} \dealcoloneq \frac{1}{2}
       *  \left[ \nabla_{0}\mathbf{u} + [\nabla_{0}\mathbf{u}]^{T} \right] \, .
       * @f]
       * where $\mathbf{u} = \mathbf{u}(\mathbf{X})$ is the displacement at
       * position
       * $\mathbf{X}$ in the referential configuration.
       * The differential operator $\nabla_{0}$ is defined as
       * $\frac{\partial}{\partial \mathbf{X}}$.
       *
       * @dealiiWriggersA{24,3.17}
       */
      template <int dim, typename Number>
      DEAL_II_ALWAYS_INLINE SymmetricTensor<2, dim, Number>
                            epsilon(const Tensor<2, dim, Number> &Grad_u);

      /**
       * Return the symmetric Almansi strain tensor,
       * as constructed from the deformation gradient tensor @p F.
       * The result is expressed as
       * @f[
       *  \mathbf{e} \dealcoloneq \frac{1}{2} \left[ \mathbf{I}
       *   - \mathbf{F}^{-T}\cdot\mathbf{F}^{-1} \right] \, .
       * @f]
       *
       * @dealiiWriggersA{30,3.35}
       * @dealiiHolzapfelA{81,2.83}
       */
      template <int dim, typename Number>
      DEAL_II_ALWAYS_INLINE SymmetricTensor<2, dim, Number>
                            e(const Tensor<2, dim, Number> &F);

      /** @} */

      /**
       * @name Strain rate tensors
       */
      /** @{ */

      /**
       * Return the spatial velocity gradient tensor,
       * as constructed from the deformation gradient tensor @p F
       * and its material time derivative @p dF_dt (the material velocity
       * gradient).
       * The result is expressed as
       * @f[
       *  \mathbf{l} \dealcoloneq \dot{\mathbf{F}}\cdot\mathbf{F}^{-1} \, .
       * @f]
       *
       * @dealiiWriggersA{32,3.47}
       * @dealiiHolzapfelA{96,2.141}
       */
      template <int dim, typename Number>
      Tensor<2, dim, Number>
      l(const Tensor<2, dim, Number> &F, const Tensor<2, dim, Number> &dF_dt);

      /**
       * Return the rate of deformation tensor (also known as the rate of strain
       * tensor), as constructed from the deformation gradient tensor @p F
       * and its material time derivative @p dF_dt (the material velocity
       * gradient).
       * The result is expressed as
       * @f[
       *  \mathbf{d} \dealcoloneq \frac{1}{2}
       *  \left[ \mathbf{l} + \mathbf{l}^{T} \right]
       * @f]
       * where
       * @f[
       *  \mathbf{l} = \dot{\mathbf{F}}\cdot\mathbf{F}^{-1}
       * @f]
       * is the spatial velocity gradient tensor.
       *
       * @dealiiWriggersA{32,3.49}
       * @dealiiHolzapfelA{97,2.148}
       */
      template <int dim, typename Number>
      SymmetricTensor<2, dim, Number>
      d(const Tensor<2, dim, Number> &F, const Tensor<2, dim, Number> &dF_dt);

      /**
       * Return the rate of rotation tensor (also known as the vorticity
       * tensor), as constructed from the deformation gradient tensor @p F
       * and its material time derivative @p dF_dt (the material velocity
       * gradient).
       * The result is expressed as
       * @f[
       *  \mathbf{w} \dealcoloneq \frac{1}{2}
       *  \left[ \mathbf{l} - \mathbf{l}^{T} \right]
       * @f]
       * where
       * @f[
       *  \mathbf{l} = \dot{\mathbf{F}}\cdot\mathbf{F}^{-1}
       * @f]
       * is the spatial velocity gradient tensor.
       *
       * @dealiiHolzapfelA{97,2.149}
       */
      template <int dim, typename Number>
      Tensor<2, dim, Number>
      w(const Tensor<2, dim, Number> &F, const Tensor<2, dim, Number> &dF_dt);

      /** @} */
    } // namespace Kinematics
  }   // namespace Elasticity
} // namespace Physics



#ifndef DOXYGEN

// ------------------------- inline functions ------------------------



template <int dim, typename Number>
inline Tensor<2, dim, Number>
Physics::Elasticity::Kinematics::F(const Tensor<2, dim, Number> &Grad_u)
{
  return StandardTensors<dim>::I + Grad_u;
}



template <int dim, typename Number>
inline Tensor<2, dim, Number>
Physics::Elasticity::Kinematics::F_iso(const Tensor<2, dim, Number> &F)
{
  // Make things work with AD types
  using std::pow;
  return internal::NumberType<Number>::value(pow(determinant(F), -1.0 / dim)) *
         F;
}



template <int dim, typename Number>
inline SymmetricTensor<2, dim, Number>
Physics::Elasticity::Kinematics::F_vol(const Tensor<2, dim, Number> &F)
{
  // Make things work with AD types
  using std::pow;
  return internal::NumberType<Number>::value(pow(determinant(F), 1.0 / dim)) *
         static_cast<SymmetricTensor<2, dim, Number>>(
           unit_symmetric_tensor<dim>());
}



template <int dim, typename Number>
inline SymmetricTensor<2, dim, Number>
Physics::Elasticity::Kinematics::C(const Tensor<2, dim, Number> &F)
{
  return symmetrize(transpose(F) * F);
}



template <int dim, typename Number>
inline SymmetricTensor<2, dim, Number>
Physics::Elasticity::Kinematics::b(const Tensor<2, dim, Number> &F)
{
  return symmetrize(F * transpose(F));
}



template <int dim, typename Number>
inline SymmetricTensor<2, dim, Number>
Physics::Elasticity::Kinematics::E(const Tensor<2, dim, Number> &F)
{
  return internal::NumberType<Number>::value(0.5) *
         (C(F) - static_cast<SymmetricTensor<2, dim, Number>>(
                   StandardTensors<dim>::I));
}



template <int dim, typename Number>
inline SymmetricTensor<2, dim, Number>
Physics::Elasticity::Kinematics::epsilon(const Tensor<2, dim, Number> &Grad_u)
{
  // This is the equivalent to 0.5*symmetrize(Grad_u + transpose(Grad_u));
  return symmetrize(Grad_u);
}



template <int dim, typename Number>
inline SymmetricTensor<2, dim, Number>
Physics::Elasticity::Kinematics::e(const Tensor<2, dim, Number> &F)
{
  const Tensor<2, dim, Number> F_inv = invert(F);
  return internal::NumberType<Number>::value(0.5) *
         symmetrize(static_cast<SymmetricTensor<2, dim, Number>>(
                      StandardTensors<dim>::I) -
                    transpose(F_inv) * F_inv);
}



template <int dim, typename Number>
inline Tensor<2, dim, Number>
Physics::Elasticity::Kinematics::l(const Tensor<2, dim, Number> &F,
                                   const Tensor<2, dim, Number> &dF_dt)
{
  return dF_dt * invert(F);
}



template <int dim, typename Number>
inline SymmetricTensor<2, dim, Number>
Physics::Elasticity::Kinematics::d(const Tensor<2, dim, Number> &F,
                                   const Tensor<2, dim, Number> &dF_dt)
{
  return symmetrize(l(F, dF_dt));
}



template <int dim, typename Number>
inline Tensor<2, dim, Number>
Physics::Elasticity::Kinematics::w(const Tensor<2, dim, Number> &F,
                                   const Tensor<2, dim, Number> &dF_dt)
{
  // This could be implemented as w = l-d, but that would mean computing "l"
  // a second time.
  const auto grad_v = l(F, dF_dt);
  return internal::NumberType<Number>::value(0.5) *
         (grad_v - transpose(grad_v));
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

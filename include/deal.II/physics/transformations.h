// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii__transformations_h
#define dealii__transformations_h

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/symmetric_tensor.h>

DEAL_II_NAMESPACE_OPEN


namespace Physics
{

  namespace Transformations
  {

    /**
     * Transformation functions and tensors that are defined in terms of
     * rotation angles and axes of rotation.
     *
     * @author Jean-Paul Pelteret, 2017
    */
    namespace Rotations
    {
      /**
       * @name Rotation matrices
       */
//@{

      /**
      * Returns the rotation matrix for 2-d Euclidean space, namely
      * @f[
      *  \mathbf{R} := \left[ \begin{array}{cc}
      *  cos(\theta) & sin(\theta) \\
      *  -sin(\theta) & cos(\theta)
      * \end{array}\right]
      * @f]
      * where $\theta$ is the rotation angle given in radians.
      * In particular, this describes the counter-clockwise rotation of a vector
      * relative to a
      * <a href="http://mathworld.wolfram.com/RotationMatrix.html">fixed set of right-handed axes</a>.
      *
      * @param[in] angle The rotation angle (about the z-axis) in radians
      */
      template<typename Number>
      Tensor<2,2,Number>
      rotation_matrix_2d (const Number &angle);


      /**
      * Returns the rotation matrix for 3-d Euclidean space.
      * Most concisely stated using the Rodrigues' rotation formula, this
      * function returns the equivalent of
      * @f[
      *  \mathbf{R} := cos(\theta)\mathbf{I} + sin(\theta)\mathbf{W}
      *              + (1-cos(\theta))\mathbf{u}\otimes\mathbf{u}
      * @f]
      * where $\mathbf{u}$ is the axial vector (an axial vector) and $\theta$
      * is the rotation angle given in radians, $\mathbf{I}$ is the identity
      * tensor and $\mathbf{W}$ is the skew symmetric tensor of $\mathbf{u}$.
      *
      * @dealiiWriggersA{374,9.194}
      * This presents Rodrigues' rotation formula, but the implementation used
      * in this function is described in this
      * <a href="https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle">wikipedia link</a>.
      * In particular, this describes the counter-clockwise rotation of a vector
      * <a href="http://mathworld.wolfram.com/RotationMatrix.html">in a plane with its normal</a>.
      * defined by the @p axis of rotation.
      * An alternative implementation is discussed at
      * <a href="https://www.gamedev.net/resources/_/technical/math-and-physics/do-we-really-need-quaternions-r1199">this link</a>,
      * but is inconsistent (sign-wise) with the Rodrigues' rotation formula as
      * it describes the rotation of a coordinate system.
      *
      * @param[in] axis  A unit vector that defines the axis of rotation
      * @param[in] angle The rotation angle in radians
      */
      template<typename Number>
      Tensor<2,3,Number>
      rotation_matrix_3d (const Point<3,Number> &axis,
                          const Number          &angle);

//@}

    }

    /**
     * Transformation of tensors that are defined in terms of a set of
     * contravariant bases. Rank-1 and rank-2 contravariant tensors
     * $\left(\bullet\right)^{\sharp} = \mathbf{T}$ (and its spatial counterpart
     * $\mathbf{t}$) typically satisfy the relation
     * @f[
     *    \int_{V_{0}} \nabla_{0} \cdot \mathbf{T} \; dV
     *      = \int_{\partial V_{0}} \mathbf{T} \cdot \mathbf{N} \; dA
     *      = \int_{\partial V_{t}} \mathbf{T} \cdot \mathbf{n} \; da
     *      = \int_{V_{t}} \nabla \cdot \mathbf{t} \; dv
     * @f]
     * where $V_{0}$ and $V_{t}$ are respectively control volumes in the
     * reference and spatial configurations, and their surfaces $\partial V_{0}$
     * and $\partial V_{t}$ have the outwards facing normals $\mathbf{N}$ and
     * $\mathbf{n}$.
     *
     * @author Jean-Paul Pelteret, Andrew McBride, 2016
    */
    namespace Contravariant
    {

      /**
       * @name Push forward operations
       */
//@{

      /**
      * Returns the result of the push forward transformation on a
      * contravariant vector, i.e.
      * @f[
      *  \chi\left(\bullet\right)^{\sharp}
      *    := \mathbf{F} \cdot \left(\bullet\right)^{\sharp}
      * @f]
      *
      * @param[in] V The (referential) vector to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\chi\left( \mathbf{V} \right)$
      */
      template <int dim, typename Number>
      Tensor<1,dim,Number>
      push_forward (const Tensor<1,dim,Number> &V,
                    const Tensor<2,dim,Number> &F);

      /**
      * Returns the result of the push forward transformation on a rank-2
      * contravariant tensor, i.e.
      * @f[
      *  \chi\left(\bullet\right)^{\sharp}
      *    := \mathbf{F} \cdot \left(\bullet\right)^{\sharp} \cdot \mathbf{F}^{T}
      * @f]
      *
      * @param[in] T The (referential) rank-2 tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\chi\left( \mathbf{T} \right)$
      */
      template <int dim, typename Number>
      Tensor<2,dim,Number>
      push_forward (const Tensor<2,dim,Number> &T,
                    const Tensor<2,dim,Number> &F);

      /**
      * Returns the result of the push forward transformation on a rank-2
      * contravariant symmetric tensor, i.e.
      * @f[
      *  \chi\left(\bullet\right)^{\sharp}
      *    := \mathbf{F} \cdot \left(\bullet\right)^{\sharp} \cdot \mathbf{F}^{T}
      * @f]
      *
      * @param[in] T The (referential) rank-2 symmetric tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\chi\left( \mathbf{T} \right)$
      */
      template <int dim, typename Number>
      SymmetricTensor<2,dim,Number>
      push_forward (const SymmetricTensor<2,dim,Number> &T,
                    const Tensor<2,dim,Number>          &F);

      /**
      * Returns the result of the push forward transformation on a rank-4
      * contravariant tensor, i.e. (in index notation)
      * @f[
      *  \left[ \chi\left(\bullet\right)^{\sharp} \right]_{ijkl}
      *    := F_{iI} F_{jJ} \left(\bullet\right)^{\sharp}_{IJKL} F_{kK} F_{lL}
      * @f]
      *
      * @param[in] H The (referential) rank-4 tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\chi\left( \mathbf{H} \right)$
      */
      template <int dim, typename Number>
      Tensor<4,dim,Number>
      push_forward (const Tensor<4,dim,Number> &H,
                    const Tensor<2,dim,Number> &F);

      /**
      * Returns the result of the push forward transformation on a rank-4
      * contravariant symmetric tensor, i.e. (in index notation)
      * @f[
      *  \left[ \chi\left(\bullet\right)^{\sharp} \right]_{ijkl}
      *    := F_{iI} F_{jJ} \left(\bullet\right)^{\sharp}_{IJKL} F_{kK} F_{lL}
      * @f]
      *
      * @param[in] H The (referential) rank-4 symmetric tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\chi\left( \mathbf{H} \right)$
      */
      template <int dim, typename Number>
      SymmetricTensor<4,dim,Number>
      push_forward (const SymmetricTensor<4,dim,Number> &H,
                    const Tensor<2,dim,Number>          &F);

//@}

      /**
       * @name Pull back operations
       */
//@{

      /**
      * Returns the result of the pull back transformation on a
      * contravariant vector, i.e.
      * @f[
      *  \chi^{-1}\left(\bullet\right)^{\sharp}
      *    := \mathbf{F}^{-1} \cdot \left(\bullet\right)^{\sharp}
      * @f]
      *
      * @param[in] v The (spatial) vector to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\chi^{-1}\left( \mathbf{v} \right)$
      */
      template <int dim, typename Number>
      Tensor<1,dim,Number>
      pull_back (const Tensor<1,dim,Number> &v,
                 const Tensor<2,dim,Number> &F);

      /**
      * Returns the result of the pull back transformation on a rank-2
      * contravariant tensor, i.e.
      * @f[
      *  \chi^{-1}\left(\bullet\right)^{\sharp}
      *    := \mathbf{F}^{-1} \cdot \left(\bullet\right)^{\sharp} \cdot \mathbf{F}^{-T}
      * @f]
      *
      * @param[in] t The (spatial) tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\chi^{-1}\left( \mathbf{t} \right)$
      */
      template <int dim, typename Number>
      Tensor<2,dim,Number>
      pull_back (const Tensor<2,dim,Number> &t,
                 const Tensor<2,dim,Number> &F);

      /**
      * Returns the result of the pull back transformation on a rank-2
      * contravariant symmetric tensor, i.e.
      * @f[
      *  \chi^{-1}\left(\bullet\right)^{\sharp}
      *    := \mathbf{F}^{-1} \cdot \left(\bullet\right)^{\sharp} \cdot \mathbf{F}^{-T}
      * @f]
      *
      * @param[in] t The (spatial) symmetric tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\chi^{-1}\left( \mathbf{t} \right)$
      */
      template <int dim, typename Number>
      SymmetricTensor<2,dim,Number>
      pull_back (const SymmetricTensor<2,dim,Number> &t,
                 const Tensor<2,dim,Number>          &F);

      /**
      * Returns the result of the pull back transformation on a rank-4
      * contravariant tensor, i.e. (in index notation)
      * @f[
      *  \left[ \chi^{-1}\left(\bullet\right)^{\sharp} \right]_{IJKL}
      *    := F^{-1}_{Ii} F^{-1}_{Jj} \left(\bullet\right)^{\sharp}_{ijkl} F^{-1}_{Kk} F^{-1}_{Ll}
      * @f]
      *
      * @param[in] h The (spatial) tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\chi^{-1}\left( \mathbf{h} \right)$
      */
      template <int dim, typename Number>
      Tensor<4,dim,Number>
      pull_back (const Tensor<4,dim,Number> &h,
                 const Tensor<2,dim,Number> &F);

      /**
      * Returns the result of the pull back transformation on a rank-4
      * contravariant symmetric tensor, i.e. (in index notation)
      * @f[
      *  \left[ \chi^{-1}\left(\bullet\right)^{\sharp} \right]_{IJKL}
      *    := F^{-1}_{Ii} F^{-1}_{Jj} \left(\bullet\right)^{\sharp}_{ijkl} F^{-1}_{Kk} F^{-1}_{Ll}
      * @f]
      *
      * @param[in] h The (spatial) symmetric tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\chi^{-1}\left( \mathbf{h} \right)$
      */
      template <int dim, typename Number>
      SymmetricTensor<4,dim,Number>
      pull_back (const SymmetricTensor<4,dim,Number> &h,
                 const Tensor<2,dim,Number>          &F);

//@}
    }

    /**
     * Transformation of tensors that are defined in terms of a set of
     * covariant basis vectors. Rank-1 and rank-2 covariant tensors
     * $\left(\bullet\right)^{\flat} = \mathbf{T}$ (and its spatial counterpart
     * $\mathbf{t}$) typically satisfy the relation
     * @f[
     *    \int_{\partial V_{0}} \left[ \nabla_{0} \times \mathbf{T} \right] \cdot \mathbf{N} \; dA
     *      = \oint_{\partial A_{0}} \mathbf{T} \cdot \mathbf{L} \; dL
     *      = \oint_{\partial A_{t}} \mathbf{t} \cdot \mathbf{l} \; dl
     *      = \int_{\partial V_{t}} \left[ \nabla \times \mathbf{t} \right] \cdot \mathbf{n} \; da
     * @f]
     * where the control surfaces $\partial V_{0}$ and $\partial V_{t}$ with
     * outwards facing normals $\mathbf{N}$ and $\mathbf{n}$ are bounded by the
     * curves $\partial A_{0}$ and $\partial A_{0}$ that are, respectively,
     * associated with the line directors $\mathbf{L}$ and $\mathbf{l}$.
     *
     * @author Jean-Paul Pelteret, Andrew McBride, 2016
    */
    namespace Covariant
    {

      /**
       * @name Push forward operations
       */
//@{

      /**
      * Returns the result of the push forward transformation on a covariant
      * vector, i.e.
      * @f[
      *  \chi\left(\bullet\right)^{\flat}
      *    := \mathbf{F}^{-T} \cdot \left(\bullet\right)^{\flat}
      * @f]
      *
      * @param[in] V The (referential) vector to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\chi\left( \mathbf{V} \right)$
      */
      template <int dim, typename Number>
      Tensor<1,dim,Number>
      push_forward (const Tensor<1,dim,Number> &V,
                    const Tensor<2,dim,Number> &F);

      /**
      * Returns the result of the push forward transformation on a rank-2
      * covariant tensor, i.e.
      * @f[
      *  \chi\left(\bullet\right)^{\flat}
      *    := \mathbf{F}^{-T} \cdot \left(\bullet\right)^{\flat} \cdot \mathbf{F}^{-1}
      * @f]
      *
      * @param[in] T The (referential) rank-2 tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\chi\left( \mathbf{T} \right)$
      */
      template <int dim, typename Number>
      Tensor<2,dim,Number>
      push_forward (const Tensor<2,dim,Number> &T,
                    const Tensor<2,dim,Number> &F);

      /**
      * Returns the result of the push forward transformation on a rank-2
      * covariant symmetric tensor, i.e.
      * @f[
      *  \chi\left(\bullet\right)^{\flat}
      *    := \mathbf{F}^{-T} \cdot \left(\bullet\right)^{\flat} \cdot \mathbf{F}^{-1}
      * @f]
      *
      * @param[in] T The (referential) rank-2 symmetric tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\chi\left( \mathbf{T} \right)$
      */
      template <int dim, typename Number>
      SymmetricTensor<2,dim,Number>
      push_forward (const SymmetricTensor<2,dim,Number> &T,
                    const Tensor<2,dim,Number>          &F);

      /**
      * Returns the result of the push forward transformation on a rank-4
      * covariant tensor, i.e. (in index notation)
      * @f[
      *  \left[ \chi\left(\bullet\right)^{\flat} \right]_{ijkl}
      *    := F^{-T}_{iI} F^{-T}_{jJ} \left(\bullet\right)^{\flat}_{IJKL} F^{-T}_{kK} F^{-T}_{lL}
      * @f]
      *
      * @param[in] H The (referential) rank-4 tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\chi\left( \mathbf{H} \right)$
      */
      template <int dim, typename Number>
      Tensor<4,dim,Number>
      push_forward (const Tensor<4,dim,Number> &H,
                    const Tensor<2,dim,Number> &F);

      /**
      * Returns the result of the push forward transformation on a rank-4
      * covariant symmetric tensor, i.e. (in index notation)
      * @f[
      *  \left[ \chi\left(\bullet\right)^{\flat} \right]_{ijkl}
      *    := F^{-T}_{iI} F^{-T}_{jJ} \left(\bullet\right)^{\flat}_{IJKL} F^{-T}_{kK} F^{-T}_{lL}
      * @f]
      *
      * @param[in] H The (referential) rank-4 symmetric tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\chi\left( \mathbf{H} \right)$
      */
      template <int dim, typename Number>
      SymmetricTensor<4,dim,Number>
      push_forward (const SymmetricTensor<4,dim,Number> &H,
                    const Tensor<2,dim,Number>          &F);

//@}

      /**
       * @name Pull back operations
       */
//@{

      /**
      * Returns the result of the pull back transformation on a
      * covariant vector, i.e.
      * @f[
      *  \chi^{-1}\left(\bullet\right)^{\flat}
      *    := \mathbf{F}^{T} \cdot \left(\bullet\right)^{\flat}
      * @f]
      *
      * @param[in] v The (spatial) vector to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\chi^{-1}\left( \mathbf{v} \right)$
      */
      template <int dim, typename Number>
      Tensor<1,dim,Number>
      pull_back (const Tensor<1,dim,Number> &v,
                 const Tensor<2,dim,Number> &F);

      /**
      * Returns the result of the pull back transformation on a rank-2
      * covariant tensor, i.e.
      * @f[
      *  \chi^{-1}\left(\bullet\right)^{\flat}
      *    := \mathbf{F}^{T} \cdot \left(\bullet\right)^{\flat} \cdot \mathbf{F}
      * @f]
      *
      * @param[in] t The (spatial) tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\chi^{-1}\left( \mathbf{t} \right)$
      */
      template <int dim, typename Number>
      Tensor<2,dim,Number>
      pull_back (const Tensor<2,dim,Number> &t,
                 const Tensor<2,dim,Number> &F);

      /**
      * Returns the result of the pull back transformation on a rank-2
      * covariant symmetric tensor, i.e.
      * @f[
      *  \chi^{-1}\left(\bullet\right)^{\flat}
      *    := \mathbf{F}^{T} \cdot \left(\bullet\right)^{\flat} \cdot \mathbf{F}
      * @f]
      *
      * @param[in] t The (spatial) symmetric tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\chi^{-1}\left( \mathbf{t} \right)$
      */
      template <int dim, typename Number>
      SymmetricTensor<2,dim,Number>
      pull_back (const SymmetricTensor<2,dim,Number> &t,
                 const Tensor<2,dim,Number>          &F);

      /**
      * Returns the result of the pull back transformation on a rank-4
      * contravariant tensor, i.e. (in index notation)
      * @f[
      *  \left[ \chi^{-1}\left(\bullet\right)^{\flat} \right]_{IJKL}
      *    := F^{T}_{Ii} F^{T}_{Jj} \left(\bullet\right)^{\flat}_{ijkl} F^{T}_{Kk} F^{T}_{Ll}
      * @f]
      *
      * @param[in] h The (spatial) tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\chi^{-1}\left( \mathbf{h} \right)$
      */
      template <int dim, typename Number>
      Tensor<4,dim,Number>
      pull_back (const Tensor<4,dim,Number> &h,
                 const Tensor<2,dim,Number> &F);

      /**
      * Returns the result of the pull back transformation on a rank-4
      * contravariant symmetric tensor, i.e. (in index notation)
      * @f[
      *  \left[ \chi^{-1}\left(\bullet\right)^{\flat} \right]_{IJKL}
      *    := F^{T}_{Ii} F^{T}_{Jj} \left(\bullet\right)^{\flat}_{ijkl} F^{T}_{Kk} F^{T}_{Ll}
      * @f]
      *
      * @param[in] h The (spatial) symmetric tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\chi^{-1}\left( \mathbf{h} \right)$
      */
      template <int dim, typename Number>
      SymmetricTensor<4,dim,Number>
      pull_back (const SymmetricTensor<4,dim,Number> &h,
                 const Tensor<2,dim,Number>          &F);

//@}
    }

    /**
     * Transformation of tensors that are defined in terms of a set of
     * contravariant basis vectors and scale with the inverse of the volume
     * change associated with the mapping.
     *
     * @author Jean-Paul Pelteret, Andrew McBride, 2016
    */
    namespace Piola
    {

      /**
       * @name Push forward operations
       */
//@{

      /**
      * Returns the result of the push forward transformation on a
      * contravariant vector, i.e.
      * @f[
      *  \textrm{det} \mathbf{F}^{-1} \; \chi\left(\bullet\right)^{\sharp}
      *    := \frac{1}{\textrm{det} \mathbf{F}} \; \mathbf{F} \cdot \left(\bullet\right)^{\sharp}
      * @f]
      *
      * @param[in] V The (referential) vector to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\frac{1}{\textrm{det} \mathbf{F}} \; \chi\left( \mathbf{V} \right)$
      */
      template <int dim, typename Number>
      Tensor<1,dim,Number>
      push_forward (const Tensor<1,dim,Number> &V,
                    const Tensor<2,dim,Number> &F);

      /**
      * Returns the result of the push forward transformation on a rank-2
      * contravariant tensor, i.e.
      * @f[
      *  \textrm{det} \mathbf{F}^{-1} \; \chi\left(\bullet\right)^{\sharp}
      *    := \frac{1}{\textrm{det} \mathbf{F}} \; \mathbf{F} \cdot \left(\bullet\right)^{\sharp} \cdot \mathbf{F}^{T}
      * @f]
      *
      * @param[in] T The (referential) rank-2 tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\frac{1}{\textrm{det} \mathbf{F}} \; \chi\left( \mathbf{T} \right)$
      */
      template <int dim, typename Number>
      Tensor<2,dim,Number>
      push_forward (const Tensor<2,dim,Number> &T,
                    const Tensor<2,dim,Number> &F);

      /**
      * Returns the result of the push forward transformation on a rank-2
      * contravariant symmetric tensor, i.e.
      * @f[
      *  \textrm{det} \mathbf{F}^{-1} \; \chi\left(\bullet\right)^{\sharp}
      *    := \frac{1}{\textrm{det} \mathbf{F}} \; \mathbf{F} \cdot \left(\bullet\right)^{\sharp} \cdot \mathbf{F}^{T}
      * @f]
      *
      * @param[in] T The (referential) rank-2 symmetric tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\frac{1}{\textrm{det} \mathbf{F}} \; \chi\left( \mathbf{T} \right)$
      */
      template <int dim, typename Number>
      SymmetricTensor<2,dim,Number>
      push_forward (const SymmetricTensor<2,dim,Number> &T,
                    const Tensor<2,dim,Number>          &F);

      /**
      * Returns the result of the push forward transformation on a rank-4
      * contravariant tensor, i.e. (in index notation)
      * @f[
      *  \textrm{det} \mathbf{F}^{-1} \; \left[ \chi\left(\bullet\right)^{\sharp} \right]_{ijkl}
      *    := \frac{1}{\textrm{det} \mathbf{F}} \; F_{iI} F_{jJ} \left(\bullet\right)^{\sharp}_{IJKL} F_{kK} F_{lL}
      * @f]
      *
      * @param[in] H The (referential) rank-4 tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\frac{1}{\textrm{det} \mathbf{F}} \; \chi\left( \mathbf{H} \right)$
      */
      template <int dim, typename Number>
      Tensor<4,dim,Number>
      push_forward (const Tensor<4,dim,Number> &H,
                    const Tensor<2,dim,Number> &F);

      /**
      * Returns the result of the push forward transformation on a rank-4
      * contravariant symmetric tensor, i.e. (in index notation)
      * @f[
      *  \textrm{det} \mathbf{F}^{-1} \; \left[ \chi\left(\bullet\right)^{\sharp} \right]_{ijkl}
      *    := \frac{1}{\textrm{det} \mathbf{F}} \; F_{iI} F_{jJ} \left(\bullet\right)^{\sharp}_{IJKL} F_{kK} F_{lL}
      * @f]
      *
      * @param[in] H The (referential) rank-4 symmetric tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\frac{1}{\textrm{det} \mathbf{F}} \; \chi\left( \mathbf{H} \right)$
      */
      template <int dim, typename Number>
      SymmetricTensor<4,dim,Number>
      push_forward (const SymmetricTensor<4,dim,Number> &H,
                    const Tensor<2,dim,Number>          &F);

//@}

      /**
       * @name Pull back operations
       */
//@{

      /**
      * Returns the result of the pull back transformation on a
      * contravariant vector, i.e.
      * @f[
      *  \textrm{det} \mathbf{F} \; \chi^{-1}\left(\bullet\right)^{\sharp}
      *    := \textrm{det} \mathbf{F} \; \mathbf{F}^{-1} \cdot \left(\bullet\right)^{\sharp}
      * @f]
      *
      * @param[in] v The (spatial) vector to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\textrm{det} \mathbf{F} \; \chi^{-1}\left( \mathbf{v} \right)$
      */
      template <int dim, typename Number>
      Tensor<1,dim,Number>
      pull_back (const Tensor<1,dim,Number> &v,
                 const Tensor<2,dim,Number> &F);

      /**
      * Returns the result of the pull back transformation on a rank-2
      * contravariant tensor, i.e.
      * @f[
      *  \textrm{det} \mathbf{F} \; \chi^{-1}\left(\bullet\right)^{\sharp}
      *    := \textrm{det} \mathbf{F} \; \mathbf{F}^{-1} \cdot \left(\bullet\right)^{\sharp} \cdot \mathbf{F}^{-T}
      * @f]
      *
      * @param[in] t The (spatial) tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\textrm{det} \mathbf{F} \; \chi^{-1}\left( \mathbf{t} \right)$
      */
      template <int dim, typename Number>
      Tensor<2,dim,Number>
      pull_back (const Tensor<2,dim,Number> &t,
                 const Tensor<2,dim,Number> &F);

      /**
      * Returns the result of the pull back transformation on a rank-2
      * contravariant symmetric tensor, i.e.
      * @f[
      *  \textrm{det} \mathbf{F} \; \chi^{-1}\left(\bullet\right)^{\sharp}
      *    := \textrm{det} \mathbf{F} \; \mathbf{F}^{-1} \cdot \left(\bullet\right)^{\sharp} \cdot \mathbf{F}^{-T}
      * @f]
      *
      * @param[in] t The (spatial) symmetric tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\textrm{det} \mathbf{F} \; \chi^{-1}\left( \mathbf{t} \right)$
      */
      template <int dim, typename Number>
      SymmetricTensor<2,dim,Number>
      pull_back (const SymmetricTensor<2,dim,Number> &t,
                 const Tensor<2,dim,Number>          &F);

      /**
      * Returns the result of the pull back transformation on a rank-4
      * contravariant tensor, i.e. (in index notation)
      * @f[
      *  \textrm{det} \mathbf{F} \; \left[ \chi^{-1}\left(\bullet\right)^{\sharp} \right]_{IJKL}
      *    := \textrm{det} \mathbf{F} \; F^{-1}_{Ii} F^{-1}_{Jj} \left(\bullet\right)^{\sharp}_{ijkl} F^{-1}_{Kk} F^{-1}_{Ll}
      * @f]
      *
      * @param[in] h The (spatial) tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\textrm{det} \mathbf{F} \; \chi^{-1}\left( \mathbf{h} \right)$
      */
      template <int dim, typename Number>
      Tensor<4,dim,Number>
      pull_back (const Tensor<4,dim,Number> &h,
                 const Tensor<2,dim,Number> &F);

      /**
      * Returns the result of the pull back transformation on a rank-4
      * contravariant symmetric tensor, i.e. (in index notation)
      * @f[
      *  \textrm{det} \mathbf{F} \; \left[ \chi^{-1}\left(\bullet\right)^{\sharp} \right]_{IJKL}
      *    := \textrm{det} \mathbf{F} \; F^{-1}_{Ii} F^{-1}_{Jj} \left(\bullet\right)^{\sharp}_{ijkl} F^{-1}_{Kk} F^{-1}_{Ll}
      * @f]
      *
      * @param[in] h The (spatial) symmetric tensor to be operated on
      * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
      * @return      $\textrm{det} \mathbf{F} \; \chi^{-1}\left( \mathbf{h} \right)$
      */
      template <int dim, typename Number>
      SymmetricTensor<4,dim,Number>
      pull_back (const SymmetricTensor<4,dim,Number> &h,
                 const Tensor<2,dim,Number>          &F);

//@}
    }

    /**
     * @name Special operations
     */
//@{

    /**
    * Returns the result of applying Nanson's formula for the transformation of
    * the material surface area element $d\mathbf{A}$ to the current surfaces
    * area element $d\mathbf{a}$ under the nonlinear transformation map
    * $\mathbf{x} = \boldsymbol{\varphi} \left( \mathbf{X} \right)$.
    *
    * The returned result is the spatial normal scaled by the ratio of areas
    * between the reference and spatial surface elements, i.e.
    * @f[
    *  \mathbf{n} \frac{da}{dA}
    *    := \textrm{det} \mathbf{F} \, \mathbf{F}^{-T} \cdot \mathbf{N}
    *     = \textrm{cof} \mathbf{F} \cdot \mathbf{N} \, .
    * @f]
    *
    * @param[in] N The referential normal unit vector $\mathbf{N}$
    * @param[in] F The deformation gradient tensor $\mathbf{F} \left( \mathbf{X} \right)$
    * @return      The scaled spatial normal vector $\mathbf{n} \frac{da}{dA}$
    *
    * @dealiiHolzapfelA{75,2.55}
    * @dealiiWriggersA{23,3.11}
    */
    template<int dim, typename Number>
    Tensor<1,dim,Number>
    nansons_formula (const Tensor<1,dim,Number> &N,
                     const Tensor<2,dim,Number> &F);

//@}
  }
}



#ifndef DOXYGEN

// ------------------------- inline functions ------------------------

namespace internal
{
  namespace Physics
  {
    namespace
    {
      template <int dim, typename Number>
      inline
      Tensor<1,dim,Number>
      transformation_contraction (const Tensor<1,dim,Number> &V,
                                  const Tensor<2,dim,Number> &F)
      {
        return contract<1,0>(F, V);
      }



      template <int dim, typename Number>
      inline
      Tensor<2,dim,Number>
      transformation_contraction (const Tensor<2,dim,Number> &T,
                                  const Tensor<2,dim,Number> &F)
      {
        return contract<1,1>(F,contract<1,0>(F, T));
      }



      template <int dim, typename Number>
      inline
      dealii::SymmetricTensor<2,dim,Number>
      transformation_contraction (const dealii::SymmetricTensor<2,dim,Number> &T,
                                  const Tensor<2,dim,Number>                  &F)
      {
        Tensor<2,dim,Number> tmp_1;
        for (unsigned int i=0; i<dim; ++i)
          for (unsigned int J=0; J<dim; ++J)
            for (unsigned int I=0; I<dim; ++I)
              tmp_1[i][J] += F[i][I] * T[I][J];

        dealii::SymmetricTensor<2,dim,Number> out;
        for (unsigned int i=0; i<dim; ++i)
          for (unsigned int j=i; j<dim; ++j)
            for (unsigned int J=0; J<dim; ++J)
              out[i][j] += F[j][J] * tmp_1[i][J];

        return out;
      }



      template <int dim, typename Number>
      inline
      Tensor<4,dim,Number>
      transformation_contraction (const Tensor<4,dim,Number> &H,
                                  const Tensor<2,dim,Number> &F)
      {
        // Its significantly quicker (in 3d) to push forward
        // each index individually
        return contract<1,3>(F,contract<1,2>(F,contract<1,1>(F,contract<1,0>(F, H))));
      }



      template <int dim, typename Number>
      inline
      dealii::SymmetricTensor<4,dim,Number>
      transformation_contraction (const dealii::SymmetricTensor<4,dim,Number> &H,
                                  const Tensor<2,dim,Number>                  &F)
      {
        // Its significantly quicker (in 3d) to push forward
        // each index individually

        Tensor<4,dim,Number> tmp_1;
        for (unsigned int i=0; i<dim; ++i)
          for (unsigned int J=0; J<dim; ++J)
            for (unsigned int K=0; K<dim; ++K)
              for (unsigned int L=0; L<dim; ++L)
                for (unsigned int I=0; I<dim; ++I)
                  tmp_1[i][J][K][L] += F[i][I] * H[I][J][K][L];

        Tensor<4,dim,Number> tmp_2;
        for (unsigned int i=0; i<dim; ++i)
          for (unsigned int j=0; j<dim; ++j)
            for (unsigned int K=0; K<dim; ++K)
              for (unsigned int L=0; L<dim; ++L)
                for (unsigned int J=0; J<dim; ++J)
                  tmp_2[i][j][K][L] += F[j][J] * tmp_1[i][J][K][L];

        tmp_1 = 0.0;
        for (unsigned int i=0; i<dim; ++i)
          for (unsigned int j=0; j<dim; ++j)
            for (unsigned int k=0; k<dim; ++k)
              for (unsigned int L=0; L<dim; ++L)
                for (unsigned int K=0; K<dim; ++K)
                  tmp_1[i][j][k][L] += F[k][K] * tmp_2[i][j][K][L];

        dealii::SymmetricTensor<4,dim,Number> out;
        for (unsigned int i=0; i<dim; ++i)
          for (unsigned int j=i; j<dim; ++j)
            for (unsigned int k=0; k<dim; ++k)
              for (unsigned int l=k; l<dim; ++l)
                for (unsigned int L=0; L<dim; ++L)
                  out[i][j][k][l] += F[l][L] * tmp_1[i][j][k][L];

        return out;
      }
    }
  }
}



template<typename Number>
Tensor<2,2,Number>
Physics::Transformations::Rotations::rotation_matrix_2d (const Number &angle)
{
  const double rotation[2][2]
  = {{
      std::cos(angle) , -std::sin(angle)
    },
    {
      std::sin(angle) , std::cos(angle)
    }
  };
  return Tensor<2,2> (rotation);
}



template<typename Number>
Tensor<2,3,Number>
Physics::Transformations::Rotations::rotation_matrix_3d (const Point<3,Number> &axis,
                                                         const Number          &angle)
{
  Assert(std::abs(axis.norm() - 1.0) < 1e-9,
         ExcMessage("The supplied axial vector is not a unit vector."));
  const Number c = std::cos(angle);
  const Number s = std::sin(angle);
  const Number t = 1.-c;
  const double rotation[3][3]
  = {{
      t *axis[0] *axis[0] + c,
      t *axis[0] *axis[1] - s *axis[2],
      t *axis[0] *axis[2] + s *axis[1]
    },
    {
      t *axis[0] *axis[1] + s *axis[2],
      t *axis[1] *axis[1] + c,
      t *axis[1] *axis[2] - s *axis[0]
    },
    {
      t *axis[0] *axis[2] - s *axis[1],
      t *axis[1] *axis[2] + s *axis[0],
      t *axis[2] *axis[2] + c
    }
  };
  return Tensor<2,3,Number>(rotation);
}



template <int dim, typename Number>
inline
Tensor<1,dim,Number>
Physics::Transformations::Contravariant::push_forward (const Tensor<1,dim,Number> &V,
                                                       const Tensor<2,dim,Number> &F)
{
  return internal::Physics::transformation_contraction(V,F);
}



template <int dim, typename Number>
inline
Tensor<2,dim,Number>
Physics::Transformations::Contravariant::push_forward (const Tensor<2,dim,Number> &T,
                                                       const Tensor<2,dim,Number> &F)
{
  return internal::Physics::transformation_contraction(T,F);
}



template <int dim, typename Number>
inline
SymmetricTensor<2,dim,Number>
Physics::Transformations::Contravariant::push_forward (const SymmetricTensor<2,dim,Number> &T,
                                                       const Tensor<2,dim,Number>          &F)
{
  return internal::Physics::transformation_contraction(T,F);
}



template <int dim, typename Number>
inline
Tensor<4,dim,Number>
Physics::Transformations::Contravariant::push_forward (const Tensor<4,dim,Number> &H,
                                                       const Tensor<2,dim,Number> &F)
{
  return internal::Physics::transformation_contraction(H,F);
}



template <int dim, typename Number>
inline
SymmetricTensor<4,dim,Number>
Physics::Transformations::Contravariant::push_forward (const SymmetricTensor<4,dim,Number> &H,
                                                       const Tensor<2,dim,Number>          &F)
{
  return internal::Physics::transformation_contraction(H,F);
}



template <int dim, typename Number>
inline
Tensor<1,dim,Number>
Physics::Transformations::Contravariant::pull_back (const Tensor<1,dim,Number> &v,
                                                    const Tensor<2,dim,Number> &F)
{
  return internal::Physics::transformation_contraction(v,invert(F));
}



template <int dim, typename Number>
inline
Tensor<2,dim,Number>
Physics::Transformations::Contravariant::pull_back (const Tensor<2,dim,Number> &t,
                                                    const Tensor<2,dim,Number> &F)
{
  return internal::Physics::transformation_contraction(t,invert(F));
}



template <int dim, typename Number>
inline
SymmetricTensor<2,dim,Number>
Physics::Transformations::Contravariant::pull_back (const SymmetricTensor<2,dim,Number> &t,
                                                    const Tensor<2,dim,Number>          &F)
{
  return internal::Physics::transformation_contraction(t,invert(F));
}



template <int dim, typename Number>
inline
Tensor<4,dim,Number>
Physics::Transformations::Contravariant::pull_back (const Tensor<4,dim,Number> &h,
                                                    const Tensor<2,dim,Number> &F)
{
  return internal::Physics::transformation_contraction(h,invert(F));
}



template <int dim, typename Number>
inline
SymmetricTensor<4,dim,Number>
Physics::Transformations::Contravariant::pull_back (const SymmetricTensor<4,dim,Number> &h,
                                                    const Tensor<2,dim,Number>          &F)
{
  return internal::Physics::transformation_contraction(h,invert(F));
}



template <int dim, typename Number>
inline
Tensor<1,dim,Number>
Physics::Transformations::Covariant::push_forward (const Tensor<1,dim,Number> &V,
                                                   const Tensor<2,dim,Number> &F)
{
  return internal::Physics::transformation_contraction(V,transpose(invert(F)));
}



template <int dim, typename Number>
inline
Tensor<2,dim,Number>
Physics::Transformations::Covariant::push_forward (const Tensor<2,dim,Number> &T,
                                                   const Tensor<2,dim,Number> &F)
{
  return internal::Physics::transformation_contraction(T,transpose(invert(F)));
}



template <int dim, typename Number>
inline
SymmetricTensor<2,dim,Number>
Physics::Transformations::Covariant::push_forward (const SymmetricTensor<2,dim,Number> &T,
                                                   const Tensor<2,dim,Number>          &F)
{
  return internal::Physics::transformation_contraction(T,transpose(invert(F)));
}



template <int dim, typename Number>
inline
Tensor<4,dim,Number>
Physics::Transformations::Covariant::push_forward (const Tensor<4,dim,Number> &H,
                                                   const Tensor<2,dim,Number> &F)
{
  return internal::Physics::transformation_contraction(H,transpose(invert(F)));
}



template <int dim, typename Number>
inline
SymmetricTensor<4,dim,Number>
Physics::Transformations::Covariant::push_forward (const SymmetricTensor<4,dim,Number> &H,
                                                   const Tensor<2,dim,Number>          &F)
{
  return internal::Physics::transformation_contraction(H,transpose(invert(F)));
}



template <int dim, typename Number>
inline
Tensor<1,dim,Number>
Physics::Transformations::Covariant::pull_back (const Tensor<1,dim,Number> &v,
                                                const Tensor<2,dim,Number> &F)
{
  return internal::Physics::transformation_contraction(v,transpose(F));
}



template <int dim, typename Number>
inline
Tensor<2,dim,Number>
Physics::Transformations::Covariant::pull_back (const Tensor<2,dim,Number> &t,
                                                const Tensor<2,dim,Number> &F)
{
  return internal::Physics::transformation_contraction(t,transpose(F));
}



template <int dim, typename Number>
inline
SymmetricTensor<2,dim,Number>
Physics::Transformations::Covariant::pull_back (const SymmetricTensor<2,dim,Number> &t,
                                                const Tensor<2,dim,Number>          &F)
{
  return internal::Physics::transformation_contraction(t,transpose(F));
}



template <int dim, typename Number>
inline
Tensor<4,dim,Number>
Physics::Transformations::Covariant::pull_back (const Tensor<4,dim,Number> &h,
                                                const Tensor<2,dim,Number> &F)
{
  return internal::Physics::transformation_contraction(h,transpose(F));
}



template <int dim, typename Number>
inline
SymmetricTensor<4,dim,Number>
Physics::Transformations::Covariant::pull_back (const SymmetricTensor<4,dim,Number> &h,
                                                const Tensor<2,dim,Number>          &F)
{
  return internal::Physics::transformation_contraction(h,transpose(F));
}



template <int dim, typename Number>
inline
Tensor<1,dim,Number>
Physics::Transformations::Piola::push_forward (const Tensor<1,dim,Number> &V,
                                               const Tensor<2,dim,Number> &F)
{
  return Number(1.0/determinant(F))*Contravariant::push_forward(V,F);
}



template <int dim, typename Number>
inline
Tensor<2,dim,Number>
Physics::Transformations::Piola::push_forward (const Tensor<2,dim,Number> &T,
                                               const Tensor<2,dim,Number> &F)
{
  return Number(1.0/determinant(F))*Contravariant::push_forward(T,F);
}



template <int dim, typename Number>
inline
SymmetricTensor<2,dim,Number>
Physics::Transformations::Piola::push_forward (const SymmetricTensor<2,dim,Number> &T,
                                               const Tensor<2,dim,Number>          &F)
{
  return Number(1.0/determinant(F))*Contravariant::push_forward(T,F);
}



template <int dim, typename Number>
inline
Tensor<4,dim,Number>
Physics::Transformations::Piola::push_forward (const Tensor<4,dim,Number> &H,
                                               const Tensor<2,dim,Number> &F)
{
  return Number(1.0/determinant(F))*Contravariant::push_forward(H,F);
}



template <int dim, typename Number>
inline
SymmetricTensor<4,dim,Number>
Physics::Transformations::Piola::push_forward (const SymmetricTensor<4,dim,Number> &H,
                                               const Tensor<2,dim,Number>          &F)
{
  return Number(1.0/determinant(F))*Contravariant::push_forward(H,F);
}



template <int dim, typename Number>
inline
Tensor<1,dim,Number>
Physics::Transformations::Piola::pull_back (const Tensor<1,dim,Number> &v,
                                            const Tensor<2,dim,Number> &F)
{
  return Number(determinant(F))*Contravariant::pull_back(v,F);
}



template <int dim, typename Number>
inline
Tensor<2,dim,Number>
Physics::Transformations::Piola::pull_back (const Tensor<2,dim,Number> &t,
                                            const Tensor<2,dim,Number> &F)
{
  return Number(determinant(F))*Contravariant::pull_back(t,F);
}



template <int dim, typename Number>
inline
SymmetricTensor<2,dim,Number>
Physics::Transformations::Piola::pull_back (const SymmetricTensor<2,dim,Number> &t,
                                            const Tensor<2,dim,Number>          &F)
{
  return Number(determinant(F))*Contravariant::pull_back(t,F);
}



template <int dim, typename Number>
inline
Tensor<4,dim,Number>
Physics::Transformations::Piola::pull_back (const Tensor<4,dim,Number> &h,
                                            const Tensor<2,dim,Number> &F)
{
  return Number(determinant(F))*Contravariant::pull_back(h,F);
}



template <int dim, typename Number>
inline
SymmetricTensor<4,dim,Number>
Physics::Transformations::Piola::pull_back (const SymmetricTensor<4,dim,Number> &h,
                                            const Tensor<2,dim,Number>          &F)
{
  return Number(determinant(F))*Contravariant::pull_back(h,F);
}



template<int dim, typename Number>
inline Tensor<1,dim,Number>
Physics::Transformations::nansons_formula (const Tensor<1,dim,Number> &N,
                                           const Tensor<2,dim,Number> &F)
{
  return cofactor(F)*N;
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

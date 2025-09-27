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

/**
 * @defgroup physics Physics
 *
 * @brief A group dedicated to the implementation of functions and
 * classes that relate to continuum physics, physical fields and materials.
 */

/**
 * A collection of namespaces and utilities to assist in the
 * definition, construction and manipulation of data related to
 * physical fields and materials.
 */
namespace Physics
{

  /**
   * Notations that reduce the order of tensors, effectively storing them in some
   * sort of consistent compressed storage pattern. An example is storing the
   * 6 independent components of $3\times 3$ symmetric tensors of rank 2 as a
   * vector with 6 components, and then representing the 36 independent elements
   * of symmetric $3\times 3 \times 3\times 3$ tensors of rank 4 (which when
   * applied to a symmetric rank-2 tensor yields another symmetric rank-2 tensor)
   * as a $6 \times 6$ matrix.
   *
   * Although this method of representing tensors is most regularly associated with
   * the efficient storage of the fourth-order elasticity tensor, with its
   * generalization it has wider applicability. This representation is also common
   * in the physics, material science and FEM literature.
   *
   * There are several variations of tensor notation, each a slightly different
   * structure. The primary difference between the various forms of tensor notation
   * is the weighting prescribed to the various elements of the compressed tensors.
   * This <a
   * href="https://en.wikipedia.org/wiki/Voigt_notation">wikipedia article</a> has
   * some further general insights on this topic.
   *
   * @ingroup physics
   */
  namespace Notation
  {

  }

  /**
  * A collection of operations to assist in the transformation of tensor
  * quantities from the reference to spatial configuration, and vice versa.
  * These types of transformation are typically used to re-express quantities
  * measured or computed in one configuration in terms of a second configuration.
  *
  * <h3>Notation</h3>
  *
  * We will use the same notation for the coordinates $\mathbf{X}, \mathbf{x}$,
  * transformations $\varphi$, differential operator $\nabla_{0}$ and deformation
  * gradient $\mathbf{F}$ as discussed for namespace Physics::Elasticity.
  *
  * As a further point on notation, we will follow Holzapfel (2007) and denote
  * the push forward transformation as $\chi\left(\bullet\right)$ and
  * the pull back transformation as $\chi^{-1}\left(\bullet\right)$.
  * We will also use the annotation $\left(\bullet\right)^{\sharp}$ to indicate
  * that a tensor $\left(\bullet\right)$ is a contravariant tensor,
  * and $\left(\bullet\right)^{\flat}$ that it is covariant. In other
  * words, these indices do not actually change the tensor, they just indicate
  * the <i>kind</i> of object a particular tensor is.
  *
  * @note For these transformations, unless otherwise stated, we will strictly
  * assume that all indices of the transformed tensors derive from one coordinate
  * system; that is to say that they are not multi-point tensors (such as the
  * Piola stress in elasticity).
  *
  * @ingroup physics
  */
  namespace Transformations
  {
  }

  /**
   * This namespace provides a collection of definitions that
   * conform to standard notation used in (nonlinear) elasticity.
   *
   * <h3>Notation</h3>
   *
   * References for this notation include:
   * @code{.bib}
       @Book{Holzapfel2007a,
          title =     {Nonlinear solid mechanics. A Continuum Approach for Engineering},
          publisher = {John Wiley \& Sons Ltd.},
          year =      {2007},
          author =    {Holzapfel, G. A.},
          address =   {West Sussex, England},
          note =      {ISBN: 0-471-82304-X}
        }
        @Book{Wriggers2008a,
          title =     {Nonlinear finite element methods},
          publisher = {Springer Berlin Heidelberg},
          year =      {2008},
          author =    {Wriggers, P.},
          volume =    {4},
          address =   {Berlin, Germany},
          note =      {ISBN: 978-3-540-71000-4},
          doi =       {10.1007/978-3-540-71001-1}
        }
   * @endcode
   *
   * For convenience we will predefine some commonly referenced tensors and
   * operations.
   * Considering the position vector $\mathbf{X}$ in the referential (material)
   * configuration, points $\mathbf{X}$ are transformed to points $\mathbf{x}$
   * in the current (spatial) configuration through the nonlinear map
   * @f[
   *  \mathbf{x}
   *   \dealcoloneq \boldsymbol{\varphi} \left( \mathbf{X} \right)
   *    = \mathbf{X} + \mathbf{u}(\mathbf{X}) \, ,
   * @f]
   * where the $\mathbf{u}(\mathbf{X})$ represents the displacement vector.
   * From this we can compute the deformation gradient tensor as
   * @f[
   *  \mathbf{F} \dealcoloneq \mathbf{I} + \nabla_{0}\mathbf{u} \, ,
   * @f]
   * wherein the differential operator $\nabla_{0}$ is defined as
   * $\frac{\partial}{\partial \mathbf{X}}$ and $\mathbf{I}$ is the identity
   * tensor.
   *
   * Finally, two common tensor operators are represented by $\cdot$ and $:$
   * operators. These respectively represent a single and double contraction over
   * the inner tensor indices.
   * Vectors and second-order tensors are highlighted by bold font, while
   * fourth-order tensors are denoted by calliagraphic font.
   *
   * One can think of fourth-order tensors as linear operators mapping second-order
   * tensors (matrices) onto themselves in much the same way as matrices map
   * vectors onto vectors.
   * To provide some context to the implemented class members and functions,
   * consider the following fundamental operations performed on tensors with special
   * properties:
   *
   * If we represent a general second-order tensor as $\mathbf{A}$, then the general
   * fourth-order unit tensors $\mathcal{I}$ and $\overline{\mathcal{I}}$ are
   * defined by
   * @f[
   *  \mathbf{A} = \mathcal{I}:\mathbf{A}
   *        \qquad \text{and} \qquad
   *      \mathbf{A}^T = \overline{\mathcal{I}}:\mathbf{A} \, ,
   * @f]
   * or, in indicial notation,
   * @f[
   *   I_{ijkl} = \delta_{ik}\delta_{jl}
   *        \qquad \text{and} \qquad
   *     \overline I_{ijkl} = \delta_{il}\delta_{jk}
   * @f]
   * with the Kronecker deltas taking their common definition.
   * Note that $\mathcal{I} \neq \overline{\mathcal{I}}^T$.
   *
   * We then define the symmetric and skew-symmetric fourth-order unit tensors by
   * @f[
   * \mathcal{S} \dealcoloneq
   * \dfrac{1}{2}[\mathcal{I} + \overline{\mathcal{I}}]
   * \qquad \text{and} \qquad
   * \mathcal{W} \dealcoloneq
   * \dfrac{1}{2}[\mathcal{I} - \overline{\mathcal{I}}] \, ,
   * @f]
   * such that
   * @f[
   *      \mathcal{S}:\mathbf{A} = \dfrac{1}{2}[\mathbf{A} + \mathbf{A}^T]
   *    \qquad \text{and} \qquad
   *      \mathcal{W}:\mathbf{A} = \dfrac{1}{2}[\mathbf{A} - \mathbf{A}^T] \, .
   * @f]
   * The fourth-order symmetric tensor returned by identity_tensor() is
   * $\mathcal{S}$.
   *
   * @ingroup physics
   */
  namespace Elasticity
  {
  }

}

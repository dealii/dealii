// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

#ifndef dealii_particles_utilities
#define dealii_particles_utilities

#include <deal.II/base/config.h>

#include <deal.II/base/index_set.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/particles/particle_handler.h>


DEAL_II_NAMESPACE_OPEN


namespace Particles
{
  /**
   * A namespace for functions offering tools to handle ParticleHandlers and
   * their coupling with DoFHandlers
   */
  namespace Utilities
  {
    /**
     * Create an interpolation sparsity pattern for particles.
     *
     * Given a triangulation representing the domain $\Omega$, a particle
     * handler of particles in $\Omega$, and a scalar finite element space
     * $V(\Omega) = \text{span}\{v_j\}_{j=0}^n$, compute the sparsity pattern
     * that would be necessary to assemble the matrix
     * \f[
     * M_{i,j} \dealcoloneq v_j(x_i) ,
     * \f]
     * where $V(\Omega)$ is the finite element space associated with the
     * `space_dh`, and the index `i` is given by the particle id whose position
     * is `x_i`.
     *
     * In the case of vector valued finite element spaces, the components on
     * which interpolation must be performed can be selected using a component
     * mask. Only primitive finite element spaces are supported.
     *
     * When selecting more than one component, the resulting sparsity will have
     * dimension equal to `particle_handler.n_global_particles() *
     * mask.n_selected_components()` times `space_dh.n_dofs()`, and the
     * corresponding matrix entries are given by
     * \f[
     *  M_{(i*n_comps+k),j} \dealcoloneq v_j(x_i) \cdot e_{comp_j},
     * \f]
     * where `comp_j` is the only non zero component of the vector valued basis
     * function `v_j` (equal to `fe.system_to_component_index(j).first`), and
     * `k` corresponds to its index within the selected components of the mask.
     *
     * The `sparsity` is filled by locating the position of the particle with
     * index `i` within the particle handler with respect to the embedding
     * triangulation $\Omega$, and coupling it with all the local degrees of
     * freedom specified in the component mask @p space_comps, following the
     * ordering in which they are selected in the mask @p space_comps.
     *
     * If a particle does not fall within $\Omega$, it is ignored, and the
     * corresponding rows of the sparsity will be empty.
     *
     * Constraints of the form supported by the AffineConstraints class may be
     * supplied with the @p constraints argument. The method
     * AffineConstraints::add_entries_local_to_global() is used to fill the
     * final sparsity pattern.
     *
     * @author Bruno Blais, Luca Heltai, 2019
     */
    template <int dim,
              int spacedim,
              typename SparsityType,
              typename number = double>
    void
    create_interpolation_sparsity_pattern(
      const DoFHandler<dim, spacedim> &                space_dh,
      const Particles::ParticleHandler<dim, spacedim> &particle_handler,
      SparsityType &                                   sparsity,
      const AffineConstraints<number> &                constraints =
        AffineConstraints<number>(),
      const ComponentMask &space_comps = ComponentMask());

    /**
     * Create an interpolation matrix for particles.
     *
     * Given a triangulation representing the domains $\Omega$, a particle
     * handler of particles in $\Omega$, and a scalar finite element space
     * $V(\Omega) = \text{span}\{v_j\}_{j=0}^n$, compute the matrix
     * \f[
     * M_{ij} \dealcoloneq v_j(x_i) ,
     * \f]
     * where $V(\Omega)$ is the finite element space associated with the
     * `space_dh`, and the index `i` is given by the particle id whose position
     * is `x_i`.
     *
     * In the case of vector valued finite element spaces, the components on
     * which interpolation must be performed can be selected using a component
     * mask. Only primitive finite element spaces are supported.
     *
     * When selecting more than one component, the resulting sparsity will have
     * dimension equal to `particle_handler.n_global_particles() *
     * mask.n_selected_components()` times `space_dh.n_dofs()`, and the
     * corresponding matrix entries are given by
     * \f[
     *  M_{(i*n_comps+k),j} \dealcoloneq v_j(x_i) \cdot e_{comp_j},
     * \f]
     * where `comp_j` is the only non zero component of the vector valued basis
     * function `v_j` (equal to `fe.system_to_component_index(j).first`), and
     * `k` corresponds to its index within the selected components of the mask.
     *
     * The matrix is filled by locating the position of the particle with
     * index `i` within the particle handler with respect to the embedding
     * triangulation $\Omega$, and coupling it with all the local degrees of
     * freedom specified in the component mask @p space_comps, following the
     * ordering in which they are selected in the mask @p space_comps.
     *
     * If a particle does not fall within $\Omega$, it is ignored, and the
     * corresponding rows of the matrix will be zero.
     *
     * Constraints of the form supported by the AffineConstraints class may be
     * supplied with the @p constraints argument. The method
     * AffineConstraints::distribute_local_to_global() is used to distribute
     * the entries of the matrix to respect the given constraints.
     *
     * @author Bruno Blais, Luca Heltai, 2019
     */
    template <int dim, int spacedim, typename MatrixType>
    void
    create_interpolation_matrix(
      const DoFHandler<dim, spacedim> &                space_dh,
      const Particles::ParticleHandler<dim, spacedim> &particle_handler,
      MatrixType &                                     matrix,
      const AffineConstraints<typename MatrixType::value_type> &constraints =
        AffineConstraints<typename MatrixType::value_type>(),
      const ComponentMask &space_comps = ComponentMask());
  } // namespace Utilities
} // namespace Particles
DEAL_II_NAMESPACE_CLOSE

#endif

// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

#ifndef dealii_non_matching_coupling
#define dealii_non_matching_coupling

#include <deal.II/base/config.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/constraint_matrix.h>

DEAL_II_NAMESPACE_OPEN

/**
 * A namespace for functions offering tools to handle two meshes with no
 * alignment requirements, but where one of the meshes is embedded
 * inside the other in the real-space.
 *
 * Typically these functions allow for computations on the real-space
 * intersection between the two meshes e.g. surface integrals and
 * construction of mass matrices.
 */
namespace NonMatching
{
  /**
   * Create a coupling sparsity pattern for non-matching, overlapping grids.
   *
   * Given two non-matching triangulations, representing the domains $\Omega$
   * and $B$, with $B \subseteq \Omega$, and two finite element spaces
   * $V(\Omega) = \text{span}\{v_i\}_{i=0}^n$ and $Q(B) =
   * \text{span}\{w_j\}_{j=0}^m$, compute the sparsity pattern that would be
   * necessary to assemble the matrix
   * \f[
   * M_{ij} := \int_{B} v_i(x) w_j(x) dx, \quad i \in [0,n), j \in [0,m),
   * \f]
   * where $V(\Omega)$ is the finite element space associated with the
   * `space_dh` passed to this function (or part of it, if specified in
   * `space_comps`), while $Q(B)$ is the finite element space associated with
   * the `immersed_dh` passed to this function (or part of it, if specified in
   * `immersed_comps`).
   *
   * The `sparsity` is filled by locating the position of quadrature points
   * (obtained by the reference quadrature `quad`) defined on elements of $B$
   * with respect to the embedding triangulation $\Omega$. For each overlapping
   * cell, the entries corresponding to `space_comps` in `space_dh` and
   * `immersed_comps` in `immersed_dh` are added to the sparsity pattern.
   *
   * The `space_comps` and `immersed_comps` masks are assumed to be ordered in
   * the same way: the first component of `space_comps` will couple with the
   * first component of `immersed_comps`, the second with the second, and so
   * on. If one of the two masks has more non-zero than the other, then the
   * excess components will be ignored.
   *
   * If the domain $B$ does not fall within $\Omega$, an exception will be
   * thrown by the algorithm that computes the quadrature point locations. In
   * particular, notice that this function only makes sens for `dim1` lower or
   * equal than `dim0`. A static assert guards that this is actually the case.
   *
   * For both spaces, it is possible to specify a custom Mapping, which
   * defaults to StaticMappingQ1 for both.
   *
   * This function will also work in parallel, provided that the immersed
   * triangulation is of type parallel::shared::Triangulation<dim1,spacedim>.
   * An exception is thrown if you use an immersed
   * parallel::distributed::Triangulation<dim1,spacedim>.
   *
   * See the tutorial program step-60 for an example on how to use this
   * function.
   *
   * @author Luca Heltai, 2018
   */
  template <int dim0, int dim1, int spacedim, typename Sparsity>
  void
  create_coupling_sparsity_pattern(
    const DoFHandler<dim0, spacedim>& space_dh,
    const DoFHandler<dim1, spacedim>& immersed_dh,
    const Quadrature<dim1>&           quad,
    Sparsity&                         sparsity,
    const ConstraintMatrix&           constraints    = ConstraintMatrix(),
    const ComponentMask&              space_comps    = ComponentMask(),
    const ComponentMask&              immersed_comps = ComponentMask(),
    const Mapping<dim0, spacedim>&    space_mapping
    = StaticMappingQ1<dim0, spacedim>::mapping,
    const Mapping<dim1, spacedim>& immersed_mapping
    = StaticMappingQ1<dim1, spacedim>::mapping);

  /**
   * Same as above, but takes an additional GridTools::Cache object, instead of
   * creating one internally. In this version of the function, the parameter @p
   * space_mapping cannot be specified, since it is taken from the @p cache
   * parameter.
   *
   * @author Luca Heltai, 2018
   */
  template <int dim0, int dim1, int spacedim, typename Sparsity>
  void
  create_coupling_sparsity_pattern(
    const GridTools::Cache<dim0, spacedim>& cache,
    const DoFHandler<dim0, spacedim>&       space_dh,
    const DoFHandler<dim1, spacedim>&       immersed_dh,
    const Quadrature<dim1>&                 quad,
    Sparsity&                               sparsity,
    const ConstraintMatrix&                 constraints    = ConstraintMatrix(),
    const ComponentMask&                    space_comps    = ComponentMask(),
    const ComponentMask&                    immersed_comps = ComponentMask(),
    const Mapping<dim1, spacedim>&          immersed_mapping
    = StaticMappingQ1<dim1, spacedim>::mapping);

  /**
   * Create a coupling mass matrix for non-matching, overlapping grids.
   *
   * Given two non-matching triangulations, representing the domains $\Omega$
   * and $B$, with $B \subseteq \Omega$, and two finite element spaces
   * $V(\Omega) = \text{span}\{v_i\}_{i=0}^n$ and $Q(B) =
   * \text{span}\{w_j\}_{j=0}^m$, compute the coupling matrix
   * \f[
   * M_{ij} := \int_{B} v_i(x) w_j(x) dx, \quad i \in [0,n), j \in [0,m),
   * \f]
   * where $V(\Omega)$ is the finite element space associated with the
   * `space_dh` passed to this function (or part of it, if specified in
   * `space_comps`), while $Q(B)$ is the finite element space associated with
   * the `immersed_dh` passed to this function (or part of it, if specified in
   * `immersed_comps`).
   *
   * The corresponding sparsity patterns can be computed by calling the
   * make_coupling_sparsity_pattern function. The elements of the matrix are
   * computed by locating the position of quadrature points defined on elements
   * of $B$ with respect to the embedding triangulation $\Omega$.
   *
   * The `space_comps` and `immersed_comps` masks are assumed to be ordered in
   * the same way: the first component of `space_comps` will couple with the
   * first component of `immersed_comps`, the second with the second, and so
   * on. If one of the two masks has more non-zero entries non-zero than the
   * other, then the excess components will be ignored.
   *
   * If the domain $B$ does not fall within $\Omega$, an exception will be
   * thrown by the algorithm that computes the quadrature point locations. In
   * particular, notice that this function only makes sense for `dim1` lower or
   * equal than `dim0`. A static assert guards that this is actually the case.
   *
   * For both spaces, it is possible to specify a custom Mapping, which
   * defaults to StaticMappingQ1 for both.
   *
   * This function will also work in parallel, provided that the immersed
   * triangulation is of type parallel::shared::Triangulation<dim1,spacedim>.
   * An exception is thrown if you use an immersed
   * parallel::distributed::Triangulation<dim1,spacedim>.
   *
   * See the tutorial program step-60 for an example on how to use this
   * function.
   *
   * @author Luca Heltai, 2018
   */
  template <int dim0, int dim1, int spacedim, typename Matrix>
  void
  create_coupling_mass_matrix(
    const DoFHandler<dim0, spacedim>& space_dh,
    const DoFHandler<dim1, spacedim>& immersed_dh,
    const Quadrature<dim1>&           quad,
    Matrix&                           matrix,
    const ConstraintMatrix&           constraints    = ConstraintMatrix(),
    const ComponentMask&              space_comps    = ComponentMask(),
    const ComponentMask&              immersed_comps = ComponentMask(),
    const Mapping<dim0, spacedim>&    space_mapping
    = StaticMappingQ1<dim0, spacedim>::mapping,
    const Mapping<dim1, spacedim>& immersed_mapping
    = StaticMappingQ1<dim1, spacedim>::mapping);

  /**
   * Same as above, but takes an additional GridTools::Cache object, instead of
   * creating one internally. In this version of the function, the parameter @p
   * space_mapping cannot specified, since it is taken from the @p cache
   * parameter.
   *
   * @author Luca Heltai, 2018
   */
  template <int dim0, int dim1, int spacedim, typename Matrix>
  void
  create_coupling_mass_matrix(
    const GridTools::Cache<dim0, spacedim>& cache,
    const DoFHandler<dim0, spacedim>&       space_dh,
    const DoFHandler<dim1, spacedim>&       immersed_dh,
    const Quadrature<dim1>&                 quad,
    Matrix&                                 matrix,
    const ConstraintMatrix&                 constraints    = ConstraintMatrix(),
    const ComponentMask&                    space_comps    = ComponentMask(),
    const ComponentMask&                    immersed_comps = ComponentMask(),
    const Mapping<dim1, spacedim>&          immersed_mapping
    = StaticMappingQ1<dim1, spacedim>::mapping);
} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE

#endif

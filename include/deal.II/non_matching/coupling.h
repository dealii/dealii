// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2022 by the deal.II authors
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

#ifndef dealii_non_matching_coupling
#define dealii_non_matching_coupling

#include <deal.II/base/config.h>

#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/lac/affine_constraints.h>

#ifdef DEAL_II_WITH_CGAL
#  include <deal.II/cgal/utilities.h>
#endif

DEAL_II_NAMESPACE_OPEN

/**
 * A namespace for functions offering tools to handle two meshes with no
 * alignment requirements.
 *
 * Typically these functions allow for computations on the real-space
 * intersection between the two meshes e.g. surface integrals and
 * construction of coupling matrices.
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
   * M_{ij} \dealcoloneq \int_{B} v_i(x) w_j(x) dx,
   *                     \quad i \in [0,n), j \in [0,m),
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
   */
  template <int dim0,
            int dim1,
            int spacedim,
            typename Sparsity,
            typename number = double>
  void
  create_coupling_sparsity_pattern(
    const DoFHandler<dim0, spacedim> &space_dh,
    const DoFHandler<dim1, spacedim> &immersed_dh,
    const Quadrature<dim1> &          quad,
    Sparsity &                        sparsity,
    const AffineConstraints<number> & constraints = AffineConstraints<number>(),
    const ComponentMask &             space_comps = ComponentMask(),
    const ComponentMask &             immersed_comps = ComponentMask(),
    const Mapping<dim0, spacedim> &   space_mapping =
      StaticMappingQ1<dim0, spacedim>::mapping,
    const Mapping<dim1, spacedim> &immersed_mapping =
      StaticMappingQ1<dim1, spacedim>::mapping,
    const AffineConstraints<number> &immersed_constraints =
      AffineConstraints<number>());

  /**
   * Same as above, but takes an additional GridTools::Cache object, instead of
   * creating one internally. In this version of the function, the parameter @p
   * space_mapping cannot be specified, since it is taken from the @p cache
   * parameter.
   */
  template <int dim0,
            int dim1,
            int spacedim,
            typename Sparsity,
            typename number = double>
  void
  create_coupling_sparsity_pattern(
    const GridTools::Cache<dim0, spacedim> &cache,
    const DoFHandler<dim0, spacedim> &      space_dh,
    const DoFHandler<dim1, spacedim> &      immersed_dh,
    const Quadrature<dim1> &                quad,
    Sparsity &                              sparsity,
    const AffineConstraints<number> &constraints = AffineConstraints<number>(),
    const ComponentMask &            space_comps = ComponentMask(),
    const ComponentMask &            immersed_comps = ComponentMask(),
    const Mapping<dim1, spacedim> &  immersed_mapping =
      StaticMappingQ1<dim1, spacedim>::mapping,
    const AffineConstraints<number> &immersed_constraints =
      AffineConstraints<number>());


  /**
   * Create a coupling mass matrix for non-matching, overlapping grids.
   *
   * Given two non-matching triangulations, representing the domains $\Omega$
   * and $B$, with $B \subseteq \Omega$, and two finite element spaces
   * $V(\Omega) = \text{span}\{v_i\}_{i=0}^n$ and $Q(B) =
   * \text{span}\{w_j\}_{j=0}^m$, compute the coupling matrix
   * \f[
   * M_{ij} \dealcoloneq \int_{B} v_i(x) w_j(x) dx,
   *                     \quad i \in [0,n), j \in [0,m),
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
   */
  template <int dim0, int dim1, int spacedim, typename Matrix>
  void
  create_coupling_mass_matrix(
    const DoFHandler<dim0, spacedim> &                    space_dh,
    const DoFHandler<dim1, spacedim> &                    immersed_dh,
    const Quadrature<dim1> &                              quad,
    Matrix &                                              matrix,
    const AffineConstraints<typename Matrix::value_type> &constraints =
      AffineConstraints<typename Matrix::value_type>(),
    const ComponentMask &          space_comps    = ComponentMask(),
    const ComponentMask &          immersed_comps = ComponentMask(),
    const Mapping<dim0, spacedim> &space_mapping =
      StaticMappingQ1<dim0, spacedim>::mapping,
    const Mapping<dim1, spacedim> &immersed_mapping =
      StaticMappingQ1<dim1, spacedim>::mapping,
    const AffineConstraints<typename Matrix::value_type> &immersed_constraints =
      AffineConstraints<typename Matrix::value_type>());

  /**
   * Same as above, but takes an additional GridTools::Cache object, instead of
   * creating one internally. In this version of the function, the parameter @p
   * space_mapping cannot specified, since it is taken from the @p cache
   * parameter.
   */
  template <int dim0, int dim1, int spacedim, typename Matrix>
  void
  create_coupling_mass_matrix(
    const GridTools::Cache<dim0, spacedim> &              cache,
    const DoFHandler<dim0, spacedim> &                    space_dh,
    const DoFHandler<dim1, spacedim> &                    immersed_dh,
    const Quadrature<dim1> &                              quad,
    Matrix &                                              matrix,
    const AffineConstraints<typename Matrix::value_type> &constraints =
      AffineConstraints<typename Matrix::value_type>(),
    const ComponentMask &          space_comps    = ComponentMask(),
    const ComponentMask &          immersed_comps = ComponentMask(),
    const Mapping<dim1, spacedim> &immersed_mapping =
      StaticMappingQ1<dim1, spacedim>::mapping,
    const AffineConstraints<typename Matrix::value_type> &immersed_constraints =
      AffineConstraints<typename Matrix::value_type>());

  /**
   * Create a coupling sparsity pattern for non-matching independent grids,
   * using a convolution kernel with compact support of radius epsilon.
   *
   * Given two non-matching triangulations, representing the domains $\Omega^0$
   * and $\Omega^1$, both embedded in $\mathbb{R}^d$, and two finite element
   * spaces $V^0(\Omega^0) = \text{span}\{v_i\}_{i=0}^n$ and $V^1(\Omega^1) =
   * \text{span}\{w_\alpha\}_{\alpha=0}^m$, compute the sparsity pattern that
   * would be necessary to assemble the matrix
   *
   * \f[
   * M_{i\alpha} \dealcoloneq \int_{\Omega^0} \int_{\Omega^1}
   * v_i(x) K^{\epsilon}(x-y) w_\alpha(y) dx \ dy,
   * \quad i \in [0,n), \alpha \in [0,m),
   * \f]
   *
   * where $V^0(\Omega^0)$ is the finite element space associated with the
   * @p dh0 passed to this function (or part of it, if specified in
   * @p comps0), while $V^1(\Omega^1)$ is the finite element space associated
   * with the @p dh1 passed to this function (or part of it, if specified
   * in @p comps1), and $K^\epsilon$ is a function derived from
   * CutOffFunctionBase with compact support included in a ball of radius
   * $\epsilon$.
   *
   * The @p comps0 and @p comps1 masks are assumed to be ordered in
   * the same way: the first component of @p comps0 will couple with the
   * first component of @p comps1, the second with the second, and so
   * on. If one of the two masks has more active components than the other, then
   * the excess components will be ignored.
   *
   * For both spaces, it is possible to specify a custom Mapping, which
   * defaults to StaticMappingQ1 for both.
   *
   * This function will also work in parallel, provided that at least one of the
   * triangulations is of type parallel::shared::Triangulation<dim1,spacedim>.
   * An exception is thrown if both triagnulations are of type
   * parallel::distributed::Triangulation<dim1,spacedim>.
   *
   * This function assumes that the convolution has support contained in a box
   * of radius @p epsilon. If epsilon is set to zero, then we assume that the
   * kernel is the Dirac delta distribution, and the call is forwarded to the
   * method in this namespace with the same name, that does not take an epsilon
   * as input (but a quadrature formula @p quad is required). In this case, more
   * restrictive conditions are required on the two spaces. See the
   * documentation of the other create_coupling_sparsity_pattern() function.
   */
  template <int dim0,
            int dim1,
            int spacedim,
            typename Sparsity,
            typename Number = double>
  void
  create_coupling_sparsity_pattern(
    const double &                          epsilon,
    const GridTools::Cache<dim0, spacedim> &cache0,
    const GridTools::Cache<dim1, spacedim> &cache1,
    const DoFHandler<dim0, spacedim> &      dh0,
    const DoFHandler<dim1, spacedim> &      dh1,
    const Quadrature<dim1> &                quad,
    Sparsity &                              sparsity,
    const AffineConstraints<Number> &constraints0 = AffineConstraints<Number>(),
    const ComponentMask &            comps0       = ComponentMask(),
    const ComponentMask &            comps1       = ComponentMask());

  /**
   * Create a coupling mass matrix for non-matching independent grids,
   * using a convolution kernel with compact support.
   *
   * Given two non-matching triangulations, representing the domains
   * $\Omega^0$ and $\Omega^1$, both embedded in $\mathbb{R}^d$, and two finite
   * element spaces $V^0(\Omega^0) = \text{span}\{v_i\}_{i=0}^n$ and
   * $V^1(\Omega^1) = \text{span}\{w_\alpha\}_{\alpha=0}^m$, compute the matrix
   *
   * \f[
   * M_{i\alpha} \dealcoloneq \int_{\Omega^0} \int_{\Omega^1}
   * v_i(x) K^{\epsilon}(x-y) w_\alpha(y) dx \ dy,
   * \quad i \in [0,n), \alpha \in [0,m),
   * \f]
   *
   * where $V^0(\Omega^0)$ is the finite element space associated with the
   * @p dh0 passed to this function (or part of it, if specified in
   * @p comps0), while $V^1(\Omega^1)$ is the finite element space associated
   * with the @p dh1 passed to this function (or part of it, if specified
   * in @p comps1), and $K^\epsilon$ is a function derived from
   * CutOffFunctionBase with compact support included in a ball of radius
   * $\epsilon$.
   *
   * The corresponding sparsity patterns can be computed by calling the
   * make_coupling_sparsity_pattern() function.
   *
   * The @p comps0 and @p comps1 masks are assumed to be ordered in
   * the same way: the first component of @p comps0 will couple with the
   * first component of @p comps1, the second with the second, and so
   * on. If one of the two masks has more active components than the other, then
   * the excess components will be ignored.
   *
   * For both spaces, it is possible to specify a custom Mapping, which
   * defaults to StaticMappingQ1 for both.
   *
   * This function will also work in parallel, provided that one of the two
   * triangulations is of type parallel::shared::Triangulation<dim1,spacedim>.
   * An exception is thrown if both triangulations are of type
   * parallel::distributed::Triangulation<dim1,spacedim>.
   *
   * The parameter @p epsilon is used to set the size of the cut-off function
   * used to compute the convolution. If epsilon is set to zero, then we assume
   * that the kernel is the Dirac delta distribution, and the call is forwarded
   * to the method in this namespace with the same name, that does not take an
   * epsilon as input.
   */
  template <int dim0, int dim1, int spacedim, typename Matrix>
  void
  create_coupling_mass_matrix(
    Functions::CutOffFunctionBase<spacedim> &             kernel,
    const double &                                        epsilon,
    const GridTools::Cache<dim0, spacedim> &              cache0,
    const GridTools::Cache<dim1, spacedim> &              cache1,
    const DoFHandler<dim0, spacedim> &                    dh0,
    const DoFHandler<dim1, spacedim> &                    dh1,
    const Quadrature<dim0> &                              quadrature0,
    const Quadrature<dim1> &                              quadrature1,
    Matrix &                                              matrix,
    const AffineConstraints<typename Matrix::value_type> &constraints0 =
      AffineConstraints<typename Matrix::value_type>(),
    const ComponentMask &comps0 = ComponentMask(),
    const ComponentMask &comps1 = ComponentMask());

#ifdef DEAL_II_WITH_CGAL
  /**
   * Given two cached, arbitrarily overlapped grids, the following function
   * computes Quadrature rules on the intersection of the embedding grid with
   * the embedded one, of degree `degree`. The return type is a vector of tuples
   * `v`, where `v[i][0]` is an iterator to a cell of the embedding grid,
   * `v[i][1]` is an iterator to a cell of the embedded grid,
   * `v[i][2]` is a Quadrature formula to integrate over the intersection of the
   * two.
   *
   * The last parameter `tol` defaults to 1e-6, and can be used to discard small
   * intersections.
   *
   * @note This function calls compute_quadrature_on_intersection().
   *
   * @param [in] space_cache First cached triangulation.
   * @param [in] immersed_cache Second cached triangulation.
   * @param [in] degree The degree of accuracy of each quadrature formula.
   * @param [in] tol Tolerance used to discard small intersections.
   * @return std::vector<std::tuple<typename Triangulation<dim0, spacedim>::cell_iterator,
   * typename Triangulation<dim1, spacedim>::cell_iterator,
   * Quadrature<spacedim>>>.
   */
  template <int dim0, int dim1, int spacedim>
  std::vector<std::tuple<typename Triangulation<dim0, spacedim>::cell_iterator,
                         typename Triangulation<dim1, spacedim>::cell_iterator,
                         Quadrature<spacedim>>>
  compute_intersection(const GridTools::Cache<dim0, spacedim> &space_cache,
                       const GridTools::Cache<dim1, spacedim> &immersed_cache,
                       const unsigned int                      degree,
                       const double                            tol = 1e-6);

  /**
   * Create a coupling sparsity pattern for non-matching, overlapping
   * grids in an "exact" way. Here exact refers to that fact, differently
   * with respect to what is done in create_coupling_sparsity_pattern(), the
   * sparsity pattern is filled by looking at the intersection of the cells, as
   * it means that DoFs there are coupled.
   *
   * Those intersections are encoded inside @p intersection_info, which needs
   * to be computed right before calling this function.
   *
   * @param [in] intersections_info A vector of tuples where the i-th entry
   * contains two `active_cell_iterator`s to the intersected cells.
   * @param [in] space_dh `DoFHandler` object for the space grid.
   * @param [in] immersed_dh `DoFHandler` object for the embedded grid.
   * @param [out] sparsity The sparsity pattern to be filled.
   * @param [in] constraints `AffineConstraints` for the space grid.
   * @param [in] space_comps Mask for the space components of the finite
   * element.
   * @param [in] immersed_comps Mask for the embedded components of the finite
   * element.
   * @param [in] immersed_constraints `AffineConstraints` for the embedded grid.
   */
  template <int dim0,
            int dim1,
            int spacedim,
            typename Sparsity,
            typename number = double>
  void
  create_coupling_sparsity_pattern_with_exact_intersections(
    const std::vector<
      std::tuple<typename dealii::Triangulation<dim0, spacedim>::cell_iterator,
                 typename dealii::Triangulation<dim1, spacedim>::cell_iterator,
                 dealii::Quadrature<spacedim>>> &intersections_info,
    const DoFHandler<dim0, spacedim> &           space_dh,
    const DoFHandler<dim1, spacedim> &           immersed_dh,
    Sparsity &                                   sparsity,
    const AffineConstraints<number> &            constraints,
    const ComponentMask &                        space_comps,
    const ComponentMask &                        immersed_comps,
    const AffineConstraints<number> &            immersed_constraints);

  /**
   * Create the coupling mass matrix for non-matching, overlapping grids
   * in an "exact" way, i.e. by computing the local contributions
   *  \f[
   *  M_{ij} \dealcoloneq \int_{B} v_i(x) w_j(x) dx,
   *                      \quad i \in [0,n), j \in [0,m),
   *  \f]
   *   as products of cellwise smooth
   *  functions on the intersection of cells between the two grids. This
   *   information
   *  is encoded in @p intersections_info, a vector whose each element is a tuple containing the
   *  two intersected cells and a Quadrature formula on their intersection.
   *
   * @param [in] space_dh `DoFHandler` object for the space grid.
   * @param [in] immersed_dh `DoFHandler` object for the embedded grid.
   * @param [in] intersections_info std::vector<std::tuple> with Quadratures and cell_iterators to the intersected cells.
   * @param [in] matrix A reference to the mass matrix that must be filled.
   * @param [in] constraints `AffineConstraints` for the space grid.
   * @param [in] space_comps Mask for the space components of the finite element.
   * @param [in] immersed_comps Mask for the embedded components of the finite element.
   * @param [in] mapping0 Mapping object describing the space grid.
   * @param [in] mapping1 Mapping object describing the embedded grid.
   * @param [in] immersed_constraints `AffineConstraints` for the embedded grid.
   */
  template <int dim0, int dim1, int spacedim, typename Matrix>
  void
  create_coupling_mass_matrix_with_exact_intersections(
    const dealii::DoFHandler<dim0, spacedim> &space_dh,
    const dealii::DoFHandler<dim1, spacedim> &immersed_dh,
    const std::vector<
      std::tuple<typename dealii::Triangulation<dim0, spacedim>::cell_iterator,
                 typename dealii::Triangulation<dim1, spacedim>::cell_iterator,
                 dealii::Quadrature<spacedim>>> &intersections_info,
    Matrix &                                     matrix,
    const dealii::AffineConstraints<typename Matrix::value_type> &constraints,
    const dealii::ComponentMask &                                 space_comps,
    const dealii::ComponentMask &          immersed_comps,
    const dealii::Mapping<dim0, spacedim> &mapping0,
    const dealii::Mapping<dim1, spacedim> &mapping1,
    const dealii::AffineConstraints<typename Matrix::value_type>
      &immersed_constraints);
#endif
} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE

#endif

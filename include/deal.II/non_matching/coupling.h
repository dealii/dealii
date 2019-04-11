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
   *
   * @author Luca Heltai, 2018
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
      StaticMappingQ1<dim1, spacedim>::mapping);

  /**
   * Same as above, but takes an additional GridTools::Cache object, instead of
   * creating one internally. In this version of the function, the parameter @p
   * space_mapping cannot be specified, since it is taken from the @p cache
   * parameter.
   *
   * @author Luca Heltai, 2018
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
      StaticMappingQ1<dim1, spacedim>::mapping);


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
   *
   * @author Luca Heltai, 2018
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
      StaticMappingQ1<dim1, spacedim>::mapping);

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
      StaticMappingQ1<dim1, spacedim>::mapping);

  /**
   * Create a coupling sparsity pattern for non-matching independent grids,
   * using a convolution kernel with compact support of radius epsilon.
   *
   * Given two non-matching triangulations, representing the domains $\Omega^0$
   * and $\Omega^1$, both embedded in $R^d$, and two finite element spaces
   * $V^0(\Omega^0) = \text{span}\{v_i\}_{i=0}^n$ and $V^1(\Omega^1) =
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
   * `dh0` passed to this function (or part of it, if specified in
   * `comps0`), while $V^1(\Omega^1)$ is the finite element space associated
   * with the `dh1` passed to this function (or part of it, if specified
   * in `comps1`), and $K^\epsilon$ is a function with compact support included,
   * in a ball of radius $\epsilon$, derived from CutOffFunctionBase.
   *
   * The `comps0` and `comps1` masks are assumed to be ordered in
   * the same way: the first component of `comps0` will couple with the
   * first component of `comps1`, the second with the second, and so
   * on. If one of the two masks has more non-zero than the other, then the
   * excess components will be ignored.
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
   * as input.
   *
   * @author Luca Heltai, 2019.
   */
  template <int dim0, int dim1, int spacedim, typename Sparsity>
  void
  create_coupling_sparsity_pattern(
    const double &                          epsilon,
    const GridTools::Cache<dim0, spacedim> &cache0,
    const GridTools::Cache<dim1, spacedim> &cache1,
    const DoFHandler<dim0, spacedim> &      dh0,
    const DoFHandler<dim1, spacedim> &      dh1,
    Sparsity &                              sparsity,
    const AffineConstraints<double> &constraints0 = AffineConstraints<double>(),
    const ComponentMask &            comps0       = ComponentMask(),
    const ComponentMask &            comps1       = ComponentMask())
  {
    if (epsilon == 0.0)
      {
        create_coupling_sparsity_pattern(cache0,
                                         dh0,
                                         dh1,
                                         QGauss<dim1>(2 * dh1.get_fe().degree +
                                                      1),
                                         sparsity,
                                         constraints0,
                                         comps0,
                                         comps1,
                                         cache1.get_mapping());
        return;
      }
    AssertDimension(sparsity.n_rows(), dh0.n_dofs());
    AssertDimension(sparsity.n_cols(), dh1.n_dofs());
    Assert((dynamic_cast<
              const parallel::distributed::Triangulation<dim1, spacedim> *>(
              &dh1.get_triangulation()) == nullptr),
           ExcNotImplemented());

    const auto &fe0 = dh0.get_fe();
    const auto &fe1 = dh1.get_fe();

    // Dof indices
    std::vector<types::global_dof_index> dofs0(fe0.dofs_per_cell);
    std::vector<types::global_dof_index> dofs1(fe1.dofs_per_cell);

    // Take care of components
    const ComponentMask mask0 =
      (comps0.size() == 0 ? ComponentMask(fe0.n_components(), true) : comps0);

    const ComponentMask mask1 =
      (comps1.size() == 0 ? ComponentMask(fe1.n_components(), true) : comps1);

    AssertDimension(mask0.size(), fe0.n_components());
    AssertDimension(mask1.size(), fe1.n_components());

    // Global to local indices
    std::vector<unsigned int> gtl0(fe0.n_components(),
                                   numbers::invalid_unsigned_int);
    std::vector<unsigned int> gtl1(fe1.n_components(),
                                   numbers::invalid_unsigned_int);

    for (unsigned int i = 0, j = 0; i < gtl0.size(); ++i)
      if (mask0[i])
        gtl0[i] = j++;

    for (unsigned int i = 0, j = 0; i < gtl1.size(); ++i)
      if (mask1[i])
        gtl1[i] = j++;

    const auto &tree1 = cache1.get_cell_bounding_boxes_rtree();

    std::vector<
      std::pair<BoundingBox<spacedim>,
                typename Triangulation<dim1, spacedim>::active_cell_iterator>>
      intersection;

    for (const auto &cell0 : dh0.active_cell_iterators())
      if (cell0->is_locally_owned())
        {
          intersection.resize(0);
          auto box0 = cell0->bounding_box();
          box0.extend(epsilon);
          boost::geometry::index::query(tree1,
                                        boost::geometry::index::intersects(
                                          box0),
                                        std::back_inserter(intersection));
          if (!intersection.empty())
            {
              cell0->get_dof_indices(dofs0);
              for (const auto &entry : intersection)
                {
                  typename DoFHandler<dim1, spacedim>::cell_iterator cell1(
                    *entry.second, &dh1);
                  cell1->get_dof_indices(dofs1);
                  constraints0.add_entries_local_to_global(dofs0,
                                                           dofs1,
                                                           sparsity);
                }
            }
        }
  }


  template <class FE0, class FE1>
  std::pair<std::vector<unsigned int>, std::vector<unsigned int>>
  compute_components_coupling(const ComponentMask &comps0,
                              const ComponentMask &comps1,
                              const FE0 &          fe0,
                              const FE1 &          fe1)
  {
    // Take care of components
    const ComponentMask mask0 =
      (comps0.size() == 0 ? ComponentMask(fe0.n_components(), true) : comps0);

    const ComponentMask mask1 =
      (comps1.size() == 0 ? ComponentMask(fe1.n_components(), true) : comps1);

    AssertDimension(mask0.size(), fe0.n_components());
    AssertDimension(mask1.size(), fe1.n_components());

    // Global to local indices
    std::vector<unsigned int> gtl0(fe0.n_components(),
                                   numbers::invalid_unsigned_int);
    std::vector<unsigned int> gtl1(fe1.n_components(),
                                   numbers::invalid_unsigned_int);

    for (unsigned int i = 0, j = 0; i < gtl0.size(); ++i)
      if (mask0[i])
        gtl0[i] = j++;

    for (unsigned int i = 0, j = 0; i < gtl1.size(); ++i)
      if (mask1[i])
        gtl1[i] = j++;
    return {gtl0, gtl1};
  }


  /**
   * Create a coupling  mass matrix for non-matching independent grids,
   * using a convolution kernel with compact support.
   *
   * Given two non-matching triangulations, representing the domains
   * $\Omega^0$ and $\Omega^1$, both embedded in $R^d$, and two finite element
   * spaces $V^0(\Omega^0) = \text{span}\{v_i\}_{i=0}^n$ and $V^1(\Omega^1) =
   * \text{span}\{w_\alpha\}_{\alpha=0}^m$, compute the matrix
   *
   * \f[
   * M_{i\alpha} \dealcoloneq \int_{\Omega^0} \int_{\Omega^1}
   * v_i(x) K^{\epsilon}(x-y) w_\alpha(y) dx \ dy,
   * \quad i \in [0,n), \alpha \in [0,m),
   * \f]
   *
   * where $V^0(\Omega^0)$ is the finite element space associated with the
   * `dh0` passed to this function (or part of it, if specified in
   * `comps0`), while $V^1(\Omega^1)$ is the finite element space associated
   * with the `dh1` passed to this function (or part of it, if specified
   * in `comps1`), and $K^\epsilon$ is a function with compact support included,
   * in a ball of radius $\epsilon$, derived from CutOffFunctionBase.
   *
   * The corresponding sparsity patterns can be computed by calling the
   * make_coupling_sparsity_pattern() function.
   *
   * The `comps0` and `comps1` masks are assumed to be ordered in
   * the same way: the first component of `comps0` will couple with the
   * first component of `comps1`, the second with the second, and so
   * on. If one of the two masks has more non-zero than the other, then the
   * excess components will be ignored.
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
   *
   * @author Luca Heltai, 2019.
   */
  template <int dim0, int dim1, int spacedim, typename Matrix>
  void
  create_coupling_mass_matrix(
    Functions::CutOffFunctionBase<spacedim> &kernel,
    const double &                           epsilon,
    const GridTools::Cache<dim0, spacedim> & cache0,
    const GridTools::Cache<dim1, spacedim> & cache1,
    const DoFHandler<dim0, spacedim> &       dh0,
    const DoFHandler<dim1, spacedim> &       dh1,
    const Quadrature<dim0> &                 quadrature0,
    const Quadrature<dim1> &                 quadrature1,
    Matrix &                                 matrix,
    const AffineConstraints<double> &constraints0 = AffineConstraints<double>(),
    const ComponentMask &            comps0       = ComponentMask(),
    const ComponentMask &            comps1       = ComponentMask())
  {
    if (epsilon == 0)
      {
        create_coupling_mass_matrix(cache0,
                                    dh0,
                                    dh1,
                                    quadrature1,
                                    matrix,
                                    constraints0,
                                    comps0,
                                    comps1,
                                    cache1.get_mapping());
        return;
      }

    AssertDimension(matrix.m(), dh0.n_dofs());
    AssertDimension(matrix.n(), dh1.n_dofs());

    const bool zero_is_distributed =
      (dynamic_cast<const parallel::distributed::Triangulation<dim1, spacedim>
                      *>(&dh0.get_triangulation()) == nullptr);
    const bool one_is_distributed =
      (dynamic_cast<const parallel::distributed::Triangulation<dim1, spacedim>
                      *>(&dh1.get_triangulation()) == nullptr);

    // We bail out if both are distributed triangulations
    Assert(zero_is_distributed && one_is_distributed, ExcNotImplemented());

    // If we can loop on both, we decide where to make the outer loop according
    // to the size of triangulation. The reasoning is the following:
    // - Access to the tree: log(N)
    // - We compute intersection for each of the outer loop cells (M)
    // Total cost (besides the setup) is: M log(N)
    // If we can, make sure M is the smallest number
    const bool outer_loop_on_zero =
      (zero_is_distributed && !one_is_distributed) ||
      (dh1.get_triangulation().n_active_cells() >
       dh0.get_triangulation().n_active_cells());

    const auto &fe0 = dh0.get_fe();
    const auto &fe1 = dh1.get_fe();

    FEValues<dim0, spacedim> fev0(cache0.get_mapping(),
                                  fe0,
                                  quadrature0,
                                  update_values | update_JxW_values |
                                    update_quadrature_points);

    FEValues<dim1, spacedim> fev1(cache1.get_mapping(),
                                  fe1,
                                  quadrature1,
                                  update_values | update_JxW_values |
                                    update_quadrature_points);

    // Dof indices
    std::vector<types::global_dof_index> dofs0(fe0.dofs_per_cell);
    std::vector<types::global_dof_index> dofs1(fe1.dofs_per_cell);

    // Local Matrix
    FullMatrix<double> cell_matrix(fe0.dofs_per_cell, fe1.dofs_per_cell);

    // Global to local indices
    auto        p    = compute_components_coupling(comps0, comps1, fe0, fe1);
    const auto &gtl0 = p.first;
    const auto &gtl1 = p.second;

    kernel.set_radius(epsilon);
    std::vector<double> kernel_values(quadrature1.size());

    auto assemble_one_pair = [&]() {
      cell_matrix = 0;
      for (unsigned int q0 = 0; q0 < quadrature0.size(); ++q0)
        {
          kernel.set_center(fev0.quadrature_point(q0));
          kernel.value_list(fev1.get_quadrature_points(), kernel_values);
          for (unsigned int q1 = 0; q1 < quadrature1.size(); ++q1)
            for (unsigned int i = 0; i < fe0.dofs_per_cell; ++i)
              {
                const auto comp_i = fe0.system_to_component_index(i).first;
                if (gtl0[comp_i] != numbers::invalid_unsigned_int)
                  for (unsigned int j = 0; j < fe1.dofs_per_cell; ++j)
                    {
                      const auto comp_j =
                        fe1.system_to_component_index(j).first;
                      if (gtl1[comp_j] == gtl0[comp_i])
                        {
                          cell_matrix(i, j) +=
                            fev0.shape_value(i, q0) * fev1.shape_value(j, q1) *
                            kernel_values[q1] * fev0.JxW(q0) * fev1.JxW(q1);
                        }
                    }
              }
        }

      constraints0.distribute_local_to_global(cell_matrix,
                                              dofs0,
                                              dofs1,
                                              matrix);
    };

    if (outer_loop_on_zero)
      {
        Assert(one_is_distributed == false, ExcInternalError());

        const auto &tree1 = cache1.get_cell_bounding_boxes_rtree();

        std::vector<std::pair<
          BoundingBox<spacedim>,
          typename Triangulation<dim1, spacedim>::active_cell_iterator>>
          intersection;

        for (const auto &cell0 : dh0.active_cell_iterators())
          if (cell0->is_locally_owned())
            {
              intersection.resize(0);
              auto box0 = cell0->bounding_box();
              box0.extend(epsilon);
              boost::geometry::index::query(tree1,
                                            boost::geometry::index::intersects(
                                              box0),
                                            std::back_inserter(intersection));
              if (!intersection.empty())
                {
                  cell0->get_dof_indices(dofs0);
                  fev0.reinit(cell0);
                  for (const auto &entry : intersection)
                    {
                      typename DoFHandler<dim1, spacedim>::cell_iterator cell1(
                        *entry.second, &dh1);
                      cell1->get_dof_indices(dofs1);
                      fev1.reinit(cell1);
                      assemble_one_pair();
                    }
                }
            }
      }
    else
      {
        Assert(zero_is_distributed == false, ExcInternalError());
        const auto &tree0 = cache0.get_cell_bounding_boxes_rtree();

        std::vector<std::pair<
          BoundingBox<spacedim>,
          typename Triangulation<dim0, spacedim>::active_cell_iterator>>
          intersection;

        for (const auto &cell1 : dh1.active_cell_iterators())
          if (cell1->is_locally_owned())
            {
              intersection.resize(0);
              auto box1 = cell1->bounding_box();
              box1.extend(epsilon);
              boost::geometry::index::query(tree0,
                                            boost::geometry::index::intersects(
                                              box1),
                                            std::back_inserter(intersection));
              if (!intersection.empty())
                {
                  cell1->get_dof_indices(dofs1);
                  fev1.reinit(cell1);
                  for (const auto &entry : intersection)
                    {
                      typename DoFHandler<dim0, spacedim>::cell_iterator cell0(
                        *entry.second, &dh0);
                      cell0->get_dof_indices(dofs0);
                      fev0.reinit(cell0);
                      assemble_one_pair();
                    }
                }
            }
      }
  }
} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE

#endif

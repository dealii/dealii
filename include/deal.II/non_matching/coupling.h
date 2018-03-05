// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2017 by the deal.II authors
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
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>



#include <vector>

DEAL_II_NAMESPACE_OPEN
namespace NonMatching
{
  /**
   * Create a coupling sparsity pattern for non-matching, overlapping grids.
   *
   * Given two non-matching triangulations, representing the domains $\Omega$
   * and $B$, with $B \subseteq \Omega$, and two finite element spaces
   * $V(\Omega)$ and $Q(B)$, compute the sparsity pattern that would be
   * necessary to assemble the matrix
   *
   * \f[
   * M_{ij} := \int_{B} v_i(x) w_j(x) dx, \quad \foreach v_i \in V(\Omega), w_j \in Q(B)
   * \f]
   *
   * where $V(\Omega)$ is the finite element space associated with `space_dh`
   * (or part of it, if specified in `space_comps`), while $Q(B)$ is the finite
   * element space associated with `immersed_dh` (or part of it, if specified
   * in a `immersed_comps`).
   *
   * The `sparsity` is constructed by locating the position of quadrature
   * points (obtained by the reference quadrature `quad`) defined on elements
   * of $B$ with respect to the embedding triangulation $\Omega$. For each
   * overlapping cell, the entries corresponding to `space_comps` in `space_dh`
   * and `immersed_comps` in `immersed_dh` are added to the sparsity pattern.
   *
   * The `space_comps` and `immersed_comps` masks are assumed to be ordered in
   * the same way: the first component of `space_comps` will couple with the
   * first component of `immersed_comps`, the second with the second, and so
   * on. If one of the two masks has more non-zero entries w.r.t. the other,
   * then the excess components will be ignored.
   *
   * If the domain $B$ does not fall withing $\Omega$, an exception will be
   * thrown by the algorithm that computes the quadrature point locations.
   *
   * For both spaces, it is possible to specify a custom Mapping, which
   * defaults to StaticMappingQ1 for both.
   *
   * @author Luca Heltai, 2018
   */
  template<int dim0, int dim1, int spacedim, typename Sparsity>
  void create_coupling_sparsity_pattern(const DoFHandler<dim0, spacedim> &space_dh,
                                        const DoFHandler<dim1, spacedim> &immersed_dh,
                                        const Quadrature<dim1>           &quad,
                                        Sparsity                         &sparsity,
                                        const ComponentMask              &space_comps = ComponentMask(),
                                        const ComponentMask              &immersed_comps = ComponentMask(),
                                        const Mapping<dim0, spacedim>    &space_mapping = StaticMappingQ1<dim0,spacedim>::mapping,
                                        const Mapping<dim1, spacedim>    &immersed_mapping = StaticMappingQ1<dim1, spacedim>::mapping);

  /**
   * Create a coupling mass matrix for non-matching, overlapping grids.
   *
   * Given two non-matching triangulations, representing the domains $\Omega$
   * and $B$, with $B \subseteq \Omega$, and two finite element spaces
   * $V(\Omega)$ and $Q(B)$, compute the matrix
   *
   * \f[
   * M_{ij} := \int_{B} v_i(x) w_j(x) dx, \quad \foreach v_i \in V(\Omega), w_j \in Q(B)
   * \f]
   *
   * where $V(\Omega)$ is the finite element space associated with `space_dh`
   * (or part of it, if specified in `space_comps`), while $Q(B)$ is the finite
   * element space associated with `immersed_dh` (or part of it, if specified
   * in a `immersed_comps`).
   *
   * The corresponding `sparsity` can be computed by calling the
   * make_coupling_sparsity_pattern function. The elements of the matrix are
   * computed by locating the position of quadrature points defined on elements
   * of $B$ with respect to the embedding triangulation $\Omega$.
   *
   * The `space_comps` and `immersed_comps` masks are assumed to be ordered in
   * the same way: the first component of `space_comps` will couple with the
   * first component of `immersed_comps`, the second with the second, and so
   * on. If one of the two masks has more non-zero entries w.r.t. the other,
   * then the excess components will be ignored.
   *
   * If the domain $B$ does not fall withing $\Omega$, an exception will be
   * thrown by the algorithm that computes the quadrature point locations.
   *
   * For both spaces, it is possible to specify a custom Mapping, which
   * defaults to StaticMappingQ1 for both.
   *
   * @author Luca Heltai, 2018
   */
  template<int dim0, int dim1, int spacedim, typename Matrix>
  void create_coupling_mass_matrix(const DoFHandler<dim0, spacedim> &space_dh,
                                   const DoFHandler<dim1, spacedim> &immersed_dh,
                                   const Quadrature<dim1>           &quad,
                                   Matrix                           &matrix,
                                   const ConstraintMatrix           &constraints = ConstraintMatrix(),
                                   const ComponentMask              &space_comps = ComponentMask(),
                                   const ComponentMask              &immersed_comps = ComponentMask(),
                                   const Mapping<dim0, spacedim>    &space_mapping = StaticMappingQ1<dim0,spacedim>::mapping,
                                   const Mapping<dim1, spacedim>    &immersed_mapping = StaticMappingQ1<dim1, spacedim>::mapping);




// === inline and template functions ===



  template<int dim0, int dim1, int spacedim, typename Sparsity>
  void create_coupling_sparsity_pattern(const DoFHandler<dim0, spacedim> &space_dh,
                                        const DoFHandler<dim1, spacedim> &immersed_dh,
                                        const Quadrature<dim1>           &quad,
                                        Sparsity                         &sparsity,
                                        const ComponentMask              &space_comps,
                                        const ComponentMask              &immersed_comps,
                                        const Mapping<dim0, spacedim>    &space_mapping,
                                        const Mapping<dim1, spacedim>    &immersed_mapping)
  {
    AssertDimension(sparsity.n_rows(), space_dh.n_dofs());
    AssertDimension(sparsity.n_cols(), immersed_dh.n_dofs());

    const auto &space_fe = space_dh.get_fe();
    const auto &immersed_fe = immersed_dh.get_fe();

    // Now we run on ech cell, get a quadrature formula
    typename DoFHandler<dim1,spacedim>::active_cell_iterator
    cell = immersed_dh.begin_active(),
    endc = immersed_dh.end();

    // Dof indices
    std::vector<types::global_dof_index> dofs(immersed_fe.dofs_per_cell);
    std::vector<types::global_dof_index> odofs(space_fe.dofs_per_cell);

    FEValues<dim1,spacedim> fe_v(immersed_mapping, immersed_fe, quad,
                                 update_quadrature_points);

    GridTools::Cache<dim0,spacedim> cache(space_dh.get_triangulation(), space_mapping);

    // Take care of components
    ComponentMask space_c = space_comps;
    ComponentMask immersed_c = immersed_comps;
    if (space_c.size() == 0)
      space_c = ComponentMask(space_fe.n_components(), true);
    if (immersed_c.size() == 0)
      immersed_c = ComponentMask(immersed_fe.n_components(), true);

    AssertDimension(space_c.size(), space_fe.n_components());
    AssertDimension(immersed_c.size(), immersed_fe.n_components());


    std::vector<unsigned int> space_gtl(space_fe.n_components(), numbers::invalid_unsigned_int);
    std::vector<unsigned int> immersed_gtl(immersed_fe.n_components(), numbers::invalid_unsigned_int);

    for (unsigned int i=0, j=0; i<space_gtl.size(); ++i)
      if (space_c[i])
        space_gtl[i] = j++;

    for (unsigned int i=0, j=0; i<immersed_gtl.size(); ++i)
      if (immersed_c[i])
        immersed_gtl[i] = j++;

    for (; cell != endc; ++cell)
      //    if(cell->subdomain_id() == this_mpi_process)
      {
        // Reinitialize the cell and the fe_values
        fe_v.reinit(cell);
        cell->get_dof_indices(dofs);

        const std::vector<Point<spacedim> > &Xpoints = fe_v.get_quadrature_points();

        // Get a list of outer cells, qpoints and maps.
        const auto cpm = GridTools::compute_point_locations(cache, Xpoints);
        const auto &cells = std::get<0>(cpm);

        for (unsigned int c=0; c<cells.size(); ++c)
          {
            // Get the ones in the current outer cell
            typename DoFHandler<dim0,spacedim>::cell_iterator
            ocell(*cells[c], &space_dh);
            ocell->get_dof_indices(odofs);
            for (unsigned int i=0; i<odofs.size(); ++i)
              {
                const auto comp_i = space_fe.system_to_component_index(i).first;
                if (space_gtl[comp_i] != numbers::invalid_unsigned_int)
                  {
                    for (unsigned int j=0; j<dofs.size(); ++j)
                      {
                        const auto comp_j = immersed_fe.system_to_component_index(j).first;
                        if (immersed_gtl[comp_j] == space_gtl[comp_i])
                          sparsity.add(odofs[i],dofs[j]);
                      }
                  }
              }
          }
      }
  }



  template<int dim0, int dim1, int spacedim, typename Matrix>
  void create_coupling_mass_matrix(const DoFHandler<dim0, spacedim> &space_dh,
                                   const DoFHandler<dim1, spacedim> &immersed_dh,
                                   const Quadrature<dim1>           &quad,
                                   Matrix                           &matrix,
                                   const ConstraintMatrix           &constraints,
                                   const ComponentMask              &space_comps,
                                   const ComponentMask              &immersed_comps,
                                   const Mapping<dim0, spacedim>    &space_mapping,
                                   const Mapping<dim1, spacedim>    &immersed_mapping)
  {
    AssertDimension(matrix.m(), space_dh.n_dofs());
    AssertDimension(matrix.n(), immersed_dh.n_dofs());

    const auto &space_fe = space_dh.get_fe();
    const auto &immersed_fe = immersed_dh.get_fe();

    // Dof indices
    std::vector<types::global_dof_index> dofs(immersed_fe.dofs_per_cell);
    std::vector<types::global_dof_index> odofs(space_fe.dofs_per_cell);

    GridTools::Cache<dim0,spacedim> cache(space_dh.get_triangulation(), space_mapping);

    // Take care of components
    ComponentMask space_c = space_comps;
    ComponentMask immersed_c = immersed_comps;
    if (space_c.size() == 0)
      space_c = ComponentMask(space_fe.n_components(), true);
    if (immersed_c.size() == 0)
      immersed_c = ComponentMask(immersed_fe.n_components(), true);

    AssertDimension(space_c.size(), space_fe.n_components());
    AssertDimension(immersed_c.size(), immersed_fe.n_components());

    std::vector<unsigned int> space_gtl(space_fe.n_components(), numbers::invalid_unsigned_int);
    std::vector<unsigned int> immersed_gtl(immersed_fe.n_components(), numbers::invalid_unsigned_int);

    for (unsigned int i=0, j=0; i<space_gtl.size(); ++i)
      if (space_c[i])
        space_gtl[i] = j++;

    for (unsigned int i=0, j=0; i<immersed_gtl.size(); ++i)
      if (immersed_c[i])
        immersed_gtl[i] = j++;

    FullMatrix<double> cell_matrix(space_dh.get_fe().dofs_per_cell,
                                   immersed_dh.get_fe().dofs_per_cell);

    FEValues<dim1,spacedim> fe_v(immersed_mapping, immersed_dh.get_fe(), quad,
                                 update_JxW_values |
                                 update_quadrature_points |
                                 update_values);

    // Now we run on ech cell, get a quadrature formula
    typename DoFHandler<dim1,spacedim>::active_cell_iterator
    cell = immersed_dh.begin_active(),
    endc = immersed_dh.end();

    for (; cell != endc; ++cell)
      {
        // Reinitialize the cell and the fe_values
        fe_v.reinit(cell);
        cell->get_dof_indices(dofs);

        const std::vector<Point<spacedim> > &Xpoints = fe_v.get_quadrature_points();

        // Get a list of outer cells, qpoints and maps.
        const auto cpm      = GridTools::compute_point_locations(cache, Xpoints);
        const auto &cells   = std::get<0>(cpm);
        const auto &qpoints = std::get<1>(cpm);
        const auto &maps    = std::get<2>(cpm);

        for (unsigned int c=0; c<cells.size(); ++c)
          {
            // Get the ones in the current outer cell
            typename DoFHandler<dim0,spacedim>::active_cell_iterator
            ocell(*cells[c], &space_dh);
            const std::vector< Point<spacedim> > &qps = qpoints[c];
            const std::vector< unsigned int > &ids = maps[c];

            FEValues<dim0,spacedim> o_fe_v(space_mapping, space_dh.get_fe(), qps,
                                           update_values);
            o_fe_v.reinit(ocell);
            ocell->get_dof_indices(odofs);

            // Reset the matrices.
            cell_matrix = 0;

            for (unsigned int i=0; i<space_dh.get_fe().dofs_per_cell; ++i)
              {
                const auto comp_i = space_dh.get_fe().system_to_component_index(i).first;
                if (space_gtl[comp_i] != numbers::invalid_unsigned_int)
                  for (unsigned int j=0; j<immersed_dh.get_fe().dofs_per_cell; ++j)
                    {
                      const auto comp_j = immersed_dh.get_fe().system_to_component_index(j).first;
                      if (space_gtl[comp_i] == immersed_gtl[comp_j])
                        for (unsigned int oq=0; oq<o_fe_v.n_quadrature_points; ++oq)
                          {
                            // Get the corrisponding q point
                            const unsigned int q=ids[oq];

                            cell_matrix(i,j) += ( fe_v.shape_value(j,q) *
                                                  o_fe_v.shape_value(i,oq) *
                                                  fe_v.JxW(q) );
                          }
                    }
              }

            // Now assemble the matrices
            constraints.distribute_local_to_global (cell_matrix, odofs, dofs, matrix);
          }
      }
  }
}
DEAL_II_NAMESPACE_CLOSE

#endif

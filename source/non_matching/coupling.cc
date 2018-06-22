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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/petsc_parallel_block_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <deal.II/non_matching/coupling.h>

DEAL_II_NAMESPACE_OPEN
namespace NonMatching
{
  template <int dim0,
            int dim1,
            int spacedim,
            typename Sparsity,
            typename number>
  void
  create_coupling_sparsity_pattern(
    const DoFHandler<dim0, spacedim> &space_dh,
    const DoFHandler<dim1, spacedim> &immersed_dh,
    const Quadrature<dim1> &          quad,
    Sparsity &                        sparsity,
    const AffineConstraints<number> & constraints,
    const ComponentMask &             space_comps,
    const ComponentMask &             immersed_comps,
    const Mapping<dim0, spacedim> &   space_mapping,
    const Mapping<dim1, spacedim> &   immersed_mapping)
  {
    GridTools::Cache<dim0, spacedim> cache(space_dh.get_triangulation(),
                                           space_mapping);
    create_coupling_sparsity_pattern(cache,
                                     space_dh,
                                     immersed_dh,
                                     quad,
                                     sparsity,
                                     constraints,
                                     space_comps,
                                     immersed_comps,
                                     immersed_mapping);
  }



  template <int dim0,
            int dim1,
            int spacedim,
            typename Sparsity,
            typename number>
  void
  create_coupling_sparsity_pattern(
    const GridTools::Cache<dim0, spacedim> &cache,
    const DoFHandler<dim0, spacedim> &      space_dh,
    const DoFHandler<dim1, spacedim> &      immersed_dh,
    const Quadrature<dim1> &                quad,
    Sparsity &                              sparsity,
    const AffineConstraints<number> &       constraints,
    const ComponentMask &                   space_comps,
    const ComponentMask &                   immersed_comps,
    const Mapping<dim1, spacedim> &         immersed_mapping)
  {
    AssertDimension(sparsity.n_rows(), space_dh.n_dofs());
    AssertDimension(sparsity.n_cols(), immersed_dh.n_dofs());
    static_assert(dim1 <= dim0, "This function can only work if dim1 <= dim0");
    Assert((dynamic_cast<
              const parallel::distributed::Triangulation<dim1, spacedim> *>(
              &immersed_dh.get_triangulation()) == nullptr),
           ExcNotImplemented());

    const auto &space_fe    = space_dh.get_fe();
    const auto &immersed_fe = immersed_dh.get_fe();

    // Now we run on ech cell, get a quadrature formula
    typename DoFHandler<dim1, spacedim>::active_cell_iterator
      cell = immersed_dh.begin_active(),
      endc = immersed_dh.end();

    // Dof indices
    std::vector<types::global_dof_index> dofs(immersed_fe.dofs_per_cell);
    std::vector<types::global_dof_index> odofs(space_fe.dofs_per_cell);

    FEValues<dim1, spacedim> fe_v(immersed_mapping,
                                  immersed_fe,
                                  quad,
                                  update_quadrature_points);

    // Take care of components
    const ComponentMask space_c =
      (space_comps.size() == 0 ? ComponentMask(space_fe.n_components(), true) :
                                 space_comps);

    const ComponentMask immersed_c =
      (immersed_comps.size() == 0 ?
         ComponentMask(immersed_fe.n_components(), true) :
         immersed_comps);

    AssertDimension(space_c.size(), space_fe.n_components());
    AssertDimension(immersed_c.size(), immersed_fe.n_components());

    std::vector<unsigned int> space_gtl(space_fe.n_components(),
                                        numbers::invalid_unsigned_int);
    std::vector<unsigned int> immersed_gtl(immersed_fe.n_components(),
                                           numbers::invalid_unsigned_int);

    for (unsigned int i = 0, j = 0; i < space_gtl.size(); ++i)
      if (space_c[i])
        space_gtl[i] = j++;

    for (unsigned int i = 0, j = 0; i < immersed_gtl.size(); ++i)
      if (immersed_c[i])
        immersed_gtl[i] = j++;

    // [TODO]: when the add_entries_local_to_global below will implement
    // the version with the dof_mask, this should be uncommented.
    //
    // // Construct a dof_mask, used to distribute entries to the sparsity
    // able< 2, bool > dof_mask(space_fe.dofs_per_cell,
    //                          immersed_fe.dofs_per_cell);
    // of_mask.fill(false);
    // or (unsigned int i=0; i<space_fe.dofs_per_cell; ++i)
    //  {
    //    const auto comp_i = space_fe.system_to_component_index(i).first;
    //    if (space_gtl[comp_i] != numbers::invalid_unsigned_int)
    //      for (unsigned int j=0; j<immersed_fe.dofs_per_cell; ++j)
    //        {
    //          const auto comp_j = immersed_fe.system_to_component_index(j).first;
    //          if (immersed_gtl[comp_j] == space_gtl[comp_i])
    //            dof_mask(i,j) = true;
    //        }
    //  }

    for (; cell != endc; ++cell)
      {
        // Reinitialize the cell and the fe_values
        fe_v.reinit(cell);
        cell->get_dof_indices(dofs);

        const std::vector<Point<spacedim>> &Xpoints =
          fe_v.get_quadrature_points();

        // Get a list of outer cells, qpoints and maps.
        const auto  cpm   = GridTools::compute_point_locations(cache, Xpoints);
        const auto &cells = std::get<0>(cpm);

        for (unsigned int c = 0; c < cells.size(); ++c)
          {
            // Get the ones in the current outer cell
            typename DoFHandler<dim0, spacedim>::cell_iterator ocell(*cells[c],
                                                                     &space_dh);
            // Make sure we act only on locally_owned cells
            if (ocell->is_locally_owned())
              {
                ocell->get_dof_indices(odofs);
                // [TODO]: When the following function will be implemented
                // for the case of non-trivial dof_mask, we should
                // uncomment the missing part.
                constraints.add_entries_local_to_global(
                  odofs, dofs, sparsity); //, true, dof_mask);
              }
          }
      }
  }



  template <int dim0, int dim1, int spacedim, typename Matrix>
  void
  create_coupling_mass_matrix(
    const DoFHandler<dim0, spacedim> &                    space_dh,
    const DoFHandler<dim1, spacedim> &                    immersed_dh,
    const Quadrature<dim1> &                              quad,
    Matrix &                                              matrix,
    const AffineConstraints<typename Matrix::value_type> &constraints,
    const ComponentMask &                                 space_comps,
    const ComponentMask &                                 immersed_comps,
    const Mapping<dim0, spacedim> &                       space_mapping,
    const Mapping<dim1, spacedim> &                       immersed_mapping)
  {
    GridTools::Cache<dim0, spacedim> cache(space_dh.get_triangulation(),
                                           space_mapping);
    create_coupling_mass_matrix(cache,
                                space_dh,
                                immersed_dh,
                                quad,
                                matrix,
                                constraints,
                                space_comps,
                                immersed_comps,
                                immersed_mapping);
  }



  template <int dim0, int dim1, int spacedim, typename Matrix>
  void
  create_coupling_mass_matrix(
    const GridTools::Cache<dim0, spacedim> &              cache,
    const DoFHandler<dim0, spacedim> &                    space_dh,
    const DoFHandler<dim1, spacedim> &                    immersed_dh,
    const Quadrature<dim1> &                              quad,
    Matrix &                                              matrix,
    const AffineConstraints<typename Matrix::value_type> &constraints,
    const ComponentMask &                                 space_comps,
    const ComponentMask &                                 immersed_comps,
    const Mapping<dim1, spacedim> &                       immersed_mapping)
  {
    AssertDimension(matrix.m(), space_dh.n_dofs());
    AssertDimension(matrix.n(), immersed_dh.n_dofs());
    static_assert(dim1 <= dim0, "This function can only work if dim1 <= dim0");
    Assert((dynamic_cast<
              const parallel::distributed::Triangulation<dim1, spacedim> *>(
              &immersed_dh.get_triangulation()) == nullptr),
           ExcNotImplemented());

    const auto &space_fe    = space_dh.get_fe();
    const auto &immersed_fe = immersed_dh.get_fe();

    // Dof indices
    std::vector<types::global_dof_index> dofs(immersed_fe.dofs_per_cell);
    std::vector<types::global_dof_index> odofs(space_fe.dofs_per_cell);

    // Take care of components
    const ComponentMask space_c =
      (space_comps.size() == 0 ? ComponentMask(space_fe.n_components(), true) :
                                 space_comps);

    const ComponentMask immersed_c =
      (immersed_comps.size() == 0 ?
         ComponentMask(immersed_fe.n_components(), true) :
         immersed_comps);

    AssertDimension(space_c.size(), space_fe.n_components());
    AssertDimension(immersed_c.size(), immersed_fe.n_components());

    std::vector<unsigned int> space_gtl(space_fe.n_components(),
                                        numbers::invalid_unsigned_int);
    std::vector<unsigned int> immersed_gtl(immersed_fe.n_components(),
                                           numbers::invalid_unsigned_int);

    for (unsigned int i = 0, j = 0; i < space_gtl.size(); ++i)
      if (space_c[i])
        space_gtl[i] = j++;

    for (unsigned int i = 0, j = 0; i < immersed_gtl.size(); ++i)
      if (immersed_c[i])
        immersed_gtl[i] = j++;

    FullMatrix<typename Matrix::value_type> cell_matrix(
      space_dh.get_fe().dofs_per_cell, immersed_dh.get_fe().dofs_per_cell);

    FEValues<dim1, spacedim> fe_v(immersed_mapping,
                                  immersed_dh.get_fe(),
                                  quad,
                                  update_JxW_values | update_quadrature_points |
                                    update_values);

    // Now we run on ech cell, get a quadrature formula
    typename DoFHandler<dim1, spacedim>::active_cell_iterator
      cell = immersed_dh.begin_active(),
      endc = immersed_dh.end();

    for (; cell != endc; ++cell)
      {
        // Reinitialize the cell and the fe_values
        fe_v.reinit(cell);
        cell->get_dof_indices(dofs);

        const std::vector<Point<spacedim>> &Xpoints =
          fe_v.get_quadrature_points();

        // Get a list of outer cells, qpoints and maps.
        const auto  cpm   = GridTools::compute_point_locations(cache, Xpoints);
        const auto &cells = std::get<0>(cpm);
        const auto &qpoints = std::get<1>(cpm);
        const auto &maps    = std::get<2>(cpm);

        for (unsigned int c = 0; c < cells.size(); ++c)
          {
            // Get the ones in the current outer cell
            typename DoFHandler<dim0, spacedim>::active_cell_iterator ocell(
              *cells[c], &space_dh);
            // Make sure we act only on locally_owned cells
            if (ocell->is_locally_owned())
              {
                const std::vector<Point<dim0>> & qps = qpoints[c];
                const std::vector<unsigned int> &ids = maps[c];

                FEValues<dim0, spacedim> o_fe_v(cache.get_mapping(),
                                                space_dh.get_fe(),
                                                qps,
                                                update_values);
                o_fe_v.reinit(ocell);
                ocell->get_dof_indices(odofs);

                // Reset the matrices.
                cell_matrix = typename Matrix::value_type();

                for (unsigned int i = 0; i < space_dh.get_fe().dofs_per_cell;
                     ++i)
                  {
                    const auto comp_i =
                      space_dh.get_fe().system_to_component_index(i).first;
                    if (space_gtl[comp_i] != numbers::invalid_unsigned_int)
                      for (unsigned int j = 0;
                           j < immersed_dh.get_fe().dofs_per_cell;
                           ++j)
                        {
                          const auto comp_j = immersed_dh.get_fe()
                                                .system_to_component_index(j)
                                                .first;
                          if (space_gtl[comp_i] == immersed_gtl[comp_j])
                            for (unsigned int oq = 0;
                                 oq < o_fe_v.n_quadrature_points;
                                 ++oq)
                              {
                                // Get the corresponding q point
                                const unsigned int q = ids[oq];

                                cell_matrix(i, j) +=
                                  (fe_v.shape_value(j, q) *
                                   o_fe_v.shape_value(i, oq) * fe_v.JxW(q));
                              }
                        }
                  }

                // Now assemble the matrices
                constraints.distribute_local_to_global(cell_matrix,
                                                       odofs,
                                                       dofs,
                                                       matrix);
              }
          }
      }
  }

#include "coupling.inst"
} // namespace NonMatching


DEAL_II_NAMESPACE_CLOSE

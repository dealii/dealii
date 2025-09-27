// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_hermite.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparsity_pattern_base.h>
#include <deal.II/lac/vector.h>

#include <algorithm>
#include <numeric>

DEAL_II_NAMESPACE_OPEN



namespace DoFTools
{
  template <int dim, int spacedim, typename number>
  void
  make_sparsity_pattern(const DoFHandler<dim, spacedim> &dof,
                        SparsityPatternBase             &sparsity,
                        const AffineConstraints<number> &constraints,
                        const bool                       keep_constrained_dofs,
                        const types::subdomain_id        subdomain_id)
  {
    const types::global_dof_index n_dofs = dof.n_dofs();
    (void)n_dofs;

    Assert(sparsity.n_rows() == n_dofs,
           ExcDimensionMismatch(sparsity.n_rows(), n_dofs));
    Assert(sparsity.n_cols() == n_dofs,
           ExcDimensionMismatch(sparsity.n_cols(), n_dofs));

    // If we have a distributed Triangulation only allow locally_owned
    // subdomain. Not setting a subdomain is also okay, because we skip
    // ghost cells in the loop below.
    if (const auto *triangulation = dynamic_cast<
          const parallel::DistributedTriangulationBase<dim, spacedim> *>(
          &dof.get_triangulation()))
      {
        Assert((subdomain_id == numbers::invalid_subdomain_id) ||
                 (subdomain_id == triangulation->locally_owned_subdomain()),
               ExcMessage(
                 "For distributed Triangulation objects and associated "
                 "DoFHandler objects, asking for any subdomain other than the "
                 "locally owned one does not make sense."));
      }

    const auto                 &fe_collection = dof.get_fe_collection();
    std::vector<Table<2, bool>> fe_dof_mask(fe_collection.size());
    for (unsigned int f = 0; f < fe_collection.size(); ++f)
      {
        fe_dof_mask[f] = fe_collection[f].get_local_dof_sparsity_pattern();
      }

    std::vector<types::global_dof_index> dofs_on_this_cell;
    dofs_on_this_cell.reserve(dof.get_fe_collection().max_dofs_per_cell());

    // In case we work with a distributed sparsity pattern of Trilinos
    // type, we only have to do the work if the current cell is owned by
    // the calling processor. Otherwise, just continue.
    for (const auto &cell : dof.active_cell_iterators())
      if (((subdomain_id == numbers::invalid_subdomain_id) ||
           (subdomain_id == cell->subdomain_id())) &&
          cell->is_locally_owned())
        {
          const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();
          dofs_on_this_cell.resize(dofs_per_cell);
          cell->get_dof_indices(dofs_on_this_cell);

          // make sparsity pattern for this cell. if no constraints pattern
          // was given, then the following call acts as if simply no
          // constraints existed
          const types::fe_index fe_index = cell->active_fe_index();
          if (fe_dof_mask[fe_index].empty())
            constraints.add_entries_local_to_global(dofs_on_this_cell,
                                                    sparsity,
                                                    keep_constrained_dofs);
          else
            constraints.add_entries_local_to_global(dofs_on_this_cell,
                                                    sparsity,
                                                    keep_constrained_dofs,
                                                    fe_dof_mask[fe_index]);
        }
  }



  template <int dim, int spacedim, typename number>
  void
  make_sparsity_pattern(const DoFHandler<dim, spacedim> &dof,
                        const Table<2, Coupling>        &couplings,
                        SparsityPatternBase             &sparsity,
                        const AffineConstraints<number> &constraints,
                        const bool                       keep_constrained_dofs,
                        const types::subdomain_id        subdomain_id)
  {
    const types::global_dof_index n_dofs = dof.n_dofs();
    (void)n_dofs;

    Assert(sparsity.n_rows() == n_dofs,
           ExcDimensionMismatch(sparsity.n_rows(), n_dofs));
    Assert(sparsity.n_cols() == n_dofs,
           ExcDimensionMismatch(sparsity.n_cols(), n_dofs));
    Assert(couplings.n_rows() == dof.get_fe(0).n_components(),
           ExcDimensionMismatch(couplings.n_rows(),
                                dof.get_fe(0).n_components()));
    Assert(couplings.n_cols() == dof.get_fe(0).n_components(),
           ExcDimensionMismatch(couplings.n_cols(),
                                dof.get_fe(0).n_components()));

    // If we have a distributed Triangulation only allow locally_owned
    // subdomain. Not setting a subdomain is also okay, because we skip
    // ghost cells in the loop below.
    if (const auto *triangulation = dynamic_cast<
          const parallel::DistributedTriangulationBase<dim, spacedim> *>(
          &dof.get_triangulation()))
      {
        Assert((subdomain_id == numbers::invalid_subdomain_id) ||
                 (subdomain_id == triangulation->locally_owned_subdomain()),
               ExcMessage(
                 "For distributed Triangulation objects and associated "
                 "DoFHandler objects, asking for any subdomain other than the "
                 "locally owned one does not make sense."));
      }

    const hp::FECollection<dim, spacedim> &fe_collection =
      dof.get_fe_collection();

    const std::vector<Table<2, Coupling>> dof_mask //(fe_collection.size())
      = dof_couplings_from_component_couplings(fe_collection, couplings);

    std::vector<Table<2, bool>> fe_dof_mask(fe_collection.size());
    for (unsigned int f = 0; f < fe_collection.size(); ++f)
      {
        fe_dof_mask[f] = fe_collection[f].get_local_dof_sparsity_pattern();
      }

    // Convert the dof_mask to bool_dof_mask so we can pass it
    // to constraints.add_entries_local_to_global()
    std::vector<Table<2, bool>> bool_dof_mask(fe_collection.size());
    for (unsigned int f = 0; f < fe_collection.size(); ++f)
      {
        bool_dof_mask[f].reinit(
          TableIndices<2>(fe_collection[f].n_dofs_per_cell(),
                          fe_collection[f].n_dofs_per_cell()));
        bool_dof_mask[f].fill(false);
        for (unsigned int i = 0; i < fe_collection[f].n_dofs_per_cell(); ++i)
          for (unsigned int j = 0; j < fe_collection[f].n_dofs_per_cell(); ++j)
            if (dof_mask[f](i, j) != none &&
                (fe_dof_mask[f].empty() || fe_dof_mask[f](i, j)))
              bool_dof_mask[f](i, j) = true;
      }

    std::vector<types::global_dof_index> dofs_on_this_cell(
      fe_collection.max_dofs_per_cell());

    // In case we work with a distributed sparsity pattern of Trilinos
    // type, we only have to do the work if the current cell is owned by
    // the calling processor. Otherwise, just continue.
    for (const auto &cell : dof.active_cell_iterators())
      if (((subdomain_id == numbers::invalid_subdomain_id) ||
           (subdomain_id == cell->subdomain_id())) &&
          cell->is_locally_owned())
        {
          const types::fe_index fe_index = cell->active_fe_index();
          const unsigned int    dofs_per_cell =
            fe_collection[fe_index].n_dofs_per_cell();

          dofs_on_this_cell.resize(dofs_per_cell);
          cell->get_dof_indices(dofs_on_this_cell);


          // make sparsity pattern for this cell. if no constraints pattern
          // was given, then the following call acts as if simply no
          // constraints existed
          constraints.add_entries_local_to_global(dofs_on_this_cell,
                                                  sparsity,
                                                  keep_constrained_dofs,
                                                  bool_dof_mask[fe_index]);
        }
  }



  template <int dim, int spacedim>
  void
  make_sparsity_pattern(const DoFHandler<dim, spacedim> &dof_row,
                        const DoFHandler<dim, spacedim> &dof_col,
                        SparsityPatternBase             &sparsity)
  {
    const types::global_dof_index n_dofs_row = dof_row.n_dofs();
    const types::global_dof_index n_dofs_col = dof_col.n_dofs();
    (void)n_dofs_row;
    (void)n_dofs_col;

    Assert(sparsity.n_rows() == n_dofs_row,
           ExcDimensionMismatch(sparsity.n_rows(), n_dofs_row));
    Assert(sparsity.n_cols() == n_dofs_col,
           ExcDimensionMismatch(sparsity.n_cols(), n_dofs_col));

    // It doesn't make sense to assemble sparsity patterns when the
    // Triangulations are both parallel (i.e., different cells are assigned to
    // different processors) and unequal: no processor will be responsible for
    // assembling coupling terms between dofs on a cell owned by one processor
    // and dofs on a cell owned by a different processor.
    if (dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
          &dof_row.get_triangulation()) != nullptr ||
        dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
          &dof_col.get_triangulation()) != nullptr)
      {
        Assert(&dof_row.get_triangulation() == &dof_col.get_triangulation(),
               ExcMessage("This function can only be used with with parallel "
                          "Triangulations when the Triangulations are equal."));
      }

    // TODO: Looks like wasteful memory management here

    using cell_iterator = typename DoFHandler<dim, spacedim>::cell_iterator;
    std::list<std::pair<cell_iterator, cell_iterator>> cell_list =
      GridTools::get_finest_common_cells(dof_row, dof_col);

#ifdef DEAL_II_WITH_MPI
    // get_finest_common_cells returns all cells (locally owned and otherwise)
    // for shared::Tria, but we don't want to assemble on cells that are not
    // locally owned so remove them
    if (dynamic_cast<const parallel::shared::Triangulation<dim, spacedim> *>(
          &dof_row.get_triangulation()) != nullptr ||
        dynamic_cast<const parallel::shared::Triangulation<dim, spacedim> *>(
          &dof_col.get_triangulation()) != nullptr)
      {
        const types::subdomain_id this_subdomain_id =
          dof_row.get_triangulation().locally_owned_subdomain();
        Assert(this_subdomain_id ==
                 dof_col.get_triangulation().locally_owned_subdomain(),
               ExcInternalError());
        cell_list.erase(
          std::remove_if(
            cell_list.begin(),
            cell_list.end(),
            [=](const std::pair<cell_iterator, cell_iterator> &pair) {
              return pair.first->subdomain_id() != this_subdomain_id ||
                     pair.second->subdomain_id() != this_subdomain_id;
            }),
          cell_list.end());
      }
#endif

    for (const auto &cell_pair : cell_list)
      {
        const cell_iterator cell_row = cell_pair.first;
        const cell_iterator cell_col = cell_pair.second;

        if (cell_row->is_active() && cell_col->is_active())
          {
            const unsigned int dofs_per_cell_row =
              cell_row->get_fe().n_dofs_per_cell();
            const unsigned int dofs_per_cell_col =
              cell_col->get_fe().n_dofs_per_cell();
            std::vector<types::global_dof_index> local_dof_indices_row(
              dofs_per_cell_row);
            std::vector<types::global_dof_index> local_dof_indices_col(
              dofs_per_cell_col);
            cell_row->get_dof_indices(local_dof_indices_row);
            cell_col->get_dof_indices(local_dof_indices_col);
            for (const auto &dof : local_dof_indices_row)
              sparsity.add_row_entries(dof,
                                       make_array_view(local_dof_indices_col));
          }
        else if (cell_row->has_children())
          {
            const std::vector<
              typename DoFHandler<dim, spacedim>::active_cell_iterator>
              child_cells =
                GridTools::get_active_child_cells<DoFHandler<dim, spacedim>>(
                  cell_row);
            for (unsigned int i = 0; i < child_cells.size(); ++i)
              {
                const typename DoFHandler<dim, spacedim>::cell_iterator
                                   cell_row_child = child_cells[i];
                const unsigned int dofs_per_cell_row =
                  cell_row_child->get_fe().n_dofs_per_cell();
                const unsigned int dofs_per_cell_col =
                  cell_col->get_fe().n_dofs_per_cell();
                std::vector<types::global_dof_index> local_dof_indices_row(
                  dofs_per_cell_row);
                std::vector<types::global_dof_index> local_dof_indices_col(
                  dofs_per_cell_col);
                cell_row_child->get_dof_indices(local_dof_indices_row);
                cell_col->get_dof_indices(local_dof_indices_col);
                for (const auto &dof : local_dof_indices_row)
                  sparsity.add_row_entries(
                    dof, make_array_view(local_dof_indices_col));
              }
          }
        else
          {
            std::vector<
              typename DoFHandler<dim, spacedim>::active_cell_iterator>
              child_cells =
                GridTools::get_active_child_cells<DoFHandler<dim, spacedim>>(
                  cell_col);
            for (unsigned int i = 0; i < child_cells.size(); ++i)
              {
                const typename DoFHandler<dim, spacedim>::active_cell_iterator
                                  &cell_col_child = child_cells[i];
                const unsigned int dofs_per_cell_row =
                  cell_row->get_fe().n_dofs_per_cell();
                const unsigned int dofs_per_cell_col =
                  cell_col_child->get_fe().n_dofs_per_cell();
                std::vector<types::global_dof_index> local_dof_indices_row(
                  dofs_per_cell_row);
                std::vector<types::global_dof_index> local_dof_indices_col(
                  dofs_per_cell_col);
                cell_row->get_dof_indices(local_dof_indices_row);
                cell_col_child->get_dof_indices(local_dof_indices_col);
                for (const auto &dof : local_dof_indices_row)
                  sparsity.add_row_entries(
                    dof, make_array_view(local_dof_indices_col));
              }
          }
      }
  }



  template <int dim, int spacedim>
  void
  make_boundary_sparsity_pattern(
    const DoFHandler<dim, spacedim>            &dof,
    const std::vector<types::global_dof_index> &dof_to_boundary_mapping,
    SparsityPatternBase                        &sparsity)
  {
    if (dim == 1)
      {
        // there are only 2 boundary indicators in 1d, so it is no
        // performance problem to call the other function
        std::map<types::boundary_id, const Function<spacedim, double> *>
          boundary_ids;
        boundary_ids[0] = nullptr;
        boundary_ids[1] = nullptr;
        make_boundary_sparsity_pattern<dim, spacedim>(dof,
                                                      boundary_ids,
                                                      dof_to_boundary_mapping,
                                                      sparsity);
        return;
      }

    const types::global_dof_index n_dofs = dof.n_dofs();
    (void)n_dofs;

    AssertDimension(dof_to_boundary_mapping.size(), n_dofs);
    AssertDimension(sparsity.n_rows(), dof.n_boundary_dofs());
    AssertDimension(sparsity.n_cols(), dof.n_boundary_dofs());
    if constexpr (running_in_debug_mode())
      {
        if (sparsity.n_rows() != 0)
          {
            types::global_dof_index max_element = 0;
            for (const types::global_dof_index index : dof_to_boundary_mapping)
              if ((index != numbers::invalid_dof_index) &&
                  (index > max_element))
                max_element = index;
            AssertDimension(max_element, sparsity.n_rows() - 1);
          }
      }

    std::vector<types::global_dof_index> dofs_on_this_face;
    dofs_on_this_face.reserve(dof.get_fe_collection().max_dofs_per_face());
    std::vector<types::global_dof_index> cols;

    // loop over all faces to check whether they are at a boundary. note
    // that we need not take special care of single lines (using
    // @p{cell->has_boundary_lines}), since we do not support boundaries of
    // dimension dim-2, and so every boundary line is also part of a
    // boundary face.
    for (const auto &cell : dof.active_cell_iterators())
      for (const unsigned int f : cell->face_indices())
        if (cell->at_boundary(f))
          {
            const unsigned int dofs_per_face =
              cell->get_fe().n_dofs_per_face(f);
            dofs_on_this_face.resize(dofs_per_face);
            cell->face(f)->get_dof_indices(dofs_on_this_face,
                                           cell->active_fe_index());

            // make sparsity pattern for this cell
            cols.clear();
            for (const auto &dof : dofs_on_this_face)
              cols.push_back(dof_to_boundary_mapping[dof]);
            // We are not guaranteed that the mapping to a second index space
            // is increasing so sort here to use the faster add_row_entries()
            // path
            std::sort(cols.begin(), cols.end());
            for (const auto &dof : dofs_on_this_face)
              sparsity.add_row_entries(dof_to_boundary_mapping[dof],
                                       make_array_view(cols),
                                       true);
          }
  }



  template <int dim, int spacedim, typename number>
  void
  make_boundary_sparsity_pattern(
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
                                               &boundary_ids,
    const std::vector<types::global_dof_index> &dof_to_boundary_mapping,
    SparsityPatternBase                        &sparsity)
  {
    if (dim == 1)
      {
        // first check left, then right boundary point
        for (unsigned int direction = 0; direction < 2; ++direction)
          {
            // if this boundary is not requested, then go on with next one
            if (boundary_ids.find(direction) == boundary_ids.end())
              continue;

            // find active cell at that boundary: first go to left/right,
            // then to children
            typename DoFHandler<dim, spacedim>::level_cell_iterator cell =
              dof.begin(0);
            while (!cell->at_boundary(direction))
              cell = cell->neighbor(direction);
            while (!cell->is_active())
              cell = cell->child(direction);

            const unsigned int dofs_per_vertex =
              cell->get_fe().n_dofs_per_vertex();
            std::vector<types::global_dof_index> boundary_dof_boundary_indices(
              dofs_per_vertex);

            // next get boundary mapped dof indices of boundary dofs
            for (unsigned int i = 0; i < dofs_per_vertex; ++i)
              boundary_dof_boundary_indices[i] =
                dof_to_boundary_mapping[cell->vertex_dof_index(direction, i)];

            std::sort(boundary_dof_boundary_indices.begin(),
                      boundary_dof_boundary_indices.end());
            for (const auto &dof : boundary_dof_boundary_indices)
              sparsity.add_row_entries(
                dof, make_array_view(boundary_dof_boundary_indices), true);
          }
        return;
      }

    const types::global_dof_index n_dofs = dof.n_dofs();
    (void)n_dofs;

    AssertDimension(dof_to_boundary_mapping.size(), n_dofs);
    Assert(boundary_ids.find(numbers::internal_face_boundary_id) ==
             boundary_ids.end(),
           (typename DoFHandler<dim, spacedim>::ExcInvalidBoundaryIndicator()));

    const bool fe_is_hermite = (dynamic_cast<const FE_Hermite<dim, spacedim> *>(
                                  &(dof.get_fe())) != nullptr);

    Assert(fe_is_hermite ||
             sparsity.n_rows() == dof.n_boundary_dofs(boundary_ids),
           ExcDimensionMismatch(sparsity.n_rows(),
                                dof.n_boundary_dofs(boundary_ids)));
    Assert(fe_is_hermite ||
             sparsity.n_cols() == dof.n_boundary_dofs(boundary_ids),
           ExcDimensionMismatch(sparsity.n_cols(),
                                dof.n_boundary_dofs(boundary_ids)));
    (void)fe_is_hermite;

    if constexpr (running_in_debug_mode())
      {
        if (sparsity.n_rows() != 0)
          {
            types::global_dof_index max_element = 0;
            for (const types::global_dof_index index : dof_to_boundary_mapping)
              if ((index != numbers::invalid_dof_index) &&
                  (index > max_element))
                max_element = index;
            AssertDimension(max_element, sparsity.n_rows() - 1);
          }
      }

    std::vector<types::global_dof_index> dofs_on_this_face;
    dofs_on_this_face.reserve(dof.get_fe_collection().max_dofs_per_face());
    std::vector<types::global_dof_index> cols;

    for (const auto &cell : dof.active_cell_iterators())
      for (const unsigned int f : cell->face_indices())
        if (boundary_ids.find(cell->face(f)->boundary_id()) !=
            boundary_ids.end())
          {
            const unsigned int dofs_per_face =
              cell->get_fe().n_dofs_per_face(f);
            dofs_on_this_face.resize(dofs_per_face);
            cell->face(f)->get_dof_indices(dofs_on_this_face,
                                           cell->active_fe_index());

            // make sparsity pattern for this cell
            cols.clear();
            for (const auto &dof : dofs_on_this_face)
              cols.push_back(dof_to_boundary_mapping[dof]);

            // Like the other one: sort once.
            std::sort(cols.begin(), cols.end());
            for (const auto &dof : dofs_on_this_face)
              sparsity.add_row_entries(dof_to_boundary_mapping[dof],
                                       make_array_view(cols),
                                       true);
          }
  }



  template <int dim, int spacedim, typename number>
  void
  make_flux_sparsity_pattern(const DoFHandler<dim, spacedim> &dof,
                             SparsityPatternBase             &sparsity,
                             const AffineConstraints<number> &constraints,
                             const bool                keep_constrained_dofs,
                             const types::subdomain_id subdomain_id)

  // TODO: QA: reduce the indentation level of this method..., Maier 2012

  {
    const types::global_dof_index n_dofs = dof.n_dofs();
    (void)n_dofs;

    AssertDimension(sparsity.n_rows(), n_dofs);
    AssertDimension(sparsity.n_cols(), n_dofs);

    // If we have a distributed Triangulation only allow locally_owned
    // subdomain. Not setting a subdomain is also okay, because we skip
    // ghost cells in the loop below.
    if (const auto *triangulation = dynamic_cast<
          const parallel::DistributedTriangulationBase<dim, spacedim> *>(
          &dof.get_triangulation()))
      {
        Assert((subdomain_id == numbers::invalid_subdomain_id) ||
                 (subdomain_id == triangulation->locally_owned_subdomain()),
               ExcMessage(
                 "For distributed Triangulation objects and associated "
                 "DoFHandler objects, asking for any subdomain other than the "
                 "locally owned one does not make sense."));
      }

    std::vector<types::global_dof_index> dofs_on_this_cell;
    std::vector<types::global_dof_index> dofs_on_other_cell;
    dofs_on_this_cell.reserve(dof.get_fe_collection().max_dofs_per_cell());
    dofs_on_other_cell.reserve(dof.get_fe_collection().max_dofs_per_cell());

    // TODO: in an old implementation, we used user flags before to tag
    // faces that were already touched. this way, we could reduce the work
    // a little bit. now, we instead add only data from one side. this
    // should be OK, but we need to actually verify it.

    // In case we work with a distributed sparsity pattern of Trilinos
    // type, we only have to do the work if the current cell is owned by
    // the calling processor. Otherwise, just continue.
    for (const auto &cell : dof.active_cell_iterators())
      if (((subdomain_id == numbers::invalid_subdomain_id) ||
           (subdomain_id == cell->subdomain_id())) &&
          cell->is_locally_owned())
        {
          const unsigned int n_dofs_on_this_cell =
            cell->get_fe().n_dofs_per_cell();
          dofs_on_this_cell.resize(n_dofs_on_this_cell);
          cell->get_dof_indices(dofs_on_this_cell);

          // make sparsity pattern for this cell. if no constraints pattern
          // was given, then the following call acts as if simply no
          // constraints existed
          constraints.add_entries_local_to_global(dofs_on_this_cell,
                                                  sparsity,
                                                  keep_constrained_dofs);

          for (const unsigned int face : cell->face_indices())
            {
              typename DoFHandler<dim, spacedim>::face_iterator cell_face =
                cell->face(face);
              const bool periodic_neighbor = cell->has_periodic_neighbor(face);
              if (!cell->at_boundary(face) || periodic_neighbor)
                {
                  typename DoFHandler<dim, spacedim>::level_cell_iterator
                    neighbor = cell->neighbor_or_periodic_neighbor(face);

                  // in 1d, we do not need to worry whether the neighbor
                  // might have children and then loop over those children.
                  // rather, we may as well go straight to the cell behind
                  // this particular cell's most terminal child
                  if (dim == 1)
                    while (neighbor->has_children())
                      neighbor = neighbor->child(face == 0 ? 1 : 0);

                  if (neighbor->has_children())
                    {
                      for (unsigned int sub_nr = 0;
                           sub_nr != cell_face->n_active_descendants();
                           ++sub_nr)
                        {
                          const typename DoFHandler<dim, spacedim>::
                            level_cell_iterator sub_neighbor =
                              periodic_neighbor ?
                                cell->periodic_neighbor_child_on_subface(
                                  face, sub_nr) :
                                cell->neighbor_child_on_subface(face, sub_nr);

                          const unsigned int n_dofs_on_neighbor =
                            sub_neighbor->get_fe().n_dofs_per_cell();
                          dofs_on_other_cell.resize(n_dofs_on_neighbor);
                          sub_neighbor->get_dof_indices(dofs_on_other_cell);

                          constraints.add_entries_local_to_global(
                            dofs_on_this_cell,
                            dofs_on_other_cell,
                            sparsity,
                            keep_constrained_dofs);
                          constraints.add_entries_local_to_global(
                            dofs_on_other_cell,
                            dofs_on_this_cell,
                            sparsity,
                            keep_constrained_dofs);
                          // only need to add this when the neighbor is not
                          // owned by the current processor, otherwise we add
                          // the entries for the neighbor there
                          if (sub_neighbor->subdomain_id() !=
                              cell->subdomain_id())
                            constraints.add_entries_local_to_global(
                              dofs_on_other_cell,
                              sparsity,
                              keep_constrained_dofs);
                        }
                    }
                  else
                    {
                      // Refinement edges are taken care of by coarser
                      // cells
                      if ((!periodic_neighbor &&
                           cell->neighbor_is_coarser(face)) ||
                          (periodic_neighbor &&
                           cell->periodic_neighbor_is_coarser(face)))
                        if (neighbor->subdomain_id() == cell->subdomain_id())
                          continue;

                      const unsigned int n_dofs_on_neighbor =
                        neighbor->get_fe().n_dofs_per_cell();
                      dofs_on_other_cell.resize(n_dofs_on_neighbor);

                      neighbor->get_dof_indices(dofs_on_other_cell);

                      constraints.add_entries_local_to_global(
                        dofs_on_this_cell,
                        dofs_on_other_cell,
                        sparsity,
                        keep_constrained_dofs);

                      // only need to add these in case the neighbor cell
                      // is not locally owned - otherwise, we touch each
                      // face twice and hence put the indices the other way
                      // around
                      if (!cell->neighbor_or_periodic_neighbor(face)
                             ->is_active() ||
                          (neighbor->subdomain_id() != cell->subdomain_id()))
                        {
                          constraints.add_entries_local_to_global(
                            dofs_on_other_cell,
                            dofs_on_this_cell,
                            sparsity,
                            keep_constrained_dofs);
                          if (neighbor->subdomain_id() != cell->subdomain_id())
                            constraints.add_entries_local_to_global(
                              dofs_on_other_cell,
                              sparsity,
                              keep_constrained_dofs);
                        }
                    }
                }
            }
        }
  }



  template <int dim, int spacedim>
  void
  make_flux_sparsity_pattern(const DoFHandler<dim, spacedim> &dof,
                             SparsityPatternBase             &sparsity)
  {
    AffineConstraints<double> dummy;
    make_flux_sparsity_pattern(dof, sparsity, dummy);
  }

  template <int dim, int spacedim>
  Table<2, Coupling>
  dof_couplings_from_component_couplings(
    const FiniteElement<dim, spacedim> &fe,
    const Table<2, Coupling>           &component_couplings)
  {
    Assert(component_couplings.n_rows() == fe.n_components(),
           ExcDimensionMismatch(component_couplings.n_rows(),
                                fe.n_components()));
    Assert(component_couplings.n_cols() == fe.n_components(),
           ExcDimensionMismatch(component_couplings.n_cols(),
                                fe.n_components()));

    const unsigned int n_dofs = fe.n_dofs_per_cell();

    Table<2, Coupling> dof_couplings(n_dofs, n_dofs);

    for (unsigned int i = 0; i < n_dofs; ++i)
      {
        const unsigned int ii =
          (fe.is_primitive(i) ?
             fe.system_to_component_index(i).first :
             fe.get_nonzero_components(i).first_selected_component());
        Assert(ii < fe.n_components(), ExcInternalError());

        for (unsigned int j = 0; j < n_dofs; ++j)
          {
            const unsigned int jj =
              (fe.is_primitive(j) ?
                 fe.system_to_component_index(j).first :
                 fe.get_nonzero_components(j).first_selected_component());
            Assert(jj < fe.n_components(), ExcInternalError());

            dof_couplings(i, j) = component_couplings(ii, jj);
          }
      }
    return dof_couplings;
  }



  template <int dim, int spacedim>
  std::vector<Table<2, Coupling>>
  dof_couplings_from_component_couplings(
    const hp::FECollection<dim, spacedim> &fe,
    const Table<2, Coupling>              &component_couplings)
  {
    std::vector<Table<2, Coupling>> return_value(fe.size());
    for (unsigned int i = 0; i < fe.size(); ++i)
      return_value[i] =
        dof_couplings_from_component_couplings(fe[i], component_couplings);

    return return_value;
  }



  namespace internal
  {
    namespace
    {
      // helper function
      template <typename Iterator, typename Iterator2>
      void
      add_cell_entries(
        const Iterator                             &cell,
        const unsigned int                          face_no,
        const Iterator2                            &neighbor,
        const unsigned int                          neighbor_face_no,
        const Table<2, Coupling>                   &flux_mask,
        const std::vector<types::global_dof_index> &dofs_on_this_cell,
        std::vector<types::global_dof_index>       &dofs_on_other_cell,
        std::vector<std::pair<SparsityPatternBase::size_type,
                              SparsityPatternBase::size_type>> &cell_entries)
      {
        dofs_on_other_cell.resize(neighbor->get_fe().n_dofs_per_cell());
        neighbor->get_dof_indices(dofs_on_other_cell);

        // Keep expensive data structures in separate vectors for inner j loop
        // in separate vectors
        boost::container::small_vector<unsigned int, 64>
          component_indices_neighbor(neighbor->get_fe().n_dofs_per_cell());
        boost::container::small_vector<bool, 64> support_on_face_i(
          neighbor->get_fe().n_dofs_per_cell());
        boost::container::small_vector<bool, 64> support_on_face_e(
          neighbor->get_fe().n_dofs_per_cell());
        for (unsigned int j = 0; j < neighbor->get_fe().n_dofs_per_cell(); ++j)
          {
            component_indices_neighbor[j] =
              (neighbor->get_fe().is_primitive(j) ?
                 neighbor->get_fe().system_to_component_index(j).first :
                 neighbor->get_fe()
                   .get_nonzero_components(j)
                   .first_selected_component());
            support_on_face_i[j] =
              neighbor->get_fe().has_support_on_face(j, face_no);
            support_on_face_e[j] =
              neighbor->get_fe().has_support_on_face(j, neighbor_face_no);
          }

        // For the parallel setting, must include also the diagonal
        // neighbor-neighbor coupling, otherwise those get included on the
        // other cell
        for (int f = 0; f < (neighbor->is_locally_owned() ? 1 : 2); ++f)
          {
            const auto &fe = (f == 0) ? cell->get_fe() : neighbor->get_fe();
            const auto &dofs_i =
              (f == 0) ? dofs_on_this_cell : dofs_on_other_cell;
            for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
              {
                const unsigned int ii =
                  (fe.is_primitive(i) ?
                     fe.system_to_component_index(i).first :
                     fe.get_nonzero_components(i).first_selected_component());

                Assert(ii < fe.n_components(), ExcInternalError());
                const bool i_non_zero_i =
                  fe.has_support_on_face(i,
                                         (f == 0 ? face_no : neighbor_face_no));

                for (unsigned int j = 0;
                     j < neighbor->get_fe().n_dofs_per_cell();
                     ++j)
                  {
                    const bool         j_non_zero_e = support_on_face_e[j];
                    const unsigned int jj = component_indices_neighbor[j];

                    Assert(jj < neighbor->get_fe().n_components(),
                           ExcInternalError());

                    if ((flux_mask(ii, jj) == always) ||
                        (flux_mask(ii, jj) == nonzero && i_non_zero_i &&
                         j_non_zero_e))
                      cell_entries.emplace_back(dofs_i[i],
                                                dofs_on_other_cell[j]);
                    if ((flux_mask(jj, ii) == always) ||
                        (flux_mask(jj, ii) == nonzero && j_non_zero_e &&
                         i_non_zero_i))
                      cell_entries.emplace_back(dofs_on_other_cell[j],
                                                dofs_i[i]);
                  }
              }
          }
      }



      // implementation of the same function in namespace DoFTools for
      // non-hp-DoFHandlers
      template <int dim, int spacedim, typename number>
      void
      make_flux_sparsity_pattern(
        const DoFHandler<dim, spacedim> &dof,
        SparsityPatternBase             &sparsity,
        const AffineConstraints<number> &constraints,
        const bool                       keep_constrained_dofs,
        const Table<2, Coupling>        &int_mask,
        const Table<2, Coupling>        &flux_mask,
        const types::subdomain_id        subdomain_id,
        const std::function<
          bool(const typename DoFHandler<dim, spacedim>::active_cell_iterator &,
               const unsigned int)> &face_has_flux_coupling)
      {
        std::vector<types::global_dof_index> rows;
        std::vector<std::pair<SparsityPatternBase::size_type,
                              SparsityPatternBase::size_type>>
          cell_entries;

        const dealii::hp::FECollection<dim, spacedim> &fe =
          dof.get_fe_collection();

        std::vector<types::global_dof_index> dofs_on_this_cell(
          dof.get_fe_collection().max_dofs_per_cell());
        std::vector<types::global_dof_index> dofs_on_other_cell(
          dof.get_fe_collection().max_dofs_per_cell());

        const unsigned int n_components = fe.n_components();
        AssertDimension(int_mask.size(0), n_components);
        AssertDimension(int_mask.size(1), n_components);
        AssertDimension(flux_mask.size(0), n_components);
        AssertDimension(flux_mask.size(1), n_components);

        // note that we also need to set the respective entries if flux_mask
        // says so. this is necessary since we need to consider all degrees
        // of freedom on a cell for interior faces.
        Table<2, Coupling> int_and_flux_mask(n_components, n_components);
        for (unsigned int c1 = 0; c1 < n_components; ++c1)
          for (unsigned int c2 = 0; c2 < n_components; ++c2)
            if (int_mask(c1, c2) != none || flux_mask(c1, c2) != none)
              int_and_flux_mask(c1, c2) = always;

        // Convert the int_dof_mask to bool_int_dof_mask so we can pass it
        // to constraints.add_entries_local_to_global()
        std::vector<Table<2, Coupling>> int_and_flux_dof_mask =
          dof_couplings_from_component_couplings(fe, int_and_flux_mask);

        // Convert int_and_flux_dof_mask to bool_int_and_flux_dof_mask so we
        // can pass it to constraints.add_entries_local_to_global()
        std::vector<Table<2, bool>> bool_int_and_flux_dof_mask(fe.size());
        for (unsigned int f = 0; f < fe.size(); ++f)
          {
            bool_int_and_flux_dof_mask[f].reinit(
              TableIndices<2>(fe[f].n_dofs_per_cell(),
                              fe[f].n_dofs_per_cell()));
            bool_int_and_flux_dof_mask[f].fill(false);
            for (unsigned int i = 0; i < fe[f].n_dofs_per_cell(); ++i)
              for (unsigned int j = 0; j < fe[f].n_dofs_per_cell(); ++j)
                if (int_and_flux_dof_mask[f](i, j) != none)
                  bool_int_and_flux_dof_mask[f](i, j) = true;
          }


        for (const auto &cell : dof.active_cell_iterators())
          if (((subdomain_id == numbers::invalid_subdomain_id) ||
               (subdomain_id == cell->subdomain_id())) &&
              cell->is_locally_owned())
            {
              dofs_on_this_cell.resize(cell->get_fe().n_dofs_per_cell());
              cell->get_dof_indices(dofs_on_this_cell);

              // make sparsity pattern for this cell also taking into
              // account the couplings due to face contributions on the same
              // cell
              constraints.add_entries_local_to_global(
                dofs_on_this_cell,
                sparsity,
                keep_constrained_dofs,
                bool_int_and_flux_dof_mask[cell->active_fe_index()]);

              // Loop over interior faces
              for (const unsigned int face : cell->face_indices())
                {
                  const bool periodic_neighbor =
                    cell->has_periodic_neighbor(face);

                  if ((!cell->at_boundary(face)) || periodic_neighbor)
                    {
                      typename DoFHandler<dim, spacedim>::level_cell_iterator
                        neighbor = cell->neighbor_or_periodic_neighbor(face);

                      // If the cells are on the same level (and both are
                      // active, locally-owned cells) then only add to the
                      // sparsity pattern if the current cell is 'greater' in
                      // the total ordering.
                      if (neighbor->level() == cell->level() &&
                          neighbor->index() > cell->index() &&
                          neighbor->is_active() && neighbor->is_locally_owned())
                        continue;

                      // If we are more refined then the neighbor, then we
                      // will automatically find the active neighbor cell when
                      // we call 'neighbor (face)' above. The opposite is not
                      // true; if the neighbor is more refined then the call
                      // 'neighbor (face)' will *not* return an active
                      // cell. Hence, only add things to the sparsity pattern
                      // if (when the levels are different) the neighbor is
                      // coarser than the current cell, except in the case
                      // when the neighbor is not locally owned.
                      if (neighbor->level() != cell->level() &&
                          ((!periodic_neighbor &&
                            !cell->neighbor_is_coarser(face)) ||
                           (periodic_neighbor &&
                            !cell->periodic_neighbor_is_coarser(face))) &&
                          neighbor->is_locally_owned())
                        continue; // (the neighbor is finer)

                      if (!face_has_flux_coupling(cell, face))
                        continue;

                      const unsigned int neighbor_face_no =
                        periodic_neighbor ?
                          cell->periodic_neighbor_face_no(face) :
                          cell->neighbor_face_no(face);

                      // In 1d, go straight to the cell behind this
                      // particular cell's most terminal cell. This makes us
                      // skip the if (neighbor->has_children()) section
                      // below. We need to do this since we otherwise
                      // iterate over the children of the face, which are
                      // always 0 in 1d.
                      if (dim == 1)
                        while (neighbor->has_children())
                          neighbor = neighbor->child(face == 0 ? 1 : 0);

                      if (neighbor->has_children())
                        {
                          for (unsigned int sub_nr = 0;
                               sub_nr != cell->face(face)->n_children();
                               ++sub_nr)
                            {
                              const typename DoFHandler<dim, spacedim>::
                                level_cell_iterator sub_neighbor =
                                  periodic_neighbor ?
                                    cell->periodic_neighbor_child_on_subface(
                                      face, sub_nr) :
                                    cell->neighbor_child_on_subface(face,
                                                                    sub_nr);
                              add_cell_entries(cell,
                                               face,
                                               sub_neighbor,
                                               neighbor_face_no,
                                               flux_mask,
                                               dofs_on_this_cell,
                                               dofs_on_other_cell,
                                               cell_entries);
                            }
                        }
                      else
                        add_cell_entries(cell,
                                         face,
                                         neighbor,
                                         neighbor_face_no,
                                         flux_mask,
                                         dofs_on_this_cell,
                                         dofs_on_other_cell,
                                         cell_entries);
                    }
                }
              sparsity.add_entries(make_array_view(cell_entries));
              cell_entries.clear();
            }
      }
    } // namespace

  } // namespace internal



  template <int dim, int spacedim>
  void
  make_flux_sparsity_pattern(const DoFHandler<dim, spacedim> &dof,
                             SparsityPatternBase             &sparsity,
                             const Table<2, Coupling>        &int_mask,
                             const Table<2, Coupling>        &flux_mask,
                             const types::subdomain_id        subdomain_id)
  {
    AffineConstraints<double> dummy;

    const bool keep_constrained_dofs = true;

    make_flux_sparsity_pattern(dof,
                               sparsity,
                               dummy,
                               keep_constrained_dofs,
                               int_mask,
                               flux_mask,
                               subdomain_id,
                               internal::always_couple_on_faces<dim, spacedim>);
  }



  template <int dim, int spacedim, typename number>
  void
  make_flux_sparsity_pattern(
    const DoFHandler<dim, spacedim> &dof,
    SparsityPatternBase             &sparsity,
    const AffineConstraints<number> &constraints,
    const bool                       keep_constrained_dofs,
    const Table<2, Coupling>        &int_mask,
    const Table<2, Coupling>        &flux_mask,
    const types::subdomain_id        subdomain_id,
    const std::function<
      bool(const typename DoFHandler<dim, spacedim>::active_cell_iterator &,
           const unsigned int)> &face_has_flux_coupling)
  {
    // do the error checking and frame code here, and then pass on to more
    // specialized functions in the internal namespace
    const types::global_dof_index n_dofs = dof.n_dofs();
    (void)n_dofs;
    const unsigned int n_comp = dof.get_fe(0).n_components();
    (void)n_comp;

    Assert(sparsity.n_rows() == n_dofs,
           ExcDimensionMismatch(sparsity.n_rows(), n_dofs));
    Assert(sparsity.n_cols() == n_dofs,
           ExcDimensionMismatch(sparsity.n_cols(), n_dofs));
    Assert(int_mask.n_rows() == n_comp,
           ExcDimensionMismatch(int_mask.n_rows(), n_comp));
    Assert(int_mask.n_cols() == n_comp,
           ExcDimensionMismatch(int_mask.n_cols(), n_comp));
    Assert(flux_mask.n_rows() == n_comp,
           ExcDimensionMismatch(flux_mask.n_rows(), n_comp));
    Assert(flux_mask.n_cols() == n_comp,
           ExcDimensionMismatch(flux_mask.n_cols(), n_comp));

    // If we have a distributed Triangulation only allow locally_owned
    // subdomain. Not setting a subdomain is also okay, because we skip
    // ghost cells in the loop below.
    if (const auto *triangulation = dynamic_cast<
          const parallel::DistributedTriangulationBase<dim, spacedim> *>(
          &dof.get_triangulation()))
      {
        Assert((subdomain_id == numbers::invalid_subdomain_id) ||
                 (subdomain_id == triangulation->locally_owned_subdomain()),
               ExcMessage(
                 "For distributed Triangulation objects and associated "
                 "DoFHandler objects, asking for any subdomain other than the "
                 "locally owned one does not make sense."));
      }

    Assert(
      face_has_flux_coupling,
      ExcMessage(
        "The function which specifies if a flux coupling occurs over a given "
        "face is empty."));

    internal::make_flux_sparsity_pattern(dof,
                                         sparsity,
                                         constraints,
                                         keep_constrained_dofs,
                                         int_mask,
                                         flux_mask,
                                         subdomain_id,
                                         face_has_flux_coupling);
  }

} // end of namespace DoFTools


// --------------------------------------------------- explicit instantiations

#include "dofs/dof_tools_sparsity.inst"



DEAL_II_NAMESPACE_CLOSE

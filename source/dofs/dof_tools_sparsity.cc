// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
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
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <algorithm>
#include <numeric>

DEAL_II_NAMESPACE_OPEN



namespace DoFTools
{
  template <int dim,
            int spacedim,
            typename SparsityPatternType,
            typename number>
  void
  make_sparsity_pattern(const DoFHandler<dim, spacedim> &dof,
                        SparsityPatternType &            sparsity,
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

    // If we have a distributed::Triangulation only allow locally_owned
    // subdomain. Not setting a subdomain is also okay, because we skip
    // ghost cells in the loop below.
    Assert((dof.get_triangulation().locally_owned_subdomain() ==
            numbers::invalid_subdomain_id) ||
             (subdomain_id == numbers::invalid_subdomain_id) ||
             (subdomain_id ==
              dof.get_triangulation().locally_owned_subdomain()),
           ExcMessage(
             "For parallel::distributed::Triangulation objects and "
             "associated DoF handler objects, asking for any subdomain other "
             "than the locally owned one does not make sense."));

    std::vector<types::global_dof_index> dofs_on_this_cell;
    dofs_on_this_cell.reserve(dof.get_fe_collection().max_dofs_per_cell());
    typename DoFHandler<dim, spacedim>::active_cell_iterator
      cell = dof.begin_active(),
      endc = dof.end();

    // In case we work with a distributed sparsity pattern of Trilinos
    // type, we only have to do the work if the current cell is owned by
    // the calling processor. Otherwise, just continue.
    for (; cell != endc; ++cell)
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
          constraints.add_entries_local_to_global(dofs_on_this_cell,
                                                  sparsity,
                                                  keep_constrained_dofs);
        }
  }



  template <int dim,
            int spacedim,
            typename SparsityPatternType,
            typename number>
  void
  make_sparsity_pattern(const DoFHandler<dim, spacedim> &dof,
                        const Table<2, Coupling> &       couplings,
                        SparsityPatternType &            sparsity,
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

    // If we have a distributed::Triangulation only allow locally_owned
    // subdomain. Not setting a subdomain is also okay, because we skip
    // ghost cells in the loop below.
    Assert((dof.get_triangulation().locally_owned_subdomain() ==
            numbers::invalid_subdomain_id) ||
             (subdomain_id == numbers::invalid_subdomain_id) ||
             (subdomain_id ==
              dof.get_triangulation().locally_owned_subdomain()),
           ExcMessage(
             "For parallel::distributed::Triangulation objects and "
             "associated DoF handler objects, asking for any subdomain other "
             "than the locally owned one does not make sense."));

    const hp::FECollection<dim, spacedim> &fe_collection =
      dof.get_fe_collection();

    const std::vector<Table<2, Coupling>> dof_mask //(fe_collection.size())
      = dof_couplings_from_component_couplings(fe_collection, couplings);

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
            if (dof_mask[f](i, j) != none)
              bool_dof_mask[f](i, j) = true;
      }

    std::vector<types::global_dof_index> dofs_on_this_cell(
      fe_collection.max_dofs_per_cell());
    typename DoFHandler<dim, spacedim>::active_cell_iterator
      cell = dof.begin_active(),
      endc = dof.end();

    // In case we work with a distributed sparsity pattern of Trilinos
    // type, we only have to do the work if the current cell is owned by
    // the calling processor. Otherwise, just continue.
    for (; cell != endc; ++cell)
      if (((subdomain_id == numbers::invalid_subdomain_id) ||
           (subdomain_id == cell->subdomain_id())) &&
          cell->is_locally_owned())
        {
          const unsigned int fe_index = cell->active_fe_index();
          const unsigned int dofs_per_cell =
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



  template <int dim, int spacedim, typename SparsityPatternType>
  void
  make_sparsity_pattern(const DoFHandler<dim, spacedim> &dof_row,
                        const DoFHandler<dim, spacedim> &dof_col,
                        SparsityPatternType &            sparsity)
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
            for (unsigned int i = 0; i < dofs_per_cell_row; ++i)
              sparsity.add_entries(local_dof_indices_row[i],
                                   local_dof_indices_col.begin(),
                                   local_dof_indices_col.end());
          }
        else if (cell_row->has_children())
          {
            const std::vector<
              typename DoFHandler<dim, spacedim>::active_cell_iterator>
              child_cells =
                GridTools::get_active_child_cells<DoFHandler<dim, spacedim>>(
                  cell_row);
            for (unsigned int i = 0; i < child_cells.size(); i++)
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
                for (unsigned int r = 0; r < dofs_per_cell_row; ++r)
                  sparsity.add_entries(local_dof_indices_row[r],
                                       local_dof_indices_col.begin(),
                                       local_dof_indices_col.end());
              }
          }
        else
          {
            std::vector<
              typename DoFHandler<dim, spacedim>::active_cell_iterator>
              child_cells =
                GridTools::get_active_child_cells<DoFHandler<dim, spacedim>>(
                  cell_col);
            for (unsigned int i = 0; i < child_cells.size(); i++)
              {
                const typename DoFHandler<dim, spacedim>::active_cell_iterator
                                   cell_col_child = child_cells[i];
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
                for (unsigned int r = 0; r < dofs_per_cell_row; ++r)
                  sparsity.add_entries(local_dof_indices_row[r],
                                       local_dof_indices_col.begin(),
                                       local_dof_indices_col.end());
              }
          }
      }
  }



  template <int dim, int spacedim, typename SparsityPatternType>
  void
  make_boundary_sparsity_pattern(
    const DoFHandler<dim, spacedim> &           dof,
    const std::vector<types::global_dof_index> &dof_to_boundary_mapping,
    SparsityPatternType &                       sparsity)
  {
    if (dim == 1)
      {
        // there are only 2 boundary indicators in 1d, so it is no
        // performance problem to call the other function
        std::map<types::boundary_id, const Function<spacedim, double> *>
          boundary_ids;
        boundary_ids[0] = nullptr;
        boundary_ids[1] = nullptr;
        make_boundary_sparsity_pattern<dim, spacedim, SparsityPatternType>(
          dof, boundary_ids, dof_to_boundary_mapping, sparsity);
        return;
      }

    const types::global_dof_index n_dofs = dof.n_dofs();
    (void)n_dofs;

    AssertDimension(dof_to_boundary_mapping.size(), n_dofs);
    AssertDimension(sparsity.n_rows(), dof.n_boundary_dofs());
    AssertDimension(sparsity.n_cols(), dof.n_boundary_dofs());
#ifdef DEBUG
    if (sparsity.n_rows() != 0)
      {
        types::global_dof_index max_element = 0;
        for (const types::global_dof_index index : dof_to_boundary_mapping)
          if ((index != numbers::invalid_dof_index) && (index > max_element))
            max_element = index;
        AssertDimension(max_element, sparsity.n_rows() - 1);
      }
#endif

    std::vector<types::global_dof_index> dofs_on_this_face;
    dofs_on_this_face.reserve(dof.get_fe_collection().max_dofs_per_face());

    // loop over all faces to check whether they are at a boundary. note
    // that we need not take special care of single lines (using
    // @p{cell->has_boundary_lines}), since we do not support boundaries of
    // dimension dim-2, and so every boundary line is also part of a
    // boundary face.
    typename DoFHandler<dim, spacedim>::active_cell_iterator
      cell = dof.begin_active(),
      endc = dof.end();
    for (; cell != endc; ++cell)
      for (const unsigned int f : cell->face_indices())
        if (cell->at_boundary(f))
          {
            const unsigned int dofs_per_face =
              cell->get_fe().n_dofs_per_face(f);
            dofs_on_this_face.resize(dofs_per_face);
            cell->face(f)->get_dof_indices(dofs_on_this_face,
                                           cell->active_fe_index());

            // make sparsity pattern for this cell
            for (unsigned int i = 0; i < dofs_per_face; ++i)
              for (unsigned int j = 0; j < dofs_per_face; ++j)
                sparsity.add(dof_to_boundary_mapping[dofs_on_this_face[i]],
                             dof_to_boundary_mapping[dofs_on_this_face[j]]);
          }
  }



  template <int dim,
            int spacedim,
            typename SparsityPatternType,
            typename number>
  void
  make_boundary_sparsity_pattern(
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                         boundary_ids,
    const std::vector<types::global_dof_index> &dof_to_boundary_mapping,
    SparsityPatternType &                       sparsity)
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

            for (unsigned int i = 0; i < dofs_per_vertex; ++i)
              sparsity.add_entries(boundary_dof_boundary_indices[i],
                                   boundary_dof_boundary_indices.begin(),
                                   boundary_dof_boundary_indices.end());
          }
        return;
      }

    const types::global_dof_index n_dofs = dof.n_dofs();
    (void)n_dofs;

    AssertDimension(dof_to_boundary_mapping.size(), n_dofs);
    Assert(boundary_ids.find(numbers::internal_face_boundary_id) ==
             boundary_ids.end(),
           (typename DoFHandler<dim, spacedim>::ExcInvalidBoundaryIndicator()));
    Assert(sparsity.n_rows() == dof.n_boundary_dofs(boundary_ids),
           ExcDimensionMismatch(sparsity.n_rows(),
                                dof.n_boundary_dofs(boundary_ids)));
    Assert(sparsity.n_cols() == dof.n_boundary_dofs(boundary_ids),
           ExcDimensionMismatch(sparsity.n_cols(),
                                dof.n_boundary_dofs(boundary_ids)));
#ifdef DEBUG
    if (sparsity.n_rows() != 0)
      {
        types::global_dof_index max_element = 0;
        for (const types::global_dof_index index : dof_to_boundary_mapping)
          if ((index != numbers::invalid_dof_index) && (index > max_element))
            max_element = index;
        AssertDimension(max_element, sparsity.n_rows() - 1);
      }
#endif

    std::vector<types::global_dof_index> dofs_on_this_face;
    dofs_on_this_face.reserve(dof.get_fe_collection().max_dofs_per_face());
    typename DoFHandler<dim, spacedim>::active_cell_iterator
      cell = dof.begin_active(),
      endc = dof.end();
    for (; cell != endc; ++cell)
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
            for (unsigned int i = 0; i < dofs_per_face; ++i)
              for (unsigned int j = 0; j < dofs_per_face; ++j)
                sparsity.add(dof_to_boundary_mapping[dofs_on_this_face[i]],
                             dof_to_boundary_mapping[dofs_on_this_face[j]]);
          }
  }



  template <int dim,
            int spacedim,
            typename SparsityPatternType,
            typename number>
  void
  make_flux_sparsity_pattern(const DoFHandler<dim, spacedim> &dof,
                             SparsityPatternType &            sparsity,
                             const AffineConstraints<number> &constraints,
                             const bool                keep_constrained_dofs,
                             const types::subdomain_id subdomain_id)

  // TODO: QA: reduce the indentation level of this method..., Maier 2012

  {
    const types::global_dof_index n_dofs = dof.n_dofs();
    (void)n_dofs;

    AssertDimension(sparsity.n_rows(), n_dofs);
    AssertDimension(sparsity.n_cols(), n_dofs);

    // If we have a distributed::Triangulation only allow locally_owned
    // subdomain. Not setting a subdomain is also okay, because we skip
    // ghost cells in the loop below.
    Assert((dof.get_triangulation().locally_owned_subdomain() ==
            numbers::invalid_subdomain_id) ||
             (subdomain_id == numbers::invalid_subdomain_id) ||
             (subdomain_id ==
              dof.get_triangulation().locally_owned_subdomain()),
           ExcMessage(
             "For parallel::distributed::Triangulation objects and "
             "associated DoF handler objects, asking for any subdomain other "
             "than the locally owned one does not make sense."));

    std::vector<types::global_dof_index> dofs_on_this_cell;
    std::vector<types::global_dof_index> dofs_on_other_cell;
    dofs_on_this_cell.reserve(dof.get_fe_collection().max_dofs_per_cell());
    dofs_on_other_cell.reserve(dof.get_fe_collection().max_dofs_per_cell());
    typename DoFHandler<dim, spacedim>::active_cell_iterator
      cell = dof.begin_active(),
      endc = dof.end();

    // TODO: in an old implementation, we used user flags before to tag
    // faces that were already touched. this way, we could reduce the work
    // a little bit. now, we instead add only data from one side. this
    // should be OK, but we need to actually verify it.

    // In case we work with a distributed sparsity pattern of Trilinos
    // type, we only have to do the work if the current cell is owned by
    // the calling processor. Otherwise, just continue.
    for (; cell != endc; ++cell)
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



  template <int dim, int spacedim, typename SparsityPatternType>
  void
  make_flux_sparsity_pattern(const DoFHandler<dim, spacedim> &dof,
                             SparsityPatternType &            sparsity)
  {
    AffineConstraints<double> dummy;
    make_flux_sparsity_pattern(dof, sparsity, dummy);
  }

  template <int dim, int spacedim>
  Table<2, Coupling>
  dof_couplings_from_component_couplings(
    const FiniteElement<dim, spacedim> &fe,
    const Table<2, Coupling> &          component_couplings)
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
    const Table<2, Coupling> &             component_couplings)
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
      // implementation of the same function in namespace DoFTools for
      // non-hp-DoFHandlers
      template <int dim,
                int spacedim,
                typename SparsityPatternType,
                typename number>
      void
      make_flux_sparsity_pattern(
        const DoFHandler<dim, spacedim> &dof,
        SparsityPatternType &            sparsity,
        const AffineConstraints<number> &constraints,
        const bool                       keep_constrained_dofs,
        const Table<2, Coupling> &       int_mask,
        const Table<2, Coupling> &       flux_mask,
        const types::subdomain_id        subdomain_id,
        const std::function<
          bool(const typename DoFHandler<dim, spacedim>::active_cell_iterator &,
               const unsigned int)> &face_has_flux_coupling)
      {
        if (dof.has_hp_capabilities() == false)
          {
            const FiniteElement<dim, spacedim> &fe = dof.get_fe();

            std::vector<types::global_dof_index> dofs_on_this_cell(
              fe.n_dofs_per_cell());
            std::vector<types::global_dof_index> dofs_on_other_cell(
              fe.n_dofs_per_cell());

            const Table<2, Coupling>
              int_dof_mask =
                dof_couplings_from_component_couplings(fe, int_mask),
              flux_dof_mask =
                dof_couplings_from_component_couplings(fe, flux_mask);

            Table<2, bool> support_on_face(fe.n_dofs_per_cell(),
                                           GeometryInfo<dim>::faces_per_cell);
            for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
              for (const unsigned int f : GeometryInfo<dim>::face_indices())
                support_on_face(i, f) = fe.has_support_on_face(i, f);

            // Convert the int_dof_mask to bool_int_dof_mask so we can pass it
            // to constraints.add_entries_local_to_global()
            Table<2, bool> bool_int_dof_mask(fe.n_dofs_per_cell(),
                                             fe.n_dofs_per_cell());
            bool_int_dof_mask.fill(false);
            for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
              for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
                if (int_dof_mask(i, j) != none)
                  bool_int_dof_mask(i, j) = true;

            typename DoFHandler<dim, spacedim>::active_cell_iterator
              cell = dof.begin_active(),
              endc = dof.end();
            for (; cell != endc; ++cell)
              if (((subdomain_id == numbers::invalid_subdomain_id) ||
                   (subdomain_id == cell->subdomain_id())) &&
                  cell->is_locally_owned())
                {
                  cell->get_dof_indices(dofs_on_this_cell);
                  // make sparsity pattern for this cell
                  constraints.add_entries_local_to_global(dofs_on_this_cell,
                                                          sparsity,
                                                          keep_constrained_dofs,
                                                          bool_int_dof_mask);
                  // Loop over all interior neighbors
                  for (const unsigned int face_n : cell->face_indices())
                    {
                      const typename DoFHandler<dim, spacedim>::face_iterator
                        cell_face = cell->face(face_n);

                      const bool periodic_neighbor =
                        cell->has_periodic_neighbor(face_n);

                      if (cell->at_boundary(face_n) && (!periodic_neighbor))
                        {
                          for (unsigned int i = 0; i < fe.n_dofs_per_cell();
                               ++i)
                            {
                              const bool i_non_zero_i =
                                support_on_face(i, face_n);
                              for (unsigned int j = 0; j < fe.n_dofs_per_cell();
                                   ++j)
                                {
                                  const bool j_non_zero_i =
                                    support_on_face(j, face_n);

                                  if (flux_dof_mask(i, j) == always ||
                                      (flux_dof_mask(i, j) == nonzero &&
                                       i_non_zero_i && j_non_zero_i))
                                    sparsity.add(dofs_on_this_cell[i],
                                                 dofs_on_this_cell[j]);
                                }
                            }
                        }
                      else
                        {
                          if (!face_has_flux_coupling(cell, face_n))
                            continue;

                          typename DoFHandler<dim,
                                              spacedim>::level_cell_iterator
                            neighbor =
                              cell->neighbor_or_periodic_neighbor(face_n);
                          // If the cells are on the same level (and both are
                          // active, locally-owned cells) then only add to the
                          // sparsity pattern if the current cell is 'greater'
                          // in the total ordering.
                          if (neighbor->level() == cell->level() &&
                              neighbor->index() > cell->index() &&
                              neighbor->is_active() &&
                              neighbor->is_locally_owned())
                            continue;
                          // If we are more refined then the neighbor, then we
                          // will automatically find the active neighbor cell
                          // when we call 'neighbor (face_n)' above. The
                          // opposite is not true; if the neighbor is more
                          // refined then the call 'neighbor (face_n)' will
                          // *not* return an active cell. Hence, only add things
                          // to the sparsity pattern if (when the levels are
                          // different) the neighbor is coarser than the current
                          // cell.
                          //
                          // Like above, do not use this optimization if the
                          // neighbor is not locally owned.
                          if (neighbor->level() != cell->level() &&
                              ((!periodic_neighbor &&
                                !cell->neighbor_is_coarser(face_n)) ||
                               (periodic_neighbor &&
                                !cell->periodic_neighbor_is_coarser(face_n))) &&
                              neighbor->is_locally_owned())
                            continue; // (the neighbor is finer)

                          const unsigned int neighbor_face_n =
                            periodic_neighbor ?
                              cell->periodic_neighbor_face_no(face_n) :
                              cell->neighbor_face_no(face_n);


                          // In 1D, go straight to the cell behind this
                          // particular cell's most terminal cell. This makes us
                          // skip the if (neighbor->has_children()) section
                          // below. We need to do this since we otherwise
                          // iterate over the children of the face, which are
                          // always 0 in 1D.
                          if (dim == 1)
                            while (neighbor->has_children())
                              neighbor = neighbor->child(face_n == 0 ? 1 : 0);

                          if (neighbor->has_children())
                            {
                              for (unsigned int sub_nr = 0;
                                   sub_nr != cell_face->n_children();
                                   ++sub_nr)
                                {
                                  const typename DoFHandler<dim, spacedim>::
                                    level_cell_iterator sub_neighbor =
                                      periodic_neighbor ?
                                        cell
                                          ->periodic_neighbor_child_on_subface(
                                            face_n, sub_nr) :
                                        cell->neighbor_child_on_subface(face_n,
                                                                        sub_nr);

                                  sub_neighbor->get_dof_indices(
                                    dofs_on_other_cell);
                                  for (unsigned int i = 0;
                                       i < fe.n_dofs_per_cell();
                                       ++i)
                                    {
                                      const bool i_non_zero_i =
                                        support_on_face(i, face_n);
                                      const bool i_non_zero_e =
                                        support_on_face(i, neighbor_face_n);
                                      for (unsigned int j = 0;
                                           j < fe.n_dofs_per_cell();
                                           ++j)
                                        {
                                          const bool j_non_zero_i =
                                            support_on_face(j, face_n);
                                          const bool j_non_zero_e =
                                            support_on_face(j, neighbor_face_n);

                                          if (flux_dof_mask(i, j) == always)
                                            {
                                              sparsity.add(
                                                dofs_on_this_cell[i],
                                                dofs_on_other_cell[j]);
                                              sparsity.add(
                                                dofs_on_other_cell[i],
                                                dofs_on_this_cell[j]);
                                              sparsity.add(
                                                dofs_on_this_cell[i],
                                                dofs_on_this_cell[j]);
                                              sparsity.add(
                                                dofs_on_other_cell[i],
                                                dofs_on_other_cell[j]);
                                            }
                                          else if (flux_dof_mask(i, j) ==
                                                   nonzero)
                                            {
                                              if (i_non_zero_i && j_non_zero_e)
                                                sparsity.add(
                                                  dofs_on_this_cell[i],
                                                  dofs_on_other_cell[j]);
                                              if (i_non_zero_e && j_non_zero_i)
                                                sparsity.add(
                                                  dofs_on_other_cell[i],
                                                  dofs_on_this_cell[j]);
                                              if (i_non_zero_i && j_non_zero_i)
                                                sparsity.add(
                                                  dofs_on_this_cell[i],
                                                  dofs_on_this_cell[j]);
                                              if (i_non_zero_e && j_non_zero_e)
                                                sparsity.add(
                                                  dofs_on_other_cell[i],
                                                  dofs_on_other_cell[j]);
                                            }

                                          if (flux_dof_mask(j, i) == always)
                                            {
                                              sparsity.add(
                                                dofs_on_this_cell[j],
                                                dofs_on_other_cell[i]);
                                              sparsity.add(
                                                dofs_on_other_cell[j],
                                                dofs_on_this_cell[i]);
                                              sparsity.add(
                                                dofs_on_this_cell[j],
                                                dofs_on_this_cell[i]);
                                              sparsity.add(
                                                dofs_on_other_cell[j],
                                                dofs_on_other_cell[i]);
                                            }
                                          else if (flux_dof_mask(j, i) ==
                                                   nonzero)
                                            {
                                              if (j_non_zero_i && i_non_zero_e)
                                                sparsity.add(
                                                  dofs_on_this_cell[j],
                                                  dofs_on_other_cell[i]);
                                              if (j_non_zero_e && i_non_zero_i)
                                                sparsity.add(
                                                  dofs_on_other_cell[j],
                                                  dofs_on_this_cell[i]);
                                              if (j_non_zero_i && i_non_zero_i)
                                                sparsity.add(
                                                  dofs_on_this_cell[j],
                                                  dofs_on_this_cell[i]);
                                              if (j_non_zero_e && i_non_zero_e)
                                                sparsity.add(
                                                  dofs_on_other_cell[j],
                                                  dofs_on_other_cell[i]);
                                            }
                                        }
                                    }
                                }
                            }
                          else
                            {
                              neighbor->get_dof_indices(dofs_on_other_cell);
                              for (unsigned int i = 0; i < fe.n_dofs_per_cell();
                                   ++i)
                                {
                                  const bool i_non_zero_i =
                                    support_on_face(i, face_n);
                                  const bool i_non_zero_e =
                                    support_on_face(i, neighbor_face_n);
                                  for (unsigned int j = 0;
                                       j < fe.n_dofs_per_cell();
                                       ++j)
                                    {
                                      const bool j_non_zero_i =
                                        support_on_face(j, face_n);
                                      const bool j_non_zero_e =
                                        support_on_face(j, neighbor_face_n);
                                      if (flux_dof_mask(i, j) == always)
                                        {
                                          sparsity.add(dofs_on_this_cell[i],
                                                       dofs_on_other_cell[j]);
                                          sparsity.add(dofs_on_other_cell[i],
                                                       dofs_on_this_cell[j]);
                                          sparsity.add(dofs_on_this_cell[i],
                                                       dofs_on_this_cell[j]);
                                          sparsity.add(dofs_on_other_cell[i],
                                                       dofs_on_other_cell[j]);
                                        }
                                      if (flux_dof_mask(i, j) == nonzero)
                                        {
                                          if (i_non_zero_i && j_non_zero_e)
                                            sparsity.add(dofs_on_this_cell[i],
                                                         dofs_on_other_cell[j]);
                                          if (i_non_zero_e && j_non_zero_i)
                                            sparsity.add(dofs_on_other_cell[i],
                                                         dofs_on_this_cell[j]);
                                          if (i_non_zero_i && j_non_zero_i)
                                            sparsity.add(dofs_on_this_cell[i],
                                                         dofs_on_this_cell[j]);
                                          if (i_non_zero_e && j_non_zero_e)
                                            sparsity.add(dofs_on_other_cell[i],
                                                         dofs_on_other_cell[j]);
                                        }

                                      if (flux_dof_mask(j, i) == always)
                                        {
                                          sparsity.add(dofs_on_this_cell[j],
                                                       dofs_on_other_cell[i]);
                                          sparsity.add(dofs_on_other_cell[j],
                                                       dofs_on_this_cell[i]);
                                          sparsity.add(dofs_on_this_cell[j],
                                                       dofs_on_this_cell[i]);
                                          sparsity.add(dofs_on_other_cell[j],
                                                       dofs_on_other_cell[i]);
                                        }
                                      if (flux_dof_mask(j, i) == nonzero)
                                        {
                                          if (j_non_zero_i && i_non_zero_e)
                                            sparsity.add(dofs_on_this_cell[j],
                                                         dofs_on_other_cell[i]);
                                          if (j_non_zero_e && i_non_zero_i)
                                            sparsity.add(dofs_on_other_cell[j],
                                                         dofs_on_this_cell[i]);
                                          if (j_non_zero_i && i_non_zero_i)
                                            sparsity.add(dofs_on_this_cell[j],
                                                         dofs_on_this_cell[i]);
                                          if (j_non_zero_e && i_non_zero_e)
                                            sparsity.add(dofs_on_other_cell[j],
                                                         dofs_on_other_cell[i]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
          }
        else
          {
            // while the implementation above is quite optimized and caches a
            // lot of data (see e.g. the int/flux_dof_mask tables), this is no
            // longer practical for the hp-version since we would have to have
            // it for all combinations of elements in the hp::FECollection.
            // consequently, the implementation here is simpler and probably
            // less efficient but at least readable...

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


            typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator
              cell = dof.begin_active(),
              endc = dof.end();
            for (; cell != endc; ++cell)
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
                      const typename dealii::DoFHandler<dim,
                                                        spacedim>::face_iterator
                        cell_face = cell->face(face);

                      const bool periodic_neighbor =
                        cell->has_periodic_neighbor(face);

                      if ((!cell->at_boundary(face)) || periodic_neighbor)
                        {
                          typename dealii::DoFHandler<dim, spacedim>::
                            level_cell_iterator neighbor =
                              cell->neighbor_or_periodic_neighbor(face);

                          if (!face_has_flux_coupling(cell, face))
                            continue;

                          // Like the non-hp-case: If the cells are on the same
                          // level (and both are active, locally-owned cells)
                          // then only add to the sparsity pattern if the
                          // current cell is 'greater' in the total ordering.
                          if (neighbor->level() == cell->level() &&
                              neighbor->index() > cell->index() &&
                              neighbor->is_active() &&
                              neighbor->is_locally_owned())
                            continue;
                          // Again, like the non-hp-case: If we are more refined
                          // then the neighbor, then we will automatically find
                          // the active neighbor cell when we call 'neighbor
                          // (face)' above. The opposite is not true; if the
                          // neighbor is more refined then the call 'neighbor
                          // (face)' will *not* return an active cell. Hence,
                          // only add things to the sparsity pattern if (when
                          // the levels are different) the neighbor is coarser
                          // than the current cell.
                          //
                          // Like above, do not use this optimization if the
                          // neighbor is not locally owned.
                          if (neighbor->level() != cell->level() &&
                              ((!periodic_neighbor &&
                                !cell->neighbor_is_coarser(face)) ||
                               (periodic_neighbor &&
                                !cell->periodic_neighbor_is_coarser(face))) &&
                              neighbor->is_locally_owned())
                            continue; // (the neighbor is finer)

                          // In 1D, go straight to the cell behind this
                          // particular cell's most terminal cell. This makes us
                          // skip the if (neighbor->has_children()) section
                          // below. We need to do this since we otherwise
                          // iterate over the children of the face, which are
                          // always 0 in 1D.
                          if (dim == 1)
                            while (neighbor->has_children())
                              neighbor = neighbor->child(face == 0 ? 1 : 0);

                          if (neighbor->has_children())
                            {
                              for (unsigned int sub_nr = 0;
                                   sub_nr != cell_face->n_children();
                                   ++sub_nr)
                                {
                                  const typename dealii::DoFHandler<dim,
                                                                    spacedim>::
                                    level_cell_iterator sub_neighbor =
                                      periodic_neighbor ?
                                        cell
                                          ->periodic_neighbor_child_on_subface(
                                            face, sub_nr) :
                                        cell->neighbor_child_on_subface(face,
                                                                        sub_nr);

                                  dofs_on_other_cell.resize(
                                    sub_neighbor->get_fe().n_dofs_per_cell());
                                  sub_neighbor->get_dof_indices(
                                    dofs_on_other_cell);
                                  for (unsigned int i = 0;
                                       i < cell->get_fe().n_dofs_per_cell();
                                       ++i)
                                    {
                                      const unsigned int ii =
                                        (cell->get_fe().is_primitive(i) ?
                                           cell->get_fe()
                                             .system_to_component_index(i)
                                             .first :
                                           cell->get_fe()
                                             .get_nonzero_components(i)
                                             .first_selected_component());

                                      Assert(ii < cell->get_fe().n_components(),
                                             ExcInternalError());

                                      for (unsigned int j = 0;
                                           j < sub_neighbor->get_fe()
                                                 .n_dofs_per_cell();
                                           ++j)
                                        {
                                          const unsigned int jj =
                                            (sub_neighbor->get_fe()
                                                 .is_primitive(j) ?
                                               sub_neighbor->get_fe()
                                                 .system_to_component_index(j)
                                                 .first :
                                               sub_neighbor->get_fe()
                                                 .get_nonzero_components(j)
                                                 .first_selected_component());

                                          Assert(jj < sub_neighbor->get_fe()
                                                        .n_components(),
                                                 ExcInternalError());

                                          if ((flux_mask(ii, jj) == always) ||
                                              (flux_mask(ii, jj) == nonzero))
                                            {
                                              sparsity.add(
                                                dofs_on_this_cell[i],
                                                dofs_on_other_cell[j]);
                                            }

                                          if ((flux_mask(jj, ii) == always) ||
                                              (flux_mask(jj, ii) == nonzero))
                                            {
                                              sparsity.add(
                                                dofs_on_other_cell[j],
                                                dofs_on_this_cell[i]);
                                            }
                                        }
                                    }
                                }
                            }
                          else
                            {
                              dofs_on_other_cell.resize(
                                neighbor->get_fe().n_dofs_per_cell());
                              neighbor->get_dof_indices(dofs_on_other_cell);
                              for (unsigned int i = 0;
                                   i < cell->get_fe().n_dofs_per_cell();
                                   ++i)
                                {
                                  const unsigned int ii =
                                    (cell->get_fe().is_primitive(i) ?
                                       cell->get_fe()
                                         .system_to_component_index(i)
                                         .first :
                                       cell->get_fe()
                                         .get_nonzero_components(i)
                                         .first_selected_component());

                                  Assert(ii < cell->get_fe().n_components(),
                                         ExcInternalError());

                                  for (unsigned int j = 0;
                                       j < neighbor->get_fe().n_dofs_per_cell();
                                       ++j)
                                    {
                                      const unsigned int jj =
                                        (neighbor->get_fe().is_primitive(j) ?
                                           neighbor->get_fe()
                                             .system_to_component_index(j)
                                             .first :
                                           neighbor->get_fe()
                                             .get_nonzero_components(j)
                                             .first_selected_component());

                                      Assert(
                                        jj < neighbor->get_fe().n_components(),
                                        ExcInternalError());

                                      if ((flux_mask(ii, jj) == always) ||
                                          (flux_mask(ii, jj) == nonzero))
                                        {
                                          sparsity.add(dofs_on_this_cell[i],
                                                       dofs_on_other_cell[j]);
                                        }

                                      if ((flux_mask(jj, ii) == always) ||
                                          (flux_mask(jj, ii) == nonzero))
                                        {
                                          sparsity.add(dofs_on_other_cell[j],
                                                       dofs_on_this_cell[i]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
          }
      }
    } // namespace

  } // namespace internal



  template <int dim, int spacedim, typename SparsityPatternType>
  void
  make_flux_sparsity_pattern(const DoFHandler<dim, spacedim> &dof,
                             SparsityPatternType &            sparsity,
                             const Table<2, Coupling> &       int_mask,
                             const Table<2, Coupling> &       flux_mask,
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



  template <int dim,
            int spacedim,
            typename SparsityPatternType,
            typename number>
  void
  make_flux_sparsity_pattern(
    const DoFHandler<dim, spacedim> &dof,
    SparsityPatternType &            sparsity,
    const AffineConstraints<number> &constraints,
    const bool                       keep_constrained_dofs,
    const Table<2, Coupling> &       int_mask,
    const Table<2, Coupling> &       flux_mask,
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

    // If we have a distributed::Triangulation only allow locally_owned
    // subdomain. Not setting a subdomain is also okay, because we skip
    // ghost cells in the loop below.
    Assert((dof.get_triangulation().locally_owned_subdomain() ==
            numbers::invalid_subdomain_id) ||
             (subdomain_id == numbers::invalid_subdomain_id) ||
             (subdomain_id ==
              dof.get_triangulation().locally_owned_subdomain()),
           ExcMessage(
             "For parallel::distributed::Triangulation objects and "
             "associated DoF handler objects, asking for any subdomain other "
             "than the locally owned one does not make sense."));

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

#include "dof_tools_sparsity.inst"



DEAL_II_NAMESPACE_CLOSE

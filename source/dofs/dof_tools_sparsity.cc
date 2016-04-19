// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2015 by the deal.II authors
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

#include <deal.II/base/thread_management.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>


#include <algorithm>
#include <numeric>

DEAL_II_NAMESPACE_OPEN



namespace DoFTools
{

  template <typename DoFHandlerType, typename SparsityPatternType>
  void
  make_sparsity_pattern (const DoFHandlerType      &dof,
                         SparsityPatternType       &sparsity,
                         const ConstraintMatrix    &constraints,
                         const bool                 keep_constrained_dofs,
                         const types::subdomain_id  subdomain_id)
  {
    const types::global_dof_index n_dofs = dof.n_dofs();
    (void)n_dofs;

    Assert (sparsity.n_rows() == n_dofs,
            ExcDimensionMismatch (sparsity.n_rows(), n_dofs));
    Assert (sparsity.n_cols() == n_dofs,
            ExcDimensionMismatch (sparsity.n_cols(), n_dofs));

    // If we have a distributed::Triangulation only allow locally_owned
    // subdomain. Not setting a subdomain is also okay, because we skip
    // ghost cells in the loop below.
    Assert (
      (dof.get_triangulation().locally_owned_subdomain() == numbers::invalid_subdomain_id)
      ||
      (subdomain_id == numbers::invalid_subdomain_id)
      ||
      (subdomain_id == dof.get_triangulation().locally_owned_subdomain()),
      ExcMessage ("For parallel::distributed::Triangulation objects and "
                  "associated DoF handler objects, asking for any subdomain other "
                  "than the locally owned one does not make sense."));

    std::vector<types::global_dof_index> dofs_on_this_cell;
    dofs_on_this_cell.reserve (max_dofs_per_cell(dof));
    typename DoFHandlerType::active_cell_iterator cell = dof.begin_active(),
                                                  endc = dof.end();

    // In case we work with a distributed sparsity pattern of Trilinos
    // type, we only have to do the work if the current cell is owned by
    // the calling processor. Otherwise, just continue.
    for (; cell!=endc; ++cell)
      if (((subdomain_id == numbers::invalid_subdomain_id)
           ||
           (subdomain_id == cell->subdomain_id()))
          &&
          cell->is_locally_owned())
        {
          const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
          dofs_on_this_cell.resize (dofs_per_cell);
          cell->get_dof_indices (dofs_on_this_cell);

          // make sparsity pattern for this cell. if no constraints pattern
          // was given, then the following call acts as if simply no
          // constraints existed
          constraints.add_entries_local_to_global (dofs_on_this_cell,
                                                   sparsity,
                                                   keep_constrained_dofs);
        }
  }



  template <typename DoFHandlerType, typename SparsityPatternType>
  void
  make_sparsity_pattern (const DoFHandlerType      &dof,
                         const Table<2,Coupling>   &couplings,
                         SparsityPatternType       &sparsity,
                         const ConstraintMatrix    &constraints,
                         const bool                 keep_constrained_dofs,
                         const types::subdomain_id  subdomain_id)
  {
    const types::global_dof_index n_dofs = dof.n_dofs();
    (void)n_dofs;

    Assert (sparsity.n_rows() == n_dofs,
            ExcDimensionMismatch (sparsity.n_rows(), n_dofs));
    Assert (sparsity.n_cols() == n_dofs,
            ExcDimensionMismatch (sparsity.n_cols(), n_dofs));
    Assert (couplings.n_rows() == dof.get_fe().n_components(),
            ExcDimensionMismatch(couplings.n_rows(), dof.get_fe().n_components()));
    Assert (couplings.n_cols() == dof.get_fe().n_components(),
            ExcDimensionMismatch(couplings.n_cols(), dof.get_fe().n_components()));

    // If we have a distributed::Triangulation only allow locally_owned
    // subdomain. Not setting a subdomain is also okay, because we skip
    // ghost cells in the loop below.
    Assert (
      (dof.get_triangulation().locally_owned_subdomain() == numbers::invalid_subdomain_id)
      ||
      (subdomain_id == numbers::invalid_subdomain_id)
      ||
      (subdomain_id == dof.get_triangulation().locally_owned_subdomain()),
      ExcMessage ("For parallel::distributed::Triangulation objects and "
                  "associated DoF handler objects, asking for any subdomain other "
                  "than the locally owned one does not make sense."));

    const hp::FECollection<DoFHandlerType::dimension,DoFHandlerType::space_dimension> fe_collection (dof.get_fe());

    // first, for each finite element, build a mask for each dof, not like
    // the one given which represents components. make sure we do the right
    // thing also with respect to non-primitive shape functions, which
    // takes some additional thought
    std::vector<Table<2,bool> > dof_mask(fe_collection.size());

    // check whether the table of couplings contains only true arguments,
    // i.e., we do not exclude any index. that is the easy case, since we
    // don't have to set up the tables
    bool need_dof_mask = false;
    for (unsigned int i=0; i<couplings.n_rows(); ++i)
      for (unsigned int j=0; j<couplings.n_cols(); ++j)
        if (couplings(i,j) == none)
          need_dof_mask = true;

    if (need_dof_mask == true)
      for (unsigned int f=0; f<fe_collection.size(); ++f)
        {
          const unsigned int dofs_per_cell = fe_collection[f].dofs_per_cell;

          dof_mask[f].reinit (dofs_per_cell, dofs_per_cell);

          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              if (fe_collection[f].is_primitive(i) &&
                  fe_collection[f].is_primitive(j))
                dof_mask[f](i,j)
                  = (couplings(fe_collection[f].system_to_component_index(i).first,
                               fe_collection[f].system_to_component_index(j).first) != none);
              else
                {
                  const unsigned int first_nonzero_comp_i
                    = fe_collection[f].get_nonzero_components(i).first_selected_component();
                  const unsigned int first_nonzero_comp_j
                    = fe_collection[f].get_nonzero_components(j).first_selected_component();
                  Assert (first_nonzero_comp_i < fe_collection[f].n_components(),
                          ExcInternalError());
                  Assert (first_nonzero_comp_j < fe_collection[f].n_components(),
                          ExcInternalError());

                  dof_mask[f](i,j)
                    = (couplings(first_nonzero_comp_i,first_nonzero_comp_j) != none);
                }
        }


    std::vector<types::global_dof_index> dofs_on_this_cell(fe_collection.max_dofs_per_cell());
    typename DoFHandlerType::active_cell_iterator cell = dof.begin_active(),
                                                  endc = dof.end();

    // In case we work with a distributed sparsity pattern of Trilinos
    // type, we only have to do the work if the current cell is owned by
    // the calling processor. Otherwise, just continue.
    for (; cell!=endc; ++cell)
      if (((subdomain_id == numbers::invalid_subdomain_id)
           ||
           (subdomain_id == cell->subdomain_id()))
          &&
          cell->is_locally_owned())
        {
          const unsigned int fe_index = cell->active_fe_index();
          const unsigned int dofs_per_cell =fe_collection[fe_index].dofs_per_cell;

          dofs_on_this_cell.resize (dofs_per_cell);
          cell->get_dof_indices (dofs_on_this_cell);


          // make sparsity pattern for this cell. if no constraints pattern
          // was given, then the following call acts as if simply no
          // constraints existed
          constraints.add_entries_local_to_global (dofs_on_this_cell,
                                                   sparsity,
                                                   keep_constrained_dofs,
                                                   dof_mask[fe_index]);
        }
  }



  template <typename DoFHandlerType, typename SparsityPatternType>
  void
  make_sparsity_pattern (const DoFHandlerType &dof_row,
                         const DoFHandlerType &dof_col,
                         SparsityPatternType  &sparsity)
  {
    const types::global_dof_index n_dofs_row = dof_row.n_dofs();
    const types::global_dof_index n_dofs_col = dof_col.n_dofs();
    (void)n_dofs_row;
    (void)n_dofs_col;

    Assert (sparsity.n_rows() == n_dofs_row,
            ExcDimensionMismatch (sparsity.n_rows(), n_dofs_row));
    Assert (sparsity.n_cols() == n_dofs_col,
            ExcDimensionMismatch (sparsity.n_cols(), n_dofs_col));

//TODO: Looks like wasteful memory management here

    const std::list<std::pair<typename DoFHandlerType::cell_iterator,
          typename DoFHandlerType::cell_iterator> >
          cell_list
          = GridTools::get_finest_common_cells (dof_row, dof_col);


    typename std::list<std::pair<typename DoFHandlerType::cell_iterator,
             typename DoFHandlerType::cell_iterator> >::const_iterator
             cell_iter = cell_list.begin();

    for (; cell_iter!=cell_list.end(); ++cell_iter)
      {
        const typename DoFHandlerType::cell_iterator cell_row = cell_iter->first;
        const typename DoFHandlerType::cell_iterator cell_col = cell_iter->second;

        if (!cell_row->has_children() && !cell_col->has_children())
          {
            const unsigned int dofs_per_cell_row =
              cell_row->get_fe().dofs_per_cell;
            const unsigned int dofs_per_cell_col =
              cell_col->get_fe().dofs_per_cell;
            std::vector<types::global_dof_index>
            local_dof_indices_row(dofs_per_cell_row);
            std::vector<types::global_dof_index>
            local_dof_indices_col(dofs_per_cell_col);
            cell_row->get_dof_indices (local_dof_indices_row);
            cell_col->get_dof_indices (local_dof_indices_col);
            for (unsigned int i=0; i<dofs_per_cell_row; ++i)
              sparsity.add_entries (local_dof_indices_row[i],
                                    local_dof_indices_col.begin(),
                                    local_dof_indices_col.end());
          }
        else if (cell_row->has_children())
          {
            const std::vector<typename DoFHandlerType::active_cell_iterator >
            child_cells = GridTools::get_active_child_cells<DoFHandlerType> (cell_row);
            for (unsigned int i=0; i<child_cells.size(); i++)
              {
                const typename DoFHandlerType::cell_iterator
                cell_row_child = child_cells[i];
                const unsigned int dofs_per_cell_row =
                  cell_row_child->get_fe().dofs_per_cell;
                const unsigned int dofs_per_cell_col =
                  cell_col->get_fe().dofs_per_cell;
                std::vector<types::global_dof_index>
                local_dof_indices_row(dofs_per_cell_row);
                std::vector<types::global_dof_index>
                local_dof_indices_col(dofs_per_cell_col);
                cell_row_child->get_dof_indices (local_dof_indices_row);
                cell_col->get_dof_indices (local_dof_indices_col);
                for (unsigned int r=0; r<dofs_per_cell_row; ++r)
                  sparsity.add_entries (local_dof_indices_row[r],
                                        local_dof_indices_col.begin(),
                                        local_dof_indices_col.end());
              }
          }
        else
          {
            std::vector<typename DoFHandlerType::active_cell_iterator>
            child_cells = GridTools::get_active_child_cells<DoFHandlerType> (cell_col);
            for (unsigned int i=0; i<child_cells.size(); i++)
              {
                const typename DoFHandlerType::active_cell_iterator
                cell_col_child = child_cells[i];
                const unsigned int dofs_per_cell_row =
                  cell_row->get_fe().dofs_per_cell;
                const unsigned int dofs_per_cell_col =
                  cell_col_child->get_fe().dofs_per_cell;
                std::vector<types::global_dof_index>
                local_dof_indices_row(dofs_per_cell_row);
                std::vector<types::global_dof_index>
                local_dof_indices_col(dofs_per_cell_col);
                cell_row->get_dof_indices (local_dof_indices_row);
                cell_col_child->get_dof_indices (local_dof_indices_col);
                for (unsigned int r=0; r<dofs_per_cell_row; ++r)
                  sparsity.add_entries (local_dof_indices_row[r],
                                        local_dof_indices_col.begin(),
                                        local_dof_indices_col.end());
              }
          }
      }
  }



  template <typename DoFHandlerType, typename SparsityPatternType>
  void
  make_boundary_sparsity_pattern
  (const DoFHandlerType                       &dof,
   const std::vector<types::global_dof_index> &dof_to_boundary_mapping,
   SparsityPatternType                        &sparsity)
  {
    if (DoFHandlerType::dimension == 1)
      {
        // there are only 2 boundary indicators in 1d, so it is no
        // performance problem to call the other function
        std::map<types::boundary_id, const Function<DoFHandlerType::space_dimension,double>*> boundary_ids;
        boundary_ids[0] = 0;
        boundary_ids[1] = 0;
        make_boundary_sparsity_pattern<DoFHandlerType, SparsityPatternType>
        (dof,
         boundary_ids,
         dof_to_boundary_mapping,
         sparsity);
        return;
      }

    const types::global_dof_index n_dofs = dof.n_dofs();
    (void)n_dofs;

    AssertDimension (dof_to_boundary_mapping.size(), n_dofs);
    AssertDimension (sparsity.n_rows(), dof.n_boundary_dofs());
    AssertDimension (sparsity.n_cols(), dof.n_boundary_dofs());
#ifdef DEBUG
    if (sparsity.n_rows() != 0)
      {
        types::global_dof_index max_element = 0;
        for (std::vector<types::global_dof_index>::const_iterator i=dof_to_boundary_mapping.begin();
             i!=dof_to_boundary_mapping.end(); ++i)
          if ((*i != DoFHandlerType::invalid_dof_index) &&
              (*i > max_element))
            max_element = *i;
        AssertDimension (max_element, sparsity.n_rows()-1);
      };
#endif

    std::vector<types::global_dof_index> dofs_on_this_face;
    dofs_on_this_face.reserve (max_dofs_per_face(dof));

    // loop over all faces to check whether they are at a boundary. note
    // that we need not take special care of single lines (using
    // @p{cell->has_boundary_lines}), since we do not support boundaries of
    // dimension dim-2, and so every boundary line is also part of a
    // boundary face.
    typename DoFHandlerType::active_cell_iterator cell = dof.begin_active(),
                                                  endc = dof.end();
    for (; cell!=endc; ++cell)
      for (unsigned int f=0; f<GeometryInfo<DoFHandlerType::dimension>::faces_per_cell; ++f)
        if (cell->at_boundary(f))
          {
            const unsigned int dofs_per_face = cell->get_fe().dofs_per_face;
            dofs_on_this_face.resize (dofs_per_face);
            cell->face(f)->get_dof_indices (dofs_on_this_face,
                                            cell->active_fe_index());

            // make sparsity pattern for this cell
            for (unsigned int i=0; i<dofs_per_face; ++i)
              for (unsigned int j=0; j<dofs_per_face; ++j)
                sparsity.add (dof_to_boundary_mapping[dofs_on_this_face[i]],
                              dof_to_boundary_mapping[dofs_on_this_face[j]]);
          }
  }



  template <typename DoFHandlerType, typename SparsityPatternType, typename number>
  void make_boundary_sparsity_pattern
  (const DoFHandlerType                                              &dof,
   const std::map<types::boundary_id, const Function<DoFHandlerType::space_dimension,number>*> &boundary_ids,
   const std::vector<types::global_dof_index>                        &dof_to_boundary_mapping,
   SparsityPatternType                                               &sparsity)
  {
    if (DoFHandlerType::dimension == 1)
      {
        // first check left, then right boundary point
        for (unsigned int direction=0; direction<2; ++direction)
          {
            // if this boundary is not requested, then go on with next one
            if (boundary_ids.find(direction) ==
                boundary_ids.end())
              continue;

            // find active cell at that boundary: first go to left/right,
            // then to children
            typename DoFHandlerType::level_cell_iterator cell = dof.begin(0);
            while (!cell->at_boundary(direction))
              cell = cell->neighbor(direction);
            while (!cell->active())
              cell = cell->child(direction);

            const unsigned int dofs_per_vertex = cell->get_fe().dofs_per_vertex;
            std::vector<types::global_dof_index> boundary_dof_boundary_indices (dofs_per_vertex);

            // next get boundary mapped dof indices of boundary dofs
            for (unsigned int i=0; i<dofs_per_vertex; ++i)
              boundary_dof_boundary_indices[i]
                = dof_to_boundary_mapping[cell->vertex_dof_index(direction,i)];

            for (unsigned int i=0; i<dofs_per_vertex; ++i)
              sparsity.add_entries (boundary_dof_boundary_indices[i],
                                    boundary_dof_boundary_indices.begin(),
                                    boundary_dof_boundary_indices.end());
          };
        return;
      }

    const types::global_dof_index n_dofs = dof.n_dofs();
    (void)n_dofs;

    AssertDimension (dof_to_boundary_mapping.size(), n_dofs);
    Assert (boundary_ids.find(numbers::internal_face_boundary_id) == boundary_ids.end(),
            typename DoFHandlerType::ExcInvalidBoundaryIndicator());
    Assert (sparsity.n_rows() == dof.n_boundary_dofs (boundary_ids),
            ExcDimensionMismatch (sparsity.n_rows(), dof.n_boundary_dofs (boundary_ids)));
    Assert (sparsity.n_cols() == dof.n_boundary_dofs (boundary_ids),
            ExcDimensionMismatch (sparsity.n_cols(), dof.n_boundary_dofs (boundary_ids)));
#ifdef DEBUG
    if (sparsity.n_rows() != 0)
      {
        types::global_dof_index max_element = 0;
        for (std::vector<types::global_dof_index>::const_iterator i=dof_to_boundary_mapping.begin();
             i!=dof_to_boundary_mapping.end(); ++i)
          if ((*i != DoFHandlerType::invalid_dof_index) &&
              (*i > max_element))
            max_element = *i;
        AssertDimension (max_element, sparsity.n_rows()-1);
      };
#endif

    std::vector<types::global_dof_index> dofs_on_this_face;
    dofs_on_this_face.reserve (max_dofs_per_face(dof));
    typename DoFHandlerType::active_cell_iterator cell = dof.begin_active(),
                                                  endc = dof.end();
    for (; cell!=endc; ++cell)
      for (unsigned int f=0; f<GeometryInfo<DoFHandlerType::dimension>::faces_per_cell; ++f)
        if (boundary_ids.find(cell->face(f)->boundary_id()) !=
            boundary_ids.end())
          {
            const unsigned int dofs_per_face = cell->get_fe().dofs_per_face;
            dofs_on_this_face.resize (dofs_per_face);
            cell->face(f)->get_dof_indices (dofs_on_this_face,
                                            cell->active_fe_index());

            // make sparsity pattern for this cell
            for (unsigned int i=0; i<dofs_per_face; ++i)
              for (unsigned int j=0; j<dofs_per_face; ++j)
                sparsity.add (dof_to_boundary_mapping[dofs_on_this_face[i]],
                              dof_to_boundary_mapping[dofs_on_this_face[j]]);
          }
  }



  template <typename DoFHandlerType, typename SparsityPatternType>
  void
  make_flux_sparsity_pattern (const DoFHandlerType      &dof,
                              SparsityPatternType       &sparsity,
                              const ConstraintMatrix    &constraints,
                              const bool                 keep_constrained_dofs,
                              const types::subdomain_id  subdomain_id)

  // TODO: QA: reduce the indentation level of this method..., Maier 2012

  {
    const types::global_dof_index n_dofs = dof.n_dofs();
    (void)n_dofs;

    AssertDimension (sparsity.n_rows(), n_dofs);
    AssertDimension (sparsity.n_cols(), n_dofs);

    // If we have a distributed::Triangulation only allow locally_owned
    // subdomain. Not setting a subdomain is also okay, because we skip
    // ghost cells in the loop below.
    Assert (
      (dof.get_triangulation().locally_owned_subdomain() == numbers::invalid_subdomain_id)
      ||
      (subdomain_id == numbers::invalid_subdomain_id)
      ||
      (subdomain_id == dof.get_triangulation().locally_owned_subdomain()),
      ExcMessage ("For parallel::distributed::Triangulation objects and "
                  "associated DoF handler objects, asking for any subdomain other "
                  "than the locally owned one does not make sense."));

    std::vector<types::global_dof_index> dofs_on_this_cell;
    std::vector<types::global_dof_index> dofs_on_other_cell;
    dofs_on_this_cell.reserve (max_dofs_per_cell(dof));
    dofs_on_other_cell.reserve (max_dofs_per_cell(dof));
    typename DoFHandlerType::active_cell_iterator cell = dof.begin_active(),
                                                  endc = dof.end();

    // TODO: in an old implementation, we used user flags before to tag
    // faces that were already touched. this way, we could reduce the work
    // a little bit. now, we instead add only data from one side. this
    // should be OK, but we need to actually verify it.

    // In case we work with a distributed sparsity pattern of Trilinos
    // type, we only have to do the work if the current cell is owned by
    // the calling processor. Otherwise, just continue.
    for (; cell!=endc; ++cell)
      if (((subdomain_id == numbers::invalid_subdomain_id)
           ||
           (subdomain_id == cell->subdomain_id()))
          &&
          cell->is_locally_owned())
        {
          const unsigned int n_dofs_on_this_cell = cell->get_fe().dofs_per_cell;
          dofs_on_this_cell.resize (n_dofs_on_this_cell);
          cell->get_dof_indices (dofs_on_this_cell);

          // make sparsity pattern for this cell. if no constraints pattern
          // was given, then the following call acts as if simply no
          // constraints existed
          constraints.add_entries_local_to_global (dofs_on_this_cell,
                                                   sparsity,
                                                   keep_constrained_dofs);

          for (unsigned int face = 0;
               face < GeometryInfo<DoFHandlerType::dimension>::faces_per_cell;
               ++face)
            {
              typename DoFHandlerType::face_iterator cell_face = cell->face(face);
              const bool periodic_neighbor = cell->has_periodic_neighbor(face);
              if (! cell->at_boundary(face) || periodic_neighbor)
                {
                  typename DoFHandlerType::level_cell_iterator neighbor
                    = cell->neighbor_or_periodic_neighbor(face);

                  // in 1d, we do not need to worry whether the neighbor
                  // might have children and then loop over those children.
                  // rather, we may as well go straight to the cell behind
                  // this particular cell's most terminal child
                  if (DoFHandlerType::dimension==1)
                    while (neighbor->has_children())
                      neighbor = neighbor->child(face==0 ? 1 : 0);

                  if (neighbor->has_children())
                    {
                      for (unsigned int sub_nr = 0;
                           sub_nr != cell_face->number_of_children();
                           ++sub_nr)
                        {
                          const typename DoFHandlerType::level_cell_iterator sub_neighbor
                            = periodic_neighbor?
                              cell->periodic_neighbor_child_on_subface (face, sub_nr):
                              cell->neighbor_child_on_subface (face, sub_nr);

                          const unsigned int n_dofs_on_neighbor
                            = sub_neighbor->get_fe().dofs_per_cell;
                          dofs_on_other_cell.resize (n_dofs_on_neighbor);
                          sub_neighbor->get_dof_indices (dofs_on_other_cell);

                          constraints.add_entries_local_to_global
                          (dofs_on_this_cell, dofs_on_other_cell,
                           sparsity, keep_constrained_dofs);
                          constraints.add_entries_local_to_global
                          (dofs_on_other_cell, dofs_on_this_cell,
                           sparsity, keep_constrained_dofs);
                          // only need to add this when the neighbor is not
                          // owned by the current processor, otherwise we add
                          // the entries for the neighbor there
                          if (sub_neighbor->subdomain_id() != cell->subdomain_id())
                            constraints.add_entries_local_to_global
                            (dofs_on_other_cell, sparsity, keep_constrained_dofs);
                        }
                    }
                  else
                    {
                      // Refinement edges are taken care of by coarser
                      // cells
                      if ((!periodic_neighbor && cell->neighbor_is_coarser(face)) ||
                          (periodic_neighbor && cell->periodic_neighbor_is_coarser(face)))
                        if (neighbor->subdomain_id() == cell->subdomain_id())
                          continue;

                      const unsigned int n_dofs_on_neighbor
                        = neighbor->get_fe().dofs_per_cell;
                      dofs_on_other_cell.resize (n_dofs_on_neighbor);

                      neighbor->get_dof_indices (dofs_on_other_cell);

                      constraints.add_entries_local_to_global
                      (dofs_on_this_cell, dofs_on_other_cell,
                       sparsity, keep_constrained_dofs);

                      // only need to add these in case the neighbor cell
                      // is not locally owned - otherwise, we touch each
                      // face twice and hence put the indices the other way
                      // around
                      if (!cell->neighbor_or_periodic_neighbor(face)->active()
                          ||
                          (neighbor->subdomain_id() != cell->subdomain_id()))
                        {
                          constraints.add_entries_local_to_global
                          (dofs_on_other_cell, dofs_on_this_cell,
                           sparsity, keep_constrained_dofs);
                          if (neighbor->subdomain_id() != cell->subdomain_id())
                            constraints.add_entries_local_to_global
                            (dofs_on_other_cell, sparsity, keep_constrained_dofs);
                        }
                    }
                }
            }
        }
  }



  template <typename DoFHandlerType, typename SparsityPatternType>
  void
  make_flux_sparsity_pattern (const DoFHandlerType &dof,
                              SparsityPatternType  &sparsity)
  {
    ConstraintMatrix constraints;
    make_flux_sparsity_pattern (dof, sparsity, constraints);
  }

  template <int dim, int spacedim>
  Table<2,Coupling>
  dof_couplings_from_component_couplings (const FiniteElement<dim,spacedim> &fe,
                                          const Table<2,Coupling> &component_couplings)
  {
    Assert(component_couplings.n_rows() == fe.n_components(),
           ExcDimensionMismatch(component_couplings.n_rows(),
                                fe.n_components()));
    Assert(component_couplings.n_cols() == fe.n_components(),
           ExcDimensionMismatch(component_couplings.n_cols(),
                                fe.n_components()));

    const unsigned int n_dofs = fe.dofs_per_cell;

    Table<2,Coupling> dof_couplings (n_dofs, n_dofs);

    for (unsigned int i=0; i<n_dofs; ++i)
      {
        const unsigned int ii
          = (fe.is_primitive(i) ?
             fe.system_to_component_index(i).first
             :
             fe.get_nonzero_components(i).first_selected_component()
            );
        Assert (ii < fe.n_components(), ExcInternalError());

        for (unsigned int j=0; j<n_dofs; ++j)
          {
            const unsigned int jj
              = (fe.is_primitive(j) ?
                 fe.system_to_component_index(j).first
                 :
                 fe.get_nonzero_components(j).first_selected_component()
                );
            Assert (jj < fe.n_components(), ExcInternalError());

            dof_couplings(i,j) = component_couplings(ii,jj);
          }
      }
    return dof_couplings;
  }



  template <int dim, int spacedim>
  std::vector<Table<2,Coupling> >
  dof_couplings_from_component_couplings
  (const hp::FECollection<dim,spacedim> &fe,
   const Table<2,Coupling> &component_couplings)
  {
    std::vector<Table<2,Coupling> > return_value (fe.size());
    for (unsigned int i=0; i<fe.size(); ++i)
      return_value[i]
        = dof_couplings_from_component_couplings(fe[i], component_couplings);

    return return_value;
  }



  namespace internal
  {
    namespace
    {

      // implementation of the same function in namespace DoFTools for
      // non-hp DoFHandlers
      template <typename DoFHandlerType, typename SparsityPatternType>
      void
      make_flux_sparsity_pattern (const DoFHandlerType    &dof,
                                  SparsityPatternType     &sparsity,
                                  const Table<2,Coupling> &int_mask,
                                  const Table<2,Coupling> &flux_mask)
      {
        const FiniteElement<DoFHandlerType::dimension,DoFHandlerType::space_dimension>
        &fe = dof.get_fe();

        std::vector<types::global_dof_index> dofs_on_this_cell (fe.dofs_per_cell);
        std::vector<types::global_dof_index> dofs_on_other_cell (fe.dofs_per_cell);

        const Table<2,Coupling>
        int_dof_mask  = dof_couplings_from_component_couplings (fe, int_mask),
        flux_dof_mask = dof_couplings_from_component_couplings (fe, flux_mask);

        Table<2,bool> support_on_face
        (fe.dofs_per_cell, GeometryInfo<DoFHandlerType::dimension>::faces_per_cell);
        for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
          for (unsigned int f=0; f<GeometryInfo<DoFHandlerType::dimension>::faces_per_cell; ++f)
            support_on_face(i,f) = fe.has_support_on_face(i,f);

        typename DoFHandlerType::active_cell_iterator cell = dof.begin_active (),
                                                      endc = dof.end ();
        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned ())
            {
              cell->get_dof_indices (dofs_on_this_cell);
              // make sparsity pattern for this cell
              for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
                for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
                  if (int_dof_mask (i,j) != none)
                    sparsity.add (dofs_on_this_cell[i],
                                  dofs_on_this_cell[j]);

              // Loop over all interior neighbors
              for (unsigned int face_n = 0;
                   face_n < GeometryInfo<DoFHandlerType::dimension>::faces_per_cell;
                   ++face_n)
                {
                  const typename DoFHandlerType::face_iterator
                  cell_face = cell->face (face_n);

                  const bool periodic_neighbor = cell->has_periodic_neighbor (face_n);

                  if (cell->at_boundary (face_n) && (!periodic_neighbor))
                    {
                      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
                        {
                          const bool i_non_zero_i = support_on_face (i, face_n);
                          for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
                            {
                              const bool j_non_zero_i = support_on_face (j, face_n);

                              if (flux_dof_mask (i,j) == always
                                  ||
                                  (flux_dof_mask (i,j) == nonzero
                                   &&
                                   i_non_zero_i
                                   &&
                                   j_non_zero_i))
                                sparsity.add (dofs_on_this_cell[i],
                                              dofs_on_this_cell[j]);
                            }
                        }
                    }
                  else
                    {
                      typename DoFHandlerType::level_cell_iterator
                      neighbor = cell->neighbor_or_periodic_neighbor (face_n);
                      // If the cells are on the same level then only add to
                      // the sparsity pattern if the current cell is 'greater'
                      // in the total ordering.
                      if (neighbor->level () == cell->level ()
                          &&
                          neighbor->index () > cell->index ())
                        continue;
                      // If we are more refined then the neighbor, then we
                      // will automatically find the active neighbor cell when
                      // we call 'neighbor (face_n)' above. The opposite is
                      // not true; if the neighbor is more refined then the
                      // call 'neighbor (face_n)' will *not* return an active
                      // cell. Hence, only add things to the sparsity pattern
                      // if the neighbor is coarser than the current cell.
                      if (neighbor->level () != cell->level ()
                          &&
                          ((!periodic_neighbor && !cell->neighbor_is_coarser (face_n))
                           ||(periodic_neighbor && !cell->periodic_neighbor_is_coarser (face_n))))
                        continue; // (the neighbor is finer)

                      const unsigned int
                      neighbor_face_n = cell->neighbor_face_no (face_n);

                      if (cell_face->has_children ())
                        {
                          for (unsigned int sub_nr = 0;
                               sub_nr != cell_face->n_children ();
                               ++sub_nr)
                            {
                              const typename DoFHandlerType::level_cell_iterator
                              sub_neighbor
                                = periodic_neighbor?
                                  cell->periodic_neighbor_child_on_subface (face_n, sub_nr):
                                  cell->neighbor_child_on_subface (face_n, sub_nr);

                              sub_neighbor->get_dof_indices (dofs_on_other_cell);
                              for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
                                {
                                  const bool i_non_zero_i = support_on_face (i, face_n);
                                  const bool i_non_zero_e = support_on_face (i, neighbor_face_n);
                                  for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
                                    {
                                      const bool j_non_zero_i = support_on_face (j, face_n);
                                      const bool j_non_zero_e = support_on_face (j, neighbor_face_n);

                                      if (flux_dof_mask (i,j) == always)
                                        {
                                          sparsity.add (dofs_on_this_cell[i],
                                                        dofs_on_other_cell[j]);
                                          sparsity.add (dofs_on_other_cell[i],
                                                        dofs_on_this_cell[j]);
                                          sparsity.add (dofs_on_this_cell[i],
                                                        dofs_on_this_cell[j]);
                                          sparsity.add (dofs_on_other_cell[i],
                                                        dofs_on_other_cell[j]);
                                        }
                                      else if (flux_dof_mask (i,j) == nonzero)
                                        {
                                          if (i_non_zero_i && j_non_zero_e)
                                            sparsity.add (dofs_on_this_cell[i],
                                                          dofs_on_other_cell[j]);
                                          if (i_non_zero_e && j_non_zero_i)
                                            sparsity.add (dofs_on_other_cell[i],
                                                          dofs_on_this_cell[j]);
                                          if (i_non_zero_i && j_non_zero_i)
                                            sparsity.add (dofs_on_this_cell[i],
                                                          dofs_on_this_cell[j]);
                                          if (i_non_zero_e && j_non_zero_e)
                                            sparsity.add (dofs_on_other_cell[i],
                                                          dofs_on_other_cell[j]);
                                        }

                                      if (flux_dof_mask (j,i) == always)
                                        {
                                          sparsity.add (dofs_on_this_cell[j],
                                                        dofs_on_other_cell[i]);
                                          sparsity.add (dofs_on_other_cell[j],
                                                        dofs_on_this_cell[i]);
                                          sparsity.add (dofs_on_this_cell[j],
                                                        dofs_on_this_cell[i]);
                                          sparsity.add (dofs_on_other_cell[j],
                                                        dofs_on_other_cell[i]);
                                        }
                                      else if (flux_dof_mask (j,i) == nonzero)
                                        {
                                          if (j_non_zero_i && i_non_zero_e)
                                            sparsity.add (dofs_on_this_cell[j],
                                                          dofs_on_other_cell[i]);
                                          if (j_non_zero_e && i_non_zero_i)
                                            sparsity.add (dofs_on_other_cell[j],
                                                          dofs_on_this_cell[i]);
                                          if (j_non_zero_i && i_non_zero_i)
                                            sparsity.add (dofs_on_this_cell[j],
                                                          dofs_on_this_cell[i]);
                                          if (j_non_zero_e && i_non_zero_e)
                                            sparsity.add (dofs_on_other_cell[j],
                                                          dofs_on_other_cell[i]);
                                        }
                                    }
                                }
                            }
                        }
                      else
                        {
                          neighbor->get_dof_indices (dofs_on_other_cell);
                          for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
                            {
                              const bool i_non_zero_i = support_on_face (i, face_n);
                              const bool i_non_zero_e = support_on_face (i, neighbor_face_n);
                              for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
                                {
                                  const bool j_non_zero_i = support_on_face (j, face_n);
                                  const bool j_non_zero_e = support_on_face (j, neighbor_face_n);
                                  if (flux_dof_mask (i,j) == always)
                                    {
                                      sparsity.add (dofs_on_this_cell[i],
                                                    dofs_on_other_cell[j]);
                                      sparsity.add (dofs_on_other_cell[i],
                                                    dofs_on_this_cell[j]);
                                      sparsity.add (dofs_on_this_cell[i],
                                                    dofs_on_this_cell[j]);
                                      sparsity.add (dofs_on_other_cell[i],
                                                    dofs_on_other_cell[j]);
                                    }
                                  if (flux_dof_mask (i,j) == nonzero)
                                    {
                                      if (i_non_zero_i && j_non_zero_e)
                                        sparsity.add (dofs_on_this_cell[i],
                                                      dofs_on_other_cell[j]);
                                      if (i_non_zero_e && j_non_zero_i)
                                        sparsity.add (dofs_on_other_cell[i],
                                                      dofs_on_this_cell[j]);
                                      if (i_non_zero_i && j_non_zero_i)
                                        sparsity.add (dofs_on_this_cell[i],
                                                      dofs_on_this_cell[j]);
                                      if (i_non_zero_e && j_non_zero_e)
                                        sparsity.add (dofs_on_other_cell[i],
                                                      dofs_on_other_cell[j]);
                                    }

                                  if (flux_dof_mask (j,i) == always)
                                    {
                                      sparsity.add (dofs_on_this_cell[j],
                                                    dofs_on_other_cell[i]);
                                      sparsity.add (dofs_on_other_cell[j],
                                                    dofs_on_this_cell[i]);
                                      sparsity.add (dofs_on_this_cell[j],
                                                    dofs_on_this_cell[i]);
                                      sparsity.add (dofs_on_other_cell[j],
                                                    dofs_on_other_cell[i]);
                                    }
                                  if (flux_dof_mask (j,i) == nonzero)
                                    {
                                      if (j_non_zero_i && i_non_zero_e)
                                        sparsity.add (dofs_on_this_cell[j],
                                                      dofs_on_other_cell[i]);
                                      if (j_non_zero_e && i_non_zero_i)
                                        sparsity.add (dofs_on_other_cell[j],
                                                      dofs_on_this_cell[i]);
                                      if (j_non_zero_i && i_non_zero_i)
                                        sparsity.add (dofs_on_this_cell[j],
                                                      dofs_on_this_cell[i]);
                                      if (j_non_zero_e && i_non_zero_e)
                                        sparsity.add (dofs_on_other_cell[j],
                                                      dofs_on_other_cell[i]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
      }


      // implementation of the same function in namespace DoFTools for hp
      // DoFHandlers
      template <int dim, int spacedim, typename SparsityPatternType>
      void
      make_flux_sparsity_pattern (const dealii::hp::DoFHandler<dim,spacedim> &dof,
                                  SparsityPatternType                        &sparsity,
                                  const Table<2,Coupling>                    &int_mask,
                                  const Table<2,Coupling>                    &flux_mask)
      {
        // while the implementation above is quite optimized and caches a
        // lot of data (see e.g. the int/flux_dof_mask tables), this is no
        // longer practical for the hp version since we would have to have
        // it for all combinations of elements in the hp::FECollection.
        // consequently, the implementation here is simpler and probably
        // less efficient but at least readable...

        const dealii::hp::FECollection<dim,spacedim> &fe = dof.get_fe();

        std::vector<types::global_dof_index> dofs_on_this_cell(DoFTools::max_dofs_per_cell(dof));
        std::vector<types::global_dof_index> dofs_on_other_cell(DoFTools::max_dofs_per_cell(dof));

        const std::vector<Table<2,Coupling> >
        int_dof_mask
          = dof_couplings_from_component_couplings(fe, int_mask);

        typename dealii::hp::DoFHandler<dim,spacedim>::active_cell_iterator
        cell = dof.begin_active(),
        endc = dof.end();
        for (; cell!=endc; ++cell)
          {
            dofs_on_this_cell.resize (cell->get_fe().dofs_per_cell);
            cell->get_dof_indices (dofs_on_this_cell);

            // make sparsity pattern for this cell
            for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
              for (unsigned int j=0; j<cell->get_fe().dofs_per_cell; ++j)
                if (int_dof_mask[cell->active_fe_index()](i,j) != none)
                  sparsity.add (dofs_on_this_cell[i],
                                dofs_on_this_cell[j]);

            // Loop over all interior neighbors
            for (unsigned int face = 0;
                 face < GeometryInfo<dim>::faces_per_cell;
                 ++face)
              {
                const typename dealii::hp::DoFHandler<dim,spacedim>::face_iterator
                cell_face = cell->face(face);
                if (cell_face->user_flag_set ())
                  continue;

                const bool periodic_neighbor = cell->has_periodic_neighbor (face);

                if (cell->at_boundary (face) && (!periodic_neighbor))
                  {
                    for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
                      for (unsigned int j=0; j<cell->get_fe().dofs_per_cell; ++j)
                        if ((flux_mask(cell->get_fe().system_to_component_index(i).first,
                                       cell->get_fe().system_to_component_index(j).first)
                             == always)
                            ||
                            (flux_mask(cell->get_fe().system_to_component_index(i).first,
                                       cell->get_fe().system_to_component_index(j).first)
                             == nonzero))
                          sparsity.add (dofs_on_this_cell[i],
                                        dofs_on_this_cell[j]);
                  }
                else
                  {
                    typename dealii::hp::DoFHandler<dim,spacedim>::level_cell_iterator
                    neighbor = cell->neighbor_or_periodic_neighbor(face);

                    // Refinement edges are taken care of by coarser cells
                    if ((!periodic_neighbor && cell->neighbor_is_coarser(face)) ||
                        (periodic_neighbor && cell->periodic_neighbor_is_coarser(face)))
                      continue;

                    const unsigned int
                    neighbor_face = periodic_neighbor?
                                    cell->periodic_neighbor_of_periodic_neighbor(face):
                                    cell->neighbor_of_neighbor(face);

                    if (cell_face->has_children())
                      {
                        for (unsigned int sub_nr = 0;
                             sub_nr != cell_face->n_children();
                             ++sub_nr)
                          {
                            const typename dealii::hp::DoFHandler<dim,spacedim>::level_cell_iterator
                            sub_neighbor
                              = periodic_neighbor?
                                cell->periodic_neighbor_child_on_subface (face, sub_nr):
                                cell->neighbor_child_on_subface (face, sub_nr);

                            dofs_on_other_cell.resize (sub_neighbor->get_fe().dofs_per_cell);
                            sub_neighbor->get_dof_indices (dofs_on_other_cell);
                            for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
                              {
                                for (unsigned int j=0; j<sub_neighbor->get_fe().dofs_per_cell;
                                     ++j)
                                  {
                                    if ((flux_mask(cell->get_fe().system_to_component_index(i).first,
                                                   sub_neighbor->get_fe().system_to_component_index(j).first)
                                         == always)
                                        ||
                                        (flux_mask(cell->get_fe().system_to_component_index(i).first,
                                                   sub_neighbor->get_fe().system_to_component_index(j).first)
                                         == nonzero))
                                      {
                                        sparsity.add (dofs_on_this_cell[i],
                                                      dofs_on_other_cell[j]);
                                        sparsity.add (dofs_on_other_cell[i],
                                                      dofs_on_this_cell[j]);
                                        sparsity.add (dofs_on_this_cell[i],
                                                      dofs_on_this_cell[j]);
                                        sparsity.add (dofs_on_other_cell[i],
                                                      dofs_on_other_cell[j]);
                                      }

                                    if ((flux_mask(sub_neighbor->get_fe().system_to_component_index(j).first,
                                                   cell->get_fe().system_to_component_index(i).first)
                                         == always)
                                        ||
                                        (flux_mask(sub_neighbor->get_fe().system_to_component_index(j).first,
                                                   cell->get_fe().system_to_component_index(i).first)
                                         == nonzero))
                                      {
                                        sparsity.add (dofs_on_this_cell[j],
                                                      dofs_on_other_cell[i]);
                                        sparsity.add (dofs_on_other_cell[j],
                                                      dofs_on_this_cell[i]);
                                        sparsity.add (dofs_on_this_cell[j],
                                                      dofs_on_this_cell[i]);
                                        sparsity.add (dofs_on_other_cell[j],
                                                      dofs_on_other_cell[i]);
                                      }
                                  }
                              }
                            sub_neighbor->face(neighbor_face)->set_user_flag ();
                          }
                      }
                    else
                      {
                        dofs_on_other_cell.resize (neighbor->get_fe().dofs_per_cell);
                        neighbor->get_dof_indices (dofs_on_other_cell);
                        for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
                          {
                            for (unsigned int j=0; j<neighbor->get_fe().dofs_per_cell; ++j)
                              {
                                if ((flux_mask(cell->get_fe().system_to_component_index(i).first,
                                               neighbor->get_fe().system_to_component_index(j).first)
                                     == always)
                                    ||
                                    (flux_mask(cell->get_fe().system_to_component_index(i).first,
                                               neighbor->get_fe().system_to_component_index(j).first)
                                     == nonzero))
                                  {
                                    sparsity.add (dofs_on_this_cell[i],
                                                  dofs_on_other_cell[j]);
                                    sparsity.add (dofs_on_other_cell[i],
                                                  dofs_on_this_cell[j]);
                                    sparsity.add (dofs_on_this_cell[i],
                                                  dofs_on_this_cell[j]);
                                    sparsity.add (dofs_on_other_cell[i],
                                                  dofs_on_other_cell[j]);
                                  }

                                if ((flux_mask(neighbor->get_fe().system_to_component_index(j).first,
                                               cell->get_fe().system_to_component_index(i).first)
                                     == always)
                                    ||
                                    (flux_mask(neighbor->get_fe().system_to_component_index(j).first,
                                               cell->get_fe().system_to_component_index(i).first)
                                     == nonzero))
                                  {
                                    sparsity.add (dofs_on_this_cell[j],
                                                  dofs_on_other_cell[i]);
                                    sparsity.add (dofs_on_other_cell[j],
                                                  dofs_on_this_cell[i]);
                                    sparsity.add (dofs_on_this_cell[j],
                                                  dofs_on_this_cell[i]);
                                    sparsity.add (dofs_on_other_cell[j],
                                                  dofs_on_other_cell[i]);
                                  }
                              }
                          }
                        neighbor->face(neighbor_face)->set_user_flag ();
                      }
                  }
              }
          }
      }
    }

  }




  template <typename DoFHandlerType, typename SparsityPatternType>
  void
  make_flux_sparsity_pattern (const DoFHandlerType    &dof,
                              SparsityPatternType     &sparsity,
                              const Table<2,Coupling> &int_mask,
                              const Table<2,Coupling> &flux_mask)
  {
    // do the error checking and frame code here, and then pass on to more
    // specialized functions in the internal namespace
    const types::global_dof_index n_dofs = dof.n_dofs();
    (void)n_dofs;
    const unsigned int n_comp = dof.get_fe().n_components();
    (void)n_comp;

    Assert (sparsity.n_rows() == n_dofs,
            ExcDimensionMismatch (sparsity.n_rows(), n_dofs));
    Assert (sparsity.n_cols() == n_dofs,
            ExcDimensionMismatch (sparsity.n_cols(), n_dofs));
    Assert (int_mask.n_rows() == n_comp,
            ExcDimensionMismatch (int_mask.n_rows(), n_comp));
    Assert (int_mask.n_cols() == n_comp,
            ExcDimensionMismatch (int_mask.n_cols(), n_comp));
    Assert (flux_mask.n_rows() == n_comp,
            ExcDimensionMismatch (flux_mask.n_rows(), n_comp));
    Assert (flux_mask.n_cols() == n_comp,
            ExcDimensionMismatch (flux_mask.n_cols(), n_comp));

    // Clear user flags because we will need them. But first we save them
    // and make sure that we restore them later such that at the end of
    // this function the Triangulation will be in the same state as it was
    // at the beginning of this function.
    std::vector<bool> user_flags;
    dof.get_triangulation().save_user_flags(user_flags);
    const_cast<Triangulation<DoFHandlerType::dimension,DoFHandlerType::space_dimension> &>
    (dof.get_triangulation()).clear_user_flags ();

    internal::make_flux_sparsity_pattern (dof, sparsity,
                                          int_mask, flux_mask);

    // finally restore the user flags
    const_cast<Triangulation<DoFHandlerType::dimension,DoFHandlerType::space_dimension> &>
    (dof.get_triangulation()).load_user_flags(user_flags);
  }


} // end of namespace DoFTools


// --------------------------------------------------- explicit instantiations

#include "dof_tools_sparsity.inst"



DEAL_II_NAMESPACE_CLOSE

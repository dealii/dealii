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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/petsc_block_sparse_matrix.h>
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
  namespace internal
  {
    /**
     * Given two triangulations, the first immersed inside the other, this
     * function computes and returns the real-space quadrature points of the
     * immersed triangulation.
     *
     * For reference:
     * cache->triangulation() is the imbdedding triangulation, which contains
     * immersed_dh->get_triangulation() the embedded triangulation
     *
     * Mapping and quadrature are those of this second triangulation.
     *
     * If the triangulation inside @p cache is parallel, only points lying over
     * locally owned cells are returned. This is why a vector of unsigned int
     * is returned: it describes the indices of cells from the immersed
     * triangulation which have been used (relative to a loop over al cells). If
     * embedding triangulation is not parallel, all cells shall be used.
     */
    template <int dim0, int dim1, int spacedim>
    std::pair<std::vector<Point<spacedim>>, std::vector<unsigned int>>
    qpoints_over_locally_owned_cells(
      const GridTools::Cache<dim0, spacedim> &cache,
      const DoFHandler<dim1, spacedim> &      immersed_dh,
      const Quadrature<dim1> &                quad,
      const Mapping<dim1, spacedim> &         immersed_mapping,
      const bool                              tria_is_parallel)
    {
      const auto &                 immersed_fe = immersed_dh.get_fe();
      std::vector<Point<spacedim>> points_over_local_cells;
      // Keep track of which cells we actually used
      std::vector<unsigned int> used_cells_ids;
      {
        FEValues<dim1, spacedim> fe_v(immersed_mapping,
                                      immersed_fe,
                                      quad,
                                      update_quadrature_points);
        unsigned int             cell_id = 0;
        for (const auto &cell : immersed_dh.active_cell_iterators())
          {
            bool use_cell = false;
            if (tria_is_parallel)
              {
                const auto bbox = cell->bounding_box();
                std::vector<std::pair<
                  BoundingBox<spacedim>,
                  typename Triangulation<dim0, spacedim>::active_cell_iterator>>
                  out_vals;
                cache.get_cell_bounding_boxes_rtree().query(
                  boost::geometry::index::intersects(bbox),
                  std::back_inserter(out_vals));
                // Each bounding box corresponds to an active cell
                // of the embedding triangulation: we now check if
                // the current cell, of the embedded triangulation,
                // overlaps a locally owned cell of the embedding one
                for (const auto &bbox_it : out_vals)
                  if (bbox_it.second->is_locally_owned())
                    {
                      use_cell = true;
                      used_cells_ids.emplace_back(cell_id);
                      break;
                    }
              }
            else
              // for sequential triangulations, simply use all cells
              use_cell = true;

            if (use_cell)
              {
                // Reinitialize the cell and the fe_values
                fe_v.reinit(cell);
                const std::vector<Point<spacedim>> &x_points =
                  fe_v.get_quadrature_points();

                // Insert the points to the vector
                points_over_local_cells.insert(points_over_local_cells.end(),
                                               x_points.begin(),
                                               x_points.end());
              }
            ++cell_id;
          }
      }
      return {std::move(points_over_local_cells), std::move(used_cells_ids)};
    }


    /**
     * Given two ComponentMasks and the corresponding finite element spaces,
     * compute a pairing between the selected components of the first finite
     * element space, and the selected components of the second finite element
     * space.
     */
    template <int dim0, int dim1, int spacedim>
    std::pair<std::vector<unsigned int>, std::vector<unsigned int>>
    compute_components_coupling(const ComponentMask &                comps0,
                                const ComponentMask &                comps1,
                                const FiniteElement<dim0, spacedim> &fe0,
                                const FiniteElement<dim1, spacedim> &fe1)
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
  } // namespace internal

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
    const Mapping<dim1, spacedim> &   immersed_mapping,
    const AffineConstraints<number> & immersed_constraints)
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
                                     immersed_mapping,
                                     immersed_constraints);
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
    const Mapping<dim1, spacedim> &         immersed_mapping,
    const AffineConstraints<number> &       immersed_constraints)
  {
    AssertDimension(sparsity.n_rows(), space_dh.n_dofs());
    AssertDimension(sparsity.n_cols(), immersed_dh.n_dofs());
    Assert(dim1 <= dim0,
           ExcMessage("This function can only work if dim1 <= dim0"));
    Assert((dynamic_cast<
              const parallel::distributed::Triangulation<dim1, spacedim> *>(
              &immersed_dh.get_triangulation()) == nullptr),
           ExcNotImplemented());

    const bool tria_is_parallel =
      (dynamic_cast<const parallel::TriangulationBase<dim1, spacedim> *>(
         &space_dh.get_triangulation()) != nullptr);
    const auto &space_fe    = space_dh.get_fe();
    const auto &immersed_fe = immersed_dh.get_fe();

    // Dof indices
    std::vector<types::global_dof_index> dofs(immersed_fe.n_dofs_per_cell());
    std::vector<types::global_dof_index> odofs(space_fe.n_dofs_per_cell());

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

    // Global to local indices
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

    const unsigned int n_q_points = quad.size();
    const unsigned int n_active_c =
      immersed_dh.get_triangulation().n_active_cells();

    const auto qpoints_cells_data = internal::qpoints_over_locally_owned_cells(
      cache, immersed_dh, quad, immersed_mapping, tria_is_parallel);

    const auto &points_over_local_cells = std::get<0>(qpoints_cells_data);
    const auto &used_cells_ids          = std::get<1>(qpoints_cells_data);

    // [TODO]: when the add_entries_local_to_global below will implement
    // the version with the dof_mask, this should be uncommented.
    //
    // // Construct a dof_mask, used to distribute entries to the sparsity
    // able< 2, bool > dof_mask(space_fe.n_dofs_per_cell(),
    //                          immersed_fe.n_dofs_per_cell());
    // of_mask.fill(false);
    // or (unsigned int i=0; i<space_fe.n_dofs_per_cell(); ++i)
    //  {
    //    const auto comp_i = space_fe.system_to_component_index(i).first;
    //    if (space_gtl[comp_i] != numbers::invalid_unsigned_int)
    //      for (unsigned int j=0; j<immersed_fe.n_dofs_per_cell(); ++j)
    //        {
    //          const auto comp_j =
    //          immersed_fe.system_to_component_index(j).first; if
    //          (immersed_gtl[comp_j] == space_gtl[comp_i])
    //            dof_mask(i,j) = true;
    //        }
    //  }


    // Get a list of outer cells, qpoints and maps.
    const auto cpm =
      GridTools::compute_point_locations(cache, points_over_local_cells);
    const auto &all_cells = std::get<0>(cpm);
    const auto &maps      = std::get<2>(cpm);

    std::vector<
      std::set<typename Triangulation<dim0, spacedim>::active_cell_iterator>>
      cell_sets(n_active_c);

    for (unsigned int i = 0; i < maps.size(); ++i)
      {
        // Quadrature points should be reasonably clustered:
        // the following index keeps track of the last id
        // where the current cell was inserted
        unsigned int last_id = std::numeric_limits<unsigned int>::max();
        unsigned int cell_id;
        for (const unsigned int idx : maps[i])
          {
            // Find in which cell of immersed triangulation the point lies
            if (tria_is_parallel)
              cell_id = used_cells_ids[idx / n_q_points];
            else
              cell_id = idx / n_q_points;

            if (last_id != cell_id)
              {
                cell_sets[cell_id].insert(all_cells[i]);
                last_id = cell_id;
              }
          }
      }

    // Now we run on each cell of the immersed
    // and build the sparsity
    unsigned int i = 0;
    for (const auto &cell : immersed_dh.active_cell_iterators())
      {
        // Reinitialize the cell
        cell->get_dof_indices(dofs);

        // List of outer cells
        const auto &cells = cell_sets[i];

        for (const auto &cell_c : cells)
          {
            // Get the ones in the current outer cell
            typename DoFHandler<dim0, spacedim>::cell_iterator ocell(*cell_c,
                                                                     &space_dh);
            // Make sure we act only on locally_owned cells
            if (ocell->is_locally_owned())
              {
                ocell->get_dof_indices(odofs);
                // [TODO]: When the following function will be implemented
                // for the case of non-trivial dof_mask, we should
                // uncomment the missing part.
                constraints.add_entries_local_to_global(
                  odofs,
                  immersed_constraints,
                  dofs,
                  sparsity); //, true, dof_mask);
              }
          }
        ++i;
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
    const Mapping<dim1, spacedim> &                       immersed_mapping,
    const AffineConstraints<typename Matrix::value_type> &immersed_constraints)
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
                                immersed_mapping,
                                immersed_constraints);
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
    const Mapping<dim1, spacedim> &                       immersed_mapping,
    const AffineConstraints<typename Matrix::value_type> &immersed_constraints)
  {
    AssertDimension(matrix.m(), space_dh.n_dofs());
    AssertDimension(matrix.n(), immersed_dh.n_dofs());
    Assert(dim1 <= dim0,
           ExcMessage("This function can only work if dim1 <= dim0"));
    Assert((dynamic_cast<
              const parallel::distributed::Triangulation<dim1, spacedim> *>(
              &immersed_dh.get_triangulation()) == nullptr),
           ExcNotImplemented());

    const bool tria_is_parallel =
      (dynamic_cast<const parallel::TriangulationBase<dim1, spacedim> *>(
         &space_dh.get_triangulation()) != nullptr);

    const auto &space_fe    = space_dh.get_fe();
    const auto &immersed_fe = immersed_dh.get_fe();

    // Dof indices
    std::vector<types::global_dof_index> dofs(immersed_fe.n_dofs_per_cell());
    std::vector<types::global_dof_index> odofs(space_fe.n_dofs_per_cell());

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
      space_dh.get_fe().n_dofs_per_cell(),
      immersed_dh.get_fe().n_dofs_per_cell());

    FEValues<dim1, spacedim> fe_v(immersed_mapping,
                                  immersed_dh.get_fe(),
                                  quad,
                                  update_JxW_values | update_quadrature_points |
                                    update_values);

    const unsigned int n_q_points = quad.size();
    const unsigned int n_active_c =
      immersed_dh.get_triangulation().n_active_cells();

    const auto used_cells_data = internal::qpoints_over_locally_owned_cells(
      cache, immersed_dh, quad, immersed_mapping, tria_is_parallel);

    const auto &points_over_local_cells = std::get<0>(used_cells_data);
    const auto &used_cells_ids          = std::get<1>(used_cells_data);

    // Get a list of outer cells, qpoints and maps.
    const auto cpm =
      GridTools::compute_point_locations(cache, points_over_local_cells);
    const auto &all_cells   = std::get<0>(cpm);
    const auto &all_qpoints = std::get<1>(cpm);
    const auto &all_maps    = std::get<2>(cpm);

    std::vector<
      std::vector<typename Triangulation<dim0, spacedim>::active_cell_iterator>>
      cell_container(n_active_c);
    std::vector<std::vector<std::vector<Point<dim0>>>> qpoints_container(
      n_active_c);
    std::vector<std::vector<std::vector<unsigned int>>> maps_container(
      n_active_c);

    // Cycle over all cells of underling mesh found
    // call it omesh, elaborating the output
    for (unsigned int o = 0; o < all_cells.size(); ++o)
      {
        for (unsigned int j = 0; j < all_maps[o].size(); ++j)
          {
            // Find the index of the "owner" cell and qpoint
            // with regard to the immersed mesh
            // Find in which cell of immersed triangulation the point lies
            unsigned int cell_id;
            if (tria_is_parallel)
              cell_id = used_cells_ids[all_maps[o][j] / n_q_points];
            else
              cell_id = all_maps[o][j] / n_q_points;

            const unsigned int n_pt = all_maps[o][j] % n_q_points;

            // If there are no cells, we just add our data
            if (cell_container[cell_id].empty())
              {
                cell_container[cell_id].emplace_back(all_cells[o]);
                qpoints_container[cell_id].emplace_back(
                  std::vector<Point<dim0>>{all_qpoints[o][j]});
                maps_container[cell_id].emplace_back(
                  std::vector<unsigned int>{n_pt});
              }
            // If there are already cells, we begin by looking
            // at the last inserted cell, which is more likely:
            else if (cell_container[cell_id].back() == all_cells[o])
              {
                qpoints_container[cell_id].back().emplace_back(
                  all_qpoints[o][j]);
                maps_container[cell_id].back().emplace_back(n_pt);
              }
            else
              {
                // We don't need to check the last element
                const auto cell_p = std::find(cell_container[cell_id].begin(),
                                              cell_container[cell_id].end() - 1,
                                              all_cells[o]);

                if (cell_p == cell_container[cell_id].end() - 1)
                  {
                    cell_container[cell_id].emplace_back(all_cells[o]);
                    qpoints_container[cell_id].emplace_back(
                      std::vector<Point<dim0>>{all_qpoints[o][j]});
                    maps_container[cell_id].emplace_back(
                      std::vector<unsigned int>{n_pt});
                  }
                else
                  {
                    const unsigned int pos =
                      cell_p - cell_container[cell_id].begin();
                    qpoints_container[cell_id][pos].emplace_back(
                      all_qpoints[o][j]);
                    maps_container[cell_id][pos].emplace_back(n_pt);
                  }
              }
          }
      }

    typename DoFHandler<dim1, spacedim>::active_cell_iterator
      cell = immersed_dh.begin_active(),
      endc = immersed_dh.end();

    for (unsigned int j = 0; cell != endc; ++cell, ++j)
      {
        // Reinitialize the cell and the fe_values
        fe_v.reinit(cell);
        cell->get_dof_indices(dofs);

        // Get a list of outer cells, qpoints and maps.
        const auto &cells   = cell_container[j];
        const auto &qpoints = qpoints_container[j];
        const auto &maps    = maps_container[j];

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

                for (unsigned int i = 0;
                     i < space_dh.get_fe().n_dofs_per_cell();
                     ++i)
                  {
                    const auto comp_i =
                      space_dh.get_fe().system_to_component_index(i).first;
                    if (space_gtl[comp_i] != numbers::invalid_unsigned_int)
                      for (unsigned int j = 0;
                           j < immersed_dh.get_fe().n_dofs_per_cell();
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
                constraints.distribute_local_to_global(
                  cell_matrix, odofs, immersed_constraints, dofs, matrix);
              }
          }
      }
  }

  template <int dim0,
            int dim1,
            int spacedim,
            typename Sparsity,
            typename Number>
  void
  create_coupling_sparsity_pattern(
    const double &                          epsilon,
    const GridTools::Cache<dim0, spacedim> &cache0,
    const GridTools::Cache<dim1, spacedim> &cache1,
    const DoFHandler<dim0, spacedim> &      dh0,
    const DoFHandler<dim1, spacedim> &      dh1,
    const Quadrature<dim1> &                quad,
    Sparsity &                              sparsity,
    const AffineConstraints<Number> &       constraints0,
    const ComponentMask &                   comps0,
    const ComponentMask &                   comps1)
  {
    if (epsilon == 0.0)
      {
        Assert(dim1 <= dim0,
               ExcMessage("When epsilon is zero, you can only "
                          "call this function with dim1 <= dim0."));
        create_coupling_sparsity_pattern(cache0,
                                         dh0,
                                         dh1,
                                         quad,
                                         sparsity,
                                         constraints0,
                                         comps0,
                                         comps1,
                                         cache1.get_mapping());
        return;
      }
    AssertDimension(sparsity.n_rows(), dh0.n_dofs());
    AssertDimension(sparsity.n_cols(), dh1.n_dofs());

    const bool zero_is_distributed =
      (dynamic_cast<const parallel::distributed::Triangulation<dim1, spacedim>
                      *>(&dh0.get_triangulation()) != nullptr);
    const bool one_is_distributed =
      (dynamic_cast<const parallel::distributed::Triangulation<dim1, spacedim>
                      *>(&dh1.get_triangulation()) != nullptr);

    // We bail out if both are distributed triangulations
    Assert(!zero_is_distributed || !one_is_distributed, ExcNotImplemented());

    // If we can loop on both, we decide where to make the outer loop according
    // to the size of the triangulation. The reasoning is the following:
    // - cost for accessing the tree: log(N)
    // - cost for computing the intersection for each of the outer loop cells: M
    // Total cost (besides the setup) is: M log(N)
    // If we can, make sure M is the smaller number of the two.
    const bool outer_loop_on_zero =
      (zero_is_distributed && !one_is_distributed) ||
      (dh1.get_triangulation().n_active_cells() >
       dh0.get_triangulation().n_active_cells());

    const auto &fe0 = dh0.get_fe();
    const auto &fe1 = dh1.get_fe();

    // Dof indices
    std::vector<types::global_dof_index> dofs0(fe0.n_dofs_per_cell());
    std::vector<types::global_dof_index> dofs1(fe1.n_dofs_per_cell());

    if (outer_loop_on_zero)
      {
        Assert(one_is_distributed == false, ExcInternalError());

        const auto &tree1 = cache1.get_cell_bounding_boxes_rtree();

        std::vector<std::pair<
          BoundingBox<spacedim>,
          typename Triangulation<dim1, spacedim>::active_cell_iterator>>
          intersection;

        for (const auto &cell0 :
             dh0.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
          {
            intersection.resize(0);
            BoundingBox<spacedim> box0 =
              cache0.get_mapping().get_bounding_box(cell0);
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
    else
      {
        Assert(zero_is_distributed == false, ExcInternalError());
        const auto &tree0 = cache0.get_cell_bounding_boxes_rtree();

        std::vector<std::pair<
          BoundingBox<spacedim>,
          typename Triangulation<dim0, spacedim>::active_cell_iterator>>
          intersection;

        for (const auto &cell1 :
             dh1.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
          {
            intersection.resize(0);
            BoundingBox<spacedim> box1 =
              cache1.get_mapping().get_bounding_box(cell1);
            box1.extend(epsilon);
            boost::geometry::index::query(tree0,
                                          boost::geometry::index::intersects(
                                            box1),
                                          std::back_inserter(intersection));
            if (!intersection.empty())
              {
                cell1->get_dof_indices(dofs1);
                for (const auto &entry : intersection)
                  {
                    typename DoFHandler<dim0, spacedim>::cell_iterator cell0(
                      *entry.second, &dh0);
                    cell0->get_dof_indices(dofs0);
                    constraints0.add_entries_local_to_global(dofs0,
                                                             dofs1,
                                                             sparsity);
                  }
              }
          }
      }
  }



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
    const AffineConstraints<typename Matrix::value_type> &constraints0,
    const ComponentMask &                                 comps0,
    const ComponentMask &                                 comps1)
  {
    if (epsilon == 0)
      {
        Assert(dim1 <= dim0,
               ExcMessage("When epsilon is zero, you can only "
                          "call this function with dim1 <= dim0."));
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
                      *>(&dh0.get_triangulation()) != nullptr);
    const bool one_is_distributed =
      (dynamic_cast<const parallel::distributed::Triangulation<dim1, spacedim>
                      *>(&dh1.get_triangulation()) != nullptr);

    // We bail out if both are distributed triangulations
    Assert(!zero_is_distributed || !one_is_distributed, ExcNotImplemented());

    // If we can loop on both, we decide where to make the outer loop according
    // to the size of the triangulation. The reasoning is the following:
    // - cost for accessing the tree: log(N)
    // - cost for computing the intersection for each of the outer loop cells: M
    // Total cost (besides the setup) is: M log(N)
    // If we can, make sure M is the smaller number of the two.
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
    std::vector<types::global_dof_index> dofs0(fe0.n_dofs_per_cell());
    std::vector<types::global_dof_index> dofs1(fe1.n_dofs_per_cell());

    // Local Matrix
    FullMatrix<typename Matrix::value_type> cell_matrix(fe0.n_dofs_per_cell(),
                                                        fe1.n_dofs_per_cell());

    // Global to local indices
    const auto p =
      internal::compute_components_coupling(comps0, comps1, fe0, fe1);
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
          for (unsigned int j = 0; j < fe1.n_dofs_per_cell(); ++j)
            {
              const auto comp_j = fe1.system_to_component_index(j).first;

              // First compute the part of the integral that does not
              // depend on i
              typename Matrix::value_type sum_q1 = {};
              for (unsigned int q1 = 0; q1 < quadrature1.size(); ++q1)
                sum_q1 +=
                  fev1.shape_value(j, q1) * kernel_values[q1] * fev1.JxW(q1);
              sum_q1 *= fev0.JxW(q0);

              // Now compute the main integral with the sum over q1 already
              // completed - this gives a cubic complexity as usual rather
              // than a quartic one with naive loops
              for (unsigned int i = 0; i < fe0.n_dofs_per_cell(); ++i)
                {
                  const auto comp_i = fe0.system_to_component_index(i).first;
                  if (gtl0[comp_i] != numbers::invalid_unsigned_int &&
                      gtl1[comp_j] == gtl0[comp_i])
                    cell_matrix(i, j) += fev0.shape_value(i, q0) * sum_q1;
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

        for (const auto &cell0 :
             dh0.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
          {
            intersection.resize(0);
            BoundingBox<spacedim> box0 =
              cache0.get_mapping().get_bounding_box(cell0);
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

        for (const auto &cell1 :
             dh1.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
          {
            intersection.resize(0);
            BoundingBox<spacedim> box1 =
              cache1.get_mapping().get_bounding_box(cell1);
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



#ifdef DEAL_II_WITH_CGAL
  template <int dim0, int dim1, int spacedim>
  std::vector<std::tuple<typename Triangulation<dim0, spacedim>::cell_iterator,
                         typename Triangulation<dim1, spacedim>::cell_iterator,
                         Quadrature<spacedim>>>
  compute_intersection(const GridTools::Cache<dim0, spacedim> &space_cache,
                       const GridTools::Cache<dim1, spacedim> &immersed_cache,
                       const unsigned int                      degree,
                       const double                            tol)
  {
    AssertThrow(
      dim1 <= dim0,
      ExcMessage(
        "Intrinsic dimension of the immersed object must be smaller than dim0."));
    AssertThrow(degree > 0, ExcMessage("Invalid quadrature degree."));
    AssertThrow(dim0 == 3 && dim1 == 3 && spacedim == 3,
                ExcNotImplemented("Not implemented for non-3D objects."));
    const auto bool_op = CGALWrappers::BooleanOperation::compute_intersection;
    std::vector<
      std::tuple<typename Triangulation<dim0, spacedim>::cell_iterator,
                 typename Triangulation<dim1, spacedim>::cell_iterator,
                 Quadrature<spacedim>>>
      cells_with_quadratures;

    const auto &space_tree =
      space_cache.get_locally_owned_cell_bounding_boxes_rtree();

    // The immersed tree *must* contain all cells, also the non-locally owned
    // ones.
    const auto &immersed_tree = immersed_cache.get_cell_bounding_boxes_rtree();

    // references to triangulations' info (cp cstrs marked as delete)
    const auto &mapping0 = space_cache.get_mapping();
    const auto &mapping1 = immersed_cache.get_mapping();
    namespace bgi        = boost::geometry::index;
    // Whenever the BB space_cell intersects the BB of an embedded cell,
    // store the space_cell in the set of intersected_cells
    for (const auto &[immersed_box, immersed_cell] : immersed_tree)
      {
        for (const auto &[space_box, space_cell] :
             space_tree | bgi::adaptors::queried(bgi::intersects(immersed_box)))
          {
            const auto test_intersection = CGALWrappers::
              compute_quadrature_on_boolean_operation<dim0, dim1, spacedim>(
                space_cell, immersed_cell, degree, bool_op, mapping0, mapping1);

            const auto & weights = test_intersection.get_weights();
            const double area =
              std::accumulate(weights.begin(), weights.end(), 0.0);
            if (area > tol) // non-trivial intersection
              {
                cells_with_quadratures.push_back(std::make_tuple(
                  space_cell, immersed_cell, test_intersection));
              }
          }
      }
    return cells_with_quadratures;
  }



  template <int dim0,
            int dim1,
            int spacedim,
            typename Sparsity,
            typename number>
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
    const AffineConstraints<number> &            immersed_constraints)
  {
    AssertDimension(sparsity.n_rows(), space_dh.n_dofs());
    AssertDimension(sparsity.n_cols(), immersed_dh.n_dofs());
    Assert(dim1 <= dim0,
           ExcMessage("This function can only work if dim1 <= dim0"));
    AssertThrow(dim0 == 3 && dim1 == 3 && spacedim == 3,
                ExcNotImplemented("Not implemented for non-3D objects."));
    Assert((dynamic_cast<
              const parallel::distributed::Triangulation<dim1, spacedim> *>(
              &immersed_dh.get_triangulation()) == nullptr),
           ExcNotImplemented());



    const auto &       space_fe                 = space_dh.get_fe();
    const auto &       immersed_fe              = immersed_dh.get_fe();
    const unsigned int n_dofs_per_space_cell    = space_fe.n_dofs_per_cell();
    const unsigned int n_dofs_per_immersed_cell = immersed_fe.n_dofs_per_cell();
    const unsigned int n_space_fe_components    = space_fe.n_components();
    const unsigned int n_immersed_fe_components = immersed_fe.n_components();
    std::vector<types::global_dof_index> space_dofs(n_dofs_per_space_cell);
    std::vector<types::global_dof_index> immersed_dofs(
      n_dofs_per_immersed_cell);


    const ComponentMask space_c =
      (space_comps.size() == 0 ? ComponentMask(n_space_fe_components, true) :
                                 space_comps);


    const ComponentMask immersed_c =
      (immersed_comps.size() == 0 ?
         ComponentMask(n_immersed_fe_components, true) :
         immersed_comps);

    AssertDimension(space_c.size(), n_space_fe_components);
    AssertDimension(immersed_c.size(), n_immersed_fe_components);


    // Global 2 Local indices
    std::vector<unsigned int> space_gtl(n_space_fe_components);
    std::vector<unsigned int> immersed_gtl(n_immersed_fe_components);
    for (unsigned int i = 0, j = 0; i < n_space_fe_components; ++i)
      {
        if (space_c[i])
          space_gtl[i] = j++;
      }


    for (unsigned int i = 0, j = 0; i < n_immersed_fe_components; ++i)
      {
        if (immersed_c[i])
          immersed_gtl[i] = j++;
      }



    Table<2, bool> dof_mask(n_dofs_per_space_cell, n_dofs_per_immersed_cell);
    dof_mask.fill(false); // start off by assuming they don't couple

    for (unsigned int i = 0; i < n_dofs_per_space_cell; ++i)
      {
        const auto comp_i = space_fe.system_to_component_index(i).first;
        if (space_gtl[comp_i] != numbers::invalid_unsigned_int)
          {
            for (unsigned int j = 0; j < n_dofs_per_immersed_cell; ++j)
              {
                const auto comp_j =
                  immersed_fe.system_to_component_index(j).first;
                if (immersed_gtl[comp_j] == space_gtl[comp_i])
                  {
                    dof_mask(i, j) = true;
                  }
              }
          }
      }

    const bool dof_mask_is_active =
      dof_mask.n_rows() == n_dofs_per_space_cell &&
      dof_mask.n_cols() == n_dofs_per_immersed_cell;

    // Whenever the BB space_cell intersects the BB of an embedded cell, those
    // DoFs have to be recorded

    for (const auto &it : intersections_info)
      {
        const auto &space_cell    = std::get<0>(it);
        const auto &immersed_cell = std::get<1>(it);
        typename DoFHandler<dim0, spacedim>::cell_iterator space_cell_dh(
          *space_cell, &space_dh);
        typename DoFHandler<dim1, spacedim>::cell_iterator immersed_cell_dh(
          *immersed_cell, &immersed_dh);

        space_cell_dh->get_dof_indices(space_dofs);
        immersed_cell_dh->get_dof_indices(immersed_dofs);

        if (dof_mask_is_active)
          {
            for (unsigned int i = 0; i < n_dofs_per_space_cell; ++i)
              {
                const unsigned int comp_i =
                  space_dh.get_fe().system_to_component_index(i).first;
                if (space_gtl[comp_i] != numbers::invalid_unsigned_int)
                  {
                    for (unsigned int j = 0; j < n_dofs_per_immersed_cell; ++j)
                      {
                        const unsigned int comp_j =
                          immersed_dh.get_fe()
                            .system_to_component_index(j)
                            .first;
                        if (space_gtl[comp_i] == immersed_gtl[comp_j])
                          {
                            // local_cell_matrix(i, j) +=
                            constraints.add_entries_local_to_global(
                              {space_dofs[i]},
                              immersed_constraints,
                              {immersed_dofs[j]},
                              sparsity,
                              true);
                          }
                      }
                  }
              }
          }
        else
          {
            constraints.add_entries_local_to_global(space_dofs,
                                                    immersed_constraints,
                                                    immersed_dofs,
                                                    sparsity,
                                                    true,
                                                    dof_mask);
          }
      }
  }



  template <int dim0, int dim1, int spacedim, typename Matrix>
  void
  create_coupling_mass_matrix_with_exact_intersections(
    const DoFHandler<dim0, spacedim> &space_dh,
    const DoFHandler<dim1, spacedim> &immersed_dh,
    const std::vector<
      std::tuple<typename Triangulation<dim0, spacedim>::cell_iterator,
                 typename Triangulation<dim1, spacedim>::cell_iterator,
                 Quadrature<spacedim>>> &                 cells_and_quads,
    Matrix &                                              matrix,
    const AffineConstraints<typename Matrix::value_type> &space_constraints,
    const ComponentMask &                                 space_comps,
    const ComponentMask &                                 immersed_comps,
    const Mapping<dim0, spacedim> &                       space_mapping,
    const Mapping<dim1, spacedim> &                       immersed_mapping,
    const AffineConstraints<typename Matrix::value_type> &immersed_constraints)
  {
    AssertDimension(matrix.m(), space_dh.n_dofs());
    AssertDimension(matrix.n(), immersed_dh.n_dofs());
    Assert(dim1 <= dim0,
           ExcMessage("This function can only work if dim1<=dim0"));
    AssertThrow(dim0 == 3 && dim1 == 3 && spacedim == 3,
                ExcNotImplemented("Not implemented for non-3D objects."));
    Assert((dynamic_cast<
              const parallel::distributed::Triangulation<dim1, spacedim> *>(
              &immersed_dh.get_triangulation()) == nullptr),
           ExcMessage("The immersed triangulation can only be a "
                      "parallel::shared::triangulation"));

    const auto &space_fe    = space_dh.get_fe();
    const auto &immersed_fe = immersed_dh.get_fe();

    const unsigned int n_dofs_per_space_cell    = space_fe.n_dofs_per_cell();
    const unsigned int n_dofs_per_immersed_cell = immersed_fe.n_dofs_per_cell();

    const unsigned int n_space_fe_components    = space_fe.n_components();
    const unsigned int n_immersed_fe_components = immersed_fe.n_components();

    FullMatrix<typename Matrix::value_type> local_cell_matrix(
      n_dofs_per_space_cell, n_dofs_per_immersed_cell);
    // DoF indices
    std::vector<types::global_dof_index> local_space_dof_indices(
      n_dofs_per_space_cell);
    std::vector<types::global_dof_index> local_immersed_dof_indices(
      n_dofs_per_immersed_cell);

    const ComponentMask space_c =
      (space_comps.size() == 0 ? ComponentMask(n_space_fe_components, true) :
                                 space_comps);
    const ComponentMask immersed_c =
      (immersed_comps.size() == 0 ?
         ComponentMask(n_immersed_fe_components, true) :
         immersed_comps);

    AssertDimension(space_c.size(), n_space_fe_components);
    AssertDimension(immersed_c.size(), n_immersed_fe_components);

    std::vector<unsigned int> space_gtl(n_space_fe_components,
                                        numbers::invalid_unsigned_int);
    std::vector<unsigned int> immersed_gtl(n_immersed_fe_components,
                                           numbers::invalid_unsigned_int);
    for (unsigned int i = 0, j = 0; i < n_space_fe_components; ++i)
      {
        if (space_c[i])
          space_gtl[i] = j++;
      }

    for (unsigned int i = 0, j = 0; i < n_immersed_fe_components; ++i)
      {
        if (immersed_c[i])
          immersed_gtl[i] = j++;
      }

    Table<2, bool> dof_mask(n_dofs_per_space_cell, n_dofs_per_immersed_cell);
    dof_mask.fill(false); // start off by assuming they don't couple

    for (unsigned int i = 0; i < n_dofs_per_space_cell; ++i)
      {
        const auto comp_i = space_fe.system_to_component_index(i).first;
        if (space_gtl[comp_i] != numbers::invalid_unsigned_int)
          {
            for (unsigned int j = 0; j < n_dofs_per_immersed_cell; ++j)
              {
                const auto comp_j =
                  immersed_fe.system_to_component_index(j).first;
                if (immersed_gtl[comp_j] == space_gtl[comp_i])
                  {
                    dof_mask(i, j) = true;
                  }
              }
          }
      }

    // Loop over vector of tuples, and gather everything together
    for (const auto &infos : cells_and_quads)
      {
        const auto &[first_cell, second_cell, quad_formula] = infos;

        local_cell_matrix = typename Matrix::value_type();

        const unsigned int       n_quad_pts = quad_formula.size();
        const auto &             real_qpts  = quad_formula.get_points();
        std::vector<Point<dim0>> ref_pts_space(n_quad_pts);
        std::vector<Point<dim1>> ref_pts_immersed(n_quad_pts);

        space_mapping.transform_points_real_to_unit_cell(first_cell,
                                                         real_qpts,
                                                         ref_pts_space);
        immersed_mapping.transform_points_real_to_unit_cell(second_cell,
                                                            real_qpts,
                                                            ref_pts_immersed);
        const auto &JxW = quad_formula.get_weights();
        for (unsigned int q = 0; q < n_quad_pts; ++q)
          {
            for (unsigned int i = 0; i < n_dofs_per_space_cell; ++i)
              {
                const unsigned int comp_i =
                  space_dh.get_fe().system_to_component_index(i).first;
                if (space_gtl[comp_i] != numbers::invalid_unsigned_int)
                  {
                    for (unsigned int j = 0; j < n_dofs_per_immersed_cell; ++j)
                      {
                        const unsigned int comp_j =
                          immersed_dh.get_fe()
                            .system_to_component_index(j)
                            .first;
                        if (space_gtl[comp_i] == immersed_gtl[comp_j])
                          {
                            local_cell_matrix(i, j) +=
                              space_fe.shape_value(i, ref_pts_space[q]) *
                              immersed_fe.shape_value(j, ref_pts_immersed[q]) *
                              JxW[q];
                          }
                      }
                  }
              }
          }
        typename DoFHandler<dim0, spacedim>::cell_iterator space_cell_dh(
          *first_cell, &space_dh);
        typename DoFHandler<dim1, spacedim>::cell_iterator immersed_cell_dh(
          *second_cell, &immersed_dh);
        space_cell_dh->get_dof_indices(local_space_dof_indices);
        immersed_cell_dh->get_dof_indices(local_immersed_dof_indices);

        space_constraints.distribute_local_to_global(local_cell_matrix,
                                                     local_space_dof_indices,
                                                     immersed_constraints,
                                                     local_immersed_dof_indices,
                                                     matrix);
      }
    matrix.compress(VectorOperation::add);
  }


  template std::vector<std::tuple<typename Triangulation<3, 3>::cell_iterator,
                                  typename Triangulation<3, 3>::cell_iterator,
                                  Quadrature<3>>>
  compute_intersection(const GridTools::Cache<3, 3> &space_cache,
                       const GridTools::Cache<3, 3> &immersed_cache,
                       const unsigned int            degree,
                       const double                  tol);
#endif
#ifndef DOXYGEN
#  include "coupling.inst"
#endif
} // namespace NonMatching

DEAL_II_NAMESPACE_CLOSE

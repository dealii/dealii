// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/logstream.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe.h>

#include <deal.II/grid/cell_id_translator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern_base.h>

#include <deal.II/multigrid/mg_base.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_tools.h>

#include <deal.II/numerics/matrix_tools.h>

#include <algorithm>
#include <numeric>
#include <set>
#include <vector>


DEAL_II_NAMESPACE_OPEN


namespace MGTools
{
  // specializations for 1d
  template <>
  void
  compute_row_length_vector(const DoFHandler<1, 1> &,
                            const unsigned int,
                            std::vector<unsigned int> &,
                            const DoFTools::Coupling)
  {
    DEAL_II_NOT_IMPLEMENTED();
  }



  template <>
  void
  compute_row_length_vector(const DoFHandler<1, 1> &,
                            const unsigned int,
                            std::vector<unsigned int> &,
                            const Table<2, DoFTools::Coupling> &,
                            const Table<2, DoFTools::Coupling> &)
  {
    DEAL_II_NOT_IMPLEMENTED();
  }



  template <>
  void
  compute_row_length_vector(const DoFHandler<1, 2> &,
                            const unsigned int,
                            std::vector<unsigned int> &,
                            const DoFTools::Coupling)
  {
    DEAL_II_NOT_IMPLEMENTED();
  }


  template <>
  void
  compute_row_length_vector(const DoFHandler<1, 2> &,
                            const unsigned int,
                            std::vector<unsigned int> &,
                            const Table<2, DoFTools::Coupling> &,
                            const Table<2, DoFTools::Coupling> &)
  {
    DEAL_II_NOT_IMPLEMENTED();
  }



  // Template for 2d and 3d. For 1d see specialization above
  template <int dim, int spacedim>
  void
  compute_row_length_vector(const DoFHandler<dim, spacedim> &dofs,
                            const unsigned int               level,
                            std::vector<unsigned int>       &row_lengths,
                            const DoFTools::Coupling         flux_coupling)
  {
    Assert(row_lengths.size() == dofs.n_dofs(),
           ExcDimensionMismatch(row_lengths.size(), dofs.n_dofs()));

    // Function starts here by
    // resetting the counters.
    std::fill(row_lengths.begin(), row_lengths.end(), 0);

    std::vector<bool> face_touched(dim == 2 ?
                                     dofs.get_triangulation().n_raw_lines() :
                                     dofs.get_triangulation().n_raw_quads());

    std::vector<types::global_dof_index> cell_indices;
    std::vector<types::global_dof_index> neighbor_indices;

    // We loop over cells and go from
    // cells to lower dimensional
    // objects. This is the only way to
    // cope with the fact, that an
    // unknown number of cells may
    // share an object of dimension
    // smaller than dim-1.
    for (const auto &cell : dofs.cell_iterators_on_level(level))
      {
        const FiniteElement<dim> &fe = cell->get_fe();
        cell_indices.resize(fe.n_dofs_per_cell());
        cell->get_mg_dof_indices(cell_indices);
        unsigned int i = 0;
        // First, dofs on
        // vertices. We assume that
        // each vertex dof couples
        // with all dofs on
        // adjacent grid cells.

        // Adding all dofs of the cells
        // will add dofs of the faces
        // of the cell adjacent to the
        // vertex twice. Therefore, we
        // subtract these here and add
        // them in a loop over the
        // faces below.

        // in 1d, faces and vertices
        // are identical. Nevertheless,
        // this will only work if
        // dofs_per_face is zero and
        // n_dofs_per_vertex() is
        // arbitrary, not the other way
        // round.
        // TODO: This assumes that the dofs per face on all faces coincide!
        const unsigned int face_no = 0;

        Assert(fe.reference_cell() == ReferenceCells::get_hypercube<dim>(),
               ExcNotImplemented());

        unsigned int increment =
          fe.n_dofs_per_cell() - dim * fe.n_dofs_per_face(face_no);
        while (i < fe.get_first_line_index())
          row_lengths[cell_indices[i++]] += increment;
        // From now on, if an object is
        // a cell, its dofs only couple
        // inside the cell. Since the
        // faces are handled below, we
        // have to subtract ALL faces
        // in this case.

        // In all other cases we
        // subtract adjacent faces to be
        // added in the loop below.
        increment =
          (dim > 1) ?
            fe.n_dofs_per_cell() - (dim - 1) * fe.n_dofs_per_face(face_no) :
            fe.n_dofs_per_cell() -
              GeometryInfo<dim>::faces_per_cell * fe.n_dofs_per_face(face_no);
        while (i < fe.get_first_quad_index(face_no))
          row_lengths[cell_indices[i++]] += increment;

        // Now quads in 2d and 3d
        increment =
          (dim > 2) ?
            fe.n_dofs_per_cell() - (dim - 2) * fe.n_dofs_per_face(face_no) :
            fe.n_dofs_per_cell() -
              GeometryInfo<dim>::faces_per_cell * fe.n_dofs_per_face(face_no);
        while (i < fe.get_first_hex_index())
          row_lengths[cell_indices[i++]] += increment;
        // Finally, cells in 3d
        increment = fe.n_dofs_per_cell() - GeometryInfo<dim>::faces_per_cell *
                                             fe.n_dofs_per_face(face_no);
        while (i < fe.n_dofs_per_cell())
          row_lengths[cell_indices[i++]] += increment;

        // At this point, we have
        // counted all dofs
        // contributing from cells
        // coupled topologically to the
        // adjacent cells, but we
        // subtracted some faces.

        // Now, let's go by the faces
        // and add the missing
        // contribution as well as the
        // flux contributions.
        for (const unsigned int iface : GeometryInfo<dim>::face_indices())
          {
            bool level_boundary = cell->at_boundary(iface);
            typename DoFHandler<dim, spacedim>::cell_iterator neighbor;
            if (!level_boundary)
              {
                neighbor = cell->neighbor(iface);
                if (static_cast<unsigned int>(neighbor->level()) != level)
                  level_boundary = true;
              }

            if (level_boundary)
              {
                for (unsigned int local_dof = 0;
                     local_dof < fe.n_dofs_per_cell();
                     ++local_dof)
                  row_lengths[cell_indices[local_dof]] +=
                    fe.n_dofs_per_face(face_no);
                continue;
              }

            const FiniteElement<dim> &nfe = neighbor->get_fe();
            typename DoFHandler<dim, spacedim>::face_iterator face =
              cell->face(iface);

            // Flux couplings are
            // computed from both sides
            // for simplicity.

            // The dofs on the common face
            // will be handled below,
            // therefore, we subtract them
            // here.
            if (flux_coupling != DoFTools::none)
              {
                const unsigned int dof_increment =
                  nfe.n_dofs_per_cell() - nfe.n_dofs_per_face(face_no);
                for (unsigned int local_dof = 0;
                     local_dof < fe.n_dofs_per_cell();
                     ++local_dof)
                  row_lengths[cell_indices[local_dof]] += dof_increment;
              }

            // Do this only once per
            // face.
            if (face_touched[face->index()])
              continue;
            face_touched[face->index()] = true;

            // At this point, we assume
            // that each cell added its
            // dofs minus the face to
            // the couplings of the
            // face dofs. Since we
            // subtracted two faces, we
            // have to re-add one.

            // If one side of the face
            // is refined, all the fine
            // face dofs couple with
            // the coarse one.
            neighbor_indices.resize(nfe.n_dofs_per_cell());
            neighbor->get_mg_dof_indices(neighbor_indices);
            for (unsigned int local_dof = 0; local_dof < fe.n_dofs_per_cell();
                 ++local_dof)
              row_lengths[cell_indices[local_dof]] +=
                nfe.n_dofs_per_face(face_no);
            for (unsigned int local_dof = 0; local_dof < nfe.n_dofs_per_cell();
                 ++local_dof)
              row_lengths[neighbor_indices[local_dof]] +=
                fe.n_dofs_per_face(face_no);
          }
      }
  }


  // This is the template for 2d and 3d. See version for 1d above
  template <int dim, int spacedim>
  void
  compute_row_length_vector(const DoFHandler<dim, spacedim>    &dofs,
                            const unsigned int                  level,
                            std::vector<unsigned int>          &row_lengths,
                            const Table<2, DoFTools::Coupling> &couplings,
                            const Table<2, DoFTools::Coupling> &flux_couplings)
  {
    Assert(row_lengths.size() == dofs.n_dofs(),
           ExcDimensionMismatch(row_lengths.size(), dofs.n_dofs()));

    // Function starts here by
    // resetting the counters.
    std::fill(row_lengths.begin(), row_lengths.end(), 0);

    std::vector<bool> face_touched(dim == 2 ?
                                     dofs.get_triangulation().n_raw_lines() :
                                     dofs.get_triangulation().n_raw_quads());

    std::vector<types::global_dof_index> cell_indices;
    std::vector<types::global_dof_index> neighbor_indices;

    // We have to translate the
    // couplings from components to
    // blocks, so this works for
    // nonprimitive elements as well.
    std::vector<Table<2, DoFTools::Coupling>> couple_cell;
    std::vector<Table<2, DoFTools::Coupling>> couple_face;
    DoFTools::convert_couplings_to_blocks(dofs, couplings, couple_cell);
    DoFTools::convert_couplings_to_blocks(dofs, flux_couplings, couple_face);

    // We loop over cells and go from
    // cells to lower dimensional
    // objects. This is the only way to
    // cope with the fact, that an
    // unknown number of cells may
    // share an object of dimension
    // smaller than dim-1.
    for (const auto &cell : dofs.cell_iterators_on_level(level))
      {
        const FiniteElement<dim> &fe       = cell->get_fe();
        const unsigned int        fe_index = cell->active_fe_index();


        // TODO: This assumes that the dofs per face on all faces coincide!
        const unsigned int face_no = 0;
        Assert(fe.reference_cell() == ReferenceCells::get_hypercube<dim>(),
               ExcNotImplemented());

        Assert(couplings.n_rows() == fe.n_components(),
               ExcDimensionMismatch(couplings.n_rows(), fe.n_components()));
        Assert(couplings.n_cols() == fe.n_components(),
               ExcDimensionMismatch(couplings.n_cols(), fe.n_components()));
        Assert(flux_couplings.n_rows() == fe.n_components(),
               ExcDimensionMismatch(flux_couplings.n_rows(),
                                    fe.n_components()));
        Assert(flux_couplings.n_cols() == fe.n_components(),
               ExcDimensionMismatch(flux_couplings.n_cols(),
                                    fe.n_components()));

        cell_indices.resize(fe.n_dofs_per_cell());
        cell->get_mg_dof_indices(cell_indices);
        unsigned int i = 0;
        // First, dofs on
        // vertices. We assume that
        // each vertex dof couples
        // with all dofs on
        // adjacent grid cells.

        // Adding all dofs of the cells
        // will add dofs of the faces
        // of the cell adjacent to the
        // vertex twice. Therefore, we
        // subtract these here and add
        // them in a loop over the
        // faces below.

        // in 1d, faces and vertices
        // are identical. Nevertheless,
        // this will only work if
        // dofs_per_face is zero and
        // n_dofs_per_vertex() is
        // arbitrary, not the other way
        // round.
        unsigned int increment;
        while (i < fe.get_first_line_index())
          {
            for (unsigned int base = 0; base < fe.n_base_elements(); ++base)
              for (unsigned int mult = 0; mult < fe.element_multiplicity(base);
                   ++mult)
                if (couple_cell[fe_index](fe.system_to_block_index(i).first,
                                          fe.first_block_of_base(base) +
                                            mult) != DoFTools::none)
                  {
                    increment =
                      fe.base_element(base).n_dofs_per_cell() -
                      dim * fe.base_element(base).n_dofs_per_face(face_no);
                    row_lengths[cell_indices[i]] += increment;
                  }
            ++i;
          }
        // From now on, if an object is
        // a cell, its dofs only couple
        // inside the cell. Since the
        // faces are handled below, we
        // have to subtract ALL faces
        // in this case.

        // In all other cases we
        // subtract adjacent faces to be
        // added in the loop below.
        while (i < fe.get_first_quad_index(face_no))
          {
            for (unsigned int base = 0; base < fe.n_base_elements(); ++base)
              for (unsigned int mult = 0; mult < fe.element_multiplicity(base);
                   ++mult)
                if (couple_cell[fe_index](fe.system_to_block_index(i).first,
                                          fe.first_block_of_base(base) +
                                            mult) != DoFTools::none)
                  {
                    increment =
                      fe.base_element(base).n_dofs_per_cell() -
                      ((dim > 1) ? (dim - 1) :
                                   GeometryInfo<dim>::faces_per_cell) *
                        fe.base_element(base).n_dofs_per_face(face_no);
                    row_lengths[cell_indices[i]] += increment;
                  }
            ++i;
          }

        // Now quads in 2d and 3d
        while (i < fe.get_first_hex_index())
          {
            for (unsigned int base = 0; base < fe.n_base_elements(); ++base)
              for (unsigned int mult = 0; mult < fe.element_multiplicity(base);
                   ++mult)
                if (couple_cell[fe_index](fe.system_to_block_index(i).first,
                                          fe.first_block_of_base(base) +
                                            mult) != DoFTools::none)
                  {
                    increment =
                      fe.base_element(base).n_dofs_per_cell() -
                      ((dim > 2) ? (dim - 2) :
                                   GeometryInfo<dim>::faces_per_cell) *
                        fe.base_element(base).n_dofs_per_face(face_no);
                    row_lengths[cell_indices[i]] += increment;
                  }
            ++i;
          }

        // Finally, cells in 3d
        while (i < fe.n_dofs_per_cell())
          {
            for (unsigned int base = 0; base < fe.n_base_elements(); ++base)
              for (unsigned int mult = 0; mult < fe.element_multiplicity(base);
                   ++mult)
                if (couple_cell[fe_index](fe.system_to_block_index(i).first,
                                          fe.first_block_of_base(base) +
                                            mult) != DoFTools::none)
                  {
                    increment =
                      fe.base_element(base).n_dofs_per_cell() -
                      GeometryInfo<dim>::faces_per_cell *
                        fe.base_element(base).n_dofs_per_face(face_no);
                    row_lengths[cell_indices[i]] += increment;
                  }
            ++i;
          }

        // At this point, we have
        // counted all dofs
        // contributing from cells
        // coupled topologically to the
        // adjacent cells, but we
        // subtracted some faces.

        // Now, let's go by the faces
        // and add the missing
        // contribution as well as the
        // flux contributions.
        for (const unsigned int iface : GeometryInfo<dim>::face_indices())
          {
            bool level_boundary = cell->at_boundary(iface);
            typename DoFHandler<dim, spacedim>::cell_iterator neighbor;
            if (!level_boundary)
              {
                neighbor = cell->neighbor(iface);
                if (static_cast<unsigned int>(neighbor->level()) != level)
                  level_boundary = true;
              }

            if (level_boundary)
              {
                for (unsigned int local_dof = 0;
                     local_dof < fe.n_dofs_per_cell();
                     ++local_dof)
                  row_lengths[cell_indices[local_dof]] +=
                    fe.n_dofs_per_face(face_no);
                continue;
              }

            const FiniteElement<dim> &nfe = neighbor->get_fe();
            typename DoFHandler<dim, spacedim>::face_iterator face =
              cell->face(iface);

            // Flux couplings are
            // computed from both sides
            // for simplicity.

            // The dofs on the common face
            // will be handled below,
            // therefore, we subtract them
            // here.
            for (unsigned int base = 0; base < nfe.n_base_elements(); ++base)
              for (unsigned int mult = 0; mult < nfe.element_multiplicity(base);
                   ++mult)
                for (unsigned int local_dof = 0;
                     local_dof < fe.n_dofs_per_cell();
                     ++local_dof)
                  if (couple_face[fe_index](
                        fe.system_to_block_index(local_dof).first,
                        nfe.first_block_of_base(base) + mult) != DoFTools::none)
                    {
                      const unsigned int dof_increment =
                        nfe.base_element(base).n_dofs_per_cell() -
                        nfe.base_element(base).n_dofs_per_face(face_no);
                      row_lengths[cell_indices[local_dof]] += dof_increment;
                    }

            // Do this only once per face and not on the hanging faces.
            if (face_touched[face->index()])
              continue;
            face_touched[face->index()] = true;

            // At this point, we assume
            // that each cell added its
            // dofs minus the face to
            // the couplings of the
            // face dofs. Since we
            // subtracted two faces, we
            // have to re-add one.

            // If one side of the face
            // is refined, all the fine
            // face dofs couple with
            // the coarse one.

            // Wolfgang, do they couple
            // with each other by
            // constraints?

            // This will not work with
            // different couplings on
            // different cells.
            neighbor_indices.resize(nfe.n_dofs_per_cell());
            neighbor->get_mg_dof_indices(neighbor_indices);
            for (unsigned int base = 0; base < nfe.n_base_elements(); ++base)
              for (unsigned int mult = 0; mult < nfe.element_multiplicity(base);
                   ++mult)
                for (unsigned int local_dof = 0;
                     local_dof < fe.n_dofs_per_cell();
                     ++local_dof)
                  if (couple_cell[fe_index](
                        fe.system_to_component_index(local_dof).first,
                        nfe.first_block_of_base(base) + mult) != DoFTools::none)
                    row_lengths[cell_indices[local_dof]] +=
                      nfe.base_element(base).n_dofs_per_face(face_no);
            for (unsigned int base = 0; base < fe.n_base_elements(); ++base)
              for (unsigned int mult = 0; mult < fe.element_multiplicity(base);
                   ++mult)
                for (unsigned int local_dof = 0;
                     local_dof < nfe.n_dofs_per_cell();
                     ++local_dof)
                  if (couple_cell[fe_index](
                        nfe.system_to_component_index(local_dof).first,
                        fe.first_block_of_base(base) + mult) != DoFTools::none)
                    row_lengths[neighbor_indices[local_dof]] +=
                      fe.base_element(base).n_dofs_per_face(face_no);
          }
      }
  }



  template <int dim, int spacedim, typename number>
  void
  make_sparsity_pattern(const DoFHandler<dim, spacedim> &dof,
                        SparsityPatternBase             &sparsity,
                        const unsigned int               level,
                        const AffineConstraints<number> &constraints,
                        const bool                       keep_constrained_dofs)
  {
    const types::global_dof_index n_dofs = dof.n_dofs(level);

    Assert(sparsity.n_rows() == n_dofs,
           ExcDimensionMismatch(sparsity.n_rows(), n_dofs));
    Assert(sparsity.n_cols() == n_dofs,
           ExcDimensionMismatch(sparsity.n_cols(), n_dofs));

    const unsigned int dofs_per_cell = dof.get_fe().n_dofs_per_cell();
    std::vector<types::global_dof_index> dofs_on_this_cell(dofs_per_cell);
    for (const auto &cell : dof.cell_iterators_on_level(level))
      if (cell->is_locally_owned_on_level())
        {
          cell->get_mg_dof_indices(dofs_on_this_cell);
          constraints.add_entries_local_to_global(dofs_on_this_cell,
                                                  sparsity,
                                                  keep_constrained_dofs);
        }
  }



  template <int dim, int spacedim, typename number>
  void
  make_flux_sparsity_pattern(const DoFHandler<dim, spacedim> &dof,
                             SparsityPatternBase             &sparsity,
                             const unsigned int               level,
                             const AffineConstraints<number> &constraints,
                             const bool keep_constrained_dofs)
  {
    const types::global_dof_index n_dofs = dof.n_dofs(level);

    Assert(sparsity.n_rows() == n_dofs,
           ExcDimensionMismatch(sparsity.n_rows(), n_dofs));
    Assert(sparsity.n_cols() == n_dofs,
           ExcDimensionMismatch(sparsity.n_cols(), n_dofs));

    const unsigned int dofs_per_cell = dof.get_fe().n_dofs_per_cell();
    std::vector<types::global_dof_index> dofs_on_this_cell(dofs_per_cell);
    std::vector<types::global_dof_index> dofs_on_other_cell(dofs_per_cell);
    for (const auto &cell : dof.cell_iterators_on_level(level))
      {
        if (!cell->is_locally_owned_on_level())
          continue;

        cell->get_mg_dof_indices(dofs_on_this_cell);
        // make sparsity pattern for this cell
        constraints.add_entries_local_to_global(dofs_on_this_cell,
                                                sparsity,
                                                keep_constrained_dofs);

        // Loop over all interior neighbors
        for (const unsigned int face : GeometryInfo<dim>::face_indices())
          {
            bool use_face = false;
            if ((!cell->at_boundary(face)) &&
                (static_cast<unsigned int>(cell->neighbor_level(face)) ==
                 level))
              use_face = true;
            else if (cell->has_periodic_neighbor(face) &&
                     (static_cast<unsigned int>(
                        cell->periodic_neighbor_level(face)) == level))
              use_face = true;

            if (use_face)
              {
                typename DoFHandler<dim, spacedim>::cell_iterator neighbor =
                  cell->neighbor_or_periodic_neighbor(face);
                neighbor->get_mg_dof_indices(dofs_on_other_cell);
                // only add one direction The other is taken care of by
                // neighbor (except when the neighbor is not owned by the same
                // processor)
                constraints.add_entries_local_to_global(dofs_on_this_cell,
                                                        dofs_on_other_cell,
                                                        sparsity,
                                                        keep_constrained_dofs);

                if (neighbor->is_locally_owned_on_level() == false)
                  {
                    constraints.add_entries_local_to_global(
                      dofs_on_other_cell, sparsity, keep_constrained_dofs);

                    constraints.add_entries_local_to_global(
                      dofs_on_other_cell,
                      dofs_on_this_cell,
                      sparsity,
                      keep_constrained_dofs);
                  }
              }
          }
      }
  }



  template <int dim, int spacedim>
  void
  make_flux_sparsity_pattern_edge(const DoFHandler<dim, spacedim> &dof,
                                  SparsityPatternBase             &sparsity,
                                  const unsigned int               level)
  {
    Assert((level >= 1) && (level < dof.get_triangulation().n_global_levels()),
           ExcIndexRange(level, 1, dof.get_triangulation().n_global_levels()));

    const types::global_dof_index fine_dofs   = dof.n_dofs(level);
    const types::global_dof_index coarse_dofs = dof.n_dofs(level - 1);

    // Matrix maps from fine level to coarse level
    Assert(sparsity.n_rows() == coarse_dofs,
           ExcDimensionMismatch(sparsity.n_rows(), coarse_dofs));
    Assert(sparsity.n_cols() == fine_dofs,
           ExcDimensionMismatch(sparsity.n_cols(), fine_dofs));

    const unsigned int dofs_per_cell = dof.get_fe().n_dofs_per_cell();
    std::vector<types::global_dof_index> dofs_on_this_cell(dofs_per_cell);
    std::vector<types::global_dof_index> dofs_on_other_cell(dofs_per_cell);
    for (const auto &cell : dof.cell_iterators_on_level(level))
      {
        if (!cell->is_locally_owned_on_level())
          continue;

        cell->get_mg_dof_indices(dofs_on_this_cell);
        // Loop over all interior neighbors
        for (const unsigned int face : GeometryInfo<dim>::face_indices())
          {
            // Neighbor is coarser
            bool use_face = false;
            if ((!cell->at_boundary(face)) &&
                (static_cast<unsigned int>(cell->neighbor_level(face)) !=
                 level))
              use_face = true;
            else if (cell->has_periodic_neighbor(face) &&
                     (static_cast<unsigned int>(
                        cell->periodic_neighbor_level(face)) != level))
              use_face = true;

            if (use_face)
              {
                typename DoFHandler<dim, spacedim>::cell_iterator neighbor =
                  cell->neighbor_or_periodic_neighbor(face);
                neighbor->get_mg_dof_indices(dofs_on_other_cell);

                std::sort(dofs_on_this_cell.begin(), dofs_on_this_cell.end());
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  sparsity.add_row_entries(dofs_on_other_cell[i],
                                           make_array_view(dofs_on_this_cell),
                                           true);
              }
          }
      }
  }



  template <int dim, int spacedim>
  void
  make_flux_sparsity_pattern(const DoFHandler<dim, spacedim>    &dof,
                             SparsityPatternBase                &sparsity,
                             const unsigned int                  level,
                             const Table<2, DoFTools::Coupling> &int_mask,
                             const Table<2, DoFTools::Coupling> &flux_mask)
  {
    const FiniteElement<dim>     &fe     = dof.get_fe();
    const types::global_dof_index n_dofs = dof.n_dofs(level);
    const unsigned int            n_comp = fe.n_components();

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

    const unsigned int                   total_dofs = fe.n_dofs_per_cell();
    std::vector<types::global_dof_index> dofs_on_this_cell(total_dofs);
    std::vector<types::global_dof_index> dofs_on_other_cell(total_dofs);
    std::vector<std::pair<types::global_dof_index, types::global_dof_index>>
      cell_entries;

    Table<2, bool> support_on_face(total_dofs,
                                   GeometryInfo<dim>::faces_per_cell);

    const Table<2, DoFTools::Coupling>
      int_dof_mask =
        DoFTools::dof_couplings_from_component_couplings(fe, int_mask),
      flux_dof_mask =
        DoFTools::dof_couplings_from_component_couplings(fe, flux_mask);

    for (unsigned int i = 0; i < total_dofs; ++i)
      for (auto f : GeometryInfo<dim>::face_indices())
        support_on_face(i, f) = fe.has_support_on_face(i, f);

    std::vector<bool> face_touched(dim == 2 ?
                                     dof.get_triangulation().n_raw_lines() :
                                     dof.get_triangulation().n_raw_quads());

    for (const auto &cell : dof.cell_iterators_on_level(level))
      {
        if (!cell->is_locally_owned_on_level())
          continue;

        cell->get_mg_dof_indices(dofs_on_this_cell);
        // make sparsity pattern for this cell
        for (unsigned int i = 0; i < total_dofs; ++i)
          for (unsigned int j = 0; j < total_dofs; ++j)
            if (int_dof_mask[i][j] != DoFTools::none)
              cell_entries.emplace_back(dofs_on_this_cell[i],
                                        dofs_on_this_cell[j]);

        // Loop over all interior neighbors
        for (const unsigned int face : GeometryInfo<dim>::face_indices())
          {
            typename DoFHandler<dim, spacedim>::face_iterator cell_face =
              cell->face(face);
            if (face_touched[cell_face->index()])
              continue;

            if (cell->at_boundary(face) && !cell->has_periodic_neighbor(face))
              {
                for (unsigned int i = 0; i < total_dofs; ++i)
                  {
                    const bool i_non_zero_i = support_on_face(i, face);
                    for (unsigned int j = 0; j < total_dofs; ++j)
                      {
                        const bool j_non_zero_i = support_on_face(j, face);

                        if (flux_dof_mask(i, j) == DoFTools::always)
                          cell_entries.emplace_back(dofs_on_this_cell[i],
                                                    dofs_on_this_cell[j]);
                        if (flux_dof_mask(i, j) == DoFTools::nonzero &&
                            i_non_zero_i && j_non_zero_i)
                          cell_entries.emplace_back(dofs_on_this_cell[i],
                                                    dofs_on_this_cell[j]);
                      }
                  }
              }
            else
              {
                typename DoFHandler<dim, spacedim>::cell_iterator neighbor =
                  cell->neighbor_or_periodic_neighbor(face);

                if (neighbor->level() < cell->level())
                  continue;

                unsigned int neighbor_face =
                  cell->has_periodic_neighbor(face) ?
                    cell->periodic_neighbor_of_periodic_neighbor(face) :
                    cell->neighbor_of_neighbor(face);

                neighbor->get_mg_dof_indices(dofs_on_other_cell);
                for (unsigned int i = 0; i < total_dofs; ++i)
                  {
                    const bool i_non_zero_i = support_on_face(i, face);
                    const bool i_non_zero_e = support_on_face(i, neighbor_face);
                    for (unsigned int j = 0; j < total_dofs; ++j)
                      {
                        const bool j_non_zero_i = support_on_face(j, face);
                        const bool j_non_zero_e =
                          support_on_face(j, neighbor_face);
                        if (flux_dof_mask(i, j) == DoFTools::always)
                          {
                            cell_entries.emplace_back(dofs_on_this_cell[i],
                                                      dofs_on_other_cell[j]);
                            cell_entries.emplace_back(dofs_on_other_cell[i],
                                                      dofs_on_this_cell[j]);
                            cell_entries.emplace_back(dofs_on_this_cell[i],
                                                      dofs_on_this_cell[j]);
                            cell_entries.emplace_back(dofs_on_other_cell[i],
                                                      dofs_on_other_cell[j]);
                          }
                        if (flux_dof_mask(i, j) == DoFTools::nonzero)
                          {
                            if (i_non_zero_i && j_non_zero_e)
                              cell_entries.emplace_back(dofs_on_this_cell[i],
                                                        dofs_on_other_cell[j]);
                            if (i_non_zero_e && j_non_zero_i)
                              cell_entries.emplace_back(dofs_on_other_cell[i],
                                                        dofs_on_this_cell[j]);
                            if (i_non_zero_i && j_non_zero_i)
                              cell_entries.emplace_back(dofs_on_this_cell[i],
                                                        dofs_on_this_cell[j]);
                            if (i_non_zero_e && j_non_zero_e)
                              cell_entries.emplace_back(dofs_on_other_cell[i],
                                                        dofs_on_other_cell[j]);
                          }

                        if (flux_dof_mask(j, i) == DoFTools::always)
                          {
                            cell_entries.emplace_back(dofs_on_this_cell[j],
                                                      dofs_on_other_cell[i]);
                            cell_entries.emplace_back(dofs_on_other_cell[j],
                                                      dofs_on_this_cell[i]);
                            cell_entries.emplace_back(dofs_on_this_cell[j],
                                                      dofs_on_this_cell[i]);
                            cell_entries.emplace_back(dofs_on_other_cell[j],
                                                      dofs_on_other_cell[i]);
                          }
                        if (flux_dof_mask(j, i) == DoFTools::nonzero)
                          {
                            if (j_non_zero_i && i_non_zero_e)
                              cell_entries.emplace_back(dofs_on_this_cell[j],
                                                        dofs_on_other_cell[i]);
                            if (j_non_zero_e && i_non_zero_i)
                              cell_entries.emplace_back(dofs_on_other_cell[j],
                                                        dofs_on_this_cell[i]);
                            if (j_non_zero_i && i_non_zero_i)
                              cell_entries.emplace_back(dofs_on_this_cell[j],
                                                        dofs_on_this_cell[i]);
                            if (j_non_zero_e && i_non_zero_e)
                              cell_entries.emplace_back(dofs_on_other_cell[j],
                                                        dofs_on_other_cell[i]);
                          }
                      }
                  }
                face_touched[neighbor->face(neighbor_face)->index()] = true;
              }
          }
        sparsity.add_entries(make_array_view(cell_entries));
        cell_entries.clear();
      }
  }



  template <int dim, int spacedim>
  void
  make_flux_sparsity_pattern_edge(const DoFHandler<dim, spacedim>    &dof,
                                  SparsityPatternBase                &sparsity,
                                  const unsigned int                  level,
                                  const Table<2, DoFTools::Coupling> &flux_mask)
  {
    const FiniteElement<dim> &fe     = dof.get_fe();
    const unsigned int        n_comp = fe.n_components();
    (void)n_comp;

    Assert((level >= 1) && (level < dof.get_triangulation().n_global_levels()),
           ExcIndexRange(level, 1, dof.get_triangulation().n_global_levels()));

    const types::global_dof_index fine_dofs   = dof.n_dofs(level);
    const types::global_dof_index coarse_dofs = dof.n_dofs(level - 1);

    // Matrix maps from fine level to coarse level
    Assert(sparsity.n_rows() == coarse_dofs,
           ExcDimensionMismatch(sparsity.n_rows(), coarse_dofs));
    Assert(sparsity.n_cols() == fine_dofs,
           ExcDimensionMismatch(sparsity.n_cols(), fine_dofs));
    Assert(flux_mask.n_rows() == n_comp,
           ExcDimensionMismatch(flux_mask.n_rows(), n_comp));
    Assert(flux_mask.n_cols() == n_comp,
           ExcDimensionMismatch(flux_mask.n_cols(), n_comp));

    const unsigned int dofs_per_cell = dof.get_fe().n_dofs_per_cell();
    std::vector<types::global_dof_index> dofs_on_this_cell(dofs_per_cell);
    std::vector<types::global_dof_index> dofs_on_other_cell(dofs_per_cell);
    std::vector<std::pair<types::global_dof_index, types::global_dof_index>>
      cell_entries;

    Table<2, bool> support_on_face(dofs_per_cell,
                                   GeometryInfo<dim>::faces_per_cell);

    const Table<2, DoFTools::Coupling> flux_dof_mask =
      DoFTools::dof_couplings_from_component_couplings(fe, flux_mask);

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      for (auto f : GeometryInfo<dim>::face_indices())
        support_on_face(i, f) = fe.has_support_on_face(i, f);

    for (const auto &cell : dof.cell_iterators_on_level(level))
      {
        if (!cell->is_locally_owned_on_level())
          continue;

        cell->get_mg_dof_indices(dofs_on_this_cell);
        // Loop over all interior neighbors
        for (const unsigned int face : GeometryInfo<dim>::face_indices())
          {
            // Neighbor is coarser
            bool use_face = false;
            if ((!cell->at_boundary(face)) &&
                (static_cast<unsigned int>(cell->neighbor_level(face)) !=
                 level))
              use_face = true;
            else if (cell->has_periodic_neighbor(face) &&
                     (static_cast<unsigned int>(
                        cell->periodic_neighbor_level(face)) != level))
              use_face = true;

            if (use_face)
              {
                typename DoFHandler<dim, spacedim>::cell_iterator neighbor =
                  cell->neighbor_or_periodic_neighbor(face);
                neighbor->get_mg_dof_indices(dofs_on_other_cell);

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  {
                    for (unsigned int j = 0; j < dofs_per_cell; ++j)
                      {
                        if (flux_dof_mask(i, j) != DoFTools::none)
                          {
                            cell_entries.emplace_back(dofs_on_other_cell[i],
                                                      dofs_on_this_cell[j]);
                            cell_entries.emplace_back(dofs_on_other_cell[j],
                                                      dofs_on_this_cell[i]);
                          }
                      }
                  }
              }
          }
        sparsity.add_entries(make_array_view(cell_entries));
        cell_entries.clear();
      }
  }



  template <int dim, int spacedim>
  void
  make_interface_sparsity_pattern(const DoFHandler<dim, spacedim> &dof,
                                  const MGConstrainedDoFs &mg_constrained_dofs,
                                  SparsityPatternBase     &sparsity,
                                  const unsigned int       level)
  {
    const types::global_dof_index n_dofs = dof.n_dofs(level);
    Assert(sparsity.n_rows() == n_dofs,
           ExcDimensionMismatch(sparsity.n_rows(), n_dofs));
    Assert(sparsity.n_cols() == n_dofs,
           ExcDimensionMismatch(sparsity.n_cols(), n_dofs));

    const unsigned int dofs_per_cell = dof.get_fe().n_dofs_per_cell();
    std::vector<types::global_dof_index> dofs_on_this_cell(dofs_per_cell);
    std::vector<types::global_dof_index> cols;
    cols.reserve(dofs_per_cell);

    for (const auto &cell : dof.cell_iterators_on_level(level))
      if (cell->is_locally_owned_on_level())
        {
          cell->get_mg_dof_indices(dofs_on_this_cell);
          std::sort(dofs_on_this_cell.begin(), dofs_on_this_cell.end());
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                if (mg_constrained_dofs.is_interface_matrix_entry(
                      level, dofs_on_this_cell[i], dofs_on_this_cell[j]))
                  cols.push_back(dofs_on_this_cell[j]);
              sparsity.add_row_entries(dofs_on_this_cell[i],
                                       make_array_view(cols),
                                       true);
              cols.clear();
            }
        }
  }



  template <int dim, int spacedim>
  void
  count_dofs_per_component(
    const DoFHandler<dim, spacedim>                   &dof_handler,
    std::vector<std::vector<types::global_dof_index>> &result,
    bool                                               only_once,
    std::vector<unsigned int>                          target_component)
  {
    const FiniteElement<dim> &fe           = dof_handler.get_fe();
    const unsigned int        n_components = fe.n_components();
    const unsigned int        nlevels =
      dof_handler.get_triangulation().n_global_levels();

    Assert(result.size() == nlevels,
           ExcDimensionMismatch(result.size(), nlevels));

    if (target_component.empty())
      {
        target_component.resize(n_components);
        for (unsigned int i = 0; i < n_components; ++i)
          target_component[i] = i;
      }

    Assert(target_component.size() == n_components,
           ExcDimensionMismatch(target_component.size(), n_components));

    for (unsigned int l = 0; l < nlevels; ++l)
      {
        result[l].resize(n_components);
        std::fill(result[l].begin(), result[l].end(), 0U);

        // special case for only one
        // component. treat this first
        // since it does not require any
        // computations
        if (n_components == 1)
          {
            result[l][0] = dof_handler.n_dofs(l);
          }
        else
          {
            // otherwise determine the number
            // of dofs in each component
            // separately. do so in parallel
            std::vector<std::vector<bool>> dofs_in_component(
              n_components, std::vector<bool>(dof_handler.n_dofs(l), false));
            std::vector<ComponentMask> component_select(n_components);
            Threads::TaskGroup<>       tasks;
            for (unsigned int i = 0; i < n_components; ++i)
              {
                void (*fun_ptr)(const unsigned int level,
                                const DoFHandler<dim, spacedim> &,
                                const ComponentMask &,
                                std::vector<bool> &) =
                  &DoFTools::extract_level_dofs<dim, spacedim>;

                std::vector<bool> tmp(n_components, false);
                tmp[i]              = true;
                component_select[i] = ComponentMask(tmp);

                tasks += Threads::new_task(fun_ptr,
                                           l,
                                           dof_handler,
                                           component_select[i],
                                           dofs_in_component[i]);
              }
            tasks.join_all();

            // next count what we got
            unsigned int component = 0;
            for (unsigned int b = 0; b < fe.n_base_elements(); ++b)
              {
                const FiniteElement<dim> &base = fe.base_element(b);
                // Dimension of base element
                unsigned int d = base.n_components();

                for (unsigned int m = 0; m < fe.element_multiplicity(b); ++m)
                  {
                    for (unsigned int dd = 0; dd < d; ++dd)
                      {
                        if (base.is_primitive() || (!only_once || dd == 0))
                          result[l][target_component[component]] +=
                            std::count(dofs_in_component[component].begin(),
                                       dofs_in_component[component].end(),
                                       true);
                        ++component;
                      }
                  }
              }
            // finally sanity check
            Assert(!dof_handler.get_fe().is_primitive() ||
                     std::accumulate(result[l].begin(),
                                     result[l].end(),
                                     types::global_dof_index(0)) ==
                       dof_handler.n_dofs(l),
                   ExcInternalError());
          }
      }
  }



  template <int dim, int spacedim>
  void
  count_dofs_per_block(
    const DoFHandler<dim, spacedim>                   &dof_handler,
    std::vector<std::vector<types::global_dof_index>> &dofs_per_block,
    std::vector<unsigned int>                          target_block)
  {
    const FiniteElement<dim, spacedim> &fe       = dof_handler.get_fe();
    const unsigned int                  n_blocks = fe.n_blocks();
    const unsigned int                  n_levels =
      dof_handler.get_triangulation().n_global_levels();

    AssertDimension(dofs_per_block.size(), n_levels);

    for (unsigned int l = 0; l < n_levels; ++l)
      std::fill(dofs_per_block[l].begin(), dofs_per_block[l].end(), 0U);
    // If the empty vector was given as
    // default argument, set up this
    // vector as identity.
    if (target_block.empty())
      {
        target_block.resize(n_blocks);
        for (unsigned int i = 0; i < n_blocks; ++i)
          target_block[i] = i;
      }
    Assert(target_block.size() == n_blocks,
           ExcDimensionMismatch(target_block.size(), n_blocks));

    const unsigned int max_block =
      *std::max_element(target_block.begin(), target_block.end());
    const unsigned int n_target_blocks = max_block + 1;
    (void)n_target_blocks;

    for (unsigned int l = 0; l < n_levels; ++l)
      AssertDimension(dofs_per_block[l].size(), n_target_blocks);

    // special case for only one
    // block. treat this first
    // since it does not require any
    // computations
    if (n_blocks == 1)
      {
        for (unsigned int l = 0; l < n_levels; ++l)
          dofs_per_block[l][0] = dof_handler.n_dofs(l);
        return;
      }
    // otherwise determine the number
    // of dofs in each block
    // separately. do so in parallel
    for (unsigned int l = 0; l < n_levels; ++l)
      {
        std::vector<std::vector<bool>> dofs_in_block(
          n_blocks, std::vector<bool>(dof_handler.n_dofs(l), false));
        std::vector<BlockMask> block_select(n_blocks);
        Threads::TaskGroup<>   tasks;
        for (unsigned int i = 0; i < n_blocks; ++i)
          {
            void (*fun_ptr)(const unsigned int level,
                            const DoFHandler<dim, spacedim> &,
                            const BlockMask &,
                            std::vector<bool> &) =
              &DoFTools::extract_level_dofs<dim, spacedim>;

            std::vector<bool> tmp(n_blocks, false);
            tmp[i]          = true;
            block_select[i] = BlockMask(tmp);

            tasks += Threads::new_task(
              fun_ptr, l, dof_handler, block_select[i], dofs_in_block[i]);
          }
        tasks.join_all();

        // next count what we got
        for (unsigned int block = 0; block < fe.n_blocks(); ++block)
          dofs_per_block[l][target_block[block]] +=
            std::count(dofs_in_block[block].begin(),
                       dofs_in_block[block].end(),
                       true);
      }
  }



  template <int dim, int spacedim>
  void
  make_boundary_list(
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim> *>
                                                   &function_map,
    std::vector<std::set<types::global_dof_index>> &boundary_indices,
    const ComponentMask                            &component_mask)
  {
    Assert(boundary_indices.size() == dof.get_triangulation().n_global_levels(),
           ExcDimensionMismatch(boundary_indices.size(),
                                dof.get_triangulation().n_global_levels()));

    std::set<types::boundary_id> boundary_ids;
    for (const auto &boundary_function : function_map)
      boundary_ids.insert(boundary_function.first);

    std::vector<IndexSet> boundary_indexset;
    make_boundary_list(dof, boundary_ids, boundary_indexset, component_mask);
    for (unsigned int i = 0; i < dof.get_triangulation().n_global_levels(); ++i)
      boundary_indices[i].insert(boundary_indexset[i].begin(),
                                 boundary_indexset[i].end());
  }


  template <int dim, int spacedim>
  void
  make_boundary_list(const DoFHandler<dim, spacedim>            &dof,
                     const std::map<types::boundary_id,
                                    const Function<spacedim> *> &function_map,
                     std::vector<IndexSet> &boundary_indices,
                     const ComponentMask   &component_mask)
  {
    Assert(boundary_indices.size() == dof.get_triangulation().n_global_levels(),
           ExcDimensionMismatch(boundary_indices.size(),
                                dof.get_triangulation().n_global_levels()));

    std::set<types::boundary_id> boundary_ids;
    for (const auto &boundary_function : function_map)
      boundary_ids.insert(boundary_function.first);

    make_boundary_list(dof, boundary_ids, boundary_indices, component_mask);
  }



  template <int dim, int spacedim>
  void
  make_boundary_list(const DoFHandler<dim, spacedim>    &dof,
                     const std::set<types::boundary_id> &boundary_ids,
                     std::vector<IndexSet>              &boundary_indices,
                     const ComponentMask                &component_mask)
  {
    boundary_indices.resize(dof.get_triangulation().n_global_levels());

    // if for whatever reason we were passed an empty set, return immediately
    if (boundary_ids.empty())
      return;

    for (unsigned int i = 0; i < dof.get_triangulation().n_global_levels(); ++i)
      if (boundary_indices[i].size() == 0)
        boundary_indices[i] = IndexSet(dof.n_dofs(i));

    const unsigned int n_components = dof.get_fe_collection().n_components();
    const bool         fe_is_system = (n_components != 1);

    std::vector<types::global_dof_index> local_dofs;
    local_dofs.reserve(dof.get_fe_collection().max_dofs_per_face());
    std::fill(local_dofs.begin(), local_dofs.end(), numbers::invalid_dof_index);

    std::vector<std::vector<types::global_dof_index>> dofs_by_level(
      dof.get_triangulation().n_levels());

    // First, deal with the simpler case when we have to identify all boundary
    // dofs
    if (component_mask.n_selected_components(n_components) == n_components)
      {
        for (const auto &cell : dof.cell_iterators())
          {
            if (cell->is_artificial_on_level())
              continue;
            const FiniteElement<dim> &fe    = cell->get_fe();
            const unsigned int        level = cell->level();

            for (const unsigned int face_no : GeometryInfo<dim>::face_indices())
              if (cell->at_boundary(face_no) == true)
                {
                  const typename DoFHandler<dim, spacedim>::face_iterator face =
                    cell->face(face_no);
                  const types::boundary_id bi = face->boundary_id();
                  // Face is listed in boundary map
                  if (boundary_ids.find(bi) != boundary_ids.end())
                    {
                      local_dofs.resize(fe.n_dofs_per_face(face_no));
                      face->get_mg_dof_indices(level, local_dofs);
                      dofs_by_level[level].insert(dofs_by_level[level].end(),
                                                  local_dofs.begin(),
                                                  local_dofs.end());
                    }
                }
          }
      }
    else
      {
        Assert(component_mask.n_selected_components(n_components) > 0,
               ExcMessage(
                 "It's probably worthwhile to select at least one component."));

        for (const auto &cell : dof.cell_iterators())
          if (!cell->is_artificial_on_level())
            for (const unsigned int face_no : GeometryInfo<dim>::face_indices())
              {
                if (cell->at_boundary(face_no) == false)
                  continue;

                const FiniteElement<dim> &fe    = cell->get_fe();
                const unsigned int        level = cell->level();

                typename DoFHandler<dim, spacedim>::face_iterator face =
                  cell->face(face_no);
                const types::boundary_id boundary_component =
                  face->boundary_id();
                if (boundary_ids.find(boundary_component) != boundary_ids.end())
                  // we want to constrain this boundary
                  {
                    for (unsigned int i = 0;
                         i < cell->get_fe().n_dofs_per_cell();
                         ++i)
                      {
                        const ComponentMask &nonzero_component_array =
                          cell->get_fe().get_nonzero_components(i);
                        // if we want to constrain one of the nonzero
                        // components, we have to constrain all of them

                        bool selected = false;
                        for (unsigned int c = 0; c < n_components; ++c)
                          if (nonzero_component_array[c] == true &&
                              component_mask[c] == true)
                            {
                              selected = true;
                              break;
                            }
                        if (selected)
                          for (unsigned int c = 0; c < n_components; ++c)
                            Assert(
                              nonzero_component_array[c] == false ||
                                component_mask[c] == true,
                              ExcMessage(
                                "You are using a non-primitive FiniteElement "
                                "and try to constrain just some of its components!"));
                      }

                    // get indices, physical location and boundary values of
                    // dofs on this face
                    local_dofs.resize(fe.n_dofs_per_face(face_no));
                    face->get_mg_dof_indices(level, local_dofs);
                    if (fe_is_system)
                      {
                        for (unsigned int i = 0; i < local_dofs.size(); ++i)
                          {
                            unsigned int component =
                              numbers::invalid_unsigned_int;
                            if (fe.is_primitive())
                              component =
                                fe.face_system_to_component_index(i, face_no)
                                  .first;
                            else
                              {
                                // Just pick the first of the components
                                // We already know that either all or none
                                // of the components are selected
                                const ComponentMask &nonzero_component_array =
                                  cell->get_fe().get_nonzero_components(i);
                                for (unsigned int c = 0; c < n_components; ++c)
                                  if (nonzero_component_array[c] == true)
                                    {
                                      component = c;
                                      break;
                                    }
                              }
                            Assert(component != numbers::invalid_unsigned_int,
                                   ExcInternalError());
                            if (component_mask[component] == true)
                              dofs_by_level[level].push_back(local_dofs[i]);
                          }
                      }
                    else
                      dofs_by_level[level].insert(dofs_by_level[level].end(),
                                                  local_dofs.begin(),
                                                  local_dofs.end());
                  }
              }
      }
    for (unsigned int level = 0; level < dof.get_triangulation().n_levels();
         ++level)
      {
        std::sort(dofs_by_level[level].begin(), dofs_by_level[level].end());
        boundary_indices[level].add_indices(dofs_by_level[level].begin(),
                                            dofs_by_level[level].end());
      }
  }



  template <int dim, int spacedim>
  void
  extract_inner_interface_dofs(const DoFHandler<dim, spacedim> &mg_dof_handler,
                               std::vector<IndexSet>           &interface_dofs)
  {
    Assert(interface_dofs.size() ==
             mg_dof_handler.get_triangulation().n_global_levels(),
           ExcDimensionMismatch(
             interface_dofs.size(),
             mg_dof_handler.get_triangulation().n_global_levels()));

    std::vector<std::vector<types::global_dof_index>> tmp_interface_dofs(
      interface_dofs.size());

    const FiniteElement<dim, spacedim> &fe = mg_dof_handler.get_fe();

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<bool> cell_dofs(dofs_per_cell, false);

    for (const auto &cell : mg_dof_handler.cell_iterators())
      {
        // Do not look at artificial level cells (in a serial computation we
        // need to ignore the level_subdomain_id() because it is never set).
        if (cell->is_artificial_on_level())
          continue;

        bool has_coarser_neighbor = false;

        std::fill(cell_dofs.begin(), cell_dofs.end(), false);

        for (const unsigned int face_nr : GeometryInfo<dim>::face_indices())
          {
            const typename DoFHandler<dim, spacedim>::face_iterator face =
              cell->face(face_nr);
            if (!face->at_boundary() || cell->has_periodic_neighbor(face_nr))
              {
                // interior face
                const typename DoFHandler<dim>::cell_iterator neighbor =
                  cell->neighbor_or_periodic_neighbor(face_nr);

                // only process cell pairs if one or both of them are owned by
                // me (ignore if running in serial)
                if (neighbor->is_artificial_on_level())
                  continue;

                // Do refinement face from the coarse side
                if (neighbor->level() < cell->level())
                  {
                    for (unsigned int j = 0; j < fe.n_dofs_per_face(face_nr);
                         ++j)
                      cell_dofs[fe.face_to_cell_index(j, face_nr)] = true;

                    has_coarser_neighbor = true;
                  }
              }
          }

        if (has_coarser_neighbor == false)
          continue;

        const unsigned int level = cell->level();
        cell->get_mg_dof_indices(local_dof_indices);

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            if (cell_dofs[i])
              tmp_interface_dofs[level].push_back(local_dof_indices[i]);
          }
      }

    for (unsigned int l = 0;
         l < mg_dof_handler.get_triangulation().n_global_levels();
         ++l)
      {
        interface_dofs[l].clear();
        std::sort(tmp_interface_dofs[l].begin(), tmp_interface_dofs[l].end());
        interface_dofs[l].add_indices(tmp_interface_dofs[l].begin(),
                                      tmp_interface_dofs[l].end());
        interface_dofs[l].compress();
      }
  }



  template <int dim, int spacedim>
  unsigned int
  max_level_for_coarse_mesh(const Triangulation<dim, spacedim> &tria)
  {
    // Find minimum level for an active cell in
    // this locally owned subdomain
    // Note: with the way active cells are traversed,
    // the first locally owned cell we find will have
    // the lowest level in the particular subdomain.
    unsigned int min_level = tria.n_global_levels();
    for (const auto &cell : tria.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          min_level = cell->level();
          break;
        }

    unsigned int global_min = min_level;
    // If necessary, communicate to find minimum
    // level for an active cell over all subdomains
    if (const parallel::TriangulationBase<dim, spacedim> *tr =
          dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
            &tria))
      global_min = Utilities::MPI::min(min_level, tr->get_mpi_communicator());

    AssertIndexRange(global_min, tria.n_global_levels());

    return global_min;
  }



  namespace internal
  {
    double
    workload_imbalance(
      const std::vector<types::global_dof_index> &n_cells_on_levels,
      const MPI_Comm                              comm)
    {
      std::vector<types::global_dof_index> n_cells_on_levels_max(
        n_cells_on_levels.size());
      std::vector<types::global_dof_index> n_cells_on_levels_sum(
        n_cells_on_levels.size());

      Utilities::MPI::max(n_cells_on_levels, comm, n_cells_on_levels_max);
      Utilities::MPI::sum(n_cells_on_levels, comm, n_cells_on_levels_sum);

      const unsigned int n_proc = Utilities::MPI::n_mpi_processes(comm);

      const double ideal_work = std::accumulate(n_cells_on_levels_sum.begin(),
                                                n_cells_on_levels_sum.end(),
                                                0) /
                                static_cast<double>(n_proc);
      const double workload_imbalance =
        std::accumulate(n_cells_on_levels_max.begin(),
                        n_cells_on_levels_max.end(),
                        0) /
        ideal_work;

      return workload_imbalance;
    }



    double
    vertical_communication_efficiency(
      const std::vector<
        std::pair<types::global_dof_index, types::global_dof_index>> &cells,
      const MPI_Comm                                                  comm)
    {
      std::vector<types::global_dof_index> cells_local(cells.size());
      std::vector<types::global_dof_index> cells_remote(cells.size());

      for (unsigned int i = 0; i < cells.size(); ++i)
        {
          cells_local[i]  = cells[i].first;
          cells_remote[i] = cells[i].second;
        }

      std::vector<types::global_dof_index> cells_local_sum(cells_local.size());
      Utilities::MPI::sum(cells_local, comm, cells_local_sum);

      std::vector<types::global_dof_index> cells_remote_sum(
        cells_remote.size());
      Utilities::MPI::sum(cells_remote, comm, cells_remote_sum);

      const auto n_cells_local =
        std::accumulate(cells_local_sum.begin(), cells_local_sum.end(), 0);
      const auto n_cells_remote =
        std::accumulate(cells_remote_sum.begin(), cells_remote_sum.end(), 0);

      return static_cast<double>(n_cells_local) /
             (n_cells_local + n_cells_remote);
    }
  } // namespace internal



  template <int dim, int spacedim>
  std::vector<types::global_dof_index>
  local_workload(const Triangulation<dim, spacedim> &tria)
  {
    if (const parallel::TriangulationBase<dim, spacedim> *tr =
          dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
            &tria))
      Assert(
        tr->is_multilevel_hierarchy_constructed(),
        ExcMessage(
          "We can only compute the workload imbalance if the multilevel hierarchy has been constructed!"));

    const unsigned int n_global_levels = tria.n_global_levels();

    std::vector<types::global_dof_index> n_cells_on_levels(n_global_levels);

    for (unsigned int lvl = 0; lvl < n_global_levels; ++lvl)
      for (const auto &cell : tria.cell_iterators_on_level(lvl))
        if (cell->is_locally_owned_on_level())
          ++n_cells_on_levels[lvl];

    return n_cells_on_levels;
  }



  template <int dim, int spacedim>
  std::vector<types::global_dof_index>
  local_workload(
    const std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
      &trias)
  {
    const unsigned int n_global_levels = trias.size();

    std::vector<types::global_dof_index> n_cells_on_levels(n_global_levels);

    for (unsigned int lvl = 0; lvl < n_global_levels; ++lvl)
      for (const auto &cell : trias[lvl]->active_cell_iterators())
        if (cell->is_locally_owned())
          ++n_cells_on_levels[lvl];

    return n_cells_on_levels;
  }



  template <int dim, int spacedim>
  double
  workload_imbalance(const Triangulation<dim, spacedim> &tria)
  {
    return internal::workload_imbalance(local_workload(tria),
                                        tria.get_mpi_communicator());
  }



  template <int dim, int spacedim>
  double
  workload_imbalance(
    const std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
      &trias)
  {
    return internal::workload_imbalance(local_workload(trias),
                                        trias.back()->get_mpi_communicator());
  }



  template <int dim, int spacedim>
  std::vector<std::pair<types::global_dof_index, types::global_dof_index>>
  local_vertical_communication_cost(const Triangulation<dim, spacedim> &tria)
  {
    const unsigned int n_global_levels = tria.n_global_levels();

    std::vector<std::pair<types::global_dof_index, types::global_dof_index>>
      cells(n_global_levels);

    const MPI_Comm communicator = tria.get_mpi_communicator();

    const unsigned int my_rank = Utilities::MPI::this_mpi_process(communicator);

    for (unsigned int lvl = 0; lvl < n_global_levels - 1; ++lvl)
      for (const auto &cell : tria.cell_iterators_on_level(lvl))
        if (cell->is_locally_owned_on_level() && cell->has_children())
          for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_cell;
               ++i)
            {
              const auto level_subdomain_id =
                cell->child(i)->level_subdomain_id();
              if (level_subdomain_id == my_rank)
                ++cells[lvl + 1].first;
              else if (level_subdomain_id != numbers::invalid_unsigned_int)
                ++cells[lvl + 1].second;
              else
                AssertThrow(false, ExcNotImplemented());
            }

    return cells;
  }



  template <int dim, int spacedim>
  std::vector<std::pair<types::global_dof_index, types::global_dof_index>>
  local_vertical_communication_cost(
    const std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
      &trias)
  {
    const unsigned int n_global_levels = trias.size();

    std::vector<std::pair<types::global_dof_index, types::global_dof_index>>
      cells(n_global_levels);

    const MPI_Comm communicator = trias.back()->get_mpi_communicator();

    const unsigned int my_rank = Utilities::MPI::this_mpi_process(communicator);

    for (unsigned int lvl = 0; lvl < n_global_levels - 1; ++lvl)
      {
        const auto &tria_coarse = *trias[lvl];
        const auto &tria_fine   = *trias[lvl + 1];

        const unsigned int n_coarse_cells  = tria_fine.n_global_coarse_cells();
        const unsigned int n_global_levels = tria_fine.n_global_levels();

        const dealii::internal::CellIDTranslator<dim> cell_id_translator(
          n_coarse_cells, n_global_levels);

        IndexSet is_fine_owned(cell_id_translator.size());
        IndexSet is_fine_required(cell_id_translator.size());

        for (const auto &cell : tria_fine.active_cell_iterators())
          if (!cell->is_artificial() && cell->is_locally_owned())
            is_fine_owned.add_index(cell_id_translator.translate(cell));

        for (const auto &cell : tria_coarse.active_cell_iterators())
          if (!cell->is_artificial() && cell->is_locally_owned())
            {
              if (cell->level() + 1u == tria_fine.n_global_levels())
                continue;

              for (unsigned int i = 0;
                   i < GeometryInfo<dim>::max_children_per_cell;
                   ++i)
                is_fine_required.add_index(
                  cell_id_translator.translate(cell, i));
            }

        const std::vector<unsigned int> is_fine_required_ranks =
          Utilities::MPI::compute_index_owner(is_fine_owned,
                                              is_fine_required,
                                              communicator);

        for (unsigned i = 0; i < is_fine_required.n_elements(); ++i)
          if (is_fine_required_ranks[i] == my_rank)
            ++cells[lvl + 1].first;
          else if (is_fine_required_ranks[i] != numbers::invalid_unsigned_int)
            ++cells[lvl + 1].second;
      }

    return cells;
  }



  template <int dim, int spacedim>
  double
  vertical_communication_efficiency(const Triangulation<dim, spacedim> &tria)
  {
    return internal::vertical_communication_efficiency(
      local_vertical_communication_cost(tria), tria.get_mpi_communicator());
  }



  template <int dim, int spacedim>
  double
  vertical_communication_efficiency(
    const std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
      &trias)
  {
    return internal::vertical_communication_efficiency(
      local_vertical_communication_cost(trias),
      trias.back()->get_mpi_communicator());
  }

} // namespace MGTools


// explicit instantiations
#include "multigrid/mg_tools.inst"

DEAL_II_NAMESPACE_CLOSE

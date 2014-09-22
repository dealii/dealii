// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2014 by the deal.II authors
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

#include <deal.II/base/logstream.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_base.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/numerics/matrix_tools.h>

#include <vector>
#include <algorithm>
#include <numeric>

DEAL_II_NAMESPACE_OPEN


namespace MGTools
{


  // specializations for 1D
  template <>
  void
  compute_row_length_vector(
    const DoFHandler<1,1> &,
    const unsigned int,
    std::vector<unsigned int> &,
    const DoFTools::Coupling)
  {
    Assert(false, ExcNotImplemented());
  }



  template <>
  void
  compute_row_length_vector(
    const DoFHandler<1,1> &,
    const unsigned int,
    std::vector<unsigned int> &,
    const Table<2,DoFTools::Coupling> &,
    const Table<2,DoFTools::Coupling> &)
  {
    Assert(false, ExcNotImplemented());
  }



  template <>
  void
  compute_row_length_vector(
    const DoFHandler<1,2> &,
    const unsigned int,
    std::vector<unsigned int> &,
    const DoFTools::Coupling)
  {
    Assert(false, ExcNotImplemented());
  }


  template <>
  void
  compute_row_length_vector(
    const DoFHandler<1,2> &,
    const unsigned int,
    std::vector<unsigned int> &,
    const Table<2,DoFTools::Coupling> &,
    const Table<2,DoFTools::Coupling> &)
  {
    Assert(false, ExcNotImplemented());
  }



// Template for 2D and 3D. For 1D see specialization above
  template <int dim, int spacedim>
  void
  compute_row_length_vector(
    const DoFHandler<dim,spacedim> &dofs,
    const unsigned int level,
    std::vector<unsigned int> &row_lengths,
    const DoFTools::Coupling             flux_coupling)
  {
    Assert (row_lengths.size() == dofs.n_dofs(),
            ExcDimensionMismatch(row_lengths.size(), dofs.n_dofs()));

    // Function starts here by
    // resetting the counters.
    std::fill(row_lengths.begin(), row_lengths.end(), 0);
    // We need the user flags, so we
    // save them for later restoration
    std::vector<bool> old_flags;
    // We need a non-constant
    // triangulation for the user
    // flags. Since we restore them in
    // the end, this cast is safe.
    Triangulation<dim,spacedim> &user_flags_triangulation =
      const_cast<Triangulation<dim,spacedim>&> (dofs.get_tria());
    user_flags_triangulation.save_user_flags(old_flags);
    user_flags_triangulation.clear_user_flags();

    const typename DoFHandler<dim,spacedim>::cell_iterator end = dofs.end(level);
    typename DoFHandler<dim,spacedim>::active_cell_iterator cell;
    std::vector<types::global_dof_index> cell_indices;
    std::vector<types::global_dof_index> neighbor_indices;

    // We loop over cells and go from
    // cells to lower dimensional
    // objects. This is the only way to
    // cope with the fact, that an
    // unknown number of cells may
    // share an object of dimension
    // smaller than dim-1.
    for (cell = dofs.begin(level); cell != end; ++cell)
      {
        const FiniteElement<dim> &fe = cell->get_fe();
        cell_indices.resize(fe.dofs_per_cell);
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
        // dofs_per_vertex is
        // arbitrary, not the other way
        // round.
//TODO: This assumes that even in hp context, the dofs per face coincide!
        unsigned int increment = fe.dofs_per_cell - dim * fe.dofs_per_face;
        while (i < fe.first_line_index)
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
        increment = (dim>1)
                    ? fe.dofs_per_cell - (dim-1) * fe.dofs_per_face
                    : fe.dofs_per_cell - GeometryInfo<dim>::faces_per_cell * fe.dofs_per_face;
        while (i < fe.first_quad_index)
          row_lengths[cell_indices[i++]] += increment;

        // Now quads in 2D and 3D
        increment = (dim>2)
                    ? fe.dofs_per_cell - (dim-2) * fe.dofs_per_face
                    : fe.dofs_per_cell - GeometryInfo<dim>::faces_per_cell * fe.dofs_per_face;
        while (i < fe.first_hex_index)
          row_lengths[cell_indices[i++]] += increment;
        // Finally, cells in 3D
        increment = fe.dofs_per_cell - GeometryInfo<dim>::faces_per_cell * fe.dofs_per_face;
        while (i < fe.dofs_per_cell)
          row_lengths[cell_indices[i++]] += increment;

        // At this point, we have
        // counted all dofs
        // contributiong from cells
        // coupled topologically to the
        // adjacent cells, but we
        // subtracted some faces.

        // Now, let's go by the faces
        // and add the missing
        // contribution as well as the
        // flux contributions.
        for (unsigned int iface=0; iface<GeometryInfo<dim>::faces_per_cell; ++iface)
          {
            bool level_boundary = cell->at_boundary(iface);
            typename DoFHandler<dim,spacedim>::cell_iterator neighbor;
            if (!level_boundary)
              {
                neighbor = cell->neighbor(iface);
                if (static_cast<unsigned int>(neighbor->level()) != level)
                  level_boundary = true;
              }

            if (level_boundary)
              {
                for (unsigned int local_dof=0; local_dof<fe.dofs_per_cell; ++local_dof)
                  row_lengths[cell_indices[local_dof]] += fe.dofs_per_face;
                continue;
              }

            const FiniteElement<dim> &nfe = neighbor->get_fe();
            typename DoFHandler<dim,spacedim>::face_iterator face = cell->face(iface);

            // Flux couplings are
            // computed from both sides
            // for simplicity.

            // The dofs on the common face
            // will be handled below,
            // therefore, we subtract them
            // here.
            if (flux_coupling != DoFTools::none)
              {
                const unsigned int dof_increment = nfe.dofs_per_cell - nfe.dofs_per_face;
                for (unsigned int local_dof=0; local_dof<fe.dofs_per_cell; ++local_dof)
                  row_lengths[cell_indices[local_dof]] += dof_increment;
              }

            // Do this only once per
            // face.
            if (face->user_flag_set())
              continue;
            face->set_user_flag();
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
            neighbor_indices.resize(nfe.dofs_per_cell);
            neighbor->get_mg_dof_indices(neighbor_indices);
            for (unsigned int local_dof=0; local_dof<fe.dofs_per_cell; ++local_dof)
              row_lengths[cell_indices[local_dof]] += nfe.dofs_per_face;
            for (unsigned int local_dof=0; local_dof<nfe.dofs_per_cell; ++local_dof)
              row_lengths[neighbor_indices[local_dof]] += fe.dofs_per_face;
          }
      }
    user_flags_triangulation.load_user_flags(old_flags);
  }


// This is the template for 2D and 3D. See version for 1D above
  template <int dim, int spacedim>
  void
  compute_row_length_vector(
    const DoFHandler<dim,spacedim> &dofs,
    const unsigned int level,
    std::vector<unsigned int> &row_lengths,
    const Table<2,DoFTools::Coupling> &couplings,
    const Table<2,DoFTools::Coupling> &flux_couplings)
  {
    Assert (row_lengths.size() == dofs.n_dofs(),
            ExcDimensionMismatch(row_lengths.size(), dofs.n_dofs()));

    // Function starts here by
    // resetting the counters.
    std::fill(row_lengths.begin(), row_lengths.end(), 0);
    // We need the user flags, so we
    // save them for later restoration
    std::vector<bool> old_flags;
    // We need a non-constant
    // triangulation for the user
    // flags. Since we restore them in
    // the end, this cast is safe.
    Triangulation<dim,spacedim> &user_flags_triangulation =
      const_cast<Triangulation<dim,spacedim>&> (dofs.get_tria());
    user_flags_triangulation.save_user_flags(old_flags);
    user_flags_triangulation.clear_user_flags();

    const typename DoFHandler<dim,spacedim>::cell_iterator end = dofs.end(level);
    typename DoFHandler<dim,spacedim>::active_cell_iterator cell;
    std::vector<types::global_dof_index> cell_indices;
    std::vector<types::global_dof_index> neighbor_indices;

    // We have to translate the
    // couplings from components to
    // blocks, so this works for
    // nonprimitive elements as well.
    std::vector<Table<2, DoFTools::Coupling> > couple_cell;
    std::vector<Table<2, DoFTools::Coupling> > couple_face;
    DoFTools::convert_couplings_to_blocks(dofs, couplings, couple_cell);
    DoFTools::convert_couplings_to_blocks(dofs, flux_couplings, couple_face);

    // We loop over cells and go from
    // cells to lower dimensional
    // objects. This is the only way to
    // cope withthe fact, that an
    // unknown number of cells may
    // share an object of dimension
    // smaller than dim-1.
    for (cell = dofs.begin_active(); cell != end; ++cell)
      {
        const FiniteElement<dim> &fe = cell->get_fe();
        const unsigned int fe_index = cell->active_fe_index();

        Assert (couplings.n_rows()==fe.n_components(),
                ExcDimensionMismatch(couplings.n_rows(), fe.n_components()));
        Assert (couplings.n_cols()==fe.n_components(),
                ExcDimensionMismatch(couplings.n_cols(), fe.n_components()));
        Assert (flux_couplings.n_rows()==fe.n_components(),
                ExcDimensionMismatch(flux_couplings.n_rows(), fe.n_components()));
        Assert (flux_couplings.n_cols()==fe.n_components(),
                ExcDimensionMismatch(flux_couplings.n_cols(), fe.n_components()));

        cell_indices.resize(fe.dofs_per_cell);
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
        // dofs_per_vertex is
        // arbitrary, not the other way
        // round.
        unsigned int increment;
        while (i < fe.first_line_index)
          {
            for (unsigned int base=0; base<fe.n_base_elements(); ++base)
              for (unsigned int mult=0; mult<fe.element_multiplicity(base); ++mult)
                if (couple_cell[fe_index](fe.system_to_block_index(i).first,
                                          fe.first_block_of_base(base) + mult) != DoFTools::none)
                  {
                    increment = fe.base_element(base).dofs_per_cell
                                - dim * fe.base_element(base).dofs_per_face;
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
        while (i < fe.first_quad_index)
          {
            for (unsigned int base=0; base<fe.n_base_elements(); ++base)
              for (unsigned int mult=0; mult<fe.element_multiplicity(base); ++mult)
                if (couple_cell[fe_index](fe.system_to_block_index(i).first,
                                          fe.first_block_of_base(base) + mult) != DoFTools::none)
                  {
                    increment = fe.base_element(base).dofs_per_cell
                                - ((dim>1)
                                   ? (dim-1)
                                   : GeometryInfo<dim>::faces_per_cell)
                                * fe.base_element(base).dofs_per_face;
                    row_lengths[cell_indices[i]] += increment;
                  }
            ++i;
          }

        // Now quads in 2D and 3D
        while (i < fe.first_hex_index)
          {
            for (unsigned int base=0; base<fe.n_base_elements(); ++base)
              for (unsigned int mult=0; mult<fe.element_multiplicity(base); ++mult)
                if (couple_cell[fe_index](fe.system_to_block_index(i).first,
                                          fe.first_block_of_base(base) + mult) != DoFTools::none)
                  {
                    increment = fe.base_element(base).dofs_per_cell
                                - ((dim>2)
                                   ? (dim-2)
                                   : GeometryInfo<dim>::faces_per_cell)
                                * fe.base_element(base).dofs_per_face;
                    row_lengths[cell_indices[i]] += increment;
                  }
            ++i;
          }

        // Finally, cells in 3D
        while (i < fe.dofs_per_cell)
          {
            for (unsigned int base=0; base<fe.n_base_elements(); ++base)
              for (unsigned int mult=0; mult<fe.element_multiplicity(base); ++mult)
                if (couple_cell[fe_index](fe.system_to_block_index(i).first,
                                          fe.first_block_of_base(base) + mult) != DoFTools::none)
                  {
                    increment = fe.base_element(base).dofs_per_cell
                                - GeometryInfo<dim>::faces_per_cell
                                * fe.base_element(base).dofs_per_face;
                    row_lengths[cell_indices[i]] += increment;
                  }
            ++i;
          }

        // At this point, we have
        // counted all dofs
        // contributiong from cells
        // coupled topologically to the
        // adjacent cells, but we
        // subtracted some faces.

        // Now, let's go by the faces
        // and add the missing
        // contribution as well as the
        // flux contributions.
        for (unsigned int iface=0; iface<GeometryInfo<dim>::faces_per_cell; ++iface)
          {
            bool level_boundary = cell->at_boundary(iface);
            typename DoFHandler<dim,spacedim>::cell_iterator neighbor;
            if (!level_boundary)
              {
                neighbor = cell->neighbor(iface);
                if (static_cast<unsigned int>(neighbor->level()) != level)
                  level_boundary = true;
              }

            if (level_boundary)
              {
                for (unsigned int local_dof=0; local_dof<fe.dofs_per_cell; ++local_dof)
                  row_lengths[cell_indices[local_dof]] += fe.dofs_per_face;
                continue;
              }

            const FiniteElement<dim> &nfe = neighbor->get_fe();
            typename DoFHandler<dim,spacedim>::face_iterator face = cell->face(iface);

            // Flux couplings are
            // computed from both sides
            // for simplicity.

            // The dofs on the common face
            // will be handled below,
            // therefore, we subtract them
            // here.
            for (unsigned int base=0; base<nfe.n_base_elements(); ++base)
              for (unsigned int mult=0; mult<nfe.element_multiplicity(base); ++mult)
                for (unsigned int local_dof=0; local_dof<fe.dofs_per_cell; ++local_dof)
                  if (couple_face[fe_index](fe.system_to_block_index(local_dof).first,
                                            nfe.first_block_of_base(base) + mult) != DoFTools::none)
                    {
                      const unsigned int dof_increment = nfe.base_element(base).dofs_per_cell
                                                         - nfe.base_element(base).dofs_per_face;
                      row_lengths[cell_indices[local_dof]] += dof_increment;
                    }

            // Do this only once per
            // face and not on the
            // hanging faces.
            if (face->user_flag_set())
              continue;
            face->set_user_flag();
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
            neighbor_indices.resize(nfe.dofs_per_cell);
            neighbor->get_mg_dof_indices(neighbor_indices);
            for (unsigned int base=0; base<nfe.n_base_elements(); ++base)
              for (unsigned int mult=0; mult<nfe.element_multiplicity(base); ++mult)
                for (unsigned int local_dof=0; local_dof<fe.dofs_per_cell; ++local_dof)
                  if (couple_cell[fe_index](fe.system_to_component_index(local_dof).first,
                                            nfe.first_block_of_base(base) + mult) != DoFTools::none)
                    row_lengths[cell_indices[local_dof]]
                    += nfe.base_element(base).dofs_per_face;
            for (unsigned int base=0; base<fe.n_base_elements(); ++base)
              for (unsigned int mult=0; mult<fe.element_multiplicity(base); ++mult)
                for (unsigned int local_dof=0; local_dof<nfe.dofs_per_cell; ++local_dof)
                  if (couple_cell[fe_index](nfe.system_to_component_index(local_dof).first,
                                            fe.first_block_of_base(base) + mult) != DoFTools::none)
                    row_lengths[neighbor_indices[local_dof]]
                    += fe.base_element(base).dofs_per_face;
          }
      }
    user_flags_triangulation.load_user_flags(old_flags);
  }



  template <class DH, class SparsityPattern>
  void make_sparsity_pattern (
    const DH &dof,
    SparsityPattern         &sparsity,
    const unsigned int       level)
  {
    const types::global_dof_index n_dofs = dof.n_dofs(level);

    Assert (sparsity.n_rows() == n_dofs,
            ExcDimensionMismatch (sparsity.n_rows(), n_dofs));
    Assert (sparsity.n_cols() == n_dofs,
            ExcDimensionMismatch (sparsity.n_cols(), n_dofs));

    const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
    std::vector<types::global_dof_index> dofs_on_this_cell(dofs_per_cell);
    typename DH::cell_iterator cell = dof.begin(level),
                               endc = dof.end(level);
    for (; cell!=endc; ++cell)
      if (dof.get_tria().locally_owned_subdomain()==numbers::invalid_subdomain_id
          || cell->level_subdomain_id()==dof.get_tria().locally_owned_subdomain())
        {
          cell->get_mg_dof_indices (dofs_on_this_cell);
          // make sparsity pattern for this cell
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              sparsity.add (dofs_on_this_cell[i],
                            dofs_on_this_cell[j]);
        }
  }



  template <int dim, class SparsityPattern, int spacedim>
  void
  make_flux_sparsity_pattern (
    const DoFHandler<dim,spacedim> &dof,
    SparsityPattern       &sparsity,
    const unsigned int level)
  {
    const types::global_dof_index n_dofs = dof.n_dofs(level);

    Assert (sparsity.n_rows() == n_dofs,
            ExcDimensionMismatch (sparsity.n_rows(), n_dofs));
    Assert (sparsity.n_cols() == n_dofs,
            ExcDimensionMismatch (sparsity.n_cols(), n_dofs));

    const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
    std::vector<types::global_dof_index> dofs_on_this_cell(dofs_per_cell);
    std::vector<types::global_dof_index> dofs_on_other_cell(dofs_per_cell);
    typename DoFHandler<dim,spacedim>::cell_iterator cell = dof.begin(level),
                                                     endc = dof.end(level);
    for (; cell!=endc; ++cell)
      {
        if (!cell->is_locally_owned_on_level()) continue;

        cell->get_mg_dof_indices (dofs_on_this_cell);
        // make sparsity pattern for this cell
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            sparsity.add (dofs_on_this_cell[i],
                          dofs_on_this_cell[j]);

        // Loop over all interior neighbors
        for (unsigned int face = 0;
             face < GeometryInfo<dim>::faces_per_cell;
             ++face)
          {
            if ( (! cell->at_boundary(face)) &&
                 (static_cast<unsigned int>(cell->neighbor_level(face)) == level) )
              {
                typename DoFHandler<dim,spacedim>::cell_iterator
                neighbor = cell->neighbor(face);
                neighbor->get_mg_dof_indices (dofs_on_other_cell);
                // only add one direction The other is taken care of by
                // neighbor (except when the neighbor is not owned by the same
                // processor)
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  {
                    for (unsigned int j=0; j<dofs_per_cell; ++j)
                      {
                        sparsity.add (dofs_on_this_cell[i],
                                      dofs_on_other_cell[j]);
                      }
                  }
                if (neighbor->is_locally_owned_on_level() == false)
                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    for (unsigned int j=0; j<dofs_per_cell; ++j)
                      {
                        sparsity.add (dofs_on_other_cell[i],
                                      dofs_on_other_cell[j]);
                        sparsity.add (dofs_on_other_cell[i],
                                      dofs_on_this_cell[j]);
                      }
              }
          }
      }
  }



  template <int dim, class SparsityPattern, int spacedim>
  void
  make_flux_sparsity_pattern_edge (
    const DoFHandler<dim,spacedim> &dof,
    SparsityPattern       &sparsity,
    const unsigned int level)
  {
    Assert ((level>=1) && (level<dof.get_tria().n_global_levels()),
            ExcIndexRange(level, 1, dof.get_tria().n_global_levels()));

    const types::global_dof_index fine_dofs = dof.n_dofs(level);
    const types::global_dof_index coarse_dofs = dof.n_dofs(level-1);

    // Matrix maps from fine level to coarse level

    Assert (sparsity.n_rows() == coarse_dofs,
            ExcDimensionMismatch (sparsity.n_rows(), coarse_dofs));
    Assert (sparsity.n_cols() == fine_dofs,
            ExcDimensionMismatch (sparsity.n_cols(), fine_dofs));

    const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
    std::vector<types::global_dof_index> dofs_on_this_cell(dofs_per_cell);
    std::vector<types::global_dof_index> dofs_on_other_cell(dofs_per_cell);
    typename DoFHandler<dim,spacedim>::cell_iterator cell = dof.begin(level),
                                                     endc = dof.end(level);
    for (; cell!=endc; ++cell)
      {
        if (!cell->is_locally_owned_on_level()) continue;

        cell->get_mg_dof_indices (dofs_on_this_cell);
        // Loop over all interior neighbors
        for (unsigned int face = 0;
             face < GeometryInfo<dim>::faces_per_cell;
             ++face)
          {
            // Neighbor is coarser

            if ( (! cell->at_boundary(face)) &&
                 (static_cast<unsigned int>(cell->neighbor_level(face)) != level) )
              {
                typename DoFHandler<dim,spacedim>::cell_iterator
                neighbor = cell->neighbor(face);
                neighbor->get_mg_dof_indices (dofs_on_other_cell);

                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  {
                    for (unsigned int j=0; j<dofs_per_cell; ++j)
                      {
                        sparsity.add (dofs_on_other_cell[i],
                                      dofs_on_this_cell[j]);
                        sparsity.add (dofs_on_other_cell[j],
                                      dofs_on_this_cell[i]);
                      }
                  }
              }
          }
      }
  }



  template <int dim, class SparsityPattern, int spacedim>
  void
  make_flux_sparsity_pattern (
    const DoFHandler<dim,spacedim> &dof,
    SparsityPattern       &sparsity,
    const unsigned int level,
    const Table<2,DoFTools::Coupling> &int_mask,
    const Table<2,DoFTools::Coupling> &flux_mask)
  {
    const FiniteElement<dim> &fe = dof.get_fe();
    const types::global_dof_index n_dofs = dof.n_dofs(level);
    const unsigned int n_comp = fe.n_components();

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

    const unsigned int total_dofs = fe.dofs_per_cell;
    std::vector<types::global_dof_index> dofs_on_this_cell(total_dofs);
    std::vector<types::global_dof_index> dofs_on_other_cell(total_dofs);
    Table<2,bool> support_on_face(total_dofs, GeometryInfo<dim>::faces_per_cell);

    typename DoFHandler<dim,spacedim>::cell_iterator cell = dof.begin(level),
                                                     endc = dof.end(level);

    const Table<2,DoFTools::Coupling>
    int_dof_mask  = DoFTools::dof_couplings_from_component_couplings(fe, int_mask),
    flux_dof_mask = DoFTools::dof_couplings_from_component_couplings(fe, flux_mask);

    for (unsigned int i=0; i<total_dofs; ++i)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        support_on_face(i,f) = fe.has_support_on_face(i,f);

    // Clear user flags because we will
    // need them. But first we save
    // them and make sure that we
    // restore them later such that at
    // the end of this function the
    // Triangulation will be in the
    // same state as it was at the
    // beginning of this function.
    std::vector<bool> user_flags;
    dof.get_tria().save_user_flags(user_flags);
    const_cast<Triangulation<dim,spacedim> &>(dof.get_tria()).clear_user_flags ();

    for (; cell!=endc; ++cell)
      {
        if (!cell->is_locally_owned_on_level()) continue;

        cell->get_mg_dof_indices (dofs_on_this_cell);
        // make sparsity pattern for this cell
        for (unsigned int i=0; i<total_dofs; ++i)
          for (unsigned int j=0; j<total_dofs; ++j)
            if (int_dof_mask[i][j] != DoFTools::none)
              sparsity.add (dofs_on_this_cell[i],
                            dofs_on_this_cell[j]);

        // Loop over all interior neighbors
        for (unsigned int face = 0;
             face < GeometryInfo<dim>::faces_per_cell;
             ++face)
          {
            typename DoFHandler<dim,spacedim>::face_iterator cell_face = cell->face(face);
            if (cell_face->user_flag_set ())
              continue;

            if (cell->at_boundary (face) )
              {
                for (unsigned int i=0; i<total_dofs; ++i)
                  {
                    const bool i_non_zero_i = support_on_face (i, face);
                    for (unsigned int j=0; j<total_dofs; ++j)
                      {
                        const bool j_non_zero_i = support_on_face (j, face);

                        if (flux_dof_mask(i,j) == DoFTools::always)
                          sparsity.add (dofs_on_this_cell[i],
                                        dofs_on_this_cell[j]);
                        if (flux_dof_mask(i,j) == DoFTools::nonzero
                            && i_non_zero_i && j_non_zero_i)
                          sparsity.add (dofs_on_this_cell[i],
                                        dofs_on_this_cell[j]);
                      }
                  }
              }
            else
              {
                typename DoFHandler<dim,spacedim>::cell_iterator
                neighbor = cell->neighbor(face);

                if (neighbor->level() < cell->level())
                  continue;

                unsigned int neighbor_face = cell->neighbor_of_neighbor(face);

                neighbor->get_mg_dof_indices (dofs_on_other_cell);
                for (unsigned int i=0; i<total_dofs; ++i)
                  {
                    const bool i_non_zero_i = support_on_face (i, face);
                    const bool i_non_zero_e = support_on_face (i, neighbor_face);
                    for (unsigned int j=0; j<total_dofs; ++j)
                      {
                        const bool j_non_zero_i = support_on_face (j, face);
                        const bool j_non_zero_e = support_on_face (j, neighbor_face);
                        if (flux_dof_mask(i,j) == DoFTools::always)
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
                        if (flux_dof_mask(i,j) == DoFTools::nonzero)
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

                        if (flux_dof_mask(j,i) == DoFTools::always)
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
                        if (flux_dof_mask(j,i) == DoFTools::nonzero)
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
                neighbor->face(neighbor_face)->set_user_flag ();
              }
          }
      }

    // finally restore the user flags
    const_cast<Triangulation<dim,spacedim> &>(dof.get_tria()).load_user_flags(user_flags);
  }



  template <int dim, class SparsityPattern, int spacedim>
  void
  make_flux_sparsity_pattern_edge (
    const DoFHandler<dim,spacedim> &dof,
    SparsityPattern       &sparsity,
    const unsigned int level,
    const Table<2,DoFTools::Coupling> &flux_mask)
  {
    const FiniteElement<dim> &fe = dof.get_fe();
    const unsigned int n_comp = fe.n_components();

    Assert ((level>=1) && (level<dof.get_tria().n_global_levels()),
            ExcIndexRange(level, 1, dof.get_tria().n_global_levels()));

    const types::global_dof_index fine_dofs = dof.n_dofs(level);
    const types::global_dof_index coarse_dofs = dof.n_dofs(level-1);

    // Matrix maps from fine level to coarse level

    Assert (sparsity.n_rows() == coarse_dofs,
            ExcDimensionMismatch (sparsity.n_rows(), coarse_dofs));
    Assert (sparsity.n_cols() == fine_dofs,
            ExcDimensionMismatch (sparsity.n_cols(), fine_dofs));
    Assert (flux_mask.n_rows() == n_comp,
            ExcDimensionMismatch (flux_mask.n_rows(), n_comp));
    Assert (flux_mask.n_cols() == n_comp,
            ExcDimensionMismatch (flux_mask.n_cols(), n_comp));

    const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
    std::vector<types::global_dof_index> dofs_on_this_cell(dofs_per_cell);
    std::vector<types::global_dof_index> dofs_on_other_cell(dofs_per_cell);
    Table<2,bool> support_on_face(dofs_per_cell, GeometryInfo<dim>::faces_per_cell);

    typename DoFHandler<dim,spacedim>::cell_iterator cell = dof.begin(level),
                                                     endc = dof.end(level);

    const Table<2,DoFTools::Coupling> flux_dof_mask
      = DoFTools::dof_couplings_from_component_couplings(fe, flux_mask);

    for (unsigned int i=0; i<dofs_per_cell; ++i)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        support_on_face(i,f) = fe.has_support_on_face(i,f);

    for (; cell!=endc; ++cell)
      {
        if (!cell->is_locally_owned_on_level()) continue;

        cell->get_mg_dof_indices (dofs_on_this_cell);
        // Loop over all interior neighbors
        for (unsigned int face = 0;
             face < GeometryInfo<dim>::faces_per_cell;
             ++face)
          {
            // Neighbor is coarser

            if ( (! cell->at_boundary(face)) &&
                 (static_cast<unsigned int>(cell->neighbor_level(face)) != level) )
              {
                typename DoFHandler<dim,spacedim>::cell_iterator
                neighbor = cell->neighbor(face);
                neighbor->get_mg_dof_indices (dofs_on_other_cell);

                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  {
                    for (unsigned int j=0; j<dofs_per_cell; ++j)
                      {
                        if (flux_dof_mask(i,j) != DoFTools::none)
                          {
                            sparsity.add (dofs_on_other_cell[i],
                                          dofs_on_this_cell[j]);
                            sparsity.add (dofs_on_other_cell[j],
                                          dofs_on_this_cell[i]);
                          }
                      }
                  }
              }
          }
      }
  }



  template <int dim, int spacedim>
  void
  count_dofs_per_component (const DoFHandler<dim,spacedim> &dof_handler,
                            std::vector<std::vector<types::global_dof_index> > &result,
                            bool                              only_once,
                            std::vector<unsigned int>         target_component)
  {
    const FiniteElement<dim> &fe = dof_handler.get_fe();
    const unsigned int n_components = fe.n_components();
    const unsigned int nlevels = dof_handler.get_tria().n_global_levels();

    Assert (result.size() == nlevels,
            ExcDimensionMismatch(result.size(), nlevels));

    if (target_component.size() == 0)
      {
        target_component.resize(n_components);
        for (unsigned int i=0; i<n_components; ++i)
          target_component[i] = i;
      }

    Assert(target_component.size() == n_components,
           ExcDimensionMismatch(target_component.size(), n_components));

    for (unsigned int l=0; l<nlevels; ++l)
      {
        result[l].resize (n_components);
        std::fill (result[l].begin(),result[l].end(), 0U);

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
            std::vector<std::vector<bool> >
            dofs_in_component (n_components,
                               std::vector<bool>(dof_handler.n_dofs(l),
                                                 false));
            std::vector<ComponentMask> component_select (n_components);
            Threads::TaskGroup<> tasks;
            for (unsigned int i=0; i<n_components; ++i)
              {
                void (*fun_ptr) (const unsigned int       level,
                                 const DoFHandler<dim,spacedim> &,
                                 const ComponentMask &,
                                 std::vector<bool> &)
                  = &DoFTools::extract_level_dofs<DoFHandler<dim,spacedim> >;

                std::vector<bool> tmp(n_components, false);
                tmp[i] = true;
                component_select[i] = ComponentMask(tmp);

                tasks += Threads::new_task (fun_ptr,
                                            l, dof_handler,
                                            component_select[i],
                                            dofs_in_component[i]);
              }
            tasks.join_all();

            // next count what we got
            unsigned int component = 0;
            for (unsigned int b=0; b<fe.n_base_elements(); ++b)
              {
                const FiniteElement<dim> &base = fe.base_element(b);
                // Dimension of base element
                unsigned int d = base.n_components();

                for (unsigned int m=0; m<fe.element_multiplicity(b); ++m)
                  {
                    for (unsigned int dd=0; dd<d; ++dd)
                      {
                        if (base.is_primitive() || (!only_once || dd==0))
                          result[l][target_component[component]]
                          += std::count(dofs_in_component[component].begin(),
                                        dofs_in_component[component].end(),
                                        true);
                        ++component;
                      }
                  }
              }
            // finally sanity check
            Assert (!dof_handler.get_fe().is_primitive()
                    ||
                    std::accumulate (result[l].begin(),
                                     result[l].end(), 0U)
                    ==
                    dof_handler.n_dofs(l),
                    ExcInternalError());
          }
      }
  }



  template <int dim, int spacedim>
  void
  count_dofs_per_component (const DoFHandler<dim,spacedim>        &dof_handler,
                            std::vector<std::vector<types::global_dof_index> > &result,
                            std::vector<unsigned int>            target_component)
  {
    count_dofs_per_component (dof_handler, result,
                              false, target_component);
  }



  template <class DH>
  void
  count_dofs_per_block (
    const DH     &dof_handler,
    std::vector<std::vector<types::global_dof_index> > &dofs_per_block,
    std::vector<unsigned int>  target_block)
  {
    const FiniteElement<DH::dimension,DH::space_dimension> &fe = dof_handler.get_fe();
    const unsigned int n_blocks = fe.n_blocks();
    const unsigned int n_levels = dof_handler.get_tria().n_global_levels();

    AssertDimension (dofs_per_block.size(), n_levels);

    for (unsigned int l=0; l<n_levels; ++l)
      std::fill (dofs_per_block[l].begin(), dofs_per_block[l].end(), 0U);
    // If the empty vector was given as
    // default argument, set up this
    // vector as identity.
    if (target_block.size()==0)
      {
        target_block.resize(n_blocks);
        for (unsigned int i=0; i<n_blocks; ++i)
          target_block[i] = i;
      }
    Assert(target_block.size()==n_blocks,
           ExcDimensionMismatch(target_block.size(),n_blocks));

    const unsigned int max_block
      = *std::max_element (target_block.begin(),
                           target_block.end());
    const unsigned int n_target_blocks = max_block + 1;

    for (unsigned int l=0; l<n_levels; ++l)
      AssertDimension (dofs_per_block[l].size(), n_target_blocks);

    // special case for only one
    // block. treat this first
    // since it does not require any
    // computations
    if (n_blocks == 1)
      {
        for (unsigned int l=0; l<n_levels; ++l)
          dofs_per_block[l][0] = dof_handler.n_dofs(l);
        return;
      }
    // otherwise determine the number
    // of dofs in each block
    // separately. do so in parallel
    for (unsigned int l=0; l<n_levels; ++l)
      {
        std::vector<std::vector<bool> >
        dofs_in_block (n_blocks, std::vector<bool>(dof_handler.n_dofs(l), false));
        std::vector<BlockMask> block_select (n_blocks);
        Threads::TaskGroup<> tasks;
        for (unsigned int i=0; i<n_blocks; ++i)
          {
            void (*fun_ptr) (const unsigned int level,
                             const DH &,
                             const BlockMask &,
                             std::vector<bool> &)
              = &DoFTools::extract_level_dofs<DH>;

            std::vector<bool> tmp(n_blocks, false);
            tmp[i] = true;
            block_select[i] = tmp;

            tasks += Threads::new_task (fun_ptr,
                                        l, dof_handler, block_select[i],
                                        dofs_in_block[i]);
          }
        tasks.join_all ();

        // next count what we got
        for (unsigned int block=0; block<fe.n_blocks(); ++block)
          dofs_per_block[l][target_block[block]]
          += std::count(dofs_in_block[block].begin(),
                        dofs_in_block[block].end(),
                        true);
      }
  }


  template <>
  void
  make_boundary_list(
    const DoFHandler<1,1> &,
    const FunctionMap<1>::type &,
    std::vector<std::set<types::global_dof_index> > &,
    const ComponentMask &)
  {
    Assert(false, ExcNotImplemented());
  }



  template <>
  void
  make_boundary_list(
    const DoFHandler<1,2> &,
    const FunctionMap<1>::type &,
    std::vector<std::set<types::global_dof_index> > &,
    const ComponentMask &)
  {
    Assert(false, ExcNotImplemented());
  }



  template <int dim, int spacedim>
  void
  make_boundary_list(
    const DoFHandler<dim,spacedim> &dof,
    const typename FunctionMap<dim>::type &function_map,
    std::vector<std::set<types::global_dof_index> > &boundary_indices,
    const ComponentMask &component_mask)
  {
    // if for whatever reason we were
    // passed an empty map, return
    // immediately
    if (function_map.size() == 0)
      return;

    const unsigned int n_levels = dof.get_tria().n_global_levels();



    const unsigned int n_components = DoFTools::n_components(dof);
    const bool          fe_is_system = (n_components != 1);

    AssertDimension (boundary_indices.size(), n_levels);

    std::vector<types::global_dof_index> local_dofs;
    local_dofs.reserve (DoFTools::max_dofs_per_face(dof));
    std::fill (local_dofs.begin (),
               local_dofs.end (),
               DoFHandler<dim,spacedim>::invalid_dof_index);

    // First, deal with the simpler
    // case when we have to identify
    // all boundary dofs
    if (component_mask.n_selected_components(n_components) == n_components)
      {
        typename DoFHandler<dim,spacedim>::cell_iterator
        cell = dof.begin(),
        endc = dof.end();
        for (; cell!=endc; ++cell)
          {
            if (dof.get_tria().locally_owned_subdomain()!=numbers::invalid_subdomain_id
                && cell->level_subdomain_id()==numbers::artificial_subdomain_id)
              continue;
            const FiniteElement<dim> &fe = cell->get_fe();
            const unsigned int level = cell->level();
            local_dofs.resize(fe.dofs_per_face);

            for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell;
                 ++face_no)
              if (cell->at_boundary(face_no) == true)
                {
                  const typename DoFHandler<dim,spacedim>::face_iterator
                  face = cell->face(face_no);
                  const types::boundary_id bi = face->boundary_indicator();
                  // Face is listed in
                  // boundary map
                  if (function_map.find(bi) != function_map.end())
                    {
                      face->get_mg_dof_indices(level, local_dofs);
                      for (unsigned int i=0; i<fe.dofs_per_face; ++i)
                        boundary_indices[level].insert(local_dofs[i]);
                    }
                }
          }
      }
    else
      {
        Assert (component_mask.n_selected_components(n_components) > 0,
                ExcMessage("It's probably worthwhile to select at least one component."));

        typename DoFHandler<dim,spacedim>::cell_iterator
        cell = dof.begin(),
        endc = dof.end();
        for (; cell!=endc; ++cell)
          if (dof.get_tria().locally_owned_subdomain()==numbers::invalid_subdomain_id
              || cell->level_subdomain_id()!=numbers::artificial_subdomain_id)
            for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell;
                 ++face_no)
              {
                if (cell->at_boundary(face_no) == false)
                  continue;

                const FiniteElement<dim> &fe = cell->get_fe();
                const unsigned int level = cell->level();

                // we can presently deal only with
                // primitive elements for boundary
                // values. this does not preclude
                // us using non-primitive elements
                // in components that we aren't
                // interested in, however. make
                // sure that all shape functions
                // that are non-zero for the
                // components we are interested in,
                // are in fact primitive
                for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
                  {
                    const ComponentMask &nonzero_component_array
                      = cell->get_fe().get_nonzero_components (i);
                    for (unsigned int c=0; c<n_components; ++c)
                      if ((nonzero_component_array[c] == true)
                          &&
                          (component_mask[c] == true))
                        Assert (cell->get_fe().is_primitive (i),
                                ExcMessage ("This function can only deal with requested boundary "
                                            "values that correspond to primitive (scalar) base "
                                            "elements"));
                  }

                typename DoFHandler<dim,spacedim>::face_iterator face = cell->face(face_no);
                const types::boundary_id boundary_component = face->boundary_indicator();
                if (function_map.find(boundary_component) != function_map.end())
                  // face is of the right component
                  {
                    // get indices, physical location and
                    // boundary values of dofs on this
                    // face
                    local_dofs.resize (fe.dofs_per_face);
                    face->get_mg_dof_indices (level, local_dofs);
                    if (fe_is_system)
                      {
                        // enter those dofs
                        // into the list that
                        // match the
                        // component
                        // signature. avoid
                        // the usual
                        // complication that
                        // we can't just use
                        // *_system_to_component_index
                        // for non-primitive
                        // FEs
                        for (unsigned int i=0; i<local_dofs.size(); ++i)
                          {
                            unsigned int component;
                            if (fe.is_primitive())
                              component = fe.face_system_to_component_index(i).first;
                            else
                              {
                                // non-primitive
                                // case. make
                                // sure that
                                // this
                                // particular
                                // shape
                                // function
                                // _is_
                                // primitive,
                                // and get at
                                // it's
                                // component. use
                                // usual
                                // trick to
                                // transfer
                                // face dof
                                // index to
                                // cell dof
                                // index
                                const unsigned int cell_i
                                  = (dim == 1 ?
                                     i
                                     :
                                     (dim == 2 ?
                                      (i<2*fe.dofs_per_vertex ? i : i+2*fe.dofs_per_vertex)
                                      :
                                      (dim == 3 ?
                                       (i<4*fe.dofs_per_vertex ?
                                        i
                                        :
                                        (i<4*fe.dofs_per_vertex+4*fe.dofs_per_line ?
                                         i+4*fe.dofs_per_vertex
                                         :
                                         i+4*fe.dofs_per_vertex+8*fe.dofs_per_line))
                                       :
                                       numbers::invalid_unsigned_int)));
                                Assert (cell_i < fe.dofs_per_cell, ExcInternalError());

                                // make sure
                                // that if
                                // this is
                                // not a
                                // primitive
                                // shape function,
                                // then all
                                // the
                                // corresponding
                                // components
                                // in the
                                // mask are
                                // not set
//                         if (!fe.is_primitive(cell_i))
//                           for (unsigned int c=0; c<n_components; ++c)
//                             if (fe.get_nonzero_components(cell_i)[c])
//                               Assert (component_mask[c] == false,
//                                       ExcFENotPrimitive());

// let's pick the first of possibly more than one non-zero
// components. if shape function is non-primitive, then we will ignore
// the result in the following anyway, otherwise there's only one
// non-zero component which we will use
                                component = fe.get_nonzero_components(cell_i).first_selected_component();
                              }

                            if (component_mask[component] == true)
                              boundary_indices[level].insert(local_dofs[i]);
                          }
                      }
                    else
                      for (unsigned int i=0; i<local_dofs.size(); ++i)
                        boundary_indices[level].insert(local_dofs[i]);
                  }
              }
      }
  }


  template <int dim, int spacedim>
  void
  make_boundary_list(const DoFHandler<dim,spacedim> &dof,
                     const typename FunctionMap<dim>::type &function_map,
                     std::vector<IndexSet> &boundary_indices,
                     const ComponentMask &component_mask)
  {
    Assert (boundary_indices.size() == dof.get_tria().n_global_levels(),
            ExcDimensionMismatch (boundary_indices.size(),
                                  dof.get_tria().n_global_levels()));

    std::vector<std::set<types::global_dof_index> >
    my_boundary_indices (dof.get_tria().n_global_levels());
    make_boundary_list (dof, function_map, my_boundary_indices, component_mask);
    for (unsigned int i=0; i<dof.get_tria().n_global_levels(); ++i)
      {
        boundary_indices[i] = IndexSet (dof.n_dofs(i));
        boundary_indices[i].add_indices (my_boundary_indices[i].begin(),
                                         my_boundary_indices[i].end());
      }
  }


  template <int dim, int spacedim>
  void
  extract_inner_interface_dofs (const DoFHandler<dim,spacedim> &mg_dof_handler,
                                std::vector<std::vector<bool> >  &interface_dofs)
  {
    std::vector<IndexSet> temp;
    temp.resize(interface_dofs.size());
    for (unsigned int l=0; l<interface_dofs.size(); ++l)
      temp[l] = IndexSet(interface_dofs[l].size());

    extract_inner_interface_dofs(mg_dof_handler, temp);

    for (unsigned int l=0; l<interface_dofs.size(); ++l)
      {
        Assert (interface_dofs[l].size() == mg_dof_handler.n_dofs(l),
                ExcDimensionMismatch (interface_dofs[l].size(),
                                      mg_dof_handler.n_dofs(l)));

        temp[l].fill_binary_vector(interface_dofs[l]);
      }
  }


  template <int dim, int spacedim>
  void
  extract_non_interface_dofs (const DoFHandler<dim,spacedim> &mg_dof_handler,
                              std::vector<std::set<types::global_dof_index> >  &non_interface_dofs)
  {
    Assert (non_interface_dofs.size() == mg_dof_handler.get_tria().n_global_levels(),
            ExcDimensionMismatch (non_interface_dofs.size(),
                                  mg_dof_handler.get_tria().n_global_levels()));

    const FiniteElement<dim,spacedim> &fe = mg_dof_handler.get_fe();

    const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
    const unsigned int   dofs_per_face   = fe.dofs_per_face;

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    std::vector<bool> cell_dofs(dofs_per_cell, false);
    std::vector<bool> cell_dofs_interface(dofs_per_cell, false);

    typename DoFHandler<dim>::cell_iterator cell = mg_dof_handler.begin(),
                                            endc = mg_dof_handler.end();


    for (; cell!=endc; ++cell)
      {
        if (mg_dof_handler.get_tria().locally_owned_subdomain()!=numbers::invalid_subdomain_id
            && cell->level_subdomain_id()!=mg_dof_handler.get_tria().locally_owned_subdomain())
          continue;

        std::fill (cell_dofs.begin(), cell_dofs.end(), false);
        std::fill (cell_dofs_interface.begin(), cell_dofs_interface.end(), false);

        for (unsigned int face_nr=0; face_nr<GeometryInfo<dim>::faces_per_cell; ++face_nr)
          {
            const typename DoFHandler<dim,spacedim>::face_iterator face = cell->face(face_nr);
            if (!face->at_boundary())
              {
                //interior face
                const typename DoFHandler<dim>::cell_iterator
                neighbor = cell->neighbor(face_nr);

                if ((neighbor->level() < cell->level()))
                  {
                    for (unsigned int j=0; j<dofs_per_face; ++j)
                      cell_dofs_interface[fe.face_to_cell_index(j,face_nr)] = true;
                  }
                else
                  {
                    for (unsigned int j=0; j<dofs_per_face; ++j)
                      cell_dofs[fe.face_to_cell_index(j,face_nr)] = true;
                  }
              }
            else
              {
                //boundary face
                for (unsigned int j=0; j<dofs_per_face; ++j)
                  cell_dofs[fe.face_to_cell_index(j,face_nr)] = true;
              }
          }

        const unsigned int level = cell->level();
        cell->get_mg_dof_indices (local_dof_indices);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          if (cell_dofs[i] && !cell_dofs_interface[i])
            non_interface_dofs[level].insert(local_dof_indices[i]);
      }
  }


  template <int dim, int spacedim>
  void
  extract_inner_interface_dofs (const DoFHandler<dim,spacedim> &mg_dof_handler,
                                std::vector<IndexSet>  &interface_dofs)
  {
    Assert (interface_dofs.size() == mg_dof_handler.get_tria().n_global_levels(),
            ExcDimensionMismatch (interface_dofs.size(),
                                  mg_dof_handler.get_tria().n_global_levels()));

    std::vector<std::vector<types::global_dof_index> >
    tmp_interface_dofs(interface_dofs.size());

    const FiniteElement<dim,spacedim> &fe = mg_dof_handler.get_fe();

    const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
    const unsigned int   dofs_per_face   = fe.dofs_per_face;

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    std::vector<bool> cell_dofs(dofs_per_cell, false);

    typename DoFHandler<dim>::cell_iterator cell = mg_dof_handler.begin(),
                                            endc = mg_dof_handler.end();

    for (; cell!=endc; ++cell)
      {
        if (mg_dof_handler.get_tria().locally_owned_subdomain()!=numbers::invalid_subdomain_id
            && cell->level_subdomain_id()!=mg_dof_handler.get_tria().locally_owned_subdomain())
          continue;

        bool has_coarser_neighbor = false;

        std::fill (cell_dofs.begin(), cell_dofs.end(), false);

        for (unsigned int face_nr=0; face_nr<GeometryInfo<dim>::faces_per_cell; ++face_nr)
          {
            const typename DoFHandler<dim,spacedim>::face_iterator face = cell->face(face_nr);
            if (!face->at_boundary())
              {
                //interior face
                const typename DoFHandler<dim>::cell_iterator
                neighbor = cell->neighbor(face_nr);

                // Do refinement face
                // from the coarse side
                if (neighbor->level() < cell->level())
                  {
                    for (unsigned int j=0; j<dofs_per_face; ++j)
                      cell_dofs[fe.face_to_cell_index(j,face_nr)] = true;

                    has_coarser_neighbor = true;
                  }
              }
          }

        if (has_coarser_neighbor == false)
          continue;

        const unsigned int level = cell->level();
        cell->get_mg_dof_indices (local_dof_indices);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            if (cell_dofs[i])
              tmp_interface_dofs[level].push_back(local_dof_indices[i]);
          }
      }

    for (unsigned int l=0; l<mg_dof_handler.get_tria().n_global_levels(); ++l)
      {
        interface_dofs[l].clear();
        std::sort(tmp_interface_dofs[l].begin(), tmp_interface_dofs[l].end());
        interface_dofs[l].add_indices(tmp_interface_dofs[l].begin(),
                                      std::unique(tmp_interface_dofs[l].begin(),
                                                  tmp_interface_dofs[l].end()));
        interface_dofs[l].compress();
      }

  }


  template <typename number>
  void
  apply_boundary_values (
    const std::set<types::global_dof_index> &boundary_dofs,
    SparseMatrix<number> &matrix,
    const bool preserve_symmetry,
    const bool /*ignore_zeros*/)
  {
    // this function is not documented and not tested in the testsuite
    // so it isn't quite clear what it's supposed to do. it also isn't
    // used anywhere else in the library. in avoiding the use of
    // a deprecated function, I therefore threw away the original function
    // and replaced it by the following, which I believe should work

    std::map<types::global_dof_index, double> boundary_values;
    for (std::set<types::global_dof_index>::const_iterator p=boundary_dofs.begin();
         p != boundary_dofs.end(); ++p)
      boundary_values[*p] = 0;

    Vector<number> dummy(matrix.m());
    MatrixTools::apply_boundary_values (boundary_values,
                                        matrix, dummy, dummy,
                                        preserve_symmetry);
  }



  template <typename number>
  void
  apply_boundary_values (
    const std::set<types::global_dof_index> &boundary_dofs,
    BlockSparseMatrix<number> &matrix,
    const bool preserve_symmetry)
  {
    Assert (matrix.n_block_rows() == matrix.n_block_cols(),
            ExcNotQuadratic());
    Assert (matrix.get_sparsity_pattern().get_row_indices() ==
            matrix.get_sparsity_pattern().get_column_indices(),
            ExcNotQuadratic());

    // this function is not documented and not tested in the testsuite
    // so it isn't quite clear what it's supposed to do. it also isn't
    // used anywhere else in the library. in avoiding the use of
    // a deprecated function, I therefore threw away the original function
    // and replaced it by the following, which I believe should work

    std::map<types::global_dof_index, double> boundary_values;
    for (std::set<types::global_dof_index>::const_iterator p=boundary_dofs.begin();
         p != boundary_dofs.end(); ++p)
      boundary_values[*p] = 0;

    BlockVector<number> dummy(matrix.n_block_rows());
    for (unsigned int i=0; i<matrix.n_block_rows(); ++i)
      dummy.block(i).reinit (matrix.block(i,i).m());
    dummy.collect_sizes();

    MatrixTools::apply_boundary_values (boundary_values,
                                        matrix, dummy, dummy,
                                        preserve_symmetry);
  }
}


// explicit instantiations
#include "mg_tools.inst"

namespace MGTools
{
  template void apply_boundary_values (
    const std::set<types::global_dof_index> &,
    SparseMatrix<float> &, const bool, const bool);
  template void apply_boundary_values (
    const std::set<types::global_dof_index> &,
    SparseMatrix<double> &, const bool, const bool);
  template void apply_boundary_values (
    const std::set<types::global_dof_index> &,
    BlockSparseMatrix<float> &, const bool);
  template void apply_boundary_values (
    const std::set<types::global_dof_index> &,
    BlockSparseMatrix<double> &, const bool);
}


DEAL_II_NAMESPACE_CLOSE

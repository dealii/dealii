// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/meshworker/scratch_data.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  template <int dim, int spacedim>
  ScratchData<dim, spacedim>::ScratchData(
    const Mapping<dim, spacedim>       &mapping,
    const FiniteElement<dim, spacedim> &fe,
    const Quadrature<dim>              &quadrature,
    const UpdateFlags                  &update_flags,
    const Quadrature<dim - 1>          &face_quadrature,
    const UpdateFlags                  &face_update_flags)
    : mapping(&mapping)
    , fe(&fe)
    , cell_quadrature(quadrature)
    , face_quadrature(face_quadrature)
    , hp_capability_enabled(false)
    , cell_update_flags(update_flags)
    , neighbor_cell_update_flags(update_flags)
    , face_update_flags(face_update_flags)
    , neighbor_face_update_flags(face_update_flags)
    , local_dof_indices(fe.n_dofs_per_cell())
    , neighbor_dof_indices(fe.n_dofs_per_cell())
  {}



  template <int dim, int spacedim>
  ScratchData<dim, spacedim>::ScratchData(
    const Mapping<dim, spacedim>       &mapping,
    const FiniteElement<dim, spacedim> &fe,
    const Quadrature<dim>              &quadrature,
    const UpdateFlags                  &update_flags,
    const UpdateFlags                  &neighbor_update_flags,
    const Quadrature<dim - 1>          &face_quadrature,
    const UpdateFlags                  &face_update_flags,
    const UpdateFlags                  &neighbor_face_update_flags)
    : mapping(&mapping)
    , fe(&fe)
    , cell_quadrature(quadrature)
    , face_quadrature(face_quadrature)
    , hp_capability_enabled(false)
    , cell_update_flags(update_flags)
    , neighbor_cell_update_flags(neighbor_update_flags)
    , face_update_flags(face_update_flags)
    , neighbor_face_update_flags(neighbor_face_update_flags)
    , local_dof_indices(fe.n_dofs_per_cell())
    , neighbor_dof_indices(fe.n_dofs_per_cell())
  {}



  template <int dim, int spacedim>
  ScratchData<dim, spacedim>::ScratchData(
    const FiniteElement<dim, spacedim> &fe,
    const Quadrature<dim>              &quadrature,
    const UpdateFlags                  &update_flags,
    const Quadrature<dim - 1>          &face_quadrature,
    const UpdateFlags                  &face_update_flags)
    : ScratchData(fe.reference_cell()
                    .template get_default_linear_mapping<dim, spacedim>(),
                  fe,
                  quadrature,
                  update_flags,
                  face_quadrature,
                  face_update_flags)
  {}



  template <int dim, int spacedim>
  ScratchData<dim, spacedim>::ScratchData(
    const FiniteElement<dim, spacedim> &fe,
    const Quadrature<dim>              &quadrature,
    const UpdateFlags                  &update_flags,
    const UpdateFlags                  &neighbor_update_flags,
    const Quadrature<dim - 1>          &face_quadrature,
    const UpdateFlags                  &face_update_flags,
    const UpdateFlags                  &neighbor_face_update_flags)
    : ScratchData(fe.reference_cell()
                    .template get_default_linear_mapping<dim, spacedim>(),
                  fe,
                  quadrature,
                  update_flags,
                  neighbor_update_flags,
                  face_quadrature,
                  face_update_flags,
                  neighbor_face_update_flags)
  {}



  template <int dim, int spacedim>
  ScratchData<dim, spacedim>::ScratchData(
    const hp::MappingCollection<dim, spacedim> &mapping_collection,
    const hp::FECollection<dim, spacedim>      &fe_collection,
    const hp::QCollection<dim>                 &cell_quadrature_collection,
    const UpdateFlags                          &cell_update_flags,
    const hp::QCollection<dim - 1>             &face_quadrature_collection,
    const UpdateFlags                          &face_update_flags)
    : mapping_collection(&mapping_collection)
    , fe_collection(&fe_collection)
    , cell_quadrature_collection(cell_quadrature_collection)
    , face_quadrature_collection(face_quadrature_collection)
    , hp_capability_enabled(true)
    , cell_update_flags(cell_update_flags)
    , neighbor_cell_update_flags(cell_update_flags)
    , face_update_flags(face_update_flags)
    , neighbor_face_update_flags(face_update_flags)
  {
    local_dof_indices.reserve(fe_collection.max_dofs_per_cell());
    neighbor_dof_indices.reserve(fe_collection.max_dofs_per_cell());
  }



  template <int dim, int spacedim>
  ScratchData<dim, spacedim>::ScratchData(
    const hp::MappingCollection<dim, spacedim> &mapping_collection,
    const hp::FECollection<dim, spacedim>      &fe_collection,
    const hp::QCollection<dim>                 &cell_quadrature_collection,
    const UpdateFlags                          &cell_update_flags,
    const UpdateFlags                          &neighbor_cell_update_flags,
    const hp::QCollection<dim - 1>             &face_quadrature_collection,
    const UpdateFlags                          &face_update_flags,
    const UpdateFlags                          &neighbor_face_update_flags)
    : mapping_collection(&mapping_collection)
    , fe_collection(&fe_collection)
    , cell_quadrature_collection(cell_quadrature_collection)
    , face_quadrature_collection(face_quadrature_collection)
    , hp_capability_enabled(true)
    , cell_update_flags(cell_update_flags)
    , neighbor_cell_update_flags(neighbor_cell_update_flags)
    , face_update_flags(face_update_flags)
    , neighbor_face_update_flags(neighbor_face_update_flags)
  {
    local_dof_indices.reserve(fe_collection.max_dofs_per_cell());
    neighbor_dof_indices.reserve(fe_collection.max_dofs_per_cell());
  }



  template <int dim, int spacedim>
  ScratchData<dim, spacedim>::ScratchData(
    const hp::FECollection<dim, spacedim> &fe_collection,
    const hp::QCollection<dim>            &cell_quadrature_collection,
    const UpdateFlags                     &cell_update_flags,
    const hp::QCollection<dim - 1>        &face_quadrature_collection,
    const UpdateFlags                     &face_update_flags)
    : ScratchData(fe_collection.get_reference_cell_default_linear_mapping(),
                  fe_collection,
                  cell_quadrature_collection,
                  cell_update_flags,
                  face_quadrature_collection,
                  face_update_flags)
  {}



  template <int dim, int spacedim>
  ScratchData<dim, spacedim>::ScratchData(
    const hp::FECollection<dim, spacedim> &fe_collection,
    const hp::QCollection<dim>            &cell_quadrature_collection,
    const UpdateFlags                     &cell_update_flags,
    const UpdateFlags                     &neighbor_cell_update_flags,
    const hp::QCollection<dim - 1>        &face_quadrature_collection,
    const UpdateFlags                     &face_update_flags,
    const UpdateFlags                     &neighbor_face_update_flags)
    : ScratchData(fe_collection.get_reference_cell_default_linear_mapping(),
                  fe_collection,
                  cell_quadrature_collection,
                  cell_update_flags,
                  neighbor_cell_update_flags,
                  face_quadrature_collection,
                  face_update_flags,
                  neighbor_face_update_flags)
  {}



  template <int dim, int spacedim>
  ScratchData<dim, spacedim>::ScratchData(
    const ScratchData<dim, spacedim> &scratch)
    : mapping(scratch.mapping)
    , fe(scratch.fe)
    , cell_quadrature(scratch.cell_quadrature)
    , face_quadrature(scratch.face_quadrature)
    , mapping_collection(scratch.mapping_collection)
    , fe_collection(scratch.fe_collection)
    , cell_quadrature_collection(scratch.cell_quadrature_collection)
    , face_quadrature_collection(scratch.face_quadrature_collection)
    , hp_capability_enabled(scratch.hp_capability_enabled)
    , cell_update_flags(scratch.cell_update_flags)
    , neighbor_cell_update_flags(scratch.neighbor_cell_update_flags)
    , face_update_flags(scratch.face_update_flags)
    , neighbor_face_update_flags(scratch.neighbor_face_update_flags)
    , local_dof_indices(scratch.local_dof_indices)
    , neighbor_dof_indices(scratch.neighbor_dof_indices)
    , user_data_storage(scratch.user_data_storage)
    , internal_data_storage(scratch.internal_data_storage)
  {}



  template <int dim, int spacedim>
  const FEValues<dim, spacedim> &
  ScratchData<dim, spacedim>::reinit(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell)
  {
    if (hp_capability_enabled == false)
      {
        if (!fe_values)
          fe_values = std::make_unique<FEValues<dim, spacedim>>(
            *mapping, *fe, cell_quadrature, cell_update_flags);

        fe_values->reinit(cell);
        local_dof_indices.resize(fe_values->dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);
        current_fe_values = fe_values.get();
        return *fe_values;
      }
    else
      {
        if (!hp_fe_values)
          hp_fe_values = std::make_unique<hp::FEValues<dim, spacedim>>(
            *mapping_collection,
            *fe_collection,
            cell_quadrature_collection,
            cell_update_flags);

        hp_fe_values->reinit(cell);
        const auto &fe_values = hp_fe_values->get_present_fe_values();

        AssertDimension(
          (*fe_collection)[cell->active_fe_index()].n_dofs_per_cell(),
          fe_values.dofs_per_cell);
        local_dof_indices.resize(fe_values.dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);

        current_fe_values = &fe_values;
        return fe_values;
      }
  }



  template <int dim, int spacedim>
  const FEFaceValues<dim, spacedim> &
  ScratchData<dim, spacedim>::reinit(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    const unsigned int                                              face_no)
  {
    if (hp_capability_enabled == false)
      {
        if (!fe_face_values)
          fe_face_values = std::make_unique<FEFaceValues<dim, spacedim>>(
            *mapping, *fe, face_quadrature, face_update_flags);

        fe_face_values->reinit(cell, face_no);
        local_dof_indices.resize(fe->n_dofs_per_cell());
        cell->get_dof_indices(local_dof_indices);
        current_fe_values = fe_face_values.get();
        return *fe_face_values;
      }
    else
      {
        return reinit(cell, cell, face_no);
      }
  }



  template <int dim, int spacedim>
  const FEFaceValues<dim, spacedim> &
  ScratchData<dim, spacedim>::reinit(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    const typename DoFHandler<dim, spacedim>::active_cell_iterator
                      &neighbor_cell,
    const unsigned int face_no)
  {
    Assert(hp_capability_enabled, ExcOnlyAvailableWithHP());

    if (!hp_fe_face_values)
      hp_fe_face_values = std::make_unique<hp::FEFaceValues<dim, spacedim>>(
        *mapping_collection,
        *fe_collection,
        face_quadrature_collection,
        face_update_flags);

    if (neighbor_cell == cell)
      {
        hp_fe_face_values->reinit(cell, face_no);
      }
    else
      {
        // When we want to ensure some agreement between the cell face and its
        // neighbor on the quadrature order and mapping to use on this face,
        // then we defer to the dominance of one FE over another. This should
        // ensure that the optimal integration order and mapping order are
        // selected for this situation.
        unsigned int dominated_fe_index = fe_collection->find_dominated_fe(
          {cell->active_fe_index(), neighbor_cell->active_fe_index()});

        // TODO: find_dominated_fe returns invalid_fe_index when no dominated FE
        // has been found. We want to pass this value to FEFaceValues, but it
        // expects an invalid_unsigned_int in this case. We need to match the
        // interfaces in the future.
        if (dominated_fe_index == numbers::invalid_fe_index)
          dominated_fe_index = numbers::invalid_unsigned_int;

        hp_fe_face_values->reinit(cell,
                                  face_no,
                                  dominated_fe_index,
                                  dominated_fe_index);
      }

    const auto &fe_face_values = hp_fe_face_values->get_present_fe_values();
    const auto &fe             = (*fe_collection)[cell->active_fe_index()];

    local_dof_indices.resize(fe.n_dofs_per_cell());
    cell->get_dof_indices(local_dof_indices);

    current_fe_values = &fe_face_values;
    return fe_face_values;
  }



  template <int dim, int spacedim>
  const FEFaceValuesBase<dim, spacedim> &
  ScratchData<dim, spacedim>::reinit(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    const unsigned int                                              face_no,
    const unsigned int                                              subface_no)
  {
    if (subface_no != numbers::invalid_unsigned_int)
      {
        if (hp_capability_enabled == false)
          {
            if (!fe_subface_values)
              fe_subface_values =
                std::make_unique<FESubfaceValues<dim, spacedim>>(
                  *mapping, *fe, face_quadrature, face_update_flags);
            fe_subface_values->reinit(cell, face_no, subface_no);
            local_dof_indices.resize(fe->n_dofs_per_cell());
            cell->get_dof_indices(local_dof_indices);

            current_fe_values = fe_subface_values.get();
            return *fe_subface_values;
          }
        else
          {
            return reinit(cell, cell, face_no, subface_no);
          }
      }
    else
      return reinit(cell, face_no);
  }



  template <int dim, int spacedim>
  const FEFaceValuesBase<dim, spacedim> &
  ScratchData<dim, spacedim>::reinit(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    const typename DoFHandler<dim, spacedim>::active_cell_iterator
                      &neighbor_cell,
    const unsigned int face_no,
    const unsigned int subface_no)
  {
    Assert(hp_capability_enabled, ExcOnlyAvailableWithHP());

    if (subface_no != numbers::invalid_unsigned_int)
      {
        if (!hp_fe_subface_values)
          hp_fe_subface_values =
            std::make_unique<hp::FESubfaceValues<dim, spacedim>>(
              *mapping_collection,
              *fe_collection,
              face_quadrature_collection,
              face_update_flags);

        if (neighbor_cell == cell)
          {
            hp_fe_subface_values->reinit(cell, face_no, subface_no);
          }
        else
          {
            // When we want to ensure some agreement between the cell face and
            // its neighbor on the quadrature order and mapping to use on this
            // face, then we defer to the dominance of one FE over another. This
            // should ensure that the optimal integration order and mapping
            // order are selected for this situation.
            unsigned int dominated_fe_index = fe_collection->find_dominated_fe(
              {cell->active_fe_index(), neighbor_cell->active_fe_index()});

            // TODO: find_dominated_fe returns invalid_fe_index when no
            // dominated FE has been found. We want to pass this value to
            // FEFaceValues, but it expects an invalid_unsigned_int in this
            // case. We need to match the interfaces in the future.
            if (dominated_fe_index == numbers::invalid_fe_index)
              dominated_fe_index = numbers::invalid_unsigned_int;

            hp_fe_subface_values->reinit(cell,
                                         face_no,
                                         subface_no,
                                         dominated_fe_index,
                                         dominated_fe_index);
          }

        const auto &fe_subface_values =
          hp_fe_subface_values->get_present_fe_values();
        const auto &fe = (*fe_collection)[cell->active_fe_index()];

        local_dof_indices.resize(fe.n_dofs_per_cell());
        cell->get_dof_indices(local_dof_indices);

        current_fe_values = &fe_subface_values;
        return fe_subface_values;
      }
    else
      {
        return reinit(cell, neighbor_cell, face_no);
      }
  }



  template <int dim, int spacedim>
  const FEInterfaceValues<dim, spacedim> &
  ScratchData<dim, spacedim>::reinit(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    const unsigned int                                              face_no,
    const typename DoFHandler<dim, spacedim>::active_cell_iterator
                      &cell_neighbor,
    const unsigned int face_no_neighbor)
  {
    return reinit(cell,
                  face_no,
                  numbers::invalid_unsigned_int,
                  cell_neighbor,
                  face_no_neighbor,
                  numbers::invalid_unsigned_int);
  }



  template <int dim, int spacedim>
  const FEInterfaceValues<dim, spacedim> &
  ScratchData<dim, spacedim>::reinit(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    const unsigned int                                              face_no,
    const unsigned int                                              sub_face_no,
    const typename DoFHandler<dim, spacedim>::active_cell_iterator
                      &cell_neighbor,
    const unsigned int face_no_neighbor,
    const unsigned int sub_face_no_neighbor)
  {
    if (hp_capability_enabled == false)
      {
        if (!interface_fe_values)
          interface_fe_values =
            std::make_unique<FEInterfaceValues<dim, spacedim>>(
              *mapping, *fe, face_quadrature, face_update_flags);

        interface_fe_values->reinit(cell,
                                    face_no,
                                    sub_face_no,
                                    cell_neighbor,
                                    face_no_neighbor,
                                    sub_face_no_neighbor);
      }
    else
      {
        if (!interface_fe_values)
          interface_fe_values =
            std::make_unique<FEInterfaceValues<dim, spacedim>>(
              *mapping_collection,
              *fe_collection,
              face_quadrature_collection,
              face_update_flags);

        // When we want to ensure some agreement between the cell face and
        // its neighbor on the quadrature order and mapping to use on this
        // face, then we defer to the dominance of one FE over another. This
        // should ensure that the optimal integration order and mapping
        // order are selected for this situation.
        unsigned int dominated_fe_index = fe_collection->find_dominated_fe(
          {cell->active_fe_index(), cell_neighbor->active_fe_index()});

        // TODO: find_dominated_fe returns invalid_fe_index when no dominated FE
        // has been found. We want to pass this value to FEFaceValues, but it
        // expects an invalid_unsigned_int in this case. We need to match the
        // interfaces in the future.
        if (dominated_fe_index == numbers::invalid_fe_index)
          dominated_fe_index = numbers::invalid_unsigned_int;

        interface_fe_values->reinit(cell,
                                    face_no,
                                    sub_face_no,
                                    cell_neighbor,
                                    face_no_neighbor,
                                    sub_face_no_neighbor,
                                    dominated_fe_index,
                                    dominated_fe_index);
      }

    current_fe_values          = &interface_fe_values->get_fe_face_values(0);
    current_neighbor_fe_values = &interface_fe_values->get_fe_face_values(1);

    local_dof_indices = interface_fe_values->get_interface_dof_indices();
    return *interface_fe_values;
  }



  template <int dim, int spacedim>
  const FEValues<dim, spacedim> &
  ScratchData<dim, spacedim>::reinit_neighbor(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell)
  {
    if (hp_capability_enabled == false)
      {
        if (!neighbor_fe_values)
          neighbor_fe_values = std::make_unique<FEValues<dim, spacedim>>(
            *mapping, *fe, cell_quadrature, neighbor_cell_update_flags);

        neighbor_fe_values->reinit(cell);
        cell->get_dof_indices(neighbor_dof_indices);
        current_neighbor_fe_values = neighbor_fe_values.get();
        return *neighbor_fe_values;
      }
    else
      {
        if (!neighbor_hp_fe_values)
          neighbor_hp_fe_values = std::make_unique<hp::FEValues<dim, spacedim>>(
            *mapping_collection,
            *fe_collection,
            cell_quadrature_collection,
            neighbor_cell_update_flags);

        neighbor_hp_fe_values->reinit(cell);
        const auto &neighbor_fe_values =
          neighbor_hp_fe_values->get_present_fe_values();

        AssertDimension(
          (*fe_collection)[cell->active_fe_index()].n_dofs_per_cell(),
          neighbor_fe_values.dofs_per_cell);
        neighbor_dof_indices.resize(neighbor_fe_values.dofs_per_cell);
        cell->get_dof_indices(neighbor_dof_indices);

        current_neighbor_fe_values = &neighbor_fe_values;
        return neighbor_fe_values;
      }
  }



  template <int dim, int spacedim>
  const FEFaceValues<dim, spacedim> &
  ScratchData<dim, spacedim>::reinit_neighbor(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    const unsigned int                                              face_no)
  {
    if (hp_capability_enabled == false)
      {
        if (!neighbor_fe_face_values)
          neighbor_fe_face_values =
            std::make_unique<FEFaceValues<dim, spacedim>>(
              *mapping, *fe, face_quadrature, neighbor_face_update_flags);
        neighbor_fe_face_values->reinit(cell, face_no);
        cell->get_dof_indices(neighbor_dof_indices);
        current_neighbor_fe_values = neighbor_fe_face_values.get();
        return *neighbor_fe_face_values;
      }
    else
      {
        return reinit_neighbor(cell, cell, face_no);
      }
  }



  template <int dim, int spacedim>
  const FEFaceValues<dim, spacedim> &
  ScratchData<dim, spacedim>::reinit_neighbor(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    const typename DoFHandler<dim, spacedim>::active_cell_iterator
                      &neighbor_cell,
    const unsigned int face_no)
  {
    Assert(hp_capability_enabled, ExcOnlyAvailableWithHP());

    if (!neighbor_hp_fe_face_values)
      neighbor_hp_fe_face_values =
        std::make_unique<hp::FEFaceValues<dim, spacedim>>(
          *mapping_collection,
          *fe_collection,
          face_quadrature_collection,
          neighbor_face_update_flags);

    if (neighbor_cell == cell)
      {
        neighbor_hp_fe_face_values->reinit(neighbor_cell, face_no);
      }
    else
      {
        // When we want to ensure some agreement between the cell face and its
        // neighbor on the quadrature order and mapping to use on this face,
        // then we defer to the dominance of one FE over another. This should
        // ensure that the optimal integration order and mapping order are
        // selected for this situation.
        unsigned int dominated_fe_index = fe_collection->find_dominated_fe(
          {cell->active_fe_index(), neighbor_cell->active_fe_index()});

        // TODO: find_dominated_fe returns invalid_fe_index when no dominated FE
        // has been found. We want to pass this value to FEFaceValues, but it
        // expects an invalid_unsigned_int in this case. We need to match the
        // interfaces in the future.
        if (dominated_fe_index == numbers::invalid_fe_index)
          dominated_fe_index = numbers::invalid_unsigned_int;

        neighbor_hp_fe_face_values->reinit(neighbor_cell,
                                           face_no,
                                           dominated_fe_index,
                                           dominated_fe_index);
      }

    const auto &neighbor_fe_face_values =
      neighbor_hp_fe_face_values->get_present_fe_values();
    const auto &neighbor_fe =
      (*fe_collection)[neighbor_cell->active_fe_index()];

    neighbor_dof_indices.resize(neighbor_fe.n_dofs_per_cell());
    neighbor_cell->get_dof_indices(neighbor_dof_indices);

    current_neighbor_fe_values = &neighbor_fe_face_values;
    return neighbor_fe_face_values;
  }



  template <int dim, int spacedim>
  const FEFaceValuesBase<dim, spacedim> &
  ScratchData<dim, spacedim>::reinit_neighbor(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    const unsigned int                                              face_no,
    const unsigned int                                              subface_no)
  {
    if (subface_no != numbers::invalid_unsigned_int)
      {
        if (hp_capability_enabled == false)
          {
            if (!neighbor_fe_subface_values)
              neighbor_fe_subface_values =
                std::make_unique<FESubfaceValues<dim, spacedim>>(
                  *mapping, *fe, face_quadrature, neighbor_face_update_flags);
            neighbor_fe_subface_values->reinit(cell, face_no, subface_no);
            cell->get_dof_indices(neighbor_dof_indices);
            current_neighbor_fe_values = neighbor_fe_subface_values.get();
            return *neighbor_fe_subface_values;
          }
        else
          {
            return reinit_neighbor(cell, cell, face_no, subface_no);
          }
      }
    else
      return reinit_neighbor(cell, face_no);
  }



  template <int dim, int spacedim>
  const FEFaceValuesBase<dim, spacedim> &
  ScratchData<dim, spacedim>::reinit_neighbor(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    const typename DoFHandler<dim, spacedim>::active_cell_iterator
                      &neighbor_cell,
    const unsigned int face_no,
    const unsigned int subface_no)
  {
    Assert(hp_capability_enabled, ExcOnlyAvailableWithHP());

    if (subface_no != numbers::invalid_unsigned_int)
      {
        if (!neighbor_hp_fe_subface_values)
          neighbor_hp_fe_subface_values =
            std::make_unique<hp::FESubfaceValues<dim, spacedim>>(
              *mapping_collection,
              *fe_collection,
              face_quadrature_collection,
              neighbor_face_update_flags);

        if (neighbor_cell == cell)
          {
            neighbor_hp_fe_subface_values->reinit(neighbor_cell,
                                                  face_no,
                                                  subface_no);
          }
        else
          {
            // When we want to ensure some agreement between the cell face and
            // its neighbor on the quadrature order and mapping to use on this
            // face, then we defer to the dominance of one FE over another. This
            // should ensure that the optimal integration order and mapping
            // order are selected for this situation.
            unsigned int dominated_fe_index = fe_collection->find_dominated_fe(
              {cell->active_fe_index(), neighbor_cell->active_fe_index()});

            // TODO: find_dominated_fe returns invalid_fe_index when no
            // dominated FE has been found. We want to pass this value to
            // FEFaceValues, but it expects an invalid_unsigned_int in this
            // case. We need to match the interfaces in the future.
            if (dominated_fe_index == numbers::invalid_fe_index)
              dominated_fe_index = numbers::invalid_unsigned_int;

            neighbor_hp_fe_subface_values->reinit(neighbor_cell,
                                                  face_no,
                                                  subface_no,
                                                  dominated_fe_index,
                                                  dominated_fe_index);
          }

        const auto &neighbor_fe_subface_values =
          neighbor_hp_fe_subface_values->get_present_fe_values();
        const auto &neighbor_fe =
          (*fe_collection)[neighbor_cell->active_fe_index()];

        neighbor_dof_indices.resize(neighbor_fe.n_dofs_per_cell());
        neighbor_cell->get_dof_indices(neighbor_dof_indices);

        current_neighbor_fe_values = &neighbor_fe_subface_values;
        return neighbor_fe_subface_values;
      }
    else
      {
        return reinit_neighbor(cell, neighbor_cell, face_no);
      }
  }



  template <int dim, int spacedim>
  const FEValuesBase<dim, spacedim> &
  ScratchData<dim, spacedim>::get_current_fe_values() const
  {
    Assert(current_fe_values != nullptr,
           ExcMessage("You have to initialize the cache using one of the "
                      "reinit functions first!"));
    return *current_fe_values;
  }



  template <int dim, int spacedim>
  const FEInterfaceValues<dim, spacedim> &
  ScratchData<dim, spacedim>::get_current_interface_fe_values() const
  {
    Assert(interface_fe_values != nullptr,
           ExcMessage("You have to initialize the cache using one of the "
                      "reinit functions first!"));
    return *interface_fe_values;
  }



  template <int dim, int spacedim>
  const FEValuesBase<dim, spacedim> &
  ScratchData<dim, spacedim>::get_current_neighbor_fe_values() const
  {
    Assert(current_neighbor_fe_values != nullptr,
           ExcMessage("You have to initialize the cache using one of the "
                      "reinit functions first!"));
    return *current_neighbor_fe_values;
  }



  template <int dim, int spacedim>
  const std::vector<Point<spacedim>> &
  ScratchData<dim, spacedim>::get_quadrature_points() const
  {
    return get_current_fe_values().get_quadrature_points();
  }



  template <int dim, int spacedim>
  const std::vector<double> &
  ScratchData<dim, spacedim>::get_JxW_values() const
  {
    return get_current_fe_values().get_JxW_values();
  }



  template <int dim, int spacedim>
  const std::vector<double> &
  ScratchData<dim, spacedim>::get_neighbor_JxW_values() const
  {
    return get_current_neighbor_fe_values().get_JxW_values();
  }



  template <int dim, int spacedim>
  const std::vector<Tensor<1, spacedim>> &
  ScratchData<dim, spacedim>::get_normal_vectors() const
  {
    return get_current_fe_values().get_normal_vectors();
  }



  template <int dim, int spacedim>
  const std::vector<Tensor<1, spacedim>> &
  ScratchData<dim, spacedim>::get_neighbor_normal_vectors()
  {
    return get_current_neighbor_fe_values().get_normal_vectors();
  }



  template <int dim, int spacedim>
  const std::vector<types::global_dof_index> &
  ScratchData<dim, spacedim>::get_local_dof_indices() const
  {
    return local_dof_indices;
  }



  template <int dim, int spacedim>
  unsigned int
  ScratchData<dim, spacedim>::n_dofs_per_cell() const
  {
    return local_dof_indices.size();
  }



  template <int dim, int spacedim>
  const std::vector<types::global_dof_index> &
  ScratchData<dim, spacedim>::get_neighbor_dof_indices() const
  {
    return neighbor_dof_indices;
  }



  template <int dim, int spacedim>
  unsigned int
  ScratchData<dim, spacedim>::n_neighbor_dofs_per_cell() const
  {
    return neighbor_dof_indices.size();
  }



  template <int dim, int spacedim>
  GeneralDataStorage &
  ScratchData<dim, spacedim>::get_general_data_storage()
  {
    return user_data_storage;
  }



  template <int dim, int spacedim>
  const GeneralDataStorage &
  ScratchData<dim, spacedim>::get_general_data_storage() const
  {
    return user_data_storage;
  }



  template <int dim, int spacedim>
  const Mapping<dim, spacedim> &
  ScratchData<dim, spacedim>::get_mapping() const
  {
    Assert(hp_capability_enabled == false, ExcOnlyAvailableWithoutHP());
    Assert(mapping, ExcNotInitialized());
    return *mapping;
  }



  template <int dim, int spacedim>
  const FiniteElement<dim, spacedim> &
  ScratchData<dim, spacedim>::get_fe() const
  {
    Assert(hp_capability_enabled == false, ExcOnlyAvailableWithoutHP());
    Assert(fe, ExcNotInitialized());
    return *fe;
  }



  template <int dim, int spacedim>
  const Quadrature<dim> &
  ScratchData<dim, spacedim>::get_cell_quadrature() const
  {
    Assert(hp_capability_enabled == false, ExcOnlyAvailableWithoutHP());
    return cell_quadrature;
  }



  template <int dim, int spacedim>
  const Quadrature<dim - 1> &
  ScratchData<dim, spacedim>::get_face_quadrature() const
  {
    Assert(hp_capability_enabled == false, ExcOnlyAvailableWithoutHP());
    return face_quadrature;
  }



  template <int dim, int spacedim>
  const hp::MappingCollection<dim, spacedim> &
  ScratchData<dim, spacedim>::get_mapping_collection() const
  {
    Assert(hp_capability_enabled == true, ExcOnlyAvailableWithHP());
    Assert(mapping_collection, ExcNotInitialized());
    return *mapping_collection;
  }



  template <int dim, int spacedim>
  const hp::FECollection<dim, spacedim> &
  ScratchData<dim, spacedim>::get_fe_collection() const
  {
    Assert(hp_capability_enabled == true, ExcOnlyAvailableWithHP());
    Assert(fe_collection, ExcNotInitialized());
    return *fe_collection;
  }



  template <int dim, int spacedim>
  const hp::QCollection<dim> &
  ScratchData<dim, spacedim>::get_cell_quadrature_collection() const
  {
    Assert(hp_capability_enabled == true, ExcOnlyAvailableWithHP());
    return cell_quadrature_collection;
  }



  template <int dim, int spacedim>
  const hp::QCollection<dim - 1> &
  ScratchData<dim, spacedim>::get_face_quadrature_collection() const
  {
    Assert(hp_capability_enabled == true, ExcOnlyAvailableWithHP());
    return face_quadrature_collection;
  }



  template <int dim, int spacedim>
  bool
  ScratchData<dim, spacedim>::has_hp_capabilities() const
  {
    return hp_capability_enabled;
  }



  template <int dim, int spacedim>
  UpdateFlags
  ScratchData<dim, spacedim>::get_cell_update_flags() const
  {
    return cell_update_flags;
  }



  template <int dim, int spacedim>
  UpdateFlags
  ScratchData<dim, spacedim>::get_neighbor_cell_update_flags() const
  {
    return neighbor_cell_update_flags;
  }



  template <int dim, int spacedim>
  UpdateFlags
  ScratchData<dim, spacedim>::get_face_update_flags() const
  {
    return face_update_flags;
  }



  template <int dim, int spacedim>
  UpdateFlags
  ScratchData<dim, spacedim>::get_neighbor_face_update_flags() const
  {
    return neighbor_face_update_flags;
  }

} // namespace MeshWorker

// Explicit instantiations
namespace MeshWorker
{
#include "meshworker/scratch_data.inst"
}
DEAL_II_NAMESPACE_CLOSE

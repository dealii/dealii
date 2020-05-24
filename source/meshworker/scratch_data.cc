// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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

#include <deal.II/meshworker/scratch_data.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  template <int dim, int spacedim>
  ScratchData<dim, spacedim>::ScratchData(
    const Mapping<dim, spacedim> &      mapping,
    const FiniteElement<dim, spacedim> &fe,
    const Quadrature<dim> &             quadrature,
    const UpdateFlags &                 update_flags,
    const Quadrature<dim - 1> &         face_quadrature,
    const UpdateFlags &                 face_update_flags)
    : mapping(&mapping)
    , fe(&fe)
    , cell_quadrature(quadrature)
    , face_quadrature(face_quadrature)
    , cell_update_flags(update_flags)
    , neighbor_cell_update_flags(update_flags)
    , face_update_flags(face_update_flags)
    , neighbor_face_update_flags(face_update_flags)
    , local_dof_indices(fe.dofs_per_cell)
    , neighbor_dof_indices(fe.dofs_per_cell)
  {}



  template <int dim, int spacedim>
  ScratchData<dim, spacedim>::ScratchData(
    const Mapping<dim, spacedim> &      mapping,
    const FiniteElement<dim, spacedim> &fe,
    const Quadrature<dim> &             quadrature,
    const UpdateFlags &                 update_flags,
    const UpdateFlags &                 neighbor_update_flags,
    const Quadrature<dim - 1> &         face_quadrature,
    const UpdateFlags &                 face_update_flags,
    const UpdateFlags &                 neighbor_face_update_flags)
    : mapping(&mapping)
    , fe(&fe)
    , cell_quadrature(quadrature)
    , face_quadrature(face_quadrature)
    , cell_update_flags(update_flags)
    , neighbor_cell_update_flags(neighbor_update_flags)
    , face_update_flags(face_update_flags)
    , neighbor_face_update_flags(neighbor_face_update_flags)
    , local_dof_indices(fe.dofs_per_cell)
    , neighbor_dof_indices(fe.dofs_per_cell)
  {}



  template <int dim, int spacedim>
  ScratchData<dim, spacedim>::ScratchData(
    const FiniteElement<dim, spacedim> &fe,
    const Quadrature<dim> &             quadrature,
    const UpdateFlags &                 update_flags,
    const Quadrature<dim - 1> &         face_quadrature,
    const UpdateFlags &                 face_update_flags)
    : ScratchData(StaticMappingQ1<dim, spacedim>::mapping,
                  fe,
                  quadrature,
                  update_flags,
                  face_quadrature,
                  face_update_flags)
  {}



  template <int dim, int spacedim>
  ScratchData<dim, spacedim>::ScratchData(
    const FiniteElement<dim, spacedim> &fe,
    const Quadrature<dim> &             quadrature,
    const UpdateFlags &                 update_flags,
    const UpdateFlags &                 neighbor_update_flags,
    const Quadrature<dim - 1> &         face_quadrature,
    const UpdateFlags &                 face_update_flags,
    const UpdateFlags &                 neighbor_face_update_flags)
    : ScratchData(StaticMappingQ1<dim, spacedim>::mapping,
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
    const ScratchData<dim, spacedim> &scratch)
    : mapping(scratch.mapping)
    , fe(scratch.fe)
    , cell_quadrature(scratch.cell_quadrature)
    , face_quadrature(scratch.face_quadrature)
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
    if (!fe_values)
      fe_values = std::make_unique<FEValues<dim, spacedim>>(*mapping,
                                                            *fe,
                                                            cell_quadrature,
                                                            cell_update_flags);

    fe_values->reinit(cell);
    cell->get_dof_indices(local_dof_indices);
    current_fe_values = fe_values.get();
    return *fe_values;
  }



  template <int dim, int spacedim>
  const FEFaceValues<dim, spacedim> &
  ScratchData<dim, spacedim>::reinit(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    const unsigned int                                              face_no)
  {
    if (!fe_face_values)
      fe_face_values = std::make_unique<FEFaceValues<dim, spacedim>>(
        *mapping, *fe, face_quadrature, face_update_flags);

    fe_face_values->reinit(cell, face_no);
    cell->get_dof_indices(local_dof_indices);
    current_fe_values = fe_face_values.get();
    return *fe_face_values;
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
        if (!fe_subface_values)
          fe_subface_values = std::make_unique<FESubfaceValues<dim, spacedim>>(
            *mapping, *fe, face_quadrature, face_update_flags);
        fe_subface_values->reinit(cell, face_no, subface_no);
        cell->get_dof_indices(local_dof_indices);

        current_fe_values = fe_subface_values.get();
        return *fe_subface_values;
      }
    else
      return reinit(cell, face_no);
  }



  template <int dim, int spacedim>
  const FEValues<dim, spacedim> &
  ScratchData<dim, spacedim>::reinit_neighbor(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell)
  {
    if (!neighbor_fe_values)
      neighbor_fe_values = std::make_unique<FEValues<dim, spacedim>>(
        *mapping, *fe, cell_quadrature, neighbor_cell_update_flags);

    neighbor_fe_values->reinit(cell);
    cell->get_dof_indices(neighbor_dof_indices);
    current_neighbor_fe_values = neighbor_fe_values.get();
    return *neighbor_fe_values;
  }



  template <int dim, int spacedim>
  const FEFaceValues<dim, spacedim> &
  ScratchData<dim, spacedim>::reinit_neighbor(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    const unsigned int                                              face_no)
  {
    if (!neighbor_fe_face_values)
      neighbor_fe_face_values = std::make_unique<FEFaceValues<dim, spacedim>>(
        *mapping, *fe, face_quadrature, neighbor_face_update_flags);
    neighbor_fe_face_values->reinit(cell, face_no);
    cell->get_dof_indices(neighbor_dof_indices);
    current_neighbor_fe_values = neighbor_fe_face_values.get();
    return *neighbor_fe_face_values;
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
      return reinit_neighbor(cell, face_no);
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
  const std::vector<types::global_dof_index> &
  ScratchData<dim, spacedim>::get_neighbor_dof_indices() const
  {
    return neighbor_dof_indices;
  }



  template <int dim, int spacedim>
  GeneralDataStorage &
  ScratchData<dim, spacedim>::get_general_data_storage()
  {
    return user_data_storage;
  }



  template <int dim, int spacedim>
  const Mapping<dim, spacedim> &
  ScratchData<dim, spacedim>::get_mapping() const
  {
    return *mapping;
  }

} // namespace MeshWorker
DEAL_II_NAMESPACE_CLOSE

// Explicit instantiations
DEAL_II_NAMESPACE_OPEN
namespace MeshWorker
{
#include "scratch_data.inst"
}
DEAL_II_NAMESPACE_CLOSE

// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2022 by the deal.II authors
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
#include <deal.II/base/numbers.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_epetra_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/non_matching/fe_values.h>

DEAL_II_NAMESPACE_OPEN

namespace NonMatching
{
  RegionUpdateFlags::RegionUpdateFlags()
    : inside(update_default)
    , outside(update_default)
    , surface(update_default)
  {}



  template <int dim>
  template <class VectorType>
  FEValues<dim>::FEValues(const hp::FECollection<dim> &fe_collection,
                          const Quadrature<1> &        quadrature,
                          const RegionUpdateFlags      region_update_flags,
                          const MeshClassifier<dim> &  mesh_classifier,
                          const DoFHandler<dim> &      dof_handler,
                          const VectorType &           level_set,
                          const AdditionalData &       additional_data)
    : mapping_collection(&dealii::hp::StaticMappingQ1<dim>::mapping_collection)
    , fe_collection(&fe_collection)
    , q_collection_1D(quadrature)
    , region_update_flags(region_update_flags)
    , mesh_classifier(&mesh_classifier)
    , quadrature_generator(q_collection_1D,
                           dof_handler,
                           level_set,
                           additional_data)
  {
    // Tensor products of each quadrature in q_collection_1D. Used on the
    // non-intersected cells.
    hp::QCollection<dim> q_collection;
    for (unsigned int i = 0; i < q_collection_1D.size(); ++i)
      q_collection.push_back(Quadrature<dim>(q_collection_1D[i]));

    initialize(q_collection);
  }



  template <int dim>
  template <class VectorType>
  FEValues<dim>::FEValues(const hp::MappingCollection<dim> &mapping_collection,
                          const hp::FECollection<dim> &     fe_collection,
                          const hp::QCollection<dim> &      q_collection,
                          const hp::QCollection<1> &        q_collection_1D,
                          const RegionUpdateFlags           region_update_flags,
                          const MeshClassifier<dim> &       mesh_classifier,
                          const DoFHandler<dim> &           dof_handler,
                          const VectorType &                level_set,
                          const AdditionalData &            additional_data)
    : mapping_collection(&mapping_collection)
    , fe_collection(&fe_collection)
    , q_collection_1D(q_collection_1D)
    , region_update_flags(region_update_flags)
    , mesh_classifier(&mesh_classifier)
    , quadrature_generator(q_collection_1D,
                           dof_handler,
                           level_set,
                           additional_data)
  {
    initialize(q_collection);
  }



  template <int dim>
  void
  FEValues<dim>::initialize(const hp::QCollection<dim> &q_collection)
  {
    current_cell_location = LocationToLevelSet::unassigned;
    active_fe_index       = numbers::invalid_unsigned_int;

    Assert(fe_collection->size() > 0,
           ExcMessage("Incoming hp::FECollection can not be empty."));
    Assert(mapping_collection->size() == fe_collection->size() ||
             mapping_collection->size() == 1,
           ExcMessage("Size of hp::MappingCollection must be "
                      "the same as hp::FECollection or 1."));
    Assert(q_collection.size() == fe_collection->size() ||
             q_collection.size() == 1,
           ExcMessage("Size of hp::QCollection<dim> must be the "
                      "same as hp::FECollection or 1."));
    Assert(q_collection_1D.size() == fe_collection->size() ||
             q_collection_1D.size() == 1,
           ExcMessage("Size of hp::QCollection<1> must be the "
                      "same as hp::FECollection or 1."));

    // For each element in fe_collection, create dealii::FEValues objects to use
    // on the non-intersected cells.
    fe_values_inside_full_quadrature.resize(fe_collection->size());
    fe_values_outside_full_quadrature.resize(fe_collection->size());
    for (unsigned int fe_index = 0; fe_index < fe_collection->size();
         ++fe_index)
      {
        const unsigned int mapping_index =
          mapping_collection->size() > 1 ? fe_index : 0;
        const unsigned int q_index = q_collection.size() > 1 ? fe_index : 0;

        fe_values_inside_full_quadrature[fe_index].emplace(
          (*mapping_collection)[mapping_index],
          (*fe_collection)[fe_index],
          q_collection[q_index],
          region_update_flags.inside);
        fe_values_outside_full_quadrature[fe_index].emplace(
          (*mapping_collection)[mapping_index],
          (*fe_collection)[fe_index],
          q_collection[q_index],
          region_update_flags.outside);
      }
  }



  template <int dim>
  template <bool level_dof_access>
  void
  FEValues<dim>::reinit(
    const TriaIterator<DoFCellAccessor<dim, dim, level_dof_access>> &cell)
  {
    current_cell_location = mesh_classifier->location_to_level_set(cell);
    active_fe_index       = cell->active_fe_index();

    // These objects were created with a quadrature based on the previous cell
    // and are thus no longer valid.
    fe_values_inside.reset();
    fe_values_surface.reset();
    fe_values_outside.reset();

    switch (current_cell_location)
      {
        case LocationToLevelSet::inside:
          {
            fe_values_inside_full_quadrature.at(active_fe_index)->reinit(cell);
            break;
          }
        case LocationToLevelSet::outside:
          {
            fe_values_outside_full_quadrature.at(active_fe_index)->reinit(cell);
            break;
          }
        case LocationToLevelSet::intersected:
          {
            const unsigned int mapping_index =
              mapping_collection->size() > 1 ? active_fe_index : 0;

            const unsigned int q1D_index =
              q_collection_1D.size() > 1 ? active_fe_index : 0;
            quadrature_generator.set_1D_quadrature(q1D_index);
            quadrature_generator.generate(cell);

            const Quadrature<dim> &inside_quadrature =
              quadrature_generator.get_inside_quadrature();
            const Quadrature<dim> &outside_quadrature =
              quadrature_generator.get_outside_quadrature();
            const ImmersedSurfaceQuadrature<dim> &surface_quadrature =
              quadrature_generator.get_surface_quadrature();

            // Even if a cell is formally intersected the number of created
            // quadrature points can be 0. Avoid creating an FEValues object
            // if that is the case.
            if (inside_quadrature.size() > 0)
              {
                fe_values_inside.emplace((*mapping_collection)[mapping_index],
                                         (*fe_collection)[active_fe_index],
                                         inside_quadrature,
                                         region_update_flags.inside);

                fe_values_inside->reinit(cell);
              }

            if (outside_quadrature.size() > 0)
              {
                fe_values_outside.emplace((*mapping_collection)[mapping_index],
                                          (*fe_collection)[active_fe_index],
                                          outside_quadrature,
                                          region_update_flags.outside);

                fe_values_outside->reinit(cell);
              }

            if (surface_quadrature.size() > 0)
              {
                fe_values_surface.emplace((*mapping_collection)[mapping_index],
                                          (*fe_collection)[active_fe_index],
                                          surface_quadrature,
                                          region_update_flags.surface);
                fe_values_surface->reinit(cell);
              }

            break;
          }
        default:
          {
            Assert(false, ExcInternalError());
            break;
          }
      }
  }



  template <int dim>
  const std_cxx17::optional<dealii::FEValues<dim>> &
  FEValues<dim>::get_inside_fe_values() const
  {
    if (current_cell_location == LocationToLevelSet::inside)
      return fe_values_inside_full_quadrature.at(active_fe_index);
    else
      return fe_values_inside;
  }



  template <int dim>
  const std_cxx17::optional<dealii::FEValues<dim>> &
  FEValues<dim>::get_outside_fe_values() const
  {
    if (current_cell_location == LocationToLevelSet::outside)
      return fe_values_outside_full_quadrature.at(active_fe_index);
    else
      return fe_values_outside;
  }



  template <int dim>
  const std_cxx17::optional<FEImmersedSurfaceValues<dim>> &
  FEValues<dim>::get_surface_fe_values() const
  {
    return fe_values_surface;
  }



  template <int dim>
  template <class VectorType>
  FEInterfaceValues<dim>::FEInterfaceValues(
    const hp::FECollection<dim> &fe_collection,
    const Quadrature<1> &        quadrature,
    const RegionUpdateFlags      region_update_flags,
    const MeshClassifier<dim> &  mesh_classifier,
    const DoFHandler<dim> &      dof_handler,
    const VectorType &           level_set,
    const AdditionalData &       additional_data)
    : mapping_collection(&dealii::hp::StaticMappingQ1<dim>::mapping_collection)
    , fe_collection(&fe_collection)
    , q_collection_1D(quadrature)
    , region_update_flags(region_update_flags)
    , mesh_classifier(&mesh_classifier)
    , face_quadrature_generator(q_collection_1D,
                                dof_handler,
                                level_set,
                                additional_data)
  {
    // Tensor products of each quadrature in q_collection_1D. Used on the
    // non-intersected cells.
    hp::QCollection<dim - 1> q_collection;
    for (unsigned int i = 0; i < q_collection_1D.size(); ++i)
      q_collection.push_back(Quadrature<dim - 1>(q_collection_1D[i]));

    initialize(q_collection);
  }



  template <int dim>
  template <class VectorType>
  FEInterfaceValues<dim>::FEInterfaceValues(
    const hp::MappingCollection<dim> &mapping_collection,
    const hp::FECollection<dim> &     fe_collection,
    const hp::QCollection<dim - 1> &  q_collection,
    const hp::QCollection<1> &        q_collection_1D,
    const RegionUpdateFlags           region_update_flags,
    const MeshClassifier<dim> &       mesh_classifier,
    const DoFHandler<dim> &           dof_handler,
    const VectorType &                level_set,
    const AdditionalData &            additional_data)
    : mapping_collection(&mapping_collection)
    , fe_collection(&fe_collection)
    , q_collection_1D(q_collection_1D)
    , region_update_flags(region_update_flags)
    , mesh_classifier(&mesh_classifier)
    , face_quadrature_generator(q_collection_1D,
                                dof_handler,
                                level_set,
                                additional_data)
  {
    initialize(q_collection);
  }



  template <int dim>
  void
  FEInterfaceValues<dim>::initialize(
    const hp::QCollection<dim - 1> &q_collection)
  {
    current_face_location = LocationToLevelSet::unassigned;
    active_fe_index       = numbers::invalid_unsigned_int;

    Assert(fe_collection->size() > 0,
           ExcMessage("Incoming hp::FECollection can not be empty."));
    Assert(
      mapping_collection->size() == fe_collection->size() ||
        mapping_collection->size() == 1,
      ExcMessage(
        "Size of hp::MappingCollection must be the same as hp::FECollection or 1."));
    Assert(
      q_collection.size() == fe_collection->size() || q_collection.size() == 1,
      ExcMessage(
        "Size of hp::QCollection<dim> must be the same as hp::FECollection or 1."));
    Assert(
      q_collection_1D.size() == fe_collection->size() ||
        q_collection_1D.size() == 1,
      ExcMessage(
        "Size of hp::QCollection<1> must be the same as hp::FECollection or 1."));

    // For each element in fe_collection, create dealii::FEInterfaceValues
    // objects to use on the non-intersected cells.
    fe_values_inside_full_quadrature.resize(fe_collection->size());
    fe_values_outside_full_quadrature.resize(fe_collection->size());
    for (unsigned int fe_index = 0; fe_index < fe_collection->size();
         ++fe_index)
      {
        const unsigned int mapping_index =
          mapping_collection->size() > 1 ? fe_index : 0;
        const unsigned int q_index = q_collection.size() > 1 ? fe_index : 0;

        fe_values_inside_full_quadrature[fe_index].emplace(
          (*mapping_collection)[mapping_index],
          (*fe_collection)[fe_index],
          q_collection[q_index],
          region_update_flags.inside);
        fe_values_outside_full_quadrature[fe_index].emplace(
          (*mapping_collection)[mapping_index],
          (*fe_collection)[fe_index],
          q_collection[q_index],
          region_update_flags.outside);
      }
  }



  template <int dim>
  template <bool level_dof_access>
  void
  FEInterfaceValues<dim>::do_reinit(
    const TriaIterator<DoFCellAccessor<dim, dim, level_dof_access>> &cell,
    const unsigned int                                               face_no,
    const std::function<void(dealii::FEInterfaceValues<dim> &)> &call_reinit)
  {
    current_face_location =
      mesh_classifier->location_to_level_set(cell, face_no);
    active_fe_index = cell->active_fe_index();

    // These objects were created with a quadrature based on the previous cell
    // and are thus no longer valid.
    fe_values_inside.reset();
    fe_values_outside.reset();

    switch (current_face_location)
      {
        case LocationToLevelSet::inside:
          {
            call_reinit(*fe_values_inside_full_quadrature.at(active_fe_index));
            break;
          }
        case LocationToLevelSet::outside:
          {
            call_reinit(*fe_values_outside_full_quadrature.at(active_fe_index));
            break;
          }
        case LocationToLevelSet::intersected:
          {
            const unsigned int mapping_index =
              mapping_collection->size() > 1 ? active_fe_index : 0;
            const unsigned int q1D_index =
              q_collection_1D.size() > 1 ? active_fe_index : 0;

            face_quadrature_generator.set_1D_quadrature(q1D_index);
            face_quadrature_generator.generate(cell, face_no);

            const Quadrature<dim - 1> &inside_quadrature =
              face_quadrature_generator.get_inside_quadrature();
            const Quadrature<dim - 1> &outside_quadrature =
              face_quadrature_generator.get_outside_quadrature();

            // Even if a cell is formally intersected the number of created
            // quadrature points can be 0. Avoid creating an FEInterfaceValues
            // object if that is the case.
            if (inside_quadrature.size() > 0)
              {
                fe_values_inside.emplace((*mapping_collection)[mapping_index],
                                         (*fe_collection)[active_fe_index],
                                         inside_quadrature,
                                         region_update_flags.inside);

                call_reinit(*fe_values_inside);
              }

            if (outside_quadrature.size() > 0)
              {
                fe_values_outside.emplace((*mapping_collection)[mapping_index],
                                          (*fe_collection)[active_fe_index],
                                          outside_quadrature,
                                          region_update_flags.outside);

                call_reinit(*fe_values_outside);
              }
            break;
          }
        default:
          {
            Assert(false, ExcInternalError());
            break;
          }
      }
  }



  template <int dim>
  const std_cxx17::optional<dealii::FEInterfaceValues<dim>> &
  FEInterfaceValues<dim>::get_inside_fe_values() const
  {
    if (current_face_location == LocationToLevelSet::inside)
      return fe_values_inside_full_quadrature.at(active_fe_index);
    else
      return fe_values_inside;
  }



  template <int dim>
  const std_cxx17::optional<dealii::FEInterfaceValues<dim>> &
  FEInterfaceValues<dim>::get_outside_fe_values() const
  {
    if (current_face_location == LocationToLevelSet::outside)
      return fe_values_outside_full_quadrature.at(active_fe_index);
    else
      return fe_values_outside;
  }


#include "fe_values.inst"

} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE

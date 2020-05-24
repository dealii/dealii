// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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

#include <deal.II/base/thread_management.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/hp/fe_values.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

namespace hp
{
  // -------------------------- FEValuesBase -------------------------

  template <int dim, int q_dim, class FEValuesType>
  FEValuesBase<dim, q_dim, FEValuesType>::FEValuesBase(
    const MappingCollection<dim, FEValuesType::space_dimension>
      &                                                     mapping_collection,
    const FECollection<dim, FEValuesType::space_dimension> &fe_collection,
    const QCollection<q_dim> &                              q_collection,
    const UpdateFlags                                       update_flags)
    : fe_collection(&fe_collection)
    , mapping_collection(&mapping_collection)
    , q_collection(q_collection)
    , fe_values_table(fe_collection.size(),
                      mapping_collection.size(),
                      q_collection.size())
    , present_fe_values_index(numbers::invalid_unsigned_int,
                              numbers::invalid_unsigned_int,
                              numbers::invalid_unsigned_int)
    , update_flags(update_flags)
  {}


  template <int dim, int q_dim, class FEValuesType>
  FEValuesBase<dim, q_dim, FEValuesType>::FEValuesBase(
    const FECollection<dim, FEValuesType::space_dimension> &fe_collection,
    const QCollection<q_dim> &                              q_collection,
    const UpdateFlags                                       update_flags)
    : fe_collection(&fe_collection)
    , mapping_collection(
        &dealii::hp::StaticMappingQ1<dim, FEValuesType::space_dimension>::
          mapping_collection)
    , q_collection(q_collection)
    , fe_values_table(fe_collection.size(), 1, q_collection.size())
    , present_fe_values_index(numbers::invalid_unsigned_int,
                              numbers::invalid_unsigned_int,
                              numbers::invalid_unsigned_int)
    , update_flags(update_flags)
  {}



  template <int dim, int q_dim, class FEValuesType>
  FEValuesBase<dim, q_dim, FEValuesType>::FEValuesBase(
    const FEValuesBase<dim, q_dim, FEValuesType> &other)
    : fe_collection(other.fe_collection)
    , mapping_collection(other.mapping_collection)
    , q_collection(other.q_collection)
    , fe_values_table(fe_collection->size(),
                      mapping_collection->size(),
                      q_collection.size())
    , present_fe_values_index(other.present_fe_values_index)
    , update_flags(other.update_flags)
  {
    // We've already resized the `fe_values_table` correctly above, but right
    // now it just contains nullptrs. Create copies of the objects that
    // `other.fe_values_table` stores
    Threads::TaskGroup<> task_group;
    for (unsigned int fe_index = 0; fe_index < fe_collection->size();
         ++fe_index)
      for (unsigned int m_index = 0; m_index < mapping_collection->size();
           ++m_index)
        for (unsigned int q_index = 0; q_index < q_collection.size(); ++q_index)
          if (other.fe_values_table[fe_index][m_index][q_index].get() !=
              nullptr)
            task_group += Threads::new_task([&, fe_index, m_index, q_index]() {
              fe_values_table[fe_index][m_index][q_index] =
                std::make_unique<FEValuesType>((*mapping_collection)[m_index],
                                               (*fe_collection)[fe_index],
                                               q_collection[q_index],
                                               update_flags);
            });

    task_group.join_all();
  }



  template <int dim, int q_dim, class FEValuesType>
  FEValuesType &
  FEValuesBase<dim, q_dim, FEValuesType>::select_fe_values(
    const unsigned int fe_index,
    const unsigned int mapping_index,
    const unsigned int q_index)
  {
    AssertIndexRange(fe_index, fe_collection->size());
    AssertIndexRange(mapping_index, mapping_collection->size());
    AssertIndexRange(q_index, q_collection.size());


    // set the triple of indices
    // that we want to work with
    present_fe_values_index = TableIndices<3>(fe_index, mapping_index, q_index);

    // first check whether we
    // already have an object for
    // this particular combination
    // of indices
    if (fe_values_table(present_fe_values_index).get() == nullptr)
      fe_values_table(present_fe_values_index) =
        std::make_unique<FEValuesType>((*mapping_collection)[mapping_index],
                                       (*fe_collection)[fe_index],
                                       q_collection[q_index],
                                       update_flags);

    // now there definitely is one!
    return *fe_values_table(present_fe_values_index);
  }



  template <int dim, int q_dim, class FEValuesType>
  void
  FEValuesBase<dim, q_dim, FEValuesType>::precalculate_fe_values(
    const std::vector<unsigned int> &fe_indices,
    const std::vector<unsigned int> &mapping_indices,
    const std::vector<unsigned int> &q_indices)
  {
    AssertDimension(fe_indices.size(), mapping_indices.size());
    AssertDimension(fe_indices.size(), q_indices.size());

    Threads::TaskGroup<> task_group;
    for (unsigned int i = 0; i < fe_indices.size(); ++i)
      {
        const unsigned int fe_index      = fe_indices[i],
                           mapping_index = mapping_indices[i],
                           q_index       = q_indices[i];

        AssertIndexRange(fe_index, fe_collection->size());
        AssertIndexRange(mapping_index, mapping_collection->size());
        AssertIndexRange(q_index, q_collection.size());

        task_group +=
          Threads::new_task([&, fe_index, mapping_index, q_index]() {
            fe_values_table(TableIndices<3>(fe_index, mapping_index, q_index)) =
              std::make_unique<FEValuesType>(
                (*mapping_collection)[mapping_index],
                (*fe_collection)[fe_index],
                q_collection[q_index],
                update_flags);
          });
      }

    task_group.join_all();
  }



  template <int dim, int q_dim, class FEValuesType>
  void
  FEValuesBase<dim, q_dim, FEValuesType>::precalculate_fe_values()
  {
    const unsigned int        size = fe_collection->size();
    std::vector<unsigned int> indices(size);
    std::iota(indices.begin(), indices.end(), 0);

    precalculate_fe_values(/*fe_indices=*/indices,
                           /*mapping_indices=*/
                           (mapping_collection->size() > 1) ?
                             indices :
                             std::vector<unsigned int>(size, 0),
                           /*q_indices=*/
                           (q_collection.size() > 1) ?
                             indices :
                             std::vector<unsigned int>(size, 0));
  }
} // namespace hp


namespace hp
{
  // -------------------------- FEValues -------------------------


  template <int dim, int spacedim>
  FEValues<dim, spacedim>::FEValues(
    const MappingCollection<dim, spacedim> &mapping,
    const FECollection<dim, spacedim> &     fe_collection,
    const QCollection<dim> &                q_collection,
    const UpdateFlags                       update_flags)
    : hp::FEValuesBase<dim, dim, dealii::FEValues<dim, spacedim>>(mapping,
                                                                  fe_collection,
                                                                  q_collection,
                                                                  update_flags)
  {}


  template <int dim, int spacedim>
  FEValues<dim, spacedim>::FEValues(
    const FECollection<dim, spacedim> &fe_collection,
    const QCollection<dim> &           q_collection,
    const UpdateFlags                  update_flags)
    : hp::FEValuesBase<dim, dim, dealii::FEValues<dim, spacedim>>(fe_collection,
                                                                  q_collection,
                                                                  update_flags)
  {}


  template <int dim, int spacedim>
  template <typename DoFHandlerType, bool lda>
  void
  FEValues<dim, spacedim>::reinit(
    const TriaIterator<DoFCellAccessor<DoFHandlerType, lda>> cell,
    const unsigned int                                       q_index,
    const unsigned int                                       mapping_index,
    const unsigned int                                       fe_index)
  {
    // determine which indices we
    // should actually use
    unsigned int real_q_index = q_index, real_mapping_index = mapping_index,
                 real_fe_index = fe_index;

    if (real_q_index == numbers::invalid_unsigned_int)
      {
        if (this->q_collection.size() > 1)
          real_q_index = cell->active_fe_index();
        else
          real_q_index = 0;
      }

    if (real_mapping_index == numbers::invalid_unsigned_int)
      {
        if (this->mapping_collection->size() > 1)
          real_mapping_index = cell->active_fe_index();
        else
          real_mapping_index = 0;
      }

    if (real_fe_index == numbers::invalid_unsigned_int)
      real_fe_index = cell->active_fe_index();

    // some checks
    AssertIndexRange(real_q_index, this->q_collection.size());
    AssertIndexRange(real_mapping_index, this->mapping_collection->size());
    AssertIndexRange(real_fe_index, this->fe_collection->size());

    // now finally actually get the
    // corresponding object and
    // initialize it
    this->select_fe_values(real_fe_index, real_mapping_index, real_q_index)
      .reinit(cell);
  }



  template <int dim, int spacedim>
  void
  FEValues<dim, spacedim>::reinit(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          q_index,
    const unsigned int                                          mapping_index,
    const unsigned int                                          fe_index)
  {
    // determine which indices we
    // should actually use
    unsigned int real_q_index = q_index, real_mapping_index = mapping_index,
                 real_fe_index = fe_index;

    if (real_q_index == numbers::invalid_unsigned_int)
      real_q_index = 0;

    if (real_mapping_index == numbers::invalid_unsigned_int)
      real_mapping_index = 0;

    if (real_fe_index == numbers::invalid_unsigned_int)
      real_fe_index = 0;

    // some checks
    AssertIndexRange(real_q_index, this->q_collection.size());
    AssertIndexRange(real_mapping_index, this->mapping_collection->size());
    AssertIndexRange(real_fe_index, this->fe_collection->size());

    // now finally actually get the
    // corresponding object and
    // initialize it
    this->select_fe_values(real_fe_index, real_mapping_index, real_q_index)
      .reinit(cell);
  }


  // -------------------------- FEFaceValues -------------------------


  template <int dim, int spacedim>
  FEFaceValues<dim, spacedim>::FEFaceValues(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const hp::FECollection<dim, spacedim> &     fe_collection,
    const hp::QCollection<dim - 1> &            q_collection,
    const UpdateFlags                           update_flags)
    : hp::FEValuesBase<dim, dim - 1, dealii::FEFaceValues<dim, spacedim>>(
        mapping,
        fe_collection,
        q_collection,
        update_flags)
  {}


  template <int dim, int spacedim>
  FEFaceValues<dim, spacedim>::FEFaceValues(
    const hp::FECollection<dim, spacedim> &fe_collection,
    const hp::QCollection<dim - 1> &       q_collection,
    const UpdateFlags                      update_flags)
    : hp::FEValuesBase<dim, dim - 1, dealii::FEFaceValues<dim, spacedim>>(
        fe_collection,
        q_collection,
        update_flags)
  {}


  template <int dim, int spacedim>
  template <typename DoFHandlerType, bool lda>
  void
  FEFaceValues<dim, spacedim>::reinit(
    const TriaIterator<DoFCellAccessor<DoFHandlerType, lda>> cell,
    const unsigned int                                       face_no,
    const unsigned int                                       q_index,
    const unsigned int                                       mapping_index,
    const unsigned int                                       fe_index)
  {
    // determine which indices we
    // should actually use
    unsigned int real_q_index = q_index, real_mapping_index = mapping_index,
                 real_fe_index = fe_index;

    if (real_q_index == numbers::invalid_unsigned_int)
      {
        if (this->q_collection.size() > 1)
          real_q_index = cell->active_fe_index();
        else
          real_q_index = 0;
      }

    if (real_mapping_index == numbers::invalid_unsigned_int)
      {
        if (this->mapping_collection->size() > 1)
          real_mapping_index = cell->active_fe_index();
        else
          real_mapping_index = 0;
      }

    if (real_fe_index == numbers::invalid_unsigned_int)
      real_fe_index = cell->active_fe_index();

    // some checks
    AssertIndexRange(real_q_index, this->q_collection.size());
    AssertIndexRange(real_mapping_index, this->mapping_collection->size());
    AssertIndexRange(real_fe_index, this->fe_collection->size());

    // now finally actually get the
    // corresponding object and
    // initialize it
    this->select_fe_values(real_fe_index, real_mapping_index, real_q_index)
      .reinit(cell, face_no);
  }



  template <int dim, int spacedim>
  void
  FEFaceValues<dim, spacedim>::reinit(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          q_index,
    const unsigned int                                          mapping_index,
    const unsigned int                                          fe_index)
  {
    // determine which indices we
    // should actually use
    unsigned int real_q_index = q_index, real_mapping_index = mapping_index,
                 real_fe_index = fe_index;

    if (real_q_index == numbers::invalid_unsigned_int)
      real_q_index = 0;

    if (real_mapping_index == numbers::invalid_unsigned_int)
      real_mapping_index = 0;

    if (real_fe_index == numbers::invalid_unsigned_int)
      real_fe_index = 0;

    // some checks
    AssertIndexRange(real_q_index, this->q_collection.size());
    AssertIndexRange(real_mapping_index, this->mapping_collection->size());
    AssertIndexRange(real_fe_index, this->fe_collection->size());

    // now finally actually get the
    // corresponding object and
    // initialize it
    this->select_fe_values(real_fe_index, real_mapping_index, real_q_index)
      .reinit(cell, face_no);
  }


  // -------------------------- FESubfaceValues -------------------------


  template <int dim, int spacedim>
  FESubfaceValues<dim, spacedim>::FESubfaceValues(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const hp::FECollection<dim, spacedim> &     fe_collection,
    const hp::QCollection<dim - 1> &            q_collection,
    const UpdateFlags                           update_flags)
    : hp::FEValuesBase<dim, dim - 1, dealii::FESubfaceValues<dim, spacedim>>(
        mapping,
        fe_collection,
        q_collection,
        update_flags)
  {}


  template <int dim, int spacedim>
  FESubfaceValues<dim, spacedim>::FESubfaceValues(
    const hp::FECollection<dim, spacedim> &fe_collection,
    const hp::QCollection<dim - 1> &       q_collection,
    const UpdateFlags                      update_flags)
    : hp::FEValuesBase<dim, dim - 1, dealii::FESubfaceValues<dim, spacedim>>(
        fe_collection,
        q_collection,
        update_flags)
  {}


  template <int dim, int spacedim>
  template <typename DoFHandlerType, bool lda>
  void
  FESubfaceValues<dim, spacedim>::reinit(
    const TriaIterator<DoFCellAccessor<DoFHandlerType, lda>> cell,
    const unsigned int                                       face_no,
    const unsigned int                                       subface_no,
    const unsigned int                                       q_index,
    const unsigned int                                       mapping_index,
    const unsigned int                                       fe_index)
  {
    // determine which indices we
    // should actually use
    unsigned int real_q_index = q_index, real_mapping_index = mapping_index,
                 real_fe_index = fe_index;

    if (real_q_index == numbers::invalid_unsigned_int)
      {
        if (this->q_collection.size() > 1)
          real_q_index = cell->active_fe_index();
        else
          real_q_index = 0;
      }

    if (real_mapping_index == numbers::invalid_unsigned_int)
      {
        if (this->mapping_collection->size() > 1)
          real_mapping_index = cell->active_fe_index();
        else
          real_mapping_index = 0;
      }

    if (real_fe_index == numbers::invalid_unsigned_int)
      real_fe_index = cell->active_fe_index();

    // some checks
    AssertIndexRange(real_q_index, this->q_collection.size());
    AssertIndexRange(real_mapping_index, this->mapping_collection->size());
    AssertIndexRange(real_fe_index, this->fe_collection->size());

    // now finally actually get the
    // corresponding object and
    // initialize it
    this->select_fe_values(real_fe_index, real_mapping_index, real_q_index)
      .reinit(cell, face_no, subface_no);
  }



  template <int dim, int spacedim>
  void
  FESubfaceValues<dim, spacedim>::reinit(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          subface_no,
    const unsigned int                                          q_index,
    const unsigned int                                          mapping_index,
    const unsigned int                                          fe_index)
  {
    // determine which indices we
    // should actually use
    unsigned int real_q_index = q_index, real_mapping_index = mapping_index,
                 real_fe_index = fe_index;

    if (real_q_index == numbers::invalid_unsigned_int)
      real_q_index = 0;

    if (real_mapping_index == numbers::invalid_unsigned_int)
      real_mapping_index = 0;

    if (real_fe_index == numbers::invalid_unsigned_int)
      real_fe_index = 0;

    // some checks
    AssertIndexRange(real_q_index, this->q_collection.size());
    AssertIndexRange(real_mapping_index, this->mapping_collection->size());
    AssertIndexRange(real_fe_index, this->fe_collection->size());

    // now finally actually get the
    // corresponding object and
    // initialize it
    this->select_fe_values(real_fe_index, real_mapping_index, real_q_index)
      .reinit(cell, face_no, subface_no);
  }
} // namespace hp


// explicit instantiations
#include "fe_values.inst"


DEAL_II_NAMESPACE_CLOSE

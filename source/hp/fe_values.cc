// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2018 by the deal.II authors
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

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/hp/fe_values.h>

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
        std::make_shared<FEValuesType>((*mapping_collection)[mapping_index],
                                       (*fe_collection)[fe_index],
                                       q_collection[q_index],
                                       update_flags);

    // now there definitely is one!
    return *fe_values_table(present_fe_values_index);
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

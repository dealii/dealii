// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2021 by the deal.II authors
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


#ifndef dealii_non_matching_mapping_info_h
#define dealii_non_matching_mapping_info_h


#include <deal.II/base/config.h>

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_cartesian.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_related_data.h>

#include <memory>


DEAL_II_NAMESPACE_OPEN

namespace NonMatching
{
  /**
   * This class provides the mapping information computation and mapping data
   * storage to be used together with FEPointEvaluation.
   */
  template <int dim, int spacedim = dim>
  class MappingInfo : public Subscriptor
  {
  public:
    using MappingData =
      dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                   spacedim>;

    /**
     * Constructor.
     *
     * @param mapping The Mapping class describing the geometry of a cell.
     *
     * @param update_flags Specify the quantities to be computed by the mapping
     * during the call of reinit(). These update flags are also handed to a
     * FEEvaluation object if you construct it with this MappingInfo object.
     */
    MappingInfo(const Mapping<dim> &mapping, const UpdateFlags update_flags);

    /**
     * Compute the mapping information for the incoming cell and unit
     * points. This overload is needed to resolve ambiguity.
     */
    void
    reinit(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
           const std::vector<Point<dim>> &unit_points);

    /**
     * Compute the mapping information for the incoming cell and unit
     * points.
     */
    void
    reinit(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
           const ArrayView<const Point<dim>> &unit_points);

    /**
     * Compute the mapping information for the given cell and
     * quadrature formula. As opposed to the other `reinit` function, this
     * method allows to access a `JxW` factor at the points.
     */
    void
    reinit(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
           const Quadrature<dim> &quadrature);

    /**
     * Compute the mapping information for the incoming iterable container of
     * cell iterators and corresponding vector of unit points.
     *
     * It is possible to give an IteratorRange<FilteredIterator> to this
     * function and together with @p n_unfiltered_cells specified this object
     * can compress its storage while giving you access to the underlying data
     * with the "unfiltered" index. The default argument for
     * @p n_unfiltered_cells disables this built-in compression.
     */
    template <typename ContainerType>
    void
    reinit_cells(
      const ContainerType &                       cell_iterator_range,
      const std::vector<std::vector<Point<dim>>> &unit_points_vector,
      const unsigned int n_unfiltered_cells = numbers::invalid_unsigned_int);

    /**
     * Compute the mapping information for the incoming iterable container of
     * cell iterators and corresponding vector of quadratures. As opposed to the
     * other `reinit` function, this method allows to access a `JxW` factor at
     * the points.
     */
    template <typename ContainerType>
    void
    reinit_cells(
      const ContainerType &               cell_iterator_range,
      const std::vector<Quadrature<dim>> &quadrature_vector,
      const unsigned int n_unfiltered_cells = numbers::invalid_unsigned_int);

    /**
     * Compute the mapping information for the incoming vector of cells and
     * corresponding vector of ImmersedSurfaceQuadrature.
     */
    template <typename Iterator>
    void
    reinit_surface(
      const IteratorRange<Iterator> &                    cell_iterator_range,
      const std::vector<ImmersedSurfaceQuadrature<dim>> &quadrature_vector,
      const unsigned int n_unfiltered_cells = numbers::invalid_unsigned_int);

    /**
     * Compute the mapping information for all faces of the incoming vector
     * of cells and corresponding vector of quadratures.
     */
    template <typename Iterator>
    void
    reinit_faces(
      const IteratorRange<Iterator> &                      cell_iterator_range,
      const std::vector<std::vector<Quadrature<dim - 1>>> &quadrature_vector,
      const unsigned int n_unfiltered_cells = numbers::invalid_unsigned_int);

    /**
     * Getter function for current unit points.
     *
     * @p cell_index and @p face_number are the indices
     * into the compressed, CRS like data storage of unit points.
     *
     * If you have initialized this object with reinit_cells() you can access
     * the stored unit points of the cell with the respective @p cell_index
     * (and the default argument for @p face_number).
     *
     * If you have initialized this object with reinit_faces() you can access
     * the stored unit points of the faces on the cell with the respective @p cell_index
     * and the respective local @p face_number.
     *
     * If you have initialized this object with reinit() you can access the
     * stored unit points of a single cell with the default arguments.
     *
     * The correct state of this object is checked in this call (in debug mode).
     */
    const ArrayView<const Point<dim>>
    get_unit_points(
      const unsigned int cell_index  = numbers::invalid_unsigned_int,
      const unsigned int face_number = numbers::invalid_unsigned_int) const;

    /**
     * Getter function for computed mapping data. This function accesses
     * internal data and is therefore not a stable interface.
     */
    const MappingData &
    get_mapping_data(
      const unsigned int cell_index  = numbers::invalid_unsigned_int,
      const unsigned int face_number = numbers::invalid_unsigned_int) const;

    /**
     * Getter function for underlying mapping.
     */
    const Mapping<dim, spacedim> &
    get_mapping() const;

    /**
     * Getter function for the update flags.
     */
    UpdateFlags
    get_update_flags() const;

    /**
     * Connects to is_reinitialized().
     */
    boost::signals2::connection
    connect_is_reinitialized(const std::function<void()> &set_is_reinitialized);

  private:
    /**
     * Enum class for reinitialized states.
     */
    enum class State
    {
      invalid,
      single_cell,
      cell_vector,
      faces_on_cells_in_vector
    };

    /**
     * Enum class that stores the currently initialized state
     * upon the last call of reinit().
     */
    State state;

    /**
     * Compute the mapping related data for the given @p mapping,
     * @p cell and @p unit_points that is required by the FEPointEvaluation
     * class.
     */
    void
    compute_mapping_data_for_quadrature(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell,
      CellSimilarity::Similarity &cell_similarity,
      const Quadrature<dim> &     quadrature,
      MappingData &               mapping_data);

    /**
     * Compute the mapping related data for the given @p mapping,
     * @p cell and @p quadrature that is required by the FEPointEvaluation
     * class.
     */
    void
    compute_mapping_data_for_immersed_surface_quadrature(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell,
      const ImmersedSurfaceQuadrature<dim> &                      quadrature,
      MappingData &                                               mapping_data);

    /**
     * Compute the mapping related data for the given @p mapping, @p cell,
     * @p face_no and @p quadrature that is required by the FEPointEvaluation
     * class.
     */
    void
    compute_mapping_data_for_face_quadrature(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell,
      const unsigned int                                          face_no,
      const Quadrature<dim - 1> &                                 quadrature,
      MappingData &                                               mapping_data);

    /**
     * The reference points specified at reinit().
     */
    std::vector<Point<dim>> unit_points;

    /**
     * Offset to point to the first unit point of a cell/face
     */
    std::vector<unsigned int> unit_points_index;

    /**
     * A pointer to the internal data of the underlying mapping.
     */
    std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
      internal_mapping_data;

    /**
     * A pointer to the underlying mapping.
     */
    const SmartPointer<const Mapping<dim, spacedim>> mapping;

    /**
     * The desired update flags for the evaluation.
     */
    const UpdateFlags update_flags;

    /**
     * The update flags for the desired mapping information.
     */
    UpdateFlags update_flags_mapping;

    /**
     * The internal data container for mapping information. The implementation
     * is subject to future changes.
     */
    std::vector<MappingData> mapping_data;

    /**
     * Offset to point to the first element of a cell in internal data
     * containers.
     */
    std::vector<unsigned int> cell_index_offset;

    /**
     * A vector that converts the cell index to a compressed cell index for e.g.
     * a filtered IteratorRange.
     */
    std::vector<unsigned int> cell_index_to_compressed_cell_index;

    /**
     * A bool that determines weather cell index compression should be done.
     */
    bool do_cell_index_compression;

    /**
     * This signal is triggered right after this object is reinitialized, to let
     * dependent objects know that they need to reinitialize as well.
     */
    boost::signals2::signal<void()> is_reinitialized;
  };

  // ----------------------- template functions ----------------------


  template <int dim, int spacedim>
  MappingInfo<dim, spacedim>::MappingInfo(const Mapping<dim> &mapping,
                                          const UpdateFlags   update_flags)
    : mapping(&mapping)
    , update_flags(update_flags)
  {
    update_flags_mapping = update_default;
    // translate update flags
    if (update_flags & update_jacobians || update_flags & update_JxW_values)
      update_flags_mapping |= update_jacobians;
    if (update_flags & update_JxW_values)
      update_flags_mapping |= update_JxW_values;
    if (update_flags & update_normal_vectors)
      update_flags_mapping |= update_normal_vectors;
    if (update_flags & update_gradients ||
        update_flags & update_inverse_jacobians)
      update_flags_mapping |= update_inverse_jacobians;

    // always save quadrature points for now
    update_flags_mapping |= update_quadrature_points;

    // construct internal_mapping_data for MappingQ to be able to reuse it in
    // reinit() calls to avoid memory allocations
    if (const MappingQ<dim, spacedim> *mapping_q =
          dynamic_cast<const MappingQ<dim, spacedim> *>(&mapping))
      {
        internal_mapping_data =
          std::make_unique<typename MappingQ<dim, spacedim>::InternalData>(
            mapping_q->get_degree());
      }
  }



  template <int dim, int spacedim>
  void
  MappingInfo<dim, spacedim>::reinit(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const std::vector<Point<dim>> &                             unit_points_in)
  {
    reinit(cell, Quadrature<dim>(unit_points_in));
  }



  template <int dim, int spacedim>
  void
  MappingInfo<dim, spacedim>::reinit(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const ArrayView<const Point<dim>> &                         unit_points_in)
  {
    reinit(cell,
           std::vector<Point<dim>>(unit_points_in.begin(),
                                   unit_points_in.end()));
  }



  template <int dim, int spacedim>
  void
  MappingInfo<dim, spacedim>::reinit(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Quadrature<dim> &                                     quadrature)
  {
    unit_points = quadrature.get_points();

    mapping_data.resize(1);
    CellSimilarity::Similarity cell_similarity =
      CellSimilarity::Similarity::none;
    compute_mapping_data_for_quadrature(cell,
                                        cell_similarity,
                                        quadrature,
                                        mapping_data[0]);

    state = State::single_cell;
    is_reinitialized();
  }



  template <int dim, int spacedim>
  template <typename ContainerType>
  void
  MappingInfo<dim, spacedim>::reinit_cells(
    const ContainerType &                       cell_iterator_range,
    const std::vector<std::vector<Point<dim>>> &unit_points_vector,
    const unsigned int                          n_unfiltered_cells)
  {
    const unsigned int n_cells = unit_points_vector.size();
    AssertDimension(n_cells,
                    std::distance(cell_iterator_range.begin(),
                                  cell_iterator_range.end()));

    std::vector<Quadrature<dim>> quadrature_vector(n_cells);
    for (unsigned int cell_index = 0; cell_index < n_cells; ++cell_index)
      quadrature_vector[cell_index] =
        Quadrature<dim>(unit_points_vector[cell_index]);

    reinit_cells(cell_iterator_range, quadrature_vector, n_unfiltered_cells);
  }



  template <int dim, int spacedim>
  template <typename ContainerType>
  void
  MappingInfo<dim, spacedim>::reinit_cells(
    const ContainerType &               cell_iterator_range,
    const std::vector<Quadrature<dim>> &quadrature_vector,
    const unsigned int                  n_unfiltered_cells)
  {
    do_cell_index_compression =
      n_unfiltered_cells != numbers::invalid_unsigned_int;

    const unsigned int n_cells = quadrature_vector.size();
    AssertDimension(n_cells,
                    std::distance(cell_iterator_range.begin(),
                                  cell_iterator_range.end()));

    // fill unit points index offset vector
    unit_points_index.reserve(n_cells + 1);
    unit_points_index.push_back(0);
    for (const auto &quadrature : quadrature_vector)
      unit_points_index.push_back(unit_points_index.back() + quadrature.size());

    const unsigned int n_unit_points = unit_points_index.back();

    unit_points.resize(n_unit_points);
    mapping_data.resize(n_cells);

    if (do_cell_index_compression)
      cell_index_to_compressed_cell_index.resize(n_unfiltered_cells,
                                                 numbers::invalid_unsigned_int);
    CellSimilarity::Similarity cell_similarity =
      CellSimilarity::Similarity::none;
    unsigned int cell_index = 0;
    for (const auto &cell : cell_iterator_range)
      {
        auto it = unit_points.begin() + unit_points_index[cell_index];
        for (const auto &unit_point :
             quadrature_vector[cell_index].get_points())
          {
            *it = unit_point;
            ++it;
          }

        compute_mapping_data_for_quadrature(cell,
                                            cell_similarity,
                                            quadrature_vector[cell_index],
                                            mapping_data[cell_index]);

        if (do_cell_index_compression)
          cell_index_to_compressed_cell_index[cell->active_cell_index()] =
            cell_index;

        ++cell_index;
      }

    state = State::cell_vector;
    is_reinitialized();
  }



  template <int dim, int spacedim>
  template <typename Iterator>
  void
  MappingInfo<dim, spacedim>::reinit_surface(
    const IteratorRange<Iterator> &                    cell_iterator_range,
    const std::vector<ImmersedSurfaceQuadrature<dim>> &quadrature_vector,
    const unsigned int                                 n_unfiltered_cells)
  {
    do_cell_index_compression =
      n_unfiltered_cells != numbers::invalid_unsigned_int;

    if (update_flags_mapping & (update_JxW_values | update_normal_vectors))
      update_flags_mapping |= update_covariant_transformation;

    const unsigned int n_cells = quadrature_vector.size();
    AssertDimension(n_cells,
                    std::distance(cell_iterator_range.begin(),
                                  cell_iterator_range.end()));

    // fill unit points index offset vector
    unit_points_index.reserve(n_cells + 1);
    unit_points_index.push_back(0);
    for (const auto &quadrature : quadrature_vector)
      unit_points_index.push_back(unit_points_index.back() +
                                  quadrature.get_points().size());

    const unsigned int n_unit_points = unit_points_index.back();

    unit_points.resize(n_unit_points);
    mapping_data.resize(n_cells);

    if (do_cell_index_compression)
      cell_index_to_compressed_cell_index.resize(n_unfiltered_cells,
                                                 numbers::invalid_unsigned_int);
    unsigned int cell_index = 0;
    for (const auto &cell : cell_iterator_range)
      {
        const auto &quadrature = quadrature_vector[cell_index];

        auto it = unit_points.begin() + unit_points_index[cell_index];
        for (const auto &unit_point : quadrature.get_points())
          {
            *it = unit_point;
            ++it;
          }

        compute_mapping_data_for_immersed_surface_quadrature(
          cell, quadrature, mapping_data[cell_index]);

        if (do_cell_index_compression)
          cell_index_to_compressed_cell_index[cell->active_cell_index()] =
            cell_index;

        ++cell_index;
      }

    state = State::cell_vector;
    is_reinitialized();
  }



  template <int dim, int spacedim>
  template <typename Iterator>
  void
  MappingInfo<dim, spacedim>::reinit_faces(
    const IteratorRange<Iterator> &                      cell_iterator_range,
    const std::vector<std::vector<Quadrature<dim - 1>>> &quadrature_vector,
    const unsigned int                                   n_unfiltered_cells)
  {
    do_cell_index_compression =
      n_unfiltered_cells != numbers::invalid_unsigned_int;

    const unsigned int n_cells = quadrature_vector.size();
    AssertDimension(n_cells,
                    std::distance(cell_iterator_range.begin(),
                                  cell_iterator_range.end()));

    // fill cell index offset vector
    cell_index_offset.resize(n_cells);
    unsigned int n_faces    = 0;
    unsigned int cell_index = 0;
    for (const auto &cell : cell_iterator_range)
      {
        cell_index_offset[cell_index] = n_faces;
        n_faces += cell->n_faces();
        ++cell_index;
      }

    // fill unit points index offset vector
    unit_points_index.resize(n_faces + 1);
    cell_index                 = 0;
    unsigned int n_unit_points = 0;
    for (const auto &cell : cell_iterator_range)
      {
        for (const auto &f : cell->face_indices())
          {
            const unsigned int current_face_index =
              cell_index_offset[cell_index] + f;

            unit_points_index[current_face_index] = n_unit_points;
            n_unit_points +=
              quadrature_vector[cell_index][f].get_points().size();
          }

        ++cell_index;
      }
    unit_points_index[n_faces] = n_unit_points;

    // compress indices
    if (do_cell_index_compression)
      cell_index_to_compressed_cell_index.resize(n_unfiltered_cells,
                                                 numbers::invalid_unsigned_int);

    // fill unit points and mapping data for every face of all cells
    unit_points.resize(n_unit_points);
    mapping_data.resize(n_faces);
    cell_index = 0;
    QProjector<dim> q_projector;
    for (const auto &cell : cell_iterator_range)
      {
        const auto &quadratures_on_faces = quadrature_vector[cell_index];

        Assert(quadratures_on_faces.size() == cell->n_faces(),
               ExcDimensionMismatch(quadratures_on_faces.size(),
                                    cell->n_faces()));

        for (const auto &f : cell->face_indices())
          {
            const auto &quadrature_on_face = quadratures_on_faces[f];

            const auto quadrature_on_cell =
              q_projector.project_to_face(cell->reference_cell(),
                                          quadrature_on_face,
                                          f);

            const auto &unit_points_on_cell = quadrature_on_cell.get_points();

            const unsigned int current_face_index =
              cell_index_offset[cell_index] + f;

            auto it =
              unit_points.begin() + unit_points_index[current_face_index];
            for (const auto &unit_point : unit_points_on_cell)
              {
                *it = unit_point;
                ++it;
              }

            compute_mapping_data_for_face_quadrature(
              cell, f, quadrature_on_face, mapping_data[current_face_index]);
          }
        if (do_cell_index_compression)
          cell_index_to_compressed_cell_index[cell->active_cell_index()] =
            cell_index;

        ++cell_index;
      }

    state = State::faces_on_cells_in_vector;
    is_reinitialized();
  }



  template <int dim, int spacedim>
  inline const ArrayView<const Point<dim>>
  MappingInfo<dim, spacedim>::get_unit_points(
    const unsigned int cell_index,
    const unsigned int face_number) const
  {
    if (cell_index == numbers::invalid_unsigned_int &&
        face_number == numbers::invalid_unsigned_int)
      {
        Assert(state == State::single_cell,
               ExcMessage(
                 "This mapping info is not reinitialized for a single cell!"));
        return unit_points;
      }
    else if (face_number == numbers::invalid_unsigned_int)
      {
        Assert(state == State::cell_vector,
               ExcMessage(
                 "This mapping info is not reinitialized for a cell vector!"));
        if (do_cell_index_compression)
          {
            Assert(cell_index_to_compressed_cell_index[cell_index] !=
                     numbers::invalid_unsigned_int,
                   ExcMessage("Mapping info object was not initialized for this"
                              " active cell index!"));
            const auto it_begin =
              unit_points.begin() +
              unit_points_index
                [cell_index_to_compressed_cell_index[cell_index]];
            const auto it_end =
              unit_points.begin() +
              unit_points_index
                [cell_index_to_compressed_cell_index[cell_index] + 1];
            return make_array_view(it_begin, it_end);
          }
        else
          {
            const auto it_begin =
              unit_points.begin() + unit_points_index[cell_index];
            const auto it_end =
              unit_points.begin() + unit_points_index[cell_index + 1];
            return make_array_view(it_begin, it_end);
          }
      }
    else if (cell_index != numbers::invalid_unsigned_int)
      {
        Assert(state == State::faces_on_cells_in_vector,
               ExcMessage("This mapping info is not reinitialized for faces"
                          " on cells in a vector!"));
        if (do_cell_index_compression)
          {
            Assert(
              cell_index_to_compressed_cell_index[cell_index] !=
                numbers::invalid_unsigned_int,
              ExcMessage(
                "Mapping info object was not initialized for this active cell index"
                " and corresponding face numbers!"));
            const unsigned int current_face_index =
              cell_index_offset
                [cell_index_to_compressed_cell_index[cell_index]] +
              face_number;
            const auto it_begin =
              unit_points.begin() + unit_points_index[current_face_index];
            const auto it_end =
              unit_points.begin() + unit_points_index[current_face_index + 1];
            return make_array_view(it_begin, it_end);
          }
        else
          {
            const unsigned int current_face_index =
              cell_index_offset[cell_index] + face_number;
            const auto it_begin =
              unit_points.begin() + unit_points_index[current_face_index];
            const auto it_end =
              unit_points.begin() + unit_points_index[current_face_index + 1];
            return make_array_view(it_begin, it_end);
          }
      }
    else
      AssertThrow(
        false,
        ExcMessage(
          "cell_index has to be specified if face_number is specified!"));
  }



  template <int dim, int spacedim>
  inline const typename MappingInfo<dim, spacedim>::MappingData &
  MappingInfo<dim, spacedim>::get_mapping_data(
    const unsigned int cell_index,
    const unsigned int face_number) const
  {
    if (cell_index == numbers::invalid_unsigned_int &&
        face_number == numbers::invalid_unsigned_int)
      {
        Assert(state == State::single_cell,
               ExcMessage(
                 "This mapping info is not reinitialized for a single cell!"));
        return mapping_data[0];
      }
    else if (face_number == numbers::invalid_unsigned_int)
      {
        Assert(state == State::cell_vector,
               ExcMessage(
                 "This mapping info is not reinitialized for a cell vector!"));
        if (do_cell_index_compression)
          {
            Assert(cell_index_to_compressed_cell_index[cell_index] !=
                     numbers::invalid_unsigned_int,
                   ExcMessage("Mapping info object was not initialized for this"
                              " active cell index!"));
            return mapping_data
              [cell_index_to_compressed_cell_index[cell_index]];
          }
        else
          return mapping_data[cell_index];
      }
    else if (cell_index != numbers::invalid_unsigned_int)
      {
        Assert(state == State::faces_on_cells_in_vector,
               ExcMessage("This mapping info is not reinitialized for faces"
                          " on cells in a vector!"));
        if (do_cell_index_compression)
          {
            Assert(
              cell_index_to_compressed_cell_index[cell_index] !=
                numbers::invalid_unsigned_int,
              ExcMessage(
                "Mapping info object was not initialized for this active cell index"
                " and corresponding face numbers!"));
            return mapping_data
              [cell_index_offset
                 [cell_index_to_compressed_cell_index[cell_index]] +
               face_number];
          }
        else
          return mapping_data[cell_index_offset[cell_index] + face_number];
      }
    else
      AssertThrow(
        false,
        ExcMessage(
          "cell_index has to be specified if face_number is specified!"));
  }



  template <int dim, int spacedim>
  const Mapping<dim, spacedim> &
  MappingInfo<dim, spacedim>::get_mapping() const
  {
    return *mapping;
  }



  template <int dim, int spacedim>
  UpdateFlags
  MappingInfo<dim, spacedim>::get_update_flags() const
  {
    return update_flags;
  }



  template <int dim, int spacedim>
  boost::signals2::connection
  MappingInfo<dim, spacedim>::connect_is_reinitialized(
    const std::function<void()> &set_is_reinitialized)
  {
    return is_reinitialized.connect(set_is_reinitialized);
  }



  template <int dim, int spacedim>
  void
  MappingInfo<dim, spacedim>::compute_mapping_data_for_quadrature(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    CellSimilarity::Similarity &                                cell_similarity,
    const Quadrature<dim> &                                     quadrature,
    MappingData &                                               mapping_data)
  {
    update_flags_mapping |=
      mapping->requires_update_flags(update_flags_mapping);

    mapping_data.initialize(quadrature.size(), update_flags_mapping);

    // reuse internal_mapping_data for MappingQ to avoid memory allocations
    if (const MappingQ<dim, spacedim> *mapping_q =
          dynamic_cast<const MappingQ<dim, spacedim> *>(&(*mapping)))
      {
        (void)mapping_q;
        auto &data =
          dynamic_cast<typename MappingQ<dim, spacedim>::InternalData &>(
            *internal_mapping_data);
        data.initialize(update_flags_mapping, quadrature, quadrature.size());
      }
    else
      {
        internal_mapping_data =
          mapping->get_data(update_flags_mapping, quadrature);
      }

    cell_similarity = mapping->fill_fe_values(
      cell, cell_similarity, quadrature, *internal_mapping_data, mapping_data);
  }



  template <int dim, int spacedim>
  void
  MappingInfo<dim, spacedim>::
    compute_mapping_data_for_immersed_surface_quadrature(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell,
      const ImmersedSurfaceQuadrature<dim> &                      quadrature,
      MappingData &                                               mapping_data)
  {
    update_flags_mapping |=
      mapping->requires_update_flags(update_flags_mapping);

    mapping_data.initialize(quadrature.size(), update_flags_mapping);

    // reuse internal_mapping_data for MappingQ to avoid memory allocations
    if (const MappingQ<dim, spacedim> *mapping_q =
          dynamic_cast<const MappingQ<dim, spacedim> *>(&(*mapping)))
      {
        (void)mapping_q;
        auto &data =
          dynamic_cast<typename MappingQ<dim, spacedim>::InternalData &>(
            *internal_mapping_data);
        data.initialize(update_flags_mapping, quadrature, quadrature.size());
      }
    else
      {
        internal_mapping_data =
          mapping->get_data(update_flags_mapping, quadrature);
      }

    mapping->fill_fe_immersed_surface_values(cell,
                                             quadrature,
                                             *internal_mapping_data,
                                             mapping_data);
  }



  template <int dim, int spacedim>
  void
  MappingInfo<dim, spacedim>::compute_mapping_data_for_face_quadrature(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const Quadrature<dim - 1> &                                 quadrature,
    MappingData &                                               mapping_data)
  {
    update_flags_mapping |=
      mapping->requires_update_flags(update_flags_mapping);

    mapping_data.initialize(quadrature.size(), update_flags_mapping);

    // reuse internal_mapping_data for MappingQ to avoid memory allocations
    if (const MappingQ<dim, spacedim> *mapping_q =
          dynamic_cast<const MappingQ<dim, spacedim> *>(&(*mapping)))
      {
        auto &data =
          dynamic_cast<typename MappingQ<dim, spacedim>::InternalData &>(
            *internal_mapping_data);
        data.initialize_face(update_flags_mapping,
                             QProjector<dim>::project_to_oriented_face(
                               ReferenceCells::get_hypercube<dim>(),
                               quadrature,
                               face_no,
                               cell->face_orientation(face_no),
                               cell->face_flip(face_no),
                               cell->face_rotation(face_no)),
                             quadrature.size());

        mapping_q->fill_mapping_data_for_face_quadrature(
          cell, face_no, quadrature, *internal_mapping_data, mapping_data);
      }
    else
      {
        auto internal_mapping_data =
          mapping->get_face_data(update_flags_mapping,
                                 hp::QCollection<dim - 1>(quadrature));

        mapping->fill_fe_face_values(cell,
                                     face_no,
                                     hp::QCollection<dim - 1>(quadrature),
                                     *internal_mapping_data,
                                     mapping_data);
      }
  }
} // namespace NonMatching

DEAL_II_NAMESPACE_CLOSE

#endif

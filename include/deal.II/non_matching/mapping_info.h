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
  namespace internal
  {
    template <int dim, int spacedim = dim>
    class ComputeMappingDataHelper
    {
      using MappingData =
        dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                     spacedim>;

    public:
      static UpdateFlags
      required_update_flags(
        const SmartPointer<const Mapping<dim, spacedim>> mapping,
        const UpdateFlags &                              update_flags)
      {
        return mapping->requires_update_flags(update_flags);
      }

      static void
      compute_mapping_data_for_quadrature(
        const SmartPointer<const Mapping<dim, spacedim>> mapping,
        const UpdateFlags &                              update_flags_mapping,
        const typename Triangulation<dim, spacedim>::cell_iterator &cell,
        CellSimilarity::Similarity &cell_similarity,
        const Quadrature<dim> &     quadrature,
        std::shared_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
                     internal_mapping_data,
        MappingData &mapping_data)
      {
        mapping_data.initialize(quadrature.size(), update_flags_mapping);

        // reuse internal_mapping_data for MappingQ to avoid memory allocations
        if (const MappingQ<dim, spacedim> *mapping_q =
              dynamic_cast<const MappingQ<dim, spacedim> *>(&(*mapping)))
          {
            (void)mapping_q;
            auto &data =
              dynamic_cast<typename MappingQ<dim, spacedim>::InternalData &>(
                *internal_mapping_data);
            data.initialize(update_flags_mapping,
                            quadrature,
                            quadrature.size());
          }
        else
          {
            internal_mapping_data =
              mapping->get_data(update_flags_mapping, quadrature);
          }

        cell_similarity = mapping->fill_fe_values(cell,
                                                  cell_similarity,
                                                  quadrature,
                                                  *internal_mapping_data,
                                                  mapping_data);
      }



      static void
      compute_mapping_data_for_immersed_surface_quadrature(
        const SmartPointer<const Mapping<dim, spacedim>> mapping,
        const UpdateFlags &                              update_flags_mapping,
        const typename Triangulation<dim, spacedim>::cell_iterator &cell,
        const ImmersedSurfaceQuadrature<dim> &                      quadrature,
        std::shared_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
                     internal_mapping_data,
        MappingData &mapping_data)
      {
        mapping_data.initialize(quadrature.size(), update_flags_mapping);

        // reuse internal_mapping_data for MappingQ to avoid memory allocations
        if (dynamic_cast<const MappingQ<dim, spacedim> *>(&(*mapping)))
          {
            auto &data =
              dynamic_cast<typename MappingQ<dim, spacedim>::InternalData &>(
                *internal_mapping_data);
            data.initialize(update_flags_mapping,
                            quadrature,
                            quadrature.size());
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



      static void
      compute_mapping_data_for_face_quadrature(
        const SmartPointer<const Mapping<dim, spacedim>> mapping,
        const UpdateFlags &                              update_flags_mapping,
        const typename Triangulation<dim, spacedim>::cell_iterator &cell,
        const unsigned int                                          face_no,
        const Quadrature<dim - 1> &                                 quadrature,
        std::shared_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
                     internal_mapping_data,
        MappingData &mapping_data)
      {
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
    };
  } // namespace internal

  /**
   * This class provides the mapping information computation and mapping data
   * storage to be used together with FEPointEvaluation.
   *
   * MappingInfo is vectorized across quadrature points, which means data can be
   * provided in vectorized format.
   *
   * Two different modes are available: partially vectorized and fully
   * vectorized across quadrature points. Partially vectorized (scalar Number
   * template argument) means the computed mapping data is provided in scalar
   * format and only unit points are provided vectorized (for seamless
   * interaction with FEPointEvaluation). Fully vectorized (vectorized Number
   * template argument, e.g. VectorizedArray<double>) provides both mapping data
   * and unit points in vectorized format. The Number template parameter of
   * MappingInfo and FEPointEvaluation has to be identical.
   */
  template <int dim, int spacedim = dim, typename Number = double>
  class MappingInfo : public Subscriptor
  {
  public:
    /**
     * The VectorizedArray type the unit points are stored in inside this class.
     */
    using VectorizedArrayType = typename dealii::internal::VectorizedArrayTrait<
      Number>::vectorized_value_type;

    /**
     * Collects the options which can be used to specify the behavior during
     * reinitialization.
     */
    struct AdditionalData
    {
      /**
       * Constructor which sets the default arguments.
       */
      AdditionalData(const bool use_global_weights = false)
        : use_global_weights(use_global_weights)
      {}

      /**
       * During initialization, assume that the Quadrature object contains
       * global weights as, e.g., obtained by
       * QSimplex::compute_affine_transformation().
       */
      bool use_global_weights;
    };

    /**
     * Constructor.
     *
     * @param mapping The Mapping class describing the geometry of a cell.
     *
     * @param update_flags Specify the quantities to be computed by the mapping
     * during the call of reinit(). These update flags are also handed to a
     * FEEvaluation object if you construct it with this MappingInfo object.
     *
     * @param additional_data Additional data for the class to specify the
     * behavior during reinitialization.
     */
    MappingInfo(const Mapping<dim> & mapping,
                const UpdateFlags    update_flags,
                const AdditionalData additional_data = AdditionalData());

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
     * Return if this MappingInfo object is reinitialized for faces (by
     * reinit_faces()) or not.
     */
    bool
    is_face_state() const;

    /**
     * Getter function for unit points. The offset can be obtained with
     * compute_unit_point_index_offset().
     */
    const Point<dim, VectorizedArrayType> *
    get_unit_point(const unsigned int offset) const;

    /**
     * Getter function for unit points on faces. The offset can be obtained with
     * compute_unit_point_index_offset().
     */
    const Point<dim - 1, VectorizedArrayType> *
    get_unit_point_faces(const unsigned int offset) const;

    /**
     * Getter function for Jacobians. The offset can be obtained with
     * compute_data_index_offset().
     */
    const DerivativeForm<1, dim, spacedim, Number> *
    get_jacobian(const unsigned int offset) const;

    /**
     * Getter function for inverse Jacobians. The offset can be obtained with
     * compute_data_index_offset().
     */
    const DerivativeForm<1, spacedim, dim, Number> *
    get_inverse_jacobian(const unsigned int offset) const;

    /**
     * Getter function for normal vectors. The offset can be obtained with
     * compute_data_index_offset().
     */
    const Tensor<1, spacedim, Number> *
    get_normal_vector(const unsigned int offset) const;

    /**
     * Getter function for Jacobian times quadrature weight (JxW). The offset
     * can be obtained with compute_data_index_offset().
     */
    const Number *
    get_JxW(const unsigned int offset) const;

    /**
     * Getter function for real points. The offset can be obtained with
     * compute_data_index_offset().
     */
    const Point<dim, Number> *
    get_real_point(const unsigned int offset) const;

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
     * Getter function for the mapping update flags.
     */
    UpdateFlags
    get_update_flags_mapping() const;

    /**
     * Connects to is_reinitialized().
     */
    boost::signals2::connection
    connect_is_reinitialized(const std::function<void()> &set_is_reinitialized);

    /**
     * Compute the unit points index offset for the current cell/face.
     */
    unsigned int
    compute_unit_point_index_offset(const unsigned int cell_index,
                                    const unsigned int face_number) const;

    /**
     * Compute the data index offset for the current cell/face.
     */
    unsigned int
    compute_data_index_offset(const unsigned int cell_index,
                              const unsigned int face_number) const;

    /**
     * Get number of unvectorized quadrature points.
     */
    unsigned int
    get_n_q_points_unvectorized(const unsigned int cell_index,
                                const unsigned int face_number) const;

  private:
    using MappingData =
      dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                   spacedim>;

    /**
     * Compute number of quadrature point batches depending on NumberType.
     */
    template <typename NumberType>
    unsigned int
    compute_n_q_points(const unsigned int n_q_points_unvectorized);

    /**
     * Resize the unit_point data field.
     */
    void
    resize_unit_points(const unsigned int n_unit_point_batches);

    /**
     * Resize the unit_point_faces data field.
     */
    void
    resize_unit_points_faces(const unsigned int n_unit_point_batches);

    /**
     * Resize the mapping data fields.
     */
    void
    resize_data_fields(const unsigned int n_data_point_batches);

    /**
     * Store the unit points.
     */
    void
    store_unit_points(const unsigned int             unit_points_index_offset,
                      const unsigned int             n_q_points,
                      const unsigned int             n_q_points_unvectorized,
                      const std::vector<Point<dim>> &points);

    /**
     * Store the unit points on faces.
     */
    void
    store_unit_points_faces(const unsigned int unit_points_index_offset,
                            const unsigned int n_q_points,
                            const unsigned int n_q_points_unvectorized,
                            const std::vector<Point<dim - 1>> &points);

    /**
     * Store the requested mapping data.
     */
    void
    store_mapping_data(const unsigned int         unit_points_index_offset,
                       const unsigned int         n_q_points,
                       const unsigned int         n_q_points_unvectorized,
                       const MappingData &        mapping_data,
                       const std::vector<double> &weights);

    /**
     * Compute the compressed cell index.
     */
    unsigned int
    compute_compressed_cell_index(const unsigned int cell_index) const;

    /**
     * Compute the geometry index offset of the current cell/face.
     */
    unsigned int
    compute_geometry_index_offset(const unsigned int cell_index,
                                  const unsigned int face_number) const;

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
     * The reference points specified at reinit().
     *
     * Indexed by @p unit_points_index.
     */
    AlignedVector<Point<dim, VectorizedArrayType>> unit_points;

    /**
     * The reference points on faces specified at reinit().
     *
     * Indexed by @p unit_points_index.
     */
    AlignedVector<Point<dim - 1, VectorizedArrayType>> unit_points_faces;

    /**
     * Offset to point to the first unit point of a cell/face.
     */
    AlignedVector<unsigned int> unit_points_index;

    /**
     * A pointer to the internal data of the underlying mapping.
     */
    std::shared_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
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
     * AdditionalData for this object.
     */
    const AdditionalData additional_data;

    /**
     * Stores the index offset into the arrays @p JxW_values, @p jacobians,
     * @p inverse_jacobians and @p normal_vectors.
     */
    AlignedVector<unsigned int> data_index_offsets;

    /**
     * The storage of the Jacobian determinant times the quadrature weight on
     * quadrature points.
     *
     * Indexed by @p data_index_offsets.
     */
    AlignedVector<Number> JxW_values;

    /**
     * Stores the normal vectors.
     *
     * Indexed by @p data_index_offsets.
     */
    AlignedVector<Tensor<1, spacedim, Number>> normal_vectors;

    /**
     * The storage of contravariant transformation on quadrature points, i.e.,
     * the Jacobians of the transformation from the unit to the real cell.
     *
     * Indexed by @p data_index_offsets.
     */
    AlignedVector<DerivativeForm<1, dim, spacedim, Number>> jacobians;

    /**
     * The storage of covariant transformation on quadrature points, i.e.,
     * the inverse Jacobians of the transformation from the
     * unit to the real cell.
     *
     * Indexed by @p data_index_offsets.
     */
    AlignedVector<DerivativeForm<1, spacedim, dim, Number>> inverse_jacobians;

    /**
     * The mapped real points.
     *
     * Indexed by @p data_index_offsets.
     */
    AlignedVector<Point<spacedim, Number>> real_points;

    /**
     * Number of unvectorized unit points per geometric entity (cell/face).
     */
    std::vector<unsigned int> n_q_points_unvectorized;

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


  template <int dim, int spacedim, typename Number>
  MappingInfo<dim, spacedim, Number>::MappingInfo(
    const Mapping<dim> & mapping,
    const UpdateFlags    update_flags,
    const AdditionalData additional_data)
    : mapping(&mapping)
    , update_flags(update_flags)
    , update_flags_mapping(update_default)
    , additional_data(additional_data)
  {
    // translate update flags
    if (update_flags & update_jacobians)
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

    update_flags_mapping =
      internal::ComputeMappingDataHelper<dim, spacedim>::required_update_flags(
        this->mapping, update_flags_mapping);

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



  template <int dim, int spacedim, typename Number>
  void
  MappingInfo<dim, spacedim, Number>::reinit(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const std::vector<Point<dim>> &                             unit_points_in)
  {
    reinit(cell, Quadrature<dim>(unit_points_in));
  }



  template <int dim, int spacedim, typename Number>
  void
  MappingInfo<dim, spacedim, Number>::reinit(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const ArrayView<const Point<dim>> &                         unit_points_in)
  {
    reinit(cell,
           std::vector<Point<dim>>(unit_points_in.begin(),
                                   unit_points_in.end()));
  }



  template <int dim, int spacedim, typename Number>
  void
  MappingInfo<dim, spacedim, Number>::reinit(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Quadrature<dim> &                                     quadrature)
  {
    n_q_points_unvectorized.resize(1);
    n_q_points_unvectorized[0] = quadrature.size();

    const unsigned int n_q_points =
      compute_n_q_points<VectorizedArrayType>(n_q_points_unvectorized[0]);

    const unsigned int n_q_points_data =
      compute_n_q_points<Number>(n_q_points_unvectorized[0]);

    // resize data vectors
    resize_unit_points(n_q_points);
    resize_data_fields(n_q_points_data);

    // store unit points
    store_unit_points(0,
                      n_q_points,
                      n_q_points_unvectorized[0],
                      quadrature.get_points());

    // compute mapping data
    MappingData                mapping_data;
    CellSimilarity::Similarity cell_similarity = CellSimilarity::none;
    internal::ComputeMappingDataHelper<dim, spacedim>::
      compute_mapping_data_for_quadrature(mapping,
                                          update_flags_mapping,
                                          cell,
                                          cell_similarity,
                                          quadrature,
                                          internal_mapping_data,
                                          mapping_data);

    // store mapping data
    store_mapping_data(0,
                       n_q_points_data,
                       n_q_points_unvectorized[0],
                       mapping_data,
                       quadrature.get_weights());

    state = State::single_cell;
    is_reinitialized();
  }



  template <int dim, int spacedim, typename Number>
  template <typename ContainerType>
  void
  MappingInfo<dim, spacedim, Number>::reinit_cells(
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



  template <int dim, int spacedim, typename Number>
  template <typename ContainerType>
  void
  MappingInfo<dim, spacedim, Number>::reinit_cells(
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

    n_q_points_unvectorized.reserve(n_cells);

    // fill unit points index offset vector
    unit_points_index.reserve(n_cells + 1);
    unit_points_index.push_back(0);
    data_index_offsets.reserve(n_cells + 1);
    data_index_offsets.push_back(0);
    for (const auto &quadrature : quadrature_vector)
      {
        const unsigned int n_points = quadrature.size();
        n_q_points_unvectorized.push_back(n_points);

        const unsigned int n_q_points =
          compute_n_q_points<VectorizedArrayType>(n_points);
        unit_points_index.push_back(unit_points_index.back() + n_q_points);

        const unsigned int n_q_points_data =
          compute_n_q_points<Number>(n_points);
        data_index_offsets.push_back(data_index_offsets.back() +
                                     n_q_points_data);
      }

    const unsigned int n_unit_points = unit_points_index.back();
    const unsigned int n_data_points = data_index_offsets.back();

    // resize data vectors
    resize_unit_points(n_unit_points);
    resize_data_fields(n_data_points);

    if (do_cell_index_compression)
      cell_index_to_compressed_cell_index.resize(n_unfiltered_cells,
                                                 numbers::invalid_unsigned_int);

    MappingData                mapping_data;
    CellSimilarity::Similarity cell_similarity =
      CellSimilarity::Similarity::none;
    unsigned int cell_index = 0;
    for (const auto &cell : cell_iterator_range)
      {
        // store unit points
        const unsigned int n_q_points = compute_n_q_points<VectorizedArrayType>(
          n_q_points_unvectorized[cell_index]);
        store_unit_points(unit_points_index[cell_index],
                          n_q_points,
                          n_q_points_unvectorized[cell_index],
                          quadrature_vector[cell_index].get_points());

        // compute mapping data
        internal::ComputeMappingDataHelper<dim, spacedim>::
          compute_mapping_data_for_quadrature(mapping,
                                              update_flags_mapping,
                                              cell,
                                              cell_similarity,
                                              quadrature_vector[cell_index],
                                              internal_mapping_data,
                                              mapping_data);

        // store mapping data
        const unsigned int n_q_points_data =
          compute_n_q_points<Number>(n_q_points_unvectorized[cell_index]);
        store_mapping_data(data_index_offsets[cell_index],
                           n_q_points_data,
                           n_q_points_unvectorized[cell_index],
                           mapping_data,
                           quadrature_vector[cell_index].get_weights());

        if (do_cell_index_compression)
          cell_index_to_compressed_cell_index[cell->active_cell_index()] =
            cell_index;

        ++cell_index;
      }

    state = State::cell_vector;
    is_reinitialized();
  }



  template <int dim, int spacedim, typename Number>
  template <typename Iterator>
  void
  MappingInfo<dim, spacedim, Number>::reinit_surface(
    const IteratorRange<Iterator> &                    cell_iterator_range,
    const std::vector<ImmersedSurfaceQuadrature<dim>> &quadrature_vector,
    const unsigned int                                 n_unfiltered_cells)
  {
    Assert(
      additional_data.use_global_weights == false,
      ExcMessage(
        "There is no known use-case for AdditionalData::use_global_weights=true and reinit_surface()"));

    do_cell_index_compression =
      n_unfiltered_cells != numbers::invalid_unsigned_int;

    if (update_flags_mapping & (update_JxW_values | update_normal_vectors))
      update_flags_mapping |= update_covariant_transformation;

    const unsigned int n_cells = quadrature_vector.size();
    AssertDimension(n_cells,
                    std::distance(cell_iterator_range.begin(),
                                  cell_iterator_range.end()));

    n_q_points_unvectorized.reserve(n_cells);

    // fill unit points index offset vector
    unit_points_index.reserve(n_cells + 1);
    unit_points_index.push_back(0);
    data_index_offsets.reserve(n_cells + 1);
    data_index_offsets.push_back(0);
    for (const auto &quadrature : quadrature_vector)
      {
        const unsigned int n_points = quadrature.size();
        n_q_points_unvectorized.push_back(n_points);

        const unsigned int n_q_points =
          compute_n_q_points<VectorizedArrayType>(n_points);
        unit_points_index.push_back(unit_points_index.back() + n_q_points);

        const unsigned int n_q_points_data =
          compute_n_q_points<Number>(n_points);
        data_index_offsets.push_back(data_index_offsets.back() +
                                     n_q_points_data);
      }

    const unsigned int n_unit_points = unit_points_index.back();
    const unsigned int n_data_points = data_index_offsets.back();

    // resize data vectors
    resize_unit_points(n_unit_points);
    resize_data_fields(n_data_points);

    if (do_cell_index_compression)
      cell_index_to_compressed_cell_index.resize(n_unfiltered_cells,
                                                 numbers::invalid_unsigned_int);

    MappingData  mapping_data;
    unsigned int cell_index = 0;
    for (const auto &cell : cell_iterator_range)
      {
        const auto &quadrature = quadrature_vector[cell_index];

        // store unit points
        const unsigned int n_q_points = compute_n_q_points<VectorizedArrayType>(
          n_q_points_unvectorized[cell_index]);
        store_unit_points(unit_points_index[cell_index],
                          n_q_points,
                          n_q_points_unvectorized[cell_index],
                          quadrature_vector[cell_index].get_points());

        // compute mapping data
        internal::ComputeMappingDataHelper<dim, spacedim>::
          compute_mapping_data_for_immersed_surface_quadrature(
            mapping,
            update_flags_mapping,
            cell,
            quadrature,
            internal_mapping_data,
            mapping_data);

        // store mapping data
        const unsigned int n_q_points_data =
          compute_n_q_points<Number>(n_q_points_unvectorized[cell_index]);
        store_mapping_data(data_index_offsets[cell_index],
                           n_q_points_data,
                           n_q_points_unvectorized[cell_index],
                           mapping_data,
                           quadrature_vector[cell_index].get_weights());

        if (do_cell_index_compression)
          cell_index_to_compressed_cell_index[cell->active_cell_index()] =
            cell_index;

        ++cell_index;
      }

    state = State::cell_vector;
    is_reinitialized();
  }



  template <int dim, int spacedim, typename Number>
  template <typename Iterator>
  void
  MappingInfo<dim, spacedim, Number>::reinit_faces(
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

    n_q_points_unvectorized.reserve(n_faces);

    // fill unit points index offset vector
    unit_points_index.resize(n_faces + 1);
    data_index_offsets.resize(n_faces + 1);
    cell_index                 = 0;
    unsigned int n_unit_points = 0;
    unsigned int n_data_points = 0;
    for (const auto &cell : cell_iterator_range)
      {
        for (const auto &f : cell->face_indices())
          {
            const unsigned int current_face_index =
              cell_index_offset[cell_index] + f;

            unit_points_index[current_face_index]  = n_unit_points;
            data_index_offsets[current_face_index] = n_data_points;

            const unsigned int n_points =
              quadrature_vector[cell_index][f].size();
            n_q_points_unvectorized.push_back(n_points);

            const unsigned int n_q_points =
              compute_n_q_points<VectorizedArrayType>(n_points);
            n_unit_points += n_q_points;

            const unsigned int n_q_points_data =
              compute_n_q_points<Number>(n_points);
            n_data_points += n_q_points_data;
          }

        ++cell_index;
      }
    unit_points_index[n_faces]  = n_unit_points;
    data_index_offsets[n_faces] = n_data_points;

    // compress indices
    if (do_cell_index_compression)
      cell_index_to_compressed_cell_index.resize(n_unfiltered_cells,
                                                 numbers::invalid_unsigned_int);

    // fill unit points and mapping data for every face of all cells
    // resize data vectors
    resize_unit_points_faces(n_unit_points);
    resize_unit_points(n_unit_points);
    resize_data_fields(n_data_points);

    MappingData mapping_data;
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

            const unsigned int current_face_index =
              cell_index_offset[cell_index] + f;

            // store unit points
            const unsigned int n_q_points =
              compute_n_q_points<VectorizedArrayType>(
                n_q_points_unvectorized[current_face_index]);
            store_unit_points(unit_points_index[current_face_index],
                              n_q_points,
                              n_q_points_unvectorized[current_face_index],
                              quadrature_on_cell.get_points());

            store_unit_points_faces(unit_points_index[current_face_index],
                                    n_q_points,
                                    n_q_points_unvectorized[current_face_index],
                                    quadrature_on_face.get_points());

            internal::ComputeMappingDataHelper<dim, spacedim>::
              compute_mapping_data_for_face_quadrature(mapping,
                                                       update_flags_mapping,
                                                       cell,
                                                       f,
                                                       quadrature_on_face,
                                                       internal_mapping_data,
                                                       mapping_data);

            const unsigned int n_q_points_data = compute_n_q_points<Number>(
              n_q_points_unvectorized[current_face_index]);
            store_mapping_data(data_index_offsets[current_face_index],
                               n_q_points_data,
                               n_q_points_unvectorized[current_face_index],
                               mapping_data,
                               quadrature_on_face.get_weights());
          }
        if (do_cell_index_compression)
          cell_index_to_compressed_cell_index[cell->active_cell_index()] =
            cell_index;

        ++cell_index;
      }

    state = State::faces_on_cells_in_vector;
    is_reinitialized();
  }



  template <int dim, int spacedim, typename Number>
  bool
  MappingInfo<dim, spacedim, Number>::is_face_state() const
  {
    return state == State::faces_on_cells_in_vector;
  }



  template <int dim, int spacedim, typename Number>
  unsigned int
  MappingInfo<dim, spacedim, Number>::get_n_q_points_unvectorized(
    const unsigned int cell_index,
    const unsigned int face_number) const
  {
    if (cell_index == numbers::invalid_unsigned_int &&
        face_number == numbers::invalid_unsigned_int)
      {
        Assert(state == State::single_cell,
               ExcMessage(
                 "This mapping info is not reinitialized for a single cell!"));
        return n_q_points_unvectorized[0];
      }
    else
      {
        return n_q_points_unvectorized[compute_geometry_index_offset(
          cell_index, face_number)];
      }
  }



  template <int dim, int spacedim, typename Number>
  template <typename NumberType>
  unsigned int
  MappingInfo<dim, spacedim, Number>::compute_n_q_points(
    const unsigned int n_q_points_unvectorized)
  {
    const unsigned int n_lanes =
      dealii::internal::VectorizedArrayTrait<NumberType>::width();
    const unsigned int n_filled_lanes_last_batch =
      n_q_points_unvectorized % n_lanes;
    unsigned int n_q_points = n_q_points_unvectorized / n_lanes;
    if (n_filled_lanes_last_batch > 0)
      ++n_q_points;
    return n_q_points;
  }



  template <int dim, int spacedim, typename Number>
  unsigned int
  MappingInfo<dim, spacedim, Number>::compute_geometry_index_offset(
    const unsigned int cell_index,
    const unsigned int face_number) const
  {
    const unsigned int compressed_cell_index =
      compute_compressed_cell_index(cell_index);
    if (face_number == numbers::invalid_unsigned_int)
      {
        Assert(state == State::cell_vector,
               ExcMessage(
                 "This mapping info is not reinitialized for a cell vector!"));
        return compressed_cell_index;
      }
    else
      {
        Assert(cell_index != numbers::invalid_unsigned_int,
               ExcMessage(
                 "cell_index has to be set if face_number is specified!"));
        Assert(state == State::faces_on_cells_in_vector,
               ExcMessage("This mapping info is not reinitialized for faces"
                          " on cells in a vector!"));
        return cell_index_offset[compressed_cell_index] + face_number;
      }
  }



  template <int dim, int spacedim, typename Number>
  unsigned int
  MappingInfo<dim, spacedim, Number>::compute_compressed_cell_index(
    const unsigned int cell_index) const
  {
    if (do_cell_index_compression)
      {
        Assert(cell_index_to_compressed_cell_index[cell_index] !=
                 numbers::invalid_unsigned_int,
               ExcMessage("Mapping info object was not initialized for this"
                          " active cell index!"));
        return cell_index_to_compressed_cell_index[cell_index];
      }
    else
      return cell_index;
  }


  template <int dim, int spacedim, typename Number>
  void
  MappingInfo<dim, spacedim, Number>::store_unit_points(
    const unsigned int             unit_points_index_offset,
    const unsigned int             n_q_points,
    const unsigned int             n_q_points_unvectorized,
    const std::vector<Point<dim>> &points)
  {
    const unsigned int n_lanes =
      dealii::internal::VectorizedArrayTrait<VectorizedArrayType>::width();

    for (unsigned int q = 0; q < n_q_points; ++q)
      {
        const unsigned int offset = unit_points_index_offset + q;
        for (unsigned int v = 0;
             v < n_lanes && q * n_lanes + v < n_q_points_unvectorized;
             ++v)
          for (unsigned int d = 0; d < dim; ++d)
            dealii::internal::VectorizedArrayTrait<VectorizedArrayType>::get(
              unit_points[offset][d], v) = points[q * n_lanes + v][d];
      }
  }



  template <int dim, int spacedim, typename Number>
  void
  MappingInfo<dim, spacedim, Number>::store_unit_points_faces(
    const unsigned int                 unit_points_index_offset,
    const unsigned int                 n_q_points,
    const unsigned int                 n_q_points_unvectorized,
    const std::vector<Point<dim - 1>> &points)
  {
    const unsigned int n_lanes =
      dealii::internal::VectorizedArrayTrait<VectorizedArrayType>::width();

    for (unsigned int q = 0; q < n_q_points; ++q)
      {
        const unsigned int offset = unit_points_index_offset + q;
        for (unsigned int v = 0;
             v < n_lanes && q * n_lanes + v < n_q_points_unvectorized;
             ++v)
          for (unsigned int d = 0; d < dim - 1; ++d)
            dealii::internal::VectorizedArrayTrait<VectorizedArrayType>::get(
              unit_points_faces[offset][d], v) = points[q * n_lanes + v][d];
      }
  }



  template <int dim, int spacedim, typename Number>
  void
  MappingInfo<dim, spacedim, Number>::store_mapping_data(
    const unsigned int              unit_points_index_offset,
    const unsigned int              n_q_points,
    const unsigned int              n_q_points_unvectorized,
    const MappingInfo::MappingData &mapping_data,
    const std::vector<double> &     weights)
  {
    const unsigned int n_lanes =
      dealii::internal::VectorizedArrayTrait<Number>::width();

    for (unsigned int q = 0; q < n_q_points; ++q)
      {
        const unsigned int offset = unit_points_index_offset + q;
        for (unsigned int v = 0;
             v < n_lanes && q * n_lanes + v < n_q_points_unvectorized;
             ++v)
          {
            if (update_flags_mapping & UpdateFlags::update_jacobians)
              for (unsigned int d = 0; d < dim; ++d)
                for (unsigned int s = 0; s < spacedim; ++s)
                  dealii::internal::VectorizedArrayTrait<Number>::get(
                    jacobians[offset][d][s], v) =
                    mapping_data.jacobians[q * n_lanes + v][d][s];
            if (update_flags_mapping & UpdateFlags::update_inverse_jacobians)
              for (unsigned int d = 0; d < dim; ++d)
                for (unsigned int s = 0; s < spacedim; ++s)
                  dealii::internal::VectorizedArrayTrait<Number>::get(
                    inverse_jacobians[offset][s][d], v) =
                    mapping_data.inverse_jacobians[q * n_lanes + v][s][d];
            if (update_flags_mapping & UpdateFlags::update_JxW_values)
              {
                if (additional_data.use_global_weights)
                  {
                    dealii::internal::VectorizedArrayTrait<Number>::get(
                      JxW_values[offset], v) = weights[q * n_lanes + v];
                  }
                else
                  {
                    dealii::internal::VectorizedArrayTrait<Number>::get(
                      JxW_values[offset], v) =
                      mapping_data.JxW_values[q * n_lanes + v];
                  }
              }
            if (update_flags_mapping & UpdateFlags::update_normal_vectors)
              for (unsigned int s = 0; s < spacedim; ++s)
                dealii::internal::VectorizedArrayTrait<Number>::get(
                  normal_vectors[offset][s], v) =
                  mapping_data.normal_vectors[q * n_lanes + v][s];
            if (update_flags_mapping & UpdateFlags::update_quadrature_points)
              for (unsigned int s = 0; s < spacedim; ++s)
                dealii::internal::VectorizedArrayTrait<Number>::get(
                  real_points[offset][s], v) =
                  mapping_data.quadrature_points[q * n_lanes + v][s];
          }
      }
  }



  template <int dim, int spacedim, typename Number>
  void
  MappingInfo<dim, spacedim, Number>::resize_unit_points(
    const unsigned int n_unit_point_batches)
  {
    unit_points.resize(n_unit_point_batches);
  }



  template <int dim, int spacedim, typename Number>
  void
  MappingInfo<dim, spacedim, Number>::resize_unit_points_faces(
    const unsigned int n_unit_point_batches)
  {
    unit_points_faces.resize(n_unit_point_batches);
  }



  template <int dim, int spacedim, typename Number>
  void
  MappingInfo<dim, spacedim, Number>::resize_data_fields(
    const unsigned int n_data_point_batches)
  {
    if (update_flags_mapping & UpdateFlags::update_jacobians)
      jacobians.resize(n_data_point_batches);
    if (update_flags_mapping & UpdateFlags::update_inverse_jacobians)
      inverse_jacobians.resize(n_data_point_batches);
    if (update_flags_mapping & UpdateFlags::update_JxW_values)
      JxW_values.resize(n_data_point_batches);
    if (update_flags_mapping & UpdateFlags::update_normal_vectors)
      normal_vectors.resize(n_data_point_batches);
    if (update_flags_mapping & UpdateFlags::update_quadrature_points)
      real_points.resize(n_data_point_batches);
  }



  template <int dim, int spacedim, typename Number>
  inline const Point<
    dim,
    typename MappingInfo<dim, spacedim, Number>::VectorizedArrayType> *
  MappingInfo<dim, spacedim, Number>::get_unit_point(
    const unsigned int offset) const
  {
    return &unit_points[offset];
  }



  template <int dim, int spacedim, typename Number>
  inline const Point<
    dim - 1,
    typename MappingInfo<dim, spacedim, Number>::VectorizedArrayType> *
  MappingInfo<dim, spacedim, Number>::get_unit_point_faces(
    const unsigned int offset) const
  {
    return &unit_points_faces[offset];
  }



  template <int dim, int spacedim, typename Number>
  inline const Point<dim, Number> *
  MappingInfo<dim, spacedim, Number>::get_real_point(
    const unsigned int offset) const
  {
    return &real_points[offset];
  }



  template <int dim, int spacedim, typename Number>
  unsigned int
  MappingInfo<dim, spacedim, Number>::compute_unit_point_index_offset(
    const unsigned int cell_index,
    const unsigned int face_number) const
  {
    if (cell_index == numbers::invalid_unsigned_int &&
        face_number == numbers::invalid_unsigned_int)
      {
        Assert(state == State::single_cell,
               ExcMessage(
                 "This mapping info is not reinitialized for a single cell!"));
        return 0;
      }
    else
      {
        const unsigned int offset =
          compute_geometry_index_offset(cell_index, face_number);
        return unit_points_index[offset];
      }
  }



  template <int dim, int spacedim, typename Number>
  unsigned int
  MappingInfo<dim, spacedim, Number>::compute_data_index_offset(
    const unsigned int cell_index,
    const unsigned int face_number) const
  {
    if (cell_index == numbers::invalid_unsigned_int &&
        face_number == numbers::invalid_unsigned_int)
      {
        Assert(state == State::single_cell,
               ExcMessage(
                 "This mapping info is not reinitialized for a single cell!"));
        return 0;
      }
    else
      {
        const unsigned int offset =
          compute_geometry_index_offset(cell_index, face_number);
        return data_index_offsets[offset];
      }
  }



  template <int dim, int spacedim, typename Number>
  inline const DerivativeForm<1, dim, spacedim, Number> *
  MappingInfo<dim, spacedim, Number>::get_jacobian(
    const unsigned int offset) const
  {
    return &jacobians[offset];
  }



  template <int dim, int spacedim, typename Number>
  inline const DerivativeForm<1, spacedim, dim, Number> *
  MappingInfo<dim, spacedim, Number>::get_inverse_jacobian(
    const unsigned int offset) const
  {
    return &inverse_jacobians[offset];
  }


  template <int dim, int spacedim, typename Number>
  inline const Tensor<1, spacedim, Number> *
  MappingInfo<dim, spacedim, Number>::get_normal_vector(
    const unsigned int offset) const
  {
    return &normal_vectors[offset];
  }



  template <int dim, int spacedim, typename Number>
  inline const Number *
  MappingInfo<dim, spacedim, Number>::get_JxW(const unsigned int offset) const
  {
    return &JxW_values[offset];
  }



  template <int dim, int spacedim, typename Number>
  const Mapping<dim, spacedim> &
  MappingInfo<dim, spacedim, Number>::get_mapping() const
  {
    return *mapping;
  }



  template <int dim, int spacedim, typename Number>
  UpdateFlags
  MappingInfo<dim, spacedim, Number>::get_update_flags() const
  {
    return update_flags;
  }



  template <int dim, int spacedim, typename Number>
  UpdateFlags
  MappingInfo<dim, spacedim, Number>::get_update_flags_mapping() const
  {
    return update_flags_mapping;
  }



  template <int dim, int spacedim, typename Number>
  boost::signals2::connection
  MappingInfo<dim, spacedim, Number>::connect_is_reinitialized(
    const std::function<void()> &set_is_reinitialized)
  {
    return is_reinitialized.connect(set_is_reinitialized);
  }
} // namespace NonMatching

DEAL_II_NAMESPACE_CLOSE

#endif

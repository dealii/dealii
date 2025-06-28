// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_non_matching_mapping_info_h
#define dealii_non_matching_mapping_info_h


#include <deal.II/base/config.h>

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/floating_point_comparator.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_cartesian.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_related_data.h>

#include <deal.II/matrix_free/mapping_info_storage.h>

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
        const ObserverPointer<const Mapping<dim, spacedim>> &mapping,
        const UpdateFlags                                   &update_flags)
      {
        return mapping->requires_update_flags(update_flags);
      }

      static void
      compute_mapping_data_for_quadrature(
        const ObserverPointer<const Mapping<dim, spacedim>> &mapping,
        const UpdateFlags &update_flags_mapping,
        const typename Triangulation<dim, spacedim>::cell_iterator &cell,
        CellSimilarity::Similarity &cell_similarity,
        const Quadrature<dim>      &quadrature,
        std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
                    &internal_mapping_data,
        MappingData &mapping_data)
      {
        mapping_data.initialize(quadrature.size(), update_flags_mapping);
        internal_mapping_data->reinit(update_flags_mapping, quadrature);

        cell_similarity = mapping->fill_fe_values(cell,
                                                  cell_similarity,
                                                  quadrature,
                                                  *internal_mapping_data,
                                                  mapping_data);
      }



      static void
      compute_mapping_data_for_immersed_surface_quadrature(
        const ObserverPointer<const Mapping<dim, spacedim>> &mapping,
        const UpdateFlags &update_flags_mapping,
        const typename Triangulation<dim, spacedim>::cell_iterator &cell,
        const ImmersedSurfaceQuadrature<dim>                       &quadrature,
        std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
                    &internal_mapping_data,
        MappingData &mapping_data)
      {
        mapping_data.initialize(quadrature.size(), update_flags_mapping);

        internal_mapping_data->reinit(update_flags_mapping, quadrature);

        mapping->fill_fe_immersed_surface_values(cell,
                                                 quadrature,
                                                 *internal_mapping_data,
                                                 mapping_data);
      }



      static void
      compute_mapping_data_for_face_quadrature(
        const ObserverPointer<const Mapping<dim, spacedim>> &mapping,
        const UpdateFlags &update_flags_mapping,
        const typename Triangulation<dim, spacedim>::cell_iterator &cell,
        const unsigned int                                          face_no,
        const Quadrature<dim - 1>                                  &quadrature,
        std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
                    &internal_mapping_data,
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
                                 QProjector<dim>::project_to_face(
                                   cell->reference_cell(),
                                   quadrature,
                                   face_no,
                                   cell->combined_face_orientation(face_no)),
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

    template <int dim, int spacedim = dim>
    dealii::internal::MatrixFreeFunctions::GeometryType
    compute_geometry_type(
      const double diameter,
      const std::vector<DerivativeForm<1, spacedim, dim, double>>
        &inverse_jacobians)
    {
      const auto   jac_0 = inverse_jacobians[0];
      const double zero_tolerance_double =
        1e4 / diameter * std::numeric_limits<double>::epsilon() * 1024.;
      bool jacobian_constant = true;
      for (unsigned int q = 1; q < inverse_jacobians.size(); ++q)
        {
          const DerivativeForm<1, spacedim, dim> &jac = inverse_jacobians[q];
          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = 0; e < spacedim; ++e)
              if (std::fabs(jac_0[d][e] - jac[d][e]) > zero_tolerance_double)
                jacobian_constant = false;
          if (!jacobian_constant)
            break;
        }

      // check whether the Jacobian is diagonal to machine
      // accuracy
      bool cell_cartesian = jacobian_constant && dim == spacedim;
      for (unsigned int d = 0; d < dim; ++d)
        for (unsigned int e = 0; e < dim; ++e)
          if (d != e)
            if (std::fabs(jac_0[d][e]) > zero_tolerance_double)
              {
                cell_cartesian = false;
                break;
              }

      // return cell type
      if (cell_cartesian)
        return dealii::internal::MatrixFreeFunctions::GeometryType::cartesian;
      else if (jacobian_constant)
        return dealii::internal::MatrixFreeFunctions::GeometryType::affine;
      else
        return dealii::internal::MatrixFreeFunctions::GeometryType::general;
    }
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
   *
   * @note This class cannot be copied or copy-constructed. The intention of
   * this class is to be available as a single use, e.g. inside
   * FEPointEvaluation or as a storage container for the whole mesh. Use a
   * shared or unique pointer of this class in case several objects should be
   * held.
   */
  template <int dim, int spacedim = dim, typename Number = double>
  class MappingInfo : public EnableObserverPointer
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
      AdditionalData(const bool use_global_weights = false,
                     const bool store_cells        = false)
        : use_global_weights(use_global_weights)
        , store_cells(store_cells)
      {}

      /**
       * During initialization, assume that the Quadrature object contains
       * global weights as, e.g., obtained by
       * QSimplex::compute_affine_transformation().
       */
      bool use_global_weights;

      /**
       * During the reinit_cells() function call, cells are passed as
       * argument. In the default case, the cell is not stored,
       * since all relevant mapping related information is precomputed.
       * However, this flag enables that the cells are stored so that
       * they can be accessed later on.
       */
      bool store_cells;
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
    MappingInfo(const Mapping<dim, spacedim> &mapping,
                const UpdateFlags             update_flags,
                const AdditionalData additional_data = AdditionalData());

    /**
     * Do not allow making copies.
     */
    MappingInfo(const MappingInfo &) = delete;

    /**
     * Do not allow copy assignment of this class.
     */
    MappingInfo &
    operator=(const MappingInfo &) = delete;

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
      const ContainerType                        &cell_iterator_range,
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
      const ContainerType                &cell_iterator_range,
      const std::vector<Quadrature<dim>> &quadrature_vector,
      const unsigned int n_unfiltered_cells = numbers::invalid_unsigned_int);

    /**
     * Compute the mapping information for the incoming vector of cells and
     * corresponding vector of ImmersedSurfaceQuadrature.
     */
    template <typename ContainerType>
    void
    reinit_surface(
      const ContainerType                               &cell_iterator_range,
      const std::vector<ImmersedSurfaceQuadrature<dim>> &quadrature_vector,
      const unsigned int n_unfiltered_cells = numbers::invalid_unsigned_int);

    /**
     * Compute the mapping information for all faces of the incoming vector
     * of cells and corresponding vector of quadratures.
     */
    template <typename ContainerType>
    void
    reinit_faces(
      const ContainerType                                 &cell_iterator_range,
      const std::vector<std::vector<Quadrature<dim - 1>>> &quadrature_vector,
      const unsigned int n_unfiltered_cells = numbers::invalid_unsigned_int);

    /**
     * Compute the mapping information incoming vector of faces and
     * corresponding vector of quadratures.
     */
    template <typename CellIteratorType>
    void
    reinit_faces(const std::vector<std::pair<CellIteratorType, unsigned int>>
                   &face_iterator_range_interior,
                 const std::vector<Quadrature<dim - 1>> &quadrature_vector);

    /**
     * Return if this MappingInfo object is reinitialized for faces (by
     * reinit_faces()) or not.
     */
    bool
    is_face_state() const;

    /**
     * Returns the face number of the interior/exterior face.
     */
    unsigned int
    get_face_number(const unsigned int offset, const bool is_interior) const;

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
    get_jacobian(const unsigned int offset,
                 const bool         is_interior = true) const;

    /**
     * Getter function for inverse Jacobians. The offset can be obtained with
     * compute_data_index_offset().
     */
    const DerivativeForm<1, spacedim, dim, Number> *
    get_inverse_jacobian(const unsigned int offset,
                         const bool         is_interior = true) const;

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
    const Point<spacedim, Number> *
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
     * Compute the geometry index offset of the current cell/face.
     */
    template <bool is_face>
    unsigned int
    compute_geometry_index_offset(const unsigned int cell_index,
                                  const unsigned int face_number) const;

    /**
     * Compute the unit points index offset for the current cell/face.
     */
    unsigned int
    compute_unit_point_index_offset(const unsigned int geometry_index) const;

    /**
     * Compute the data index offset for the current cell/face.
     */
    unsigned int
    compute_data_index_offset(const unsigned int geometry_index) const;

    /**
     * Compute the data index offset for the current cell/face.
     */
    unsigned int
    compute_compressed_data_index_offset(
      const unsigned int geometry_index) const;

    /**
     * Get number of unvectorized quadrature points.
     */
    unsigned int
    get_n_q_points_unvectorized(const unsigned int geometry_index) const;

    /**
     * Get cell geometry type.
     */
    dealii::internal::MatrixFreeFunctions::GeometryType
    get_cell_type(const unsigned int geometry_index) const;

    /**
     * Return cell iterator.
     *
     * @note This call is only possible if AdditionalData::store_cells is enabled.
     */
    typename Triangulation<dim, spacedim>::cell_iterator
    get_cell_iterator(const unsigned int cell_index) const;

    /**
     * Return the memory consumption of this class in bytes.
     */
    std::size_t
    memory_consumption() const;

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
     * Clears fields to make the object reusable.
     */
    void
    clear();

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
    resize_data_fields(const unsigned int n_data_point_batches,
                       const bool         is_face_centric = false);

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
                       const MappingData         &mapping_data,
                       const std::vector<double> &weights,
                       const unsigned int compressed_unit_point_index_offset,
                       const bool         affine_cell,
                       const bool         is_interior = true);

    /**
     * Compute the compressed cell index.
     */
    unsigned int
    compute_compressed_cell_index(const unsigned int cell_index) const;

    /**
     * Compute the mapping information for cells/surface.
     */
    template <typename ContainerType, typename QuadratureType>
    void
    do_reinit_cells(
      const ContainerType               &cell_iterator_range,
      const std::vector<QuadratureType> &quadrature_vector,
      const unsigned int                 n_unfiltered_cells,
      const std::function<
        void(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
             const QuadratureType &quadrature,
             MappingData          &mapping_data)> &compute_mapping_data);

    /**
     * Enum class for reinitialized states.
     */
    enum class State
    {
      invalid,
      single_cell,
      cell_vector,
      faces_on_cells_in_vector,
      face_vector
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
    std::vector<unsigned int> unit_points_index;

    /**
     * A pointer to the internal data of the underlying mapping.
     */
    std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
      internal_mapping_data;

    /**
     * Quadrature object that is used for the reinit(ArrayView) function.
     */
    Quadrature<dim> quadrature;

    /**
     * Helper class that temporarily holds the data requested for one cell by
     * Mapping::fill_fe_values() before it is filled into appropriate data
     * structures for consumption by FEPointEvaluation.
     */
    MappingData mapping_data;

    /**
     * A pointer to the underlying mapping.
     */
    const ObserverPointer<const Mapping<dim, spacedim>> mapping;

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
     * Stores whether a cell is Cartesian (cell type 0), has constant
     * transform data (Jacobians) (cell type 1), or is general (cell type
     * 3). Type 2 is only used for faces and no cells are assigned this
     * value.
     */
    std::vector<dealii::internal::MatrixFreeFunctions::GeometryType> cell_type;

    /**
     * Stores the index offset into the arrays @p JxW_values and @p normal_vectors.
     */
    std::vector<unsigned int> data_index_offsets;

    /**
     * Stores the index offset into the arrays @p jacobians and @p inverse_jacobians.
     */
    std::vector<unsigned int> compressed_data_index_offsets;

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
     * Indexed by @p compressed_data_index_offsets.
     */
    std::array<AlignedVector<DerivativeForm<1, dim, spacedim, Number>>, 2>
      jacobians;

    /**
     * The storage of covariant transformation on quadrature points, i.e.,
     * the inverse Jacobians of the transformation from the
     * unit to the real cell.
     *
     * Indexed by @p compressed_data_index_offsets.
     */
    std::array<AlignedVector<DerivativeForm<1, spacedim, dim, Number>>, 2>
      inverse_jacobians;

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
     * Reference to the triangulation passed via the cells to the
     * reinit functions. This field is only set if
     * AdditionalData::store_cells is enabled.
     */
    ObserverPointer<const Triangulation<dim, spacedim>> triangulation;

    /**
     * Level and indices of cells passed to the reinit functions. This
     * vector is only filled if AdditionalData::store_cells is enabled.
     */
    std::vector<std::pair<int, int>> cell_level_and_indices;

    /**
     * Store the face number of interior or exterior faces as set up in the
     * initialization.
     */
    std::vector<std::pair<unsigned char, unsigned char>> face_number;
  };

  template <int dim, int spacedim, typename Number>
  inline unsigned int
  MappingInfo<dim, spacedim, Number>::get_face_number(
    const unsigned int offset,
    const bool         is_interior) const
  {
    const auto &face_pair = face_number[offset];
    return is_interior ? face_pair.first : face_pair.second;
  }

  // ----------------------- template functions ----------------------


  template <int dim, int spacedim, typename Number>
  MappingInfo<dim, spacedim, Number>::MappingInfo(
    const Mapping<dim, spacedim> &mapping,
    const UpdateFlags             update_flags,
    const AdditionalData          additional_data)
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

    // construct internal_mapping_data for mappings for reuse in reinit()
    // calls to avoid frequent memory allocations
    internal_mapping_data = mapping.get_data(update_flags, Quadrature<dim>());
  }



  template <int dim, int spacedim, typename Number>
  void
  MappingInfo<dim, spacedim, Number>::clear()
  {
    n_q_points_unvectorized.clear();
    unit_points_index.clear();
    data_index_offsets.clear();
    compressed_data_index_offsets.clear();
    cell_type.clear();
  }



  template <int dim, int spacedim, typename Number>
  void
  MappingInfo<dim, spacedim, Number>::reinit(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const std::vector<Point<dim>>                              &unit_points_in)
  {
    reinit(cell, make_array_view(unit_points_in));
  }



  template <int dim, int spacedim, typename Number>
  void
  MappingInfo<dim, spacedim, Number>::reinit(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const ArrayView<const Point<dim>>                          &unit_points_in)
  {
    quadrature.initialize(unit_points_in);
    reinit(cell, quadrature);
  }



  template <int dim, int spacedim, typename Number>
  void
  MappingInfo<dim, spacedim, Number>::reinit(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Quadrature<dim>                                      &quadrature)
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
    CellSimilarity::Similarity cell_similarity = CellSimilarity::none;
    internal::ComputeMappingDataHelper<dim, spacedim>::
      compute_mapping_data_for_quadrature(mapping,
                                          update_flags_mapping,
                                          cell,
                                          cell_similarity,
                                          quadrature,
                                          internal_mapping_data,
                                          mapping_data);

    // check for cartesian/affine cell
    if (!quadrature.empty() &&
        update_flags_mapping & UpdateFlags::update_inverse_jacobians)
      {
        cell_type.push_back(
          internal::compute_geometry_type(cell->diameter(),
                                          mapping_data.inverse_jacobians));
      }
    else
      cell_type.push_back(
        dealii::internal::MatrixFreeFunctions::GeometryType::general);

    // store mapping data
    store_mapping_data(
      0,
      n_q_points_data,
      n_q_points_unvectorized[0],
      mapping_data,
      quadrature.get_weights(),
      0,
      cell_type.back() <=
        dealii::internal::MatrixFreeFunctions::GeometryType::affine);

    unit_points_index.push_back(0);
    data_index_offsets.push_back(0);
    compressed_data_index_offsets.push_back(0);

    state = State::single_cell;
  }



  template <int dim, int spacedim, typename Number>
  template <typename ContainerType>
  void
  MappingInfo<dim, spacedim, Number>::reinit_cells(
    const ContainerType                        &cell_iterator_range,
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
  template <typename ContainerType, typename QuadratureType>
  void
  MappingInfo<dim, spacedim, Number>::do_reinit_cells(
    const ContainerType               &cell_iterator_range,
    const std::vector<QuadratureType> &quadrature_vector,
    const unsigned int                 n_unfiltered_cells,
    const std::function<
      void(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
           const QuadratureType &quadrature,
           MappingData          &mapping_data)> &compute_mapping_data)
  {
    clear();

    do_cell_index_compression =
      n_unfiltered_cells != numbers::invalid_unsigned_int;

    const unsigned int n_cells = quadrature_vector.size();
    AssertDimension(n_cells,
                    std::distance(cell_iterator_range.begin(),
                                  cell_iterator_range.end()));

    n_q_points_unvectorized.reserve(n_cells);

    cell_type.reserve(n_cells);

    if (additional_data.store_cells)
      cell_level_and_indices.resize(n_cells);

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

    MappingData  mapping_data_previous_cell;
    unsigned int size_compressed_data = 0;
    unsigned int cell_index           = 0;
    for (const auto &cell : cell_iterator_range)
      {
        if (additional_data.store_cells)
          {
            this->triangulation                = &cell->get_triangulation();
            cell_level_and_indices[cell_index] = {cell->level(), cell->index()};
          }

        const auto &quadrature = quadrature_vector[cell_index];
        const bool  empty      = quadrature.empty();

        // store unit points
        const unsigned int n_q_points = compute_n_q_points<VectorizedArrayType>(
          n_q_points_unvectorized[cell_index]);
        store_unit_points(unit_points_index[cell_index],
                          n_q_points,
                          n_q_points_unvectorized[cell_index],
                          quadrature.get_points());

        // compute mapping data
        compute_mapping_data(cell, quadrature, mapping_data);

        // store mapping data
        const unsigned int n_q_points_data =
          compute_n_q_points<Number>(n_q_points_unvectorized[cell_index]);

        // check for cartesian/affine cell
        if (!empty &&
            update_flags_mapping & UpdateFlags::update_inverse_jacobians)
          {
            cell_type.push_back(
              internal::compute_geometry_type(cell->diameter(),
                                              mapping_data.inverse_jacobians));
          }
        else
          cell_type.push_back(
            dealii::internal::MatrixFreeFunctions::GeometryType::general);

        if (cell_index > 0)
          {
            // check if current and previous cell are affine
            const bool affine_cells =
              cell_type[cell_index] <=
                dealii::internal::MatrixFreeFunctions::affine &&
              cell_type[cell_index - 1] <=
                dealii::internal::MatrixFreeFunctions::affine;

            // create a comparator to compare inverse Jacobian of current
            // and previous cell
            FloatingPointComparator<double> comparator(
              1e4 / cell->diameter() * std::numeric_limits<double>::epsilon() *
              1024.);

            // we can only compare if current and previous cell have at least
            // one quadrature point and both cells are at least affine
            const auto comparison_result =
              (!affine_cells || mapping_data.inverse_jacobians.empty() ||
               mapping_data_previous_cell.inverse_jacobians.empty()) ?
                FloatingPointComparator<double>::ComparisonResult::less :
                comparator.compare(
                  mapping_data.inverse_jacobians[0],
                  mapping_data_previous_cell.inverse_jacobians[0]);

            // we can compress the Jacobians and inverse Jacobians if
            // inverse Jacobians are equal and cells are affine
            if (affine_cells &&
                comparison_result ==
                  FloatingPointComparator<double>::ComparisonResult::equal)
              {
                compressed_data_index_offsets.push_back(
                  compressed_data_index_offsets.back());
              }
            else
              {
                const unsigned int n_compressed_data_last_cell =
                  cell_type[cell_index - 1] <=
                      dealii::internal::MatrixFreeFunctions::affine ?
                    1 :
                    compute_n_q_points<Number>(
                      n_q_points_unvectorized[cell_index - 1]);

                compressed_data_index_offsets.push_back(
                  compressed_data_index_offsets.back() +
                  n_compressed_data_last_cell);
              }
          }
        else
          compressed_data_index_offsets.push_back(0);

        // cache mapping_data from previous cell
        mapping_data_previous_cell = mapping_data;

        store_mapping_data(data_index_offsets[cell_index],
                           n_q_points_data,
                           n_q_points_unvectorized[cell_index],
                           mapping_data,
                           quadrature.get_weights(),
                           compressed_data_index_offsets[cell_index],
                           cell_type[cell_index] <=
                             dealii::internal::MatrixFreeFunctions::affine);

        // update size of compressed data depending on cell type and handle
        // empty quadratures
        if (cell_type[cell_index] <=
            dealii::internal::MatrixFreeFunctions::affine)
          size_compressed_data = compressed_data_index_offsets.back() + 1;
        else
          size_compressed_data =
            std::max(size_compressed_data,
                     compressed_data_index_offsets.back() + n_q_points_data);

        if (do_cell_index_compression)
          cell_index_to_compressed_cell_index[cell->active_cell_index()] =
            cell_index;

        ++cell_index;
      }

    if (update_flags_mapping & UpdateFlags::update_jacobians)
      {
        jacobians[0].resize(size_compressed_data);
        jacobians[0].shrink_to_fit();
      }
    if (update_flags_mapping & UpdateFlags::update_inverse_jacobians)
      {
        inverse_jacobians[0].resize(size_compressed_data);
        inverse_jacobians[0].shrink_to_fit();
      }

    state = State::cell_vector;
  }



  template <int dim, int spacedim, typename Number>
  template <typename ContainerType>
  void
  MappingInfo<dim, spacedim, Number>::reinit_cells(
    const ContainerType                &cell_iterator_range,
    const std::vector<Quadrature<dim>> &quadrature_vector,
    const unsigned int                  n_unfiltered_cells)
  {
    auto compute_mapping_data_for_cells =
      [&](const typename Triangulation<dim, spacedim>::cell_iterator &cell,
          const Quadrature<dim> &quadrature,
          MappingData           &mapping_data) {
        CellSimilarity::Similarity cell_similarity =
          CellSimilarity::Similarity::none;
        internal::ComputeMappingDataHelper<dim, spacedim>::
          compute_mapping_data_for_quadrature(mapping,
                                              update_flags_mapping,
                                              cell,
                                              cell_similarity,
                                              quadrature,
                                              internal_mapping_data,
                                              mapping_data);
      };

    do_reinit_cells<ContainerType, Quadrature<dim>>(
      cell_iterator_range,
      quadrature_vector,
      n_unfiltered_cells,
      compute_mapping_data_for_cells);
  }



  template <int dim, int spacedim, typename Number>
  template <typename ContainerType>
  void
  MappingInfo<dim, spacedim, Number>::reinit_surface(
    const ContainerType                               &cell_iterator_range,
    const std::vector<ImmersedSurfaceQuadrature<dim>> &quadrature_vector,
    const unsigned int                                 n_unfiltered_cells)
  {
    Assert(
      additional_data.use_global_weights == false,
      ExcMessage(
        "There is no known use-case for AdditionalData::use_global_weights=true and reinit_surface()"));

    Assert(additional_data.store_cells == false, ExcNotImplemented());

    if (update_flags_mapping & (update_JxW_values | update_normal_vectors))
      update_flags_mapping |= update_covariant_transformation;

    auto compute_mapping_data_for_surface =
      [&](const typename Triangulation<dim, spacedim>::cell_iterator &cell,
          const ImmersedSurfaceQuadrature<dim> &quadrature,
          MappingData                          &mapping_data) {
        internal::ComputeMappingDataHelper<dim, spacedim>::
          compute_mapping_data_for_immersed_surface_quadrature(
            mapping,
            update_flags_mapping,
            cell,
            quadrature,
            internal_mapping_data,
            mapping_data);
      };

    do_reinit_cells<ContainerType, ImmersedSurfaceQuadrature<dim>>(
      cell_iterator_range,
      quadrature_vector,
      n_unfiltered_cells,
      compute_mapping_data_for_surface);
  }



  template <int dim, int spacedim, typename Number>
  template <typename ContainerType>
  void
  MappingInfo<dim, spacedim, Number>::reinit_faces(
    const ContainerType                                 &cell_iterator_range,
    const std::vector<std::vector<Quadrature<dim - 1>>> &quadrature_vector,
    const unsigned int                                   n_unfiltered_cells)
  {
    clear();

    Assert(additional_data.store_cells == false, ExcNotImplemented());

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

    cell_type.reserve(n_faces);

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
    resize_data_fields(n_data_points);

    MappingData  mapping_data_previous_cell;
    MappingData  mapping_data_first;
    bool         first_set            = false;
    unsigned int size_compressed_data = 0;
    cell_index                        = 0;
    for (const auto &cell : cell_iterator_range)
      {
        const auto &quadratures_on_faces = quadrature_vector[cell_index];

        Assert(quadratures_on_faces.size() == cell->n_faces(),
               ExcDimensionMismatch(quadratures_on_faces.size(),
                                    cell->n_faces()));

        for (const auto &f : cell->face_indices())
          {
            const auto &quadrature_on_face = quadratures_on_faces[f];
            const bool  empty              = quadrature_on_face.empty();

            const unsigned int current_face_index =
              cell_index_offset[cell_index] + f;

            // store unit points
            const unsigned int n_q_points =
              compute_n_q_points<VectorizedArrayType>(
                n_q_points_unvectorized[current_face_index]);
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

            // check for cartesian/affine cell
            if (!empty &&
                update_flags_mapping & UpdateFlags::update_inverse_jacobians)
              {
                cell_type.push_back(internal::compute_geometry_type(
                  cell->diameter(), mapping_data.inverse_jacobians));

                if (!first_set)
                  {
                    mapping_data_first = mapping_data;
                    first_set          = true;
                  }
              }
            else
              cell_type.push_back(
                dealii::internal::MatrixFreeFunctions::GeometryType::general);

            if (current_face_index > 0)
              {
                // check if current and previous cell are affine
                const bool affine_cells =
                  cell_type[current_face_index] <=
                    dealii::internal::MatrixFreeFunctions::affine &&
                  cell_type[current_face_index - 1] <=
                    dealii::internal::MatrixFreeFunctions::affine;

                // create a comparator to compare inverse Jacobian of current
                // and previous cell
                FloatingPointComparator<double> comparator(
                  1e4 / cell->diameter() *
                  std::numeric_limits<double>::epsilon() * 1024.);

                // we can only compare if current and previous cell have at
                // least one quadrature point and both cells are at least affine
                const auto comparison_result =
                  (!affine_cells || mapping_data.inverse_jacobians.empty() ||
                   mapping_data_previous_cell.inverse_jacobians.empty()) ?
                    FloatingPointComparator<double>::ComparisonResult::less :
                    comparator.compare(
                      mapping_data.inverse_jacobians[0],
                      mapping_data_previous_cell.inverse_jacobians[0]);

                // we can compress the Jacobians and inverse Jacobians if
                // inverse Jacobians are equal and cells are affine
                if (affine_cells &&
                    comparison_result ==
                      FloatingPointComparator<double>::ComparisonResult::equal)
                  {
                    compressed_data_index_offsets.push_back(
                      compressed_data_index_offsets.back());
                  }
                else if (first_set &&
                         (cell_type[current_face_index] <=
                          dealii::internal::MatrixFreeFunctions::affine) &&
                         (comparator.compare(
                            mapping_data.inverse_jacobians[0],
                            mapping_data_first.inverse_jacobians[0]) ==
                          FloatingPointComparator<
                            double>::ComparisonResult::equal))
                  {
                    compressed_data_index_offsets.push_back(0);
                  }
                else
                  {
                    const unsigned int n_compressed_data_last_cell =
                      cell_type[current_face_index - 1] <=
                          dealii::internal::MatrixFreeFunctions::affine ?
                        1 :
                        compute_n_q_points<Number>(
                          n_q_points_unvectorized[current_face_index - 1]);

                    compressed_data_index_offsets.push_back(
                      compressed_data_index_offsets.back() +
                      n_compressed_data_last_cell);
                  }
              }
            else
              compressed_data_index_offsets.push_back(0);

            // cache mapping_data from previous cell
            mapping_data_previous_cell = mapping_data;

            const unsigned int n_q_points_data = compute_n_q_points<Number>(
              n_q_points_unvectorized[current_face_index]);
            store_mapping_data(data_index_offsets[current_face_index],
                               n_q_points_data,
                               n_q_points_unvectorized[current_face_index],
                               mapping_data,
                               quadrature_on_face.get_weights(),
                               data_index_offsets[current_face_index],
                               cell_type[current_face_index] <=
                                 dealii::internal::MatrixFreeFunctions::affine);

            // update size of compressed data depending on cell type and handle
            // empty quadratures
            if (cell_type[current_face_index] <=
                dealii::internal::MatrixFreeFunctions::affine)
              size_compressed_data = compressed_data_index_offsets.back() + 1;
            else
              size_compressed_data =
                std::max(size_compressed_data,
                         compressed_data_index_offsets.back() +
                           n_q_points_data);
          }
        if (do_cell_index_compression)
          cell_index_to_compressed_cell_index[cell->active_cell_index()] =
            cell_index;

        ++cell_index;
      }

    if (update_flags_mapping & UpdateFlags::update_jacobians)
      {
        jacobians[0].resize(size_compressed_data);
        jacobians[0].shrink_to_fit();
      }
    if (update_flags_mapping & UpdateFlags::update_inverse_jacobians)
      {
        inverse_jacobians[0].resize(size_compressed_data);
        inverse_jacobians[0].shrink_to_fit();
      }

    state = State::faces_on_cells_in_vector;
  }



  template <int dim, int spacedim, typename Number>
  template <typename CellIteratorType>
  void
  MappingInfo<dim, spacedim, Number>::reinit_faces(
    const std::vector<std::pair<CellIteratorType, unsigned int>>
                                           &face_iterator_range_interior,
    const std::vector<Quadrature<dim - 1>> &quadrature_vector)
  {
    clear();

    do_cell_index_compression = false;

    Assert(additional_data.store_cells == false, ExcNotImplemented());


    const unsigned int n_faces = quadrature_vector.size();
    AssertDimension(n_faces,
                    std::distance(face_iterator_range_interior.begin(),
                                  face_iterator_range_interior.end()));

    n_q_points_unvectorized.reserve(n_faces);

    cell_type.reserve(n_faces);
    face_number.reserve(n_faces);

    // fill unit points index offset vector
    unit_points_index.reserve(n_faces + 1);
    unit_points_index.push_back(0);
    data_index_offsets.reserve(n_faces + 1);
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
    resize_unit_points_faces(n_unit_points);
    resize_data_fields(n_data_points, true);

    std::array<MappingData, 2> mapping_data;
    std::array<MappingData, 2> mapping_data_previous_cell;
    std::array<MappingData, 2> mapping_data_first;
    bool                       first_set            = false;
    unsigned int               size_compressed_data = 0;
    unsigned int               face_index           = 0;
    for (const auto &cell_and_f : face_iterator_range_interior)
      {
        const auto &quadrature_on_face = quadrature_vector[face_index];
        const bool  empty              = quadrature_on_face.empty();

        // get interior cell and face number
        const auto &cell_m = cell_and_f.first;
        const auto  f_m    = cell_and_f.second;

        // get exterior cell and face number
        const auto &cell_p =
          cell_m->at_boundary(f_m) ? cell_m : cell_m->neighbor(f_m);
        const auto f_p =
          cell_m->at_boundary(f_m) ? f_m : cell_m->neighbor_face_no(f_m);

        Assert(
          empty || (cell_m->level() == cell_p->level()),
          ExcMessage(
            "Intersected faces with quadrature points need to have the same "
            "refinement level!"));

        face_number.emplace_back(f_m, f_p);

        Assert(
          cell_m->combined_face_orientation(f_m) ==
              numbers::default_geometric_orientation &&
            cell_p->combined_face_orientation(f_p) ==
              numbers::default_geometric_orientation,
          ExcMessage(
            "Non standard face orientation is currently not implemented."));

        // store unit points
        const unsigned int n_q_points = compute_n_q_points<VectorizedArrayType>(
          n_q_points_unvectorized[face_index]);
        store_unit_points_faces(unit_points_index[face_index],
                                n_q_points,
                                n_q_points_unvectorized[face_index],
                                quadrature_on_face.get_points());

        // compute mapping for interior face
        internal::ComputeMappingDataHelper<dim, spacedim>::
          compute_mapping_data_for_face_quadrature(mapping,
                                                   update_flags_mapping,
                                                   cell_m,
                                                   f_m,
                                                   quadrature_on_face,
                                                   internal_mapping_data,
                                                   mapping_data[0]);

        // compute mapping for exterior face
        internal::ComputeMappingDataHelper<dim, spacedim>::
          compute_mapping_data_for_face_quadrature(mapping,
                                                   update_flags_mapping,
                                                   cell_p,
                                                   f_p,
                                                   quadrature_on_face,
                                                   internal_mapping_data,
                                                   mapping_data[1]);

        // check for cartesian/affine cell
        if (!empty &&
            update_flags_mapping & UpdateFlags::update_inverse_jacobians)
          {
            // select more general type of interior and exterior cell
            cell_type.push_back(std::max(
              internal::compute_geometry_type(
                cell_m->diameter(), mapping_data[0].inverse_jacobians),
              internal::compute_geometry_type(
                cell_m->diameter(), mapping_data[1].inverse_jacobians)));

            // cache mapping data of first cell pair with non-empty quadrature
            // on the face
            if (!first_set)
              {
                mapping_data_first = mapping_data;
                first_set          = true;
              }
          }
        else
          cell_type.push_back(
            dealii::internal::MatrixFreeFunctions::GeometryType::general);

        if (face_index > 0)
          {
            // check if current and previous cell pairs are affine
            const bool affine_cells =
              cell_type[face_index] <=
                dealii::internal::MatrixFreeFunctions::affine &&
              cell_type[face_index - 1] <=
                dealii::internal::MatrixFreeFunctions::affine;

            // create a comparator to compare inverse Jacobian of current
            // and previous cell pair
            FloatingPointComparator<double> comparator(
              1e4 / cell_m->diameter() *
              std::numeric_limits<double>::epsilon() * 1024.);

            // we can only compare if current and previous cell have at
            // least one quadrature point and both cells are at least affine
            const auto comparison_result_m =
              (!affine_cells || mapping_data[0].inverse_jacobians.empty() ||
               mapping_data_previous_cell[0].inverse_jacobians.empty()) ?
                FloatingPointComparator<double>::ComparisonResult::less :
                comparator.compare(
                  mapping_data[0].inverse_jacobians[0],
                  mapping_data_previous_cell[0].inverse_jacobians[0]);

            const auto comparison_result_p =
              (!affine_cells || mapping_data[1].inverse_jacobians.empty() ||
               mapping_data_previous_cell[1].inverse_jacobians.empty()) ?
                FloatingPointComparator<double>::ComparisonResult::less :
                comparator.compare(
                  mapping_data[1].inverse_jacobians[0],
                  mapping_data_previous_cell[1].inverse_jacobians[0]);

            // we can compress the Jacobians and inverse Jacobians if
            // inverse Jacobians are equal and cells are affine
            if (affine_cells &&
                comparison_result_m ==
                  FloatingPointComparator<double>::ComparisonResult::equal &&
                comparison_result_p ==
                  FloatingPointComparator<double>::ComparisonResult::equal)
              {
                compressed_data_index_offsets.push_back(
                  compressed_data_index_offsets.back());
              }
            else if (first_set &&
                     (cell_type[face_index] <=
                      dealii::internal::MatrixFreeFunctions::affine) &&
                     (comparator.compare(
                        mapping_data[0].inverse_jacobians[0],
                        mapping_data_first[0].inverse_jacobians[0]) ==
                      FloatingPointComparator<
                        double>::ComparisonResult::equal) &&
                     (comparator.compare(
                        mapping_data[1].inverse_jacobians[0],
                        mapping_data_first[1].inverse_jacobians[0]) ==
                      FloatingPointComparator<double>::ComparisonResult::equal))
              {
                compressed_data_index_offsets.push_back(0);
              }
            else
              {
                const unsigned int n_compressed_data_last_cell =
                  cell_type[face_index - 1] <=
                      dealii::internal::MatrixFreeFunctions::affine ?
                    1 :
                    compute_n_q_points<Number>(
                      n_q_points_unvectorized[face_index - 1]);

                compressed_data_index_offsets.push_back(
                  compressed_data_index_offsets.back() +
                  n_compressed_data_last_cell);
              }
          }
        else
          compressed_data_index_offsets.push_back(0);

        // cache mapping_data from previous cell pair
        mapping_data_previous_cell = mapping_data;

        const unsigned int n_q_points_data =
          compute_n_q_points<Number>(n_q_points_unvectorized[face_index]);

        // store mapping data of interior face
        store_mapping_data(data_index_offsets[face_index],
                           n_q_points_data,
                           n_q_points_unvectorized[face_index],
                           mapping_data[0],
                           quadrature_on_face.get_weights(),
                           data_index_offsets[face_index],
                           cell_type[face_index] <=
                             dealii::internal::MatrixFreeFunctions::affine,
                           true);

        // store only necessary mapping data for exterior face (Jacobians and
        // inverse Jacobians)
        store_mapping_data(data_index_offsets[face_index],
                           n_q_points_data,
                           n_q_points_unvectorized[face_index],
                           mapping_data[1],
                           quadrature_on_face.get_weights(),
                           data_index_offsets[face_index],
                           cell_type[face_index] <=
                             dealii::internal::MatrixFreeFunctions::affine,
                           false);

        // update size of compressed data depending on cell type and handle
        // empty quadratures
        if (cell_type[face_index] <=
            dealii::internal::MatrixFreeFunctions::affine)
          size_compressed_data = compressed_data_index_offsets.back() + 1;
        else
          size_compressed_data =
            std::max(size_compressed_data,
                     compressed_data_index_offsets.back() + n_q_points_data);

        ++face_index;
      }

    if (update_flags_mapping & UpdateFlags::update_jacobians)
      {
        jacobians[0].resize(size_compressed_data);
        jacobians[0].shrink_to_fit();
        jacobians[1].resize(size_compressed_data);
        jacobians[1].shrink_to_fit();
      }
    if (update_flags_mapping & UpdateFlags::update_inverse_jacobians)
      {
        inverse_jacobians[0].resize(size_compressed_data);
        inverse_jacobians[0].shrink_to_fit();
        inverse_jacobians[1].resize(size_compressed_data);
        inverse_jacobians[1].shrink_to_fit();
      }

    state = State::face_vector;
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
    const unsigned int geometry_index) const
  {
    return n_q_points_unvectorized[geometry_index];
  }


  template <int dim, int spacedim, typename Number>
  dealii::internal::MatrixFreeFunctions::GeometryType
  MappingInfo<dim, spacedim, Number>::get_cell_type(
    const unsigned int geometry_index) const
  {
    return cell_type[geometry_index];
  }



  template <int dim, int spacedim, typename Number>
  typename Triangulation<dim, spacedim>::cell_iterator
  MappingInfo<dim, spacedim, Number>::get_cell_iterator(
    const unsigned int cell_index) const
  {
    Assert(
      additional_data.store_cells,
      ExcMessage(
        "Cells have been not stored. You can enable this by Additional::store_cells."));
    return {triangulation.get(),
            cell_level_and_indices[cell_index].first,
            cell_level_and_indices[cell_index].second};
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
  template <bool is_face>
  unsigned int
  MappingInfo<dim, spacedim, Number>::compute_geometry_index_offset(
    const unsigned int cell_index,
    const unsigned int face_number) const
  {
    if (cell_index == numbers::invalid_unsigned_int)
      return 0;

    const unsigned int compressed_cell_index =
      compute_compressed_cell_index(cell_index);
    if (!is_face)
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
        Assert(state == State::faces_on_cells_in_vector ||
                 state == State::face_vector,
               ExcMessage("This mapping info is not reinitialized for faces"
                          " on cells in a vector!"));
        if (state == State::faces_on_cells_in_vector)
          return cell_index_offset[compressed_cell_index] + face_number;
        else if (state == State::face_vector)
          return cell_index;
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
    const std::vector<double>      &weights,
    const unsigned int              compressed_unit_point_index_offset,
    const bool                      affine_cell,
    const bool                      is_interior)
  {
    const unsigned int n_lanes =
      dealii::internal::VectorizedArrayTrait<Number>::width();

    for (unsigned int q = 0; q < n_q_points; ++q)
      {
        const unsigned int offset = unit_points_index_offset + q;
        const unsigned int compressed_offset =
          compressed_unit_point_index_offset + q;
        for (unsigned int v = 0;
             v < n_lanes && q * n_lanes + v < n_q_points_unvectorized;
             ++v)
          {
            if (q == 0 || !affine_cell)
              {
                if (update_flags_mapping & UpdateFlags::update_jacobians)
                  for (unsigned int d = 0; d < dim; ++d)
                    for (unsigned int s = 0; s < spacedim; ++s)
                      dealii::internal::VectorizedArrayTrait<Number>::get(
                        jacobians[is_interior ? 0 : 1][compressed_offset][d][s],
                        v) = mapping_data.jacobians[q * n_lanes + v][d][s];
                if (update_flags_mapping &
                    UpdateFlags::update_inverse_jacobians)
                  for (unsigned int d = 0; d < dim; ++d)
                    for (unsigned int s = 0; s < spacedim; ++s)
                      dealii::internal::VectorizedArrayTrait<Number>::get(
                        inverse_jacobians[is_interior ? 0 : 1]
                                         [compressed_offset][d][s],
                        v) =
                        mapping_data.inverse_jacobians[q * n_lanes + v][d][s];
              }

            if (is_interior)
              {
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
                if (update_flags_mapping &
                    UpdateFlags::update_quadrature_points)
                  for (unsigned int s = 0; s < spacedim; ++s)
                    dealii::internal::VectorizedArrayTrait<Number>::get(
                      real_points[offset][s], v) =
                      mapping_data.quadrature_points[q * n_lanes + v][s];
              }
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
    const unsigned int n_data_point_batches,
    const bool         is_face_centric)
  {
    if (update_flags_mapping & UpdateFlags::update_jacobians)
      {
        jacobians[0].resize(n_data_point_batches);
        if (is_face_centric)
          jacobians[1].resize(n_data_point_batches);
      }
    if (update_flags_mapping & UpdateFlags::update_inverse_jacobians)
      {
        inverse_jacobians[0].resize(n_data_point_batches);
        if (is_face_centric)
          inverse_jacobians[1].resize(n_data_point_batches);
      }
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
    return unit_points.data() + offset;
  }



  template <int dim, int spacedim, typename Number>
  inline const Point<
    dim - 1,
    typename MappingInfo<dim, spacedim, Number>::VectorizedArrayType> *
  MappingInfo<dim, spacedim, Number>::get_unit_point_faces(
    const unsigned int offset) const
  {
    return unit_points_faces.data() + offset;
  }



  template <int dim, int spacedim, typename Number>
  inline const Point<spacedim, Number> *
  MappingInfo<dim, spacedim, Number>::get_real_point(
    const unsigned int offset) const
  {
    return real_points.data() + offset;
  }



  template <int dim, int spacedim, typename Number>
  unsigned int
  MappingInfo<dim, spacedim, Number>::compute_unit_point_index_offset(
    const unsigned int geometry_index) const
  {
    return unit_points_index[geometry_index];
  }



  template <int dim, int spacedim, typename Number>
  unsigned int
  MappingInfo<dim, spacedim, Number>::compute_data_index_offset(
    const unsigned int geometry_index) const
  {
    return data_index_offsets[geometry_index];
  }


  template <int dim, int spacedim, typename Number>
  unsigned int
  MappingInfo<dim, spacedim, Number>::compute_compressed_data_index_offset(
    const unsigned int geometry_index) const
  {
    return compressed_data_index_offsets[geometry_index];
  }



  template <int dim, int spacedim, typename Number>
  inline const DerivativeForm<1, dim, spacedim, Number> *
  MappingInfo<dim, spacedim, Number>::get_jacobian(const unsigned int offset,
                                                   const bool is_interior) const
  {
    return jacobians[is_interior ? 0 : 1].data() + offset;
  }



  template <int dim, int spacedim, typename Number>
  inline const DerivativeForm<1, spacedim, dim, Number> *
  MappingInfo<dim, spacedim, Number>::get_inverse_jacobian(
    const unsigned int offset,
    const bool         is_interior) const
  {
    return inverse_jacobians[is_interior ? 0 : 1].data() + offset;
  }



  template <int dim, int spacedim, typename Number>
  inline const Tensor<1, spacedim, Number> *
  MappingInfo<dim, spacedim, Number>::get_normal_vector(
    const unsigned int offset) const
  {
    return normal_vectors.data() + offset;
  }



  template <int dim, int spacedim, typename Number>
  inline const Number *
  MappingInfo<dim, spacedim, Number>::get_JxW(const unsigned int offset) const
  {
    return JxW_values.data() + offset;
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
  std::size_t
  MappingInfo<dim, spacedim, Number>::memory_consumption() const
  {
    std::size_t memory = MemoryConsumption::memory_consumption(unit_points);
    memory += MemoryConsumption::memory_consumption(unit_points_faces);
    memory += MemoryConsumption::memory_consumption(unit_points_index);
    memory += cell_type.capacity() *
              sizeof(dealii::internal::MatrixFreeFunctions::GeometryType);
    memory += MemoryConsumption::memory_consumption(data_index_offsets);
    memory +=
      MemoryConsumption::memory_consumption(compressed_data_index_offsets);
    memory += MemoryConsumption::memory_consumption(JxW_values);
    memory += MemoryConsumption::memory_consumption(normal_vectors);
    memory += MemoryConsumption::memory_consumption(jacobians);
    memory += MemoryConsumption::memory_consumption(inverse_jacobians);
    memory += MemoryConsumption::memory_consumption(real_points);
    memory += MemoryConsumption::memory_consumption(n_q_points_unvectorized);
    memory += MemoryConsumption::memory_consumption(cell_index_offset);
    memory += MemoryConsumption::memory_consumption(
      cell_index_to_compressed_cell_index);
    memory += MemoryConsumption::memory_consumption(cell_level_and_indices);
    memory += sizeof(*this);
    return memory;
  }
} // namespace NonMatching

DEAL_II_NAMESPACE_CLOSE

#endif

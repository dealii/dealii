// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2020 by the deal.II authors
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


#ifndef dealii_matrix_free_mapping_info_h
#define dealii_matrix_free_mapping_info_h


#include <deal.II/base/config.h>

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/hp/q_collection.h>

#include <deal.II/matrix_free/face_info.h>
#include <deal.II/matrix_free/helper_functions.h>

#include <memory>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace MatrixFreeFunctions
  {
    /**
     * An enum to identify various types of cells and faces. The most general
     * type is what we typically compute in the FEValues context but for many
     * geometries we can save significant storage.
     *
     * @ingroup matrixfree
     */
    enum GeometryType : unsigned char
    {
      /**
       * The cell or face is Cartesian.
       */
      cartesian = 0,

      /**
       * The cell or face can be described with an affine mapping.
       */
      affine = 1,

      /**
       * The face is flat, i.e., the normal factor on a face is the same on
       * all quadrature points. This type is not assigned for cells.
       */
      flat_faces = 2,

      /**
       * There is no special information available for compressing the
       * representation of the object under consideration.
       */
      general = 3
    };



    /**
     * Definition of a structure that stores all cached data related to the
     * evaluated geometry from the mapping. In order to support hp-adaptivity
     * and compressed storage (in particular for Jacobians, JxW values, and
     * normals), storage length can be different for different rows. Thus, it
     * allows to jump at the data of individual rows similar to compressed row
     * storage in sparse matrices. We have two different start indices for
     * fields with different sizes. The first category of offsets are the
     * indices for Jacobians of the transformation from unit to real cell (we
     * store the inverse Jacobian), second derivatives, JxW values, and normal
     * vectors. We keep separate arrays for all these data structures because
     * a user code might access only some of them. In such a case, one array
     * will be gone through in a contiguous order with access to all entries,
     * which makes it easy for the processor to prefetch data. Having all data
     * in a single array would require some strides in the access pattern,
     * which is much more complicated for the processor to predict (and indeed
     * leads to prefetching of data that does not get used on Intel processors
     * such as BroadwellEP).
     *
     * The second category of indices are the offsets for the quadrature
     * points. Quadrature points can be compressed less than the other fields
     * and thus need longer fields. Quadrature point indices are often used in
     * other contexts such as evaluation of right hand sides.
     *
     * The third component is a descriptor of data from the unit cells, called
     * QuadratureDescriptor, which contains the quadrature weights and
     * permutations of how to go through quadrature points in case of face
     * data. The latter comes in a vector for the support of hp adaptivity,
     * with several data fields for the individual quadrature formulas.
     *
     * @ingroup matrixfree
     *
     * @author Katharina Kormann, Martin Kronbichler, 2018
     */
    template <int structdim,
              int spacedim,
              typename Number,
              typename VectorizedArrayType>
    struct MappingInfoStorage
    {
      static_assert(
        std::is_same<Number, typename VectorizedArrayType::value_type>::value,
        "Type of Number and of VectorizedArrayType do not match.");

      struct QuadratureDescriptor
      {
        /**
         * Constructor. Does nothing.
         */
        QuadratureDescriptor();

        /**
         * Set up the lengths in the various members of this struct.
         */
        void
        initialize(const Quadrature<1> &quadrature_1d,
                   const UpdateFlags update_flags_inner_faces = update_default);

        /**
         * Returns the memory consumption in bytes.
         */
        std::size_t
        memory_consumption() const;

        /**
         * Number of quadrature points applied on the given cell or face.
         */
        unsigned int n_q_points;

        /**
         * Original one-dimensional quadrature formula applied on the given
         * cell or face.
         */
        Quadrature<1> quadrature_1d;

        /**
         * Quadrature formula applied on the given cell or face.
         */
        Quadrature<structdim> quadrature;

        /**
         * Quadrature weights separated by dimension for use in specific
         * situations.
         */
        std::array<AlignedVector<Number>, structdim> tensor_quadrature_weights;

        /**
         * A cached vector of quadrature weights in the given number format
         * (non-vectorized, as it is cheap to broadcast the value to all lanes
         * when it is used in a vectorized context).
         */
        AlignedVector<Number> quadrature_weights;

        /**
         * For quadrature on faces, the evaluation of basis functions is not
         * in the correct order if a face is not in the standard orientation
         * to a given element. This data structure is used to re-order the
         * data evaluated on quadrature points to represent the correct order.
         */
        dealii::Table<2, unsigned int> face_orientations;
      };

      /**
       * A class describing the layout of the sections in the @p data_storage
       * field and also includes some data that depends on the number of
       * quadrature points in the hp context such as the inner quadrature
       * formula and re-indexing for faces that are not in the standard
       * orientation.
       */
      std::vector<QuadratureDescriptor> descriptor;

      /**
       * Stores the index offset into the arrays @p jxw_values, @p jacobians,
       * @p normal_vectors and the second derivatives. Note that affine cells
       * have shorter fields of length 1, where the others have lengths equal
       * to the number of quadrature points of the given cell.
       */
      AlignedVector<unsigned int> data_index_offsets;

      /**
       * The storage of the Jacobian determinant (times the quadrature weight
       * in case the transformation is non-affine) on quadrature
       * points.
       *
       * Indexed by @p data_index_offsets.
       */
      AlignedVector<VectorizedArrayType> JxW_values;

      /**
       * Stores the normal vectors.
       *
       * Indexed by @p data_index_offsets.
       */
      AlignedVector<Tensor<1, spacedim, VectorizedArrayType>> normal_vectors;

      /**
       * The storage of covariant transformation on quadrature points, i.e.,
       * the inverse and transposed Jacobians of the transformation from the
       * unit to the real cell.
       *
       * Indexed by @p data_index_offsets.
       *
       * Contains two fields for access from both sides for interior faces,
       * but the default case (cell integrals or boundary integrals) only
       * fills the zeroth component and ignores the first one.
       */
      AlignedVector<Tensor<2, spacedim, VectorizedArrayType>> jacobians[2];

      /**
       * The storage of the gradients of the inverse Jacobian
       * transformation. Because of symmetry, only the upper diagonal and
       * diagonal part are needed. The first index runs through the
       * derivatives, starting with the diagonal and then continuing row-wise,
       * i.e., $\partial^2/\partial x_1 \partial x_2$ first, then
       * $\partial^2/\partial x_1 \partial x_3$, and so on. The second index
       * is the spatial coordinate.
       *
       * Indexed by @p data_index_offsets.
       *
       * Contains two fields for access from both sides for interior faces,
       * but the default case (cell integrals or boundary integrals) only
       * fills the zeroth component and ignores the first one.
       */
      AlignedVector<Tensor<1,
                           spacedim *(spacedim + 1) / 2,
                           Tensor<1, spacedim, VectorizedArrayType>>>
        jacobian_gradients[2];

      /**
       * Stores the Jacobian transformations times the normal vector (this
       * represents a shortcut that is accessed often and can thus get higher
       * performance).
       *
       * Indexed by @p data_index_offsets.
       */
      AlignedVector<Tensor<1, spacedim, VectorizedArrayType>>
        normals_times_jacobians[2];

      /**
       * Stores the index offset of a particular cell into the quadrature
       * points array in real coordinates. Note that Cartesian cells have
       * shorter fields (length is @p n_q_points_1d) than non-Cartesian cells
       * (length is @p n_q_points) or faces.
       */
      AlignedVector<unsigned int> quadrature_point_offsets;

      /**
       * Stores the quadrature points in real coordinates, including a
       * compression scheme for Cartesian cells where we do not need to store
       * the full data on all points.
       *
       * Indexed by @p quadrature_point_offsets.
       */
      AlignedVector<Point<spacedim, VectorizedArrayType>> quadrature_points;

      /**
       * Clears all data fields except the descriptor vector.
       */
      void
      clear_data_fields();

      /**
       * Returns the quadrature index for a given number of quadrature
       * points. If not in hp mode or if the index is not found, this
       * function always returns index 0. Hence, this function does not
       * check whether the given degree is actually present.
       */
      unsigned int
      quad_index_from_n_q_points(const unsigned int n_q_points) const;

      /**
       * Prints a detailed summary of memory consumption in the different
       * structures of this class to the given output stream.
       */
      template <typename StreamType>
      void
      print_memory_consumption(StreamType &    out,
                               const TaskInfo &task_info) const;

      /**
       * Returns the memory consumption in bytes.
       */
      std::size_t
      memory_consumption() const;
    };



    /**
     * The class that stores all geometry-dependent data related with cell
     * interiors for use in the matrix-free class.
     *
     * @ingroup matrixfree
     *
     * @author Katharina Kormann and Martin Kronbichler, 2010, 2011, 2017
     */
    template <int dim, typename Number, typename VectorizedArrayType>
    struct MappingInfo
    {
      /**
       * Compute the information in the given cells and faces. The cells are
       * specified by the level and the index within the level (as given by
       * CellIterator::level() and CellIterator::index(), in order to allow
       * for different kinds of iterators, e.g. standard DoFHandler,
       * multigrid, etc.)  on a fixed Triangulation. In addition, a mapping
       * and several quadrature formulas are given.
       */
      void
      initialize(
        const dealii::Triangulation<dim> &                        tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cells,
        const FaceInfo<VectorizedArrayType::size()> &             faces,
        const std::vector<unsigned int> &              active_fe_index,
        const Mapping<dim> &                           mapping,
        const std::vector<dealii::hp::QCollection<1>> &quad,
        const UpdateFlags                              update_flags_cells,
        const UpdateFlags update_flags_boundary_faces,
        const UpdateFlags update_flags_inner_faces,
        const UpdateFlags update_flags_faces_by_cells);

      /**
       * Update the information in the given cells and faces that is the
       * result of a change in the given `mapping` class, keeping the cells,
       * quadrature formulas and other unknowns unchanged. This call is only
       * valid if MappingInfo::initialize() has been called before.
       */
      void
      update_mapping(
        const dealii::Triangulation<dim> &                        tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cells,
        const FaceInfo<VectorizedArrayType::size()> &             faces,
        const std::vector<unsigned int> &active_fe_index,
        const Mapping<dim> &             mapping);

      /**
       * Return the type of a given cell as detected during initialization.
       */
      GeometryType
      get_cell_type(const unsigned int cell_chunk_no) const;

      /**
       * Clear all data fields in this class.
       */
      void
      clear();

      /**
       * Return the memory consumption of this class in bytes.
       */
      std::size_t
      memory_consumption() const;

      /**
       * Prints a detailed summary of memory consumption in the different
       * structures of this class to the given output stream.
       */
      template <typename StreamType>
      void
      print_memory_consumption(StreamType &    out,
                               const TaskInfo &task_info) const;

      /**
       * The given update flags for computing the geometry on the cells.
       */
      UpdateFlags update_flags_cells;

      /**
       * The given update flags for computing the geometry on the boundary
       * faces.
       */
      UpdateFlags update_flags_boundary_faces;

      /**
       * The given update flags for computing the geometry on the interior
       * faces.
       */
      UpdateFlags update_flags_inner_faces;

      /**
       * The given update flags for computing the geometry on the faces for
       * cell-centric loops.
       */
      UpdateFlags update_flags_faces_by_cells;

      /**
       * Stores whether a cell is Cartesian (cell type 0), has constant
       * transform data (Jacobians) (cell type 1), or is general (cell type
       * 3). Type 2 is only used for faces and no cells are assigned this
       * value.
       */
      std::vector<GeometryType> cell_type;

      /**
       * Stores whether a face (and both cells adjacent to the face) is
       * Cartesian (face type 0), whether it represents an affine situation
       * (face type 1), whether it is a flat face where the normal vector is
       * the same throughout the face (face type 2), or is general (face type
       * 3).
       */
      std::vector<GeometryType> face_type;

      /**
       * The data cache for the cells.
       */
      std::vector<MappingInfoStorage<dim, dim, Number, VectorizedArrayType>>
        cell_data;

      /**
       * The data cache for the faces.
       */
      std::vector<MappingInfoStorage<dim - 1, dim, Number, VectorizedArrayType>>
        face_data;

      /**
       * The data cache for the face-associated-with-cell topology, following
       * the @p cell_type variable for the cell types.
       */
      std::vector<MappingInfoStorage<dim - 1, dim, Number, VectorizedArrayType>>
        face_data_by_cells;

      /**
       * The pointer to the underlying Mapping object.
       */
      SmartPointer<const Mapping<dim>> mapping;

      /**
       * Internal function to compute the geometry for the case the mapping is
       * a MappingQ and a single quadrature formula per slot (non-hp case) is
       * used. This method computes all data from the underlying cell
       * quadrature points using the fast operator evaluation techniques from
       * the matrix-free framework itself, i.e., it uses a polynomial
       * description of the cell geometry (that is computed in a first step)
       * and then computes all Jacobians and normal vectors based on this
       * information. This optimized approach is much faster than going
       * through FEValues and FEFaceValues, especially when several different
       * quadrature formulas are involved, and consumes less memory.
       *
       * @param tria The triangulation to be used for setup
       *
       * @param cells The actual cells of the triangulation to be worked on,
       * given as a tuple of the level and index within the level as used in
       * the main initialization of the class
       *
       * @param faces The description of the connectivity from faces to cells
       * as filled in the MatrixFree class
       */
      void
      compute_mapping_q(
        const dealii::Triangulation<dim> &                        tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cells,
        const std::vector<FaceToCellTopology<VectorizedArrayType::size()>>
          &faces);

      /**
       * Computes the information in the given cells, called within
       * initialize.
       */
      void
      initialize_cells(
        const dealii::Triangulation<dim> &                        tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cells,
        const std::vector<unsigned int> &active_fe_index,
        const Mapping<dim> &             mapping);

      /**
       * Computes the information in the given faces, called within
       * initialize.
       */
      void
      initialize_faces(
        const dealii::Triangulation<dim> &                        tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cells,
        const std::vector<FaceToCellTopology<VectorizedArrayType::size()>>
          &                 faces,
        const Mapping<dim> &mapping);

      /**
       * Computes the information in the given faces, called within
       * initialize.
       */
      void
      initialize_faces_by_cells(
        const dealii::Triangulation<dim> &                        tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cells,
        const Mapping<dim> &                                      mapping);

      /**
       * Helper function to determine which update flags must be set in the
       * internal functions to initialize all data as requested by the user.
       */
      static UpdateFlags
      compute_update_flags(const UpdateFlags update_flags,
                           const std::vector<dealii::hp::QCollection<1>> &quad =
                             std::vector<dealii::hp::QCollection<1>>());
    };



    /**
     * A helper class to extract either cell or face data from mapping info
     * for use in FEEvaluationBase.
     *
     * @author Katharina Kormann, Martin Kronbichler, 2018
     */
    template <int, typename, bool, typename>
    struct MappingInfoCellsOrFaces;

    template <int dim, typename Number, typename VectorizedArrayType>
    struct MappingInfoCellsOrFaces<dim, Number, false, VectorizedArrayType>
    {
      static const MappingInfoStorage<dim, dim, Number, VectorizedArrayType> *
      get(const MappingInfo<dim, Number, VectorizedArrayType> &mapping_info,
          const unsigned int                                   quad_no)
      {
        AssertIndexRange(quad_no, mapping_info.cell_data.size());
        return &mapping_info.cell_data[quad_no];
      }
    };

    template <int dim, typename Number, typename VectorizedArrayType>
    struct MappingInfoCellsOrFaces<dim, Number, true, VectorizedArrayType>
    {
      static const MappingInfoStorage<dim - 1, dim, Number, VectorizedArrayType>
        *
        get(const MappingInfo<dim, Number, VectorizedArrayType> &mapping_info,
            const unsigned int                                   quad_no)
      {
        AssertIndexRange(quad_no, mapping_info.face_data.size());
        return &mapping_info.face_data[quad_no];
      }
    };



    /**
     * A class that is used to compare floating point arrays (e.g. std::vectors,
     * Tensor<1,dim>, etc.). The idea of this class is to consider two arrays as
     * equal if they are the same within a given tolerance. We use this
     * comparator class within a std::map<> of the given arrays. Note that this
     * comparison operator does not satisfy all the mathematical properties one
     * usually wants to have (consider e.g. the numbers a=0, b=0.1, c=0.2 with
     * tolerance 0.15; the operator gives a<c, but neither a<b? nor b<c? is
     * satisfied). This is not a problem in the use cases for this class, but be
     * careful when using it in other contexts.
     */
    template <typename Number,
              typename VectorizedArrayType = VectorizedArray<Number>>
    struct FPArrayComparator
    {
      FPArrayComparator(const Number scaling);

      /**
       * Compare two vectors of numbers (not necessarily of the same length)
       */
      bool
      operator()(const std::vector<Number> &v1,
                 const std::vector<Number> &v2) const;

      /**
       * Compare two vectorized arrays (stored as tensors to avoid alignment
       * issues).
       */
      bool
      operator()(
        const Tensor<1, VectorizedArrayType::size(), Number> &t1,
        const Tensor<1, VectorizedArrayType::size(), Number> &t2) const;

      /**
       * Compare two rank-1 tensors of vectorized arrays (stored as tensors to
       * avoid alignment issues).
       */
      template <int dim>
      bool
      operator()(
        const Tensor<1, dim, Tensor<1, VectorizedArrayType::size(), Number>>
          &t1,
        const Tensor<1, dim, Tensor<1, VectorizedArrayType::size(), Number>>
          &t2) const;

      /**
       * Compare two rank-2 tensors of vectorized arrays (stored as tensors to
       * avoid alignment issues).
       */
      template <int dim>
      bool
      operator()(
        const Tensor<2, dim, Tensor<1, VectorizedArrayType::size(), Number>>
          &t1,
        const Tensor<2, dim, Tensor<1, VectorizedArrayType::size(), Number>>
          &t2) const;

      /**
       * Compare two arrays of tensors.
       */
      template <int dim>
      bool
      operator()(const std::array<Tensor<2, dim, Number>, dim + 1> &t1,
                 const std::array<Tensor<2, dim, Number>, dim + 1> &t2) const;

      Number tolerance;
    };



    /* ------------------- inline functions ----------------------------- */

    template <int structdim,
              int spacedim,
              typename Number,
              typename VectorizedArrayType>
    inline unsigned int
    MappingInfoStorage<structdim, spacedim, Number, VectorizedArrayType>::
      quad_index_from_n_q_points(const unsigned int n_q_points) const
    {
      for (unsigned int i = 0; i < descriptor.size(); ++i)
        if (n_q_points == descriptor[i].n_q_points)
          return i;
      return 0;
    }



    template <int dim, typename Number, typename VectorizedArrayType>
    inline GeometryType
    MappingInfo<dim, Number, VectorizedArrayType>::get_cell_type(
      const unsigned int cell_no) const
    {
      AssertIndexRange(cell_no, cell_type.size());
      return cell_type[cell_no];
    }

  } // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_matrix_free_mapping_info_storage_h
#define dealii_matrix_free_mapping_info_storage_h


#include <deal.II/base/config.h>

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/fe/fe_update_flags.h>

#include <deal.II/hp/q_collection.h>

#include <memory>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace MatrixFreeFunctions
  {
    // forward declaration
    struct TaskInfo;

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
     * data. The latter comes in a vector for the support of hp-adaptivity,
     * with several data fields for the individual quadrature formulas.
     *
     * @ingroup matrixfree
     */
    template <int structdim, int spacedim, typename Number>
    struct MappingInfoStorage
    {
      struct QuadratureDescriptor
      {
        /**
         * In case this class is instantiated for VectorizedArray types, this
         * indicates the underlying scalar type for data which is the same on
         * all lanes like the quadrature weights.
         */
        using ScalarNumber = typename VectorizedArrayTrait<Number>::value_type;

        /**
         * Constructor. Does nothing.
         */
        QuadratureDescriptor();

        /**
         * Set up the lengths in the various members of this struct.
         */
        template <int dim_q>
        void
        initialize(const Quadrature<dim_q> &quadrature);

        /**
         * Set up the lengths in the various members of this struct.
         */
        void
        initialize(const Quadrature<1> &quadrature_1d);

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
        std::array<AlignedVector<ScalarNumber>, structdim>
          tensor_quadrature_weights;

        /**
         * A cached vector of quadrature weights in the given number format
         * (non-vectorized, as it is cheap to broadcast the value to all lanes
         * when it is used in a vectorized context).
         */
        AlignedVector<ScalarNumber> quadrature_weights;
      };

      /**
       * A class describing the layout of the sections in the @p data_storage
       * field and also includes some data that depends on the number of
       * quadrature points in the hp-context such as the inner quadrature
       * formula and re-indexing for faces that are not in the standard
       * orientation.
       */
      std::vector<QuadratureDescriptor> descriptor;

      /**
       * Collection of quadrature formulae applied on the given face.
       *
       * @note Only filled for faces, since faces might be quadrilateral or
       *   triangle shaped.
       */
      std::vector<dealii::hp::QCollection<structdim>> q_collection;

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
      AlignedVector<Number> JxW_values;

      /**
       * Stores the normal vectors.
       *
       * Indexed by @p data_index_offsets.
       */
      AlignedVector<Tensor<1, spacedim, Number>> normal_vectors;

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
       *
       * If the cell is Cartesian/affine then the Jacobian is stored at index 1
       * of the AlignedVector. For faces on hypercube elements, the derivatives
       * are reorder s.t the derivative orthogonal to the face is stored last,
       * i.e for dim = 3 and face_no = 0 or 1, the derivatives are ordered as
       * [dy, dz, dx], face_no = 2 or 3: [dz, dx, dy], and face_no = 5 or 6:
       * [dx, dy, dz]. If the Jacobian also is stored, the components are
       * instead reordered in the same way.
       */
      std::array<AlignedVector<Tensor<2, spacedim, Number>>, 2> jacobians;

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
      std::array<
        AlignedVector<
          Tensor<1, spacedim *(spacedim + 1) / 2, Tensor<1, spacedim, Number>>>,
        2>
        jacobian_gradients;

      /**
       * The storage of the gradients of the Jacobian transformation. Because of
       * symmetry, only the upper diagonal and diagonal part are needed. The
       * first index runs through the derivatives, starting with the diagonal
       * and then continuing row-wise, i.e., $\partial^2/\partial x_1 \partial
       * x_2$ first, then
       * $\partial^2/\partial x_1 \partial x_3$, and so on. The second index
       * is the spatial coordinate.
       *
       * Indexed by @p data_index_offsets.
       *
       * Contains two fields for access from both sides for interior faces,
       * but the default case (cell integrals or boundary integrals) only
       * fills the zeroth component and ignores the first one.
       */
      std::array<
        AlignedVector<
          Tensor<1, spacedim *(spacedim + 1) / 2, Tensor<1, spacedim, Number>>>,
        2>
        jacobian_gradients_non_inverse;

      /**
       * Stores the Jacobian transformations times the normal vector (this
       * represents a shortcut that is accessed often and can thus get higher
       * performance).
       *
       * Indexed by @p data_index_offsets.
       */
      std::array<AlignedVector<Tensor<1, spacedim, Number>>, 2>
        normals_times_jacobians;

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
      AlignedVector<Point<spacedim, Number>> quadrature_points;

      /**
       * Clears all data fields except the descriptor vector.
       */
      void
      clear_data_fields();

      /**
       * Returns the quadrature index for a given number of quadrature
       * points. If not in hp-mode or if the index is not found, this
       * function always returns index 0. Hence, this function does not
       * check whether the given degree is actually present.
       */
      unsigned int
      quad_index_from_n_q_points(const unsigned int n_q_points) const;

      /**
       * Helper function to determine which update flags must be set in the
       * internal functions to initialize all data as requested by the user.
       */
      static UpdateFlags
      compute_update_flags(
        const UpdateFlags                                     update_flags,
        const std::vector<dealii::hp::QCollection<spacedim>> &quads =
          std::vector<dealii::hp::QCollection<spacedim>>(),
        const bool piola_transform = false);

      /**
       * Prints a detailed summary of memory consumption in the different
       * structures of this class to the given output stream.
       */
      template <typename StreamType>
      void
      print_memory_consumption(StreamType     &out,
                               const TaskInfo &task_info) const;

      /**
       * Returns the memory consumption in bytes.
       */
      std::size_t
      memory_consumption() const;
    };



    /* ------------------- inline functions ----------------------------- */

    template <int structdim, int spacedim, typename Number>
    inline unsigned int
    MappingInfoStorage<structdim, spacedim, Number>::quad_index_from_n_q_points(
      const unsigned int n_q_points) const
    {
      for (unsigned int i = 0; i < descriptor.size(); ++i)
        if (n_q_points == descriptor[i].n_q_points)
          return i;
      return 0;
    }

  } // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif

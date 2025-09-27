// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_matrix_free_mapping_info_h
#define dealii_matrix_free_mapping_info_h


#include <deal.II/base/config.h>

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/reference_cell.h>

#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/matrix_free/face_info.h>
#include <deal.II/matrix_free/mapping_info_storage.h>

#include <memory>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace MatrixFreeFunctions
  {
    /**
     * The class that stores all geometry-dependent data related with cell
     * interiors for use in the matrix-free class.
     *
     * @ingroup matrixfree
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
       * and several 1d quadrature formulas are given.
       */
      void
      initialize(
        const dealii::Triangulation<dim>                         &tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cells,
        const FaceInfo<VectorizedArrayType::size()>              &face_info,
        const std::vector<unsigned int> &active_fe_index,
        const std::shared_ptr<dealii::hp::MappingCollection<dim>> &mapping,
        const std::vector<dealii::hp::QCollection<dim>>           &quad,
        const UpdateFlags update_flags_cells,
        const UpdateFlags update_flags_boundary_faces,
        const UpdateFlags update_flags_inner_faces,
        const UpdateFlags update_flags_faces_by_cells,
        const bool        piola_transform);

      /**
       * Update the information in the given cells and faces that is the
       * result of a change in the given `mapping` class, keeping the cells,
       * quadrature formulas and other unknowns unchanged. This call is only
       * valid if MappingInfo::initialize() has been called before.
       */
      void
      update_mapping(
        const dealii::Triangulation<dim>                         &tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cells,
        const FaceInfo<VectorizedArrayType::size()>              &face_info,
        const std::vector<unsigned int> &active_fe_index,
        const std::shared_ptr<dealii::hp::MappingCollection<dim>> &mapping);

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
      print_memory_consumption(StreamType     &out,
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
       * Store whether the face geometry for the face_data_by_cells data
       * structure is represented as Cartesian (cell type 0), with constant
       * transform data (Jacobians) (cell type 1), or a general type (cell
       * type 3). Note that both the interior and exterior agree on the type
       * of the data structure, using the more general of the two.
       */
      std::vector<std::array<GeometryType, GeometryInfo<dim>::faces_per_cell>>
        faces_by_cells_type;

      /**
       * The data cache for the cells.
       */
      std::vector<MappingInfoStorage<dim, dim, VectorizedArrayType>> cell_data;

      /**
       * The data cache for the faces.
       */
      std::vector<MappingInfoStorage<dim - 1, dim, VectorizedArrayType>>
        face_data;

      /**
       * The data cache for the face-associated-with-cell topology, following
       * the @p cell_type variable for the cell types.
       */
      std::vector<MappingInfoStorage<dim - 1, dim, VectorizedArrayType>>
        face_data_by_cells;

      /**
       * The pointer to the underlying hp::MappingCollection object.
       */
      std::shared_ptr<dealii::hp::MappingCollection<dim>> mapping_collection;

      /**
       * The pointer to the first entry of mapping_collection.
       */
      ObserverPointer<const Mapping<dim>> mapping;

      /**
       * Reference-cell type related to each quadrature and active quadrature
       * index.
       */
      std::vector<std::vector<ReferenceCell>> reference_cell_types;

      /**
       * Internal function to compute the geometry for the case the mapping is
       * a MappingQ and a single quadrature formula per slot (non-hp-case) is
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
       * @param face_info The description of the connectivity from faces to
       * cells as filled in the MatrixFree class
       */
      void
      compute_mapping_q(
        const dealii::Triangulation<dim>                         &tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cells,
        const FaceInfo<VectorizedArrayType::size()>              &face_info);

      /**
       * Computes the information in the given cells, called within
       * initialize.
       */
      void
      initialize_cells(
        const dealii::Triangulation<dim>                         &tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cells,
        const std::vector<unsigned int>          &active_fe_index,
        const dealii::hp::MappingCollection<dim> &mapping);

      /**
       * Computes the information in the given faces, called within
       * initialize.
       */
      void
      initialize_faces(
        const dealii::Triangulation<dim>                         &tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cells,
        const std::vector<FaceToCellTopology<VectorizedArrayType::size()>>
                                                 &faces,
        const std::vector<unsigned int>          &active_fe_index,
        const dealii::hp::MappingCollection<dim> &mapping);

      /**
       * Computes the information in the given faces, called within
       * initialize.
       */
      void
      initialize_faces_by_cells(
        const dealii::Triangulation<dim>                         &tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cells,
        const FaceInfo<VectorizedArrayType::size()>              &face_info,
        const dealii::hp::MappingCollection<dim>                 &mapping);
    };



    /**
     * A helper class to extract either cell or face data from mapping info
     * for use in FEEvaluationBase.
     */
    template <int, typename, bool, typename>
    struct MappingInfoCellsOrFaces;

    template <int dim, typename Number, typename VectorizedArrayType>
    struct MappingInfoCellsOrFaces<dim, Number, false, VectorizedArrayType>
    {
      static const MappingInfoStorage<dim, dim, VectorizedArrayType> &
      get(const MappingInfo<dim, Number, VectorizedArrayType> &mapping_info,
          const unsigned int                                   quad_no)
      {
        AssertIndexRange(quad_no, mapping_info.cell_data.size());
        return mapping_info.cell_data[quad_no];
      }
    };

    template <int dim, typename Number, typename VectorizedArrayType>
    struct MappingInfoCellsOrFaces<dim, Number, true, VectorizedArrayType>
    {
      static const MappingInfoStorage<dim - 1, dim, VectorizedArrayType> &
      get(const MappingInfo<dim, Number, VectorizedArrayType> &mapping_info,
          const unsigned int                                   quad_no)
      {
        AssertIndexRange(quad_no, mapping_info.face_data.size());
        return mapping_info.face_data[quad_no];
      }
    };


    /* ------------------- inline functions ----------------------------- */

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

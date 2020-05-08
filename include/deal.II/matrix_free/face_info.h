// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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


#ifndef dealii_matrix_free_face_info_h
#define dealii_matrix_free_face_info_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/table.h>


DEAL_II_NAMESPACE_OPEN



namespace internal
{
  namespace MatrixFreeFunctions
  {
    /**
     * Data type for information about the batches build for vectorization of
     * the face integrals. The setup of the batches for the faces is
     * independent of the cells, and thus, we must store the relation to the
     * cell indexing for accessing the degrees of freedom.
     *
     * Interior faces are stored by the two adjacent cells, which we label as
     * "interior" and "exterior" side of the face. Normal vectors stored in
     * MappingInfo are only stored once and are the outer normals to the cells
     * on the "interior" side, whereas the sign is the opposite for the
     * "exterior" side.
     *
     * This data field is stored as a vector for all faces involved in the
     * computation. In order to avoid gaps in the memory representation, the
     * four 'char' variables are put next to each other which occupies the
     * same size as the unsigned integers on most architectures.
     *
     * @author Katharina Kormann, Martin Kronbichler, 2018
     */
    template <int vectorization_width>
    struct FaceToCellTopology
    {
      /**
       * Indices of the faces in the current face batch as compared to the
       * numbers of the cells on the logical "interior" side of the face which
       * is aligned to the direction of FEEvaluation::get_normal_vector().
       */
      unsigned int cells_interior[vectorization_width];

      /**
       * Indices of the faces in the current face batch as compared to the
       * numbers of the cells on the logical "exterior" side of the face which
       * is aligned to the opposite direction of
       * FEEvaluation::get_normal_vector(). Note that the distinction into
       * interior and exterior faces is purely logical and refers to the
       * direction of the normal only. In the actual discretization of a
       * problem, the discretization typically needs to make sure that interior
       * and exterior sides are treated properly, such as with upwind fluxes.
       *
       * For boundary faces, the numbers are set to
       * `numbers::invalid_unsigned_int`.
       */
      unsigned int cells_exterior[vectorization_width];

      /**
       * Index of the face between 0 and GeometryInfo::faces_per_cell within
       * the cells on the "interior" side of the faces.
       */
      unsigned char interior_face_no;

      /**
       * Index of the face between 0 and GeometryInfo::faces_per_cell within
       * the cells on the "exterior" side of the faces.
       *
       * For a boundary face, this data field stores the boundary id.
       */
      unsigned char exterior_face_no;

      /**
       * For adaptively refined meshes, the cell on the exterior side of the
       * face might be less refined than the interior side. This index
       * indicates the possible subface index on the exterior side.
       */
      unsigned char subface_index;

      /**
       * In 3D, one of the two cells adjacent to a face might use a different
       * orientation (also called as face orientation, face flip and face
       * rotation) than the standard orientation. This variable stores the
       * value (for one of the interior or exterior side) for the present batch
       * of faces.
       */
      unsigned char face_orientation;

      /**
       * Return the memory consumption of the present data structure.
       */
      std::size_t
      memory_consumption() const
      {
        return sizeof(*this);
      }
    };



    /**
     * A data structure that holds the connectivity between the faces and the
     * cells.
     */
    template <int vectorization_width>
    struct FaceInfo
    {
      /**
       * Clear all data fields to be in a state similar to after having
       * called the default constructor.
       */
      void
      clear()
      {
        faces = std::vector<FaceToCellTopology<vectorization_width>>();
        cell_and_face_to_plain_faces.reinit(TableIndices<3>(0, 0, 0));
        cell_and_face_boundary_id.reinit(TableIndices<3>(0, 0, 0));
      }

      /**
       * Return the memory consumption of the present data structure.
       */
      std::size_t
      memory_consumption() const
      {
        return sizeof(faces) +
               cell_and_face_to_plain_faces.memory_consumption() +
               cell_and_face_boundary_id.memory_consumption();
      }

      /**
       * Vectorized storage of interior faces, linking to the two cells in the
       * vectorized cell storage.
       */
      std::vector<FaceToCellTopology<vectorization_width>> faces;

      /**
       * This table translates a triple of the macro cell number, the index of a
       * face within a cell and the index within the cell batch of vectorization
       * into the index within the @p faces array.
       */
      ::dealii::Table<3, unsigned int> cell_and_face_to_plain_faces;

      /**
       * Stores the boundary ids of the faces in vectorized format using the
       * same indexing as the cell_and_face_to_plain_faces data structure
       */
      ::dealii::Table<3, types::boundary_id> cell_and_face_boundary_id;
    };
  } // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif

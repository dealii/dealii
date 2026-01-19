// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_tria_connectivity_h
#define dealii_tria_connectivity_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/ndarray.h>

#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria_description.h>
#include <deal.II/grid/tria_objects_orientations.h>

#include <numeric>
#include <vector>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace TriangulationImplementation
  {
    /**
     * A class that represents a (compressed) array of arrays. It is used to
     * store the cell-to-vertex, cell-to-face, etc. lists.
     */
    struct ArrayOfArrays
    {
      /**
       * Default constructor.
       */
      ArrayOfArrays()
        : offsets{0} {};

      /**
       * Constructor which allows to set the internal fields directly.
       */
      ArrayOfArrays(const std::vector<std::size_t>  &offsets,
                    const std::vector<unsigned int> &columns)
        : offsets(offsets)
        , columns(columns)
      {}

      /**
       * Constructor which allows to set the internal fields directly by
       * moving data from the given arguments.
       */
      ArrayOfArrays(std::vector<std::size_t>  &&offsets,
                    std::vector<unsigned int> &&columns)
        : offsets(std::move(offsets))
        , columns(std::move(columns))
      {}

      /**
       * Return the number of (inner) arrays stored by this object.
       */
      unsigned int
      size() const
      {
        return offsets.size() - 1;
      }

      /**
       * Return a view of the `i`th array stored by this object.
       */
      ArrayView<const unsigned int>
      operator[](const unsigned int i) const
      {
        AssertIndexRange(i, offsets.size() - 1);
        return ArrayView<const unsigned int>(&columns[offsets[i]],
                                             offsets[i + 1] - offsets[i]);
      }

      /**
       * For each row in the CRS, store the corresponding offsets in the
       * 'elements' array.
       */
      std::vector<std::size_t> offsets;

      /**
       * The elements (typically column indices) stored in the CRS.
       */
      std::vector<unsigned int> columns;
    };



    /**
     * Class for storing the reduced connectivity table.
     *
     * A full connectivity table contains all possible connectivities of
     * entities of dimension d and entities of dimension d' with 0<=d,d'<=dim.
     * However, in the library we only need the following types of
     * connectivities:
     *  - dim-dimensional neighbors of dim-dimensional entities (connected via
     *    faces)
     *  - d-dimensional entity to it's (d-1)-dimension bounding entities
     *  - quad (2 - 3D), line (1 - 2d/3d) to vertices (0) to be able to process
     *    the user provided SubCellData during
     *    Triangulation::create_triangulation().
     * We call a table, which computes the corresponding entries of a full
     * connectivity table a reduced table.
     *
     * The entries of the reduced table are as follows for 1d-3d:
     *
     * 1D :    | 0 1    2d:    | 0 1 2    3d:    | 0 1 2 3
     *      ---+-----       ---+-------       ---+--------
     *       0 |             0 |               0 |
     *       1 | x n         1 | x             1 | x
     *                       2 | s x n         2 | s x
     *                                         3 |     x n
     *
     * with markers highlighting the reason for the entry x:=bounding entities;
     * n:= neighboring entities; s:=sub-cell data
     */
    struct Connectivity
    {
      Connectivity(const unsigned int                dim,
                   const std::vector<ReferenceCell> &cell_types)
        : dim(dim)
        , cell_types(cell_types)
      {}

      inline TriaObjectsOrientations &
      entity_orientations(const unsigned int structdim)
      {
        if (structdim == 1)
          return line_orientation;

        AssertDimension(structdim, 2);

        return quad_orientation;
      }

      inline const TriaObjectsOrientations &
      entity_orientations(const unsigned int structdim) const
      {
        if (structdim == 1)
          return line_orientation;

        AssertDimension(structdim, 2);

        return quad_orientation;
      }

      inline std::vector<ReferenceCell> &
      entity_types(const unsigned int structdim)
      {
        if (structdim == dim)
          return cell_types;

        // for vertices/lines the entity types are clear (0/1)
        AssertDimension(structdim, 2);
        AssertDimension(dim, 3);

        return quad_types;
      }

      inline const std::vector<ReferenceCell> &
      entity_types(const unsigned int structdim) const
      {
        if (structdim == dim)
          return cell_types;

        // for vertices/lines the entity types are clear (0/1)
        AssertDimension(structdim, 2);
        AssertDimension(dim, 3);

        return quad_types;
      }

      inline ArrayOfArrays &
      entity_to_entities(const unsigned int from, const unsigned int to)
      {
        if (from == dim && to == dim)
          return neighbors;
        else if (from == dim && to == dim - 1)
          return cell_entities;
        else if (dim == 3 && from == 2 && to == 0)
          return quad_vertices;
        else if (dim == 3 && from == 2 && to == 1)
          return quad_lines;
        else if (from == 1 && to == 0)
          return line_vertices;

        DEAL_II_NOT_IMPLEMENTED();

        return cell_entities;
      }

      inline const ArrayOfArrays &
      entity_to_entities(const unsigned int from, const unsigned int to) const
      {
        if (from == dim && to == dim)
          return neighbors;
        else if (from == dim && to == dim - 1)
          return cell_entities;
        else if (dim == 3 && from == 2 && to == 0)
          return quad_vertices;
        else if (dim == 3 && from == 2 && to == 1)
          return quad_lines;
        else if (from == 1 && to == 0)
          return line_vertices;

        DEAL_II_NOT_IMPLEMENTED();

        return cell_entities;
      }

    private:
      const unsigned int         dim;
      std::vector<ReferenceCell> cell_types;

      ArrayOfArrays line_vertices;

      TriaObjectsOrientations line_orientation;

      ArrayOfArrays quad_vertices;
      ArrayOfArrays quad_lines;

      TriaObjectsOrientations quad_orientation;

      ArrayOfArrays cell_entities;
      ArrayOfArrays neighbors;

      std::vector<ReferenceCell> quad_types;
    };



    /**
     * Determine the neighbors of all cells.
     *
     * @p con_cf connectivity cell-face
     * @return connectivity cell-cell (for each cell-face it contains the
     *   index of the neighboring cell or -1 for boundary face)
     */
    ArrayOfArrays
    determine_neighbors(const ArrayOfArrays &con_cf)
    {
      const unsigned int n_faces =
        *std::max_element(con_cf.columns.begin(), con_cf.columns.end()) + 1;

      // clear and initialize with -1 (assume that all faces are at the
      // boundary)
      std::vector<unsigned int> columns_cc(con_cf.columns.size(), -1);
      std::vector<std::size_t>  offsets_cc = con_cf.offsets;

      std::vector<std::pair<unsigned int, unsigned int>> neighbors(n_faces,
                                                                   {-1, -1});

      // loop over all cells
      unsigned int global_face_index = 0;
      for (unsigned int cell = 0; cell < con_cf.size(); ++cell)
        {
          // ... and all its faces
          const auto faces = con_cf[cell];
          for (unsigned int f = 0; f < faces.size(); ++f, ++global_face_index)
            {
              if (neighbors[faces[f]].first == static_cast<unsigned int>(-1))
                {
                  // face is visited the first time -> save the visiting cell
                  // and the face pointer
                  neighbors[faces[f]] =
                    std::pair<unsigned int, unsigned int>(cell,
                                                          global_face_index);
                }
              else
                {
                  // face is visited the second time -> now we know the cells
                  // on both sides of the face and we can determine for both
                  // cells the neighbor
                  columns_cc[global_face_index] = neighbors[faces[f]].first;
                  columns_cc[neighbors[faces[f]].second] = cell;
                }
            }
        }

      return ArrayOfArrays(std::move(offsets_cc), std::move(columns_cc));
    }



    /**
     * Build entities of dimension d (with 0<d<dim). Entities are described by
     * a set of vertices.
     *
     * Furthermore, the function determines for each cell of which d-dimensional
     * entity it consists of and its orientation relative to the cell.
     */
    template <int max_n_vertices, typename FU>
    void
    build_face_entities_templated(
      const unsigned int                face_dimensionality,
      const std::vector<ReferenceCell> &cell_types,
      const ArrayOfArrays              &cells_to_vertices,
      ArrayOfArrays                    &cells_to_faces, // result
      ArrayOfArrays                    &crs_0,          // result
      TriaObjectsOrientations          &orientations,   // result
      const FU                         &second_key_function)
    {
      const bool compatibility_mode = true;

      const std::vector<std::size_t> &cells_to_vertices_offsets =
        cells_to_vertices.offsets;
      const std::vector<unsigned int> &cells_to_vertices_indices =
        cells_to_vertices.columns;
      std::vector<std::size_t> &cells_to_faces_offsets = cells_to_faces.offsets;
      std::vector<unsigned int> &cells_to_faces_indices =
        cells_to_faces.columns;

      // note: we do not pre-allocate memory for these arrays because it turned
      // out that counting unique entities is more expensive than push_back().
      std::vector<std::size_t>  &offsets_0 = crs_0.offsets;
      std::vector<unsigned int> &columns_0 = crs_0.columns;

      // clear
      offsets_0 = {};
      columns_0 = {};

      unsigned int n_entities = 0;

      for (const auto &c : cell_types)
        {
          // Make sure that there are only two possibilities for
          // face_dimensionality that we can cover with the ?: statement below:
          Assert((face_dimensionality == c.get_dimension() - 1) ||
                   ((c.get_dimension() == 3) && (face_dimensionality == 1)),
                 ExcInternalError());
          n_entities +=
            (face_dimensionality == c.get_dimension() - 1 ? c.n_faces() :
                                                            c.n_lines());
        }

      // step 1: store each d-dimensional entity of a cell (described by their
      // vertices) into a vector and create a key for them
      //
      // note: it turned out to be more efficient to have a vector of tuples
      // than to have two vectors (sorting becomes inefficient)
      std::vector<
        std::tuple<std::array<unsigned int, max_n_vertices>, unsigned int>>
        keys; // key (sorted vertices), cell-entity index

      std::vector<std::array<unsigned int, max_n_vertices>> ad_entity_vertices;
      std::vector<ReferenceCell>                            ad_entity_types;
      std::vector<std::array<unsigned int, max_n_vertices>> ad_compatibility;

      keys.reserve(n_entities);
      ad_entity_vertices.reserve(n_entities);
      ad_entity_types.reserve(n_entities);
      ad_compatibility.reserve(n_entities);

      cells_to_faces_offsets.resize(cell_types.size() + 1);
      cells_to_faces_offsets[0] = 0;

      static const unsigned int offset = 1;

      // loop over all cells
      for (unsigned int c = 0, counter = 0; c < cell_types.size(); ++c)
        {
          const auto &cell_type = cell_types[c];

          // Make sure that there are only two possibilities for
          // face_dimensionality that we can cover with the ?: statement below:
          Assert((face_dimensionality == cell_type.get_dimension() - 1) ||
                   ((cell_type.get_dimension() == 3) &&
                    (face_dimensionality == 1)),
                 ExcInternalError());
          const unsigned int n_face_entities =
            (face_dimensionality == cell_type.get_dimension() - 1 ?
               cell_type.n_faces() :
               cell_type.n_lines());
          cells_to_faces_offsets[c + 1] =
            cells_to_faces_offsets[c] + n_face_entities;

          // ... collect vertices of cell
          const dealii::ArrayView<const unsigned int> local_vertices(
            cells_to_vertices_indices.data() + cells_to_vertices_offsets[c],
            cells_to_vertices_offsets[c + 1] - cells_to_vertices_offsets[c]);

          // ... loop over all its entities
          for (unsigned int e = 0; e < n_face_entities; ++e)
            {
              // ... determine global entity vertices
              std::array<unsigned int, max_n_vertices> entity_vertices;
              std::fill(entity_vertices.begin(), entity_vertices.end(), 0);

              // Same as above, make sure that there are only two possibilities
              // for face_dimensionality that we can cover with the ?: statement
              // below:
              Assert((face_dimensionality == cell_type.get_dimension() - 1) ||
                       ((cell_type.get_dimension() == 3) &&
                        (face_dimensionality == 1)),
                     ExcInternalError());
              for (unsigned int i = 0;
                   i < (face_dimensionality == cell_type.get_dimension() - 1 ?
                          cell_type.face_reference_cell(e).n_vertices() :
                          ReferenceCells::Line.n_vertices());
                   ++i)
                entity_vertices[i] =
                  local_vertices
                    [face_dimensionality == cell_type.get_dimension() - 1 ?
                       cell_type.face_to_cell_vertices(
                         e, i, numbers::default_geometric_orientation) :
                       cell_type.line_to_cell_vertices(e, i)] +
                  offset;

              // ... create key
              std::array<unsigned int, max_n_vertices> key = entity_vertices;
              std::sort(key.begin(), key.end());
              keys.emplace_back(key, counter++);

              ad_entity_vertices.emplace_back(entity_vertices);

              if (face_dimensionality == cell_type.get_dimension() - 1)
                ad_entity_types.emplace_back(cell_type.face_reference_cell(e));
              else if (face_dimensionality == cell_type.get_dimension() - 2)
                {
                  // Since we only deal with meshes up to dimension 3,
                  // something that has co-dimensionality -2 must either be
                  // a vertex or a line, regardless of what the object we are
                  // working on actually us:
                  if (face_dimensionality == 0)
                    ad_entity_types.emplace_back(ReferenceCells::Vertex);
                  else if (face_dimensionality == 1)
                    ad_entity_types.emplace_back(ReferenceCells::Line);
                  else
                    // But it's probably useful to be conservative if someone
                    // ever comes along to implement 4d meshes:
                    DEAL_II_ASSERT_UNREACHABLE();
                }
              else
                DEAL_II_ASSERT_UNREACHABLE();

              if (compatibility_mode)
                ad_compatibility.emplace_back(
                  second_key_function(entity_vertices, cell_type, c, e));
            }
        }

      cells_to_faces_indices.resize(keys.size());
      orientations.reinit(keys.size());

      // step 2: sort according to key so that entities with same key can be
      // merged
      std::sort(keys.begin(), keys.end());


      if (compatibility_mode)
        {
          unsigned int n_unique_entities        = 0;
          unsigned int n_unique_entity_vertices = 0;

          std::array<unsigned int, max_n_vertices> ref_key, new_key;
          std::fill(ref_key.begin(), ref_key.end(), 0);
          for (unsigned int i = 0; i < keys.size(); ++i)
            {
              const auto offset_i = std::get<1>(keys[i]);

              if (ref_key != std::get<0>(keys[i]))
                {
                  ref_key = std::get<0>(keys[i]);

                  ++n_unique_entities;
                  n_unique_entity_vertices +=
                    ad_entity_types[offset_i].n_vertices();

                  new_key = ad_compatibility[offset_i];
                }

              std::get<0>(keys[i]) = new_key;
            }

          std::sort(keys.begin(), keys.end());

          offsets_0.reserve(n_unique_entities + 1);
          columns_0.reserve(n_unique_entity_vertices);
        }


      std::array<unsigned int, max_n_vertices> ref_key;
      std::array<unsigned int, max_n_vertices> ref_indices;
      std::fill(ref_key.begin(), ref_key.end(), 0);

      unsigned int counter = dealii::numbers::invalid_unsigned_int;
      for (unsigned int i = 0; i < keys.size(); i++)
        {
          const auto offset_i = std::get<1>(keys[i]);

          if (ref_key != std::get<0>(keys[i]))
            {
              // new key: default orientation is correct
              ++counter;
              ref_key     = std::get<0>(keys[i]);
              ref_indices = ad_entity_vertices[offset_i];

              offsets_0.push_back(columns_0.size());
              for (const auto j : ad_entity_vertices[offset_i])
                if (j != 0)
                  columns_0.push_back(j - offset);
            }
          else
            {
              // previously seen key: set orientation relative to the first
              // occurrence
              orientations.set_combined_orientation(
                offset_i,
                ad_entity_types[offset_i]
                  .template get_combined_orientation<unsigned int>(
                    make_array_view(ad_entity_vertices[offset_i].begin(),
                                    ad_entity_vertices[offset_i].begin() +
                                      ad_entity_types[offset_i].n_vertices()),
                    make_array_view(ref_indices.begin(),
                                    ref_indices.begin() +
                                      ad_entity_types[offset_i].n_vertices())));
            }
          cells_to_faces_indices[offset_i] = counter;
        }
      offsets_0.push_back(columns_0.size());
    }



    /**
     * Call the right templated function to be able to use std::array instead
     * of std::vector.
     */
    template <typename FU>
    void
    build_face_entities(const unsigned int                face_dimensionality,
                        const std::vector<ReferenceCell> &cell_types,
                        const ArrayOfArrays              &cells_to_vertices,
                        ArrayOfArrays                    &cells_to_faces,
                        ArrayOfArrays                    &crs_0,
                        TriaObjectsOrientations          &orientations,
                        const FU                         &second_key_function)
    {
      unsigned int max_n_vertices = 0;

      // If we are dealing with faces of cells, figure out how many vertices
      // each face may have. Otherwise, we're in 3d and are dealing with
      // lines, for which we know the number of vertices:
      for (const auto &c : cell_types)
        if (face_dimensionality == c.get_dimension() - 1)
          {
            for (unsigned int f = 0; f < c.n_faces(); ++f)
              max_n_vertices =
                std::max(max_n_vertices, c.face_reference_cell(f).n_vertices());
          }
        else if (face_dimensionality == 1)
          max_n_vertices = std::max(max_n_vertices, 2u);
        else
          DEAL_II_ASSERT_UNREACHABLE();

      if (max_n_vertices == 2)
        build_face_entities_templated<2>(face_dimensionality,
                                         cell_types,
                                         cells_to_vertices,
                                         cells_to_faces,
                                         crs_0,
                                         orientations,
                                         second_key_function);
      else if (max_n_vertices == 3)
        build_face_entities_templated<3>(face_dimensionality,
                                         cell_types,
                                         cells_to_vertices,
                                         cells_to_faces,
                                         crs_0,
                                         orientations,
                                         second_key_function);
      else if (max_n_vertices == 4)
        build_face_entities_templated<4>(face_dimensionality,
                                         cell_types,
                                         cells_to_vertices,
                                         cells_to_faces,
                                         crs_0,
                                         orientations,
                                         second_key_function);
      else
        AssertThrow(false, dealii::StandardExceptions::ExcNotImplemented());
    }



    /**
     * Build surface lines described by:
     *  - connectivity quad -> line
     *  - orientation of line relative to the quad
     *
     * Furthermore, the type of the quad is determined.
     */
    inline void
    build_intersection(const std::vector<ReferenceCell> &cell_types,
                       const ArrayOfArrays              &cells_to_vertices,
                       const ArrayOfArrays              &cells_to_lines,
                       const ArrayOfArrays              &lines_to_vertices,
                       const ArrayOfArrays              &cells_to_quads,
                       const ArrayOfArrays              &quads_to_vertices,
                       const TriaObjectsOrientations    &ori_cq,
                       ArrayOfArrays              &quads_to_lines, // result
                       TriaObjectsOrientations    &ori_ql,         // result
                       std::vector<ReferenceCell> &quad_t_id       // result
    )
    {
      // reset output
      quads_to_lines.offsets = {};
      quads_to_lines.columns = {};

      quads_to_lines.offsets.resize(quads_to_vertices.offsets.size());
      quads_to_lines.offsets[0] = 0;

      quad_t_id.resize(quads_to_vertices.size());

      // count the number of lines of each face
      for (unsigned int c = 0; c < cells_to_quads.size(); ++c)
        {
          // loop over faces
          for (unsigned int f_ = cells_to_quads.offsets[c], f_index = 0;
               f_ < cells_to_quads.offsets[c + 1];
               ++f_, ++f_index)
            {
              const unsigned int f = cells_to_quads.columns[f_];

              quads_to_lines.offsets[f + 1] =
                cell_types[c].face_reference_cell(f_index).n_lines();
            }
        }

      // use the counts to determine the offsets -> prefix sum
      for (unsigned int i = 0; i < quads_to_lines.size(); ++i)
        quads_to_lines.offsets[i + 1] += quads_to_lines.offsets[i];

      // allocate memory
      quads_to_lines.columns.resize(quads_to_lines.offsets.back());
      ori_ql.reinit(quads_to_lines.offsets.back());

      // loop over cells
      for (unsigned int c = 0; c < cells_to_quads.size(); ++c)
        {
          const ReferenceCell cell_type = cell_types[c];

          // loop over faces
          for (unsigned int f_ = cells_to_quads.offsets[c], f_index = 0;
               f_ < cells_to_quads.offsets[c + 1];
               ++f_, ++f_index)
            {
              const unsigned int f = cells_to_quads.columns[f_];

              // only faces with default orientation have to do something
              if (ori_cq.get_combined_orientation(f_) !=
                  numbers::default_geometric_orientation)
                continue;

              // loop over the lines of this face
              quad_t_id[f] = cell_type.face_reference_cell(f_index);
              for (unsigned int l = 0; l < quad_t_id[f].n_lines(); ++l)
                {
                  // determine global index of line
                  const unsigned int local_line_index =
                    cell_type.face_to_cell_lines(
                      f_index, l, numbers::default_geometric_orientation);
                  const unsigned int global_line_index =
                    cells_to_lines
                      .columns[cells_to_lines.offsets[c] + local_line_index];
                  quads_to_lines.columns[quads_to_lines.offsets[f] + l] =
                    global_line_index;

                  // determine orientation of line
                  bool same = true;
                  for (unsigned int v = 0; v < 2; ++v)
                    if (cells_to_vertices
                          .columns[cells_to_vertices.offsets[c] +
                                   cell_type.face_and_line_to_cell_vertices(
                                     f_index,
                                     l,
                                     v,
                                     numbers::default_geometric_orientation)] !=
                        lines_to_vertices.columns
                          [lines_to_vertices.offsets[global_line_index] + v])
                      {
                        same = false;
                        break;
                      }

                  // ... comparison gives orientation
                  ori_ql.set_combined_orientation(
                    quads_to_lines.offsets[f] + l,
                    same ? numbers::default_geometric_orientation :
                           numbers::reverse_line_orientation);
                }
            }
        }
    }



    /**
     * Build the reduced connectivity table for the given dimension @p dim.
     *
     * This function is inspired by the publication Anders Logg "Efficient
     * Representation of Computational Meshes" and the FEniCS's DOLFIN mesh
     * implementation. It has been strongly adjusted to efficiently solely meet
     * our connectivity needs while sacrificing some of the flexibility there.
     */
    Connectivity
    build_connectivity(const unsigned int                dim,
                       const std::vector<ReferenceCell> &cell_types,
                       const ArrayOfArrays              &cells_to_vertices)
    {
      Connectivity connectivity(dim, cell_types);

      ArrayOfArrays temp1; // needed for 3d

      if (dim == 1)
        connectivity.entity_to_entities(1, 0) = cells_to_vertices;

      if (dim == 2 || dim == 3) // build lines
        {
          TriaObjectsOrientations dummy;

          build_face_entities(
            1,
            connectivity.entity_types(dim),
            cells_to_vertices,
            dim == 2 ? connectivity.entity_to_entities(2, 1) : temp1,
            connectivity.entity_to_entities(1, 0),
            dim == 2 ? connectivity.entity_orientations(1) : dummy,
            [](auto key, const auto &, const auto &, const auto &) {
              //  to ensure same enumeration as in deal.II
              return key;
            });
        }

      if (dim == 3) // build quads
        {
          build_face_entities(
            2,
            connectivity.entity_types(3),
            cells_to_vertices,
            connectivity.entity_to_entities(3, 2),
            connectivity.entity_to_entities(2, 0),
            connectivity.entity_orientations(2),
            [&](auto key, // of type std::array<unsigned int, max_n_vertices>
                          // but max_n_vertices is not known here
                const ReferenceCell &cell_type,
                const unsigned int  &c,
                const unsigned int  &f) {
              //  to ensure same enumeration as in deal.II
              AssertIndexRange(cell_type.face_reference_cell(f).n_lines(),
                               key.size() + 1);

              unsigned int l = 0;

              for (; l < cell_type.face_reference_cell(f).n_lines(); ++l)
                {
                  AssertIndexRange(l, key.size());
                  key[l] = temp1[c][cell_type.face_to_cell_lines(
                             f, l, numbers::default_geometric_orientation)] +
                           1 /*offset!*/;
                }

              for (; l < key.size(); ++l)
                key[l] = 0;

              return key;
            });

          // create connectivity: quad -> line
          build_intersection(connectivity.entity_types(3),
                             cells_to_vertices,
                             temp1,
                             connectivity.entity_to_entities(1, 0),
                             connectivity.entity_to_entities(3, 2),
                             connectivity.entity_to_entities(2, 0),
                             connectivity.entity_orientations(2),
                             connectivity.entity_to_entities(2, 1),
                             connectivity.entity_orientations(1),
                             connectivity.entity_types(2));
        }

      // determine neighbors
      connectivity.entity_to_entities(dim, dim) =
        determine_neighbors(connectivity.entity_to_entities(dim, dim - 1));

      return connectivity;
    }



    /**
     * Preprocessing step to remove the template argument dim.
     */
    template <int dim>
    Connectivity
    build_connectivity(const std::vector<CellData<dim>> &cells)
    {
      AssertThrow(cells.size() > 0, ExcMessage("No cells have been provided!"));

      // determine cell types and process vertices
      std::vector<unsigned int> cell_vertices;
      cell_vertices.reserve(std::accumulate(cells.begin(),
                                            cells.end(),
                                            0u,
                                            [](const unsigned int   accumulator,
                                               const CellData<dim> &cell) {
                                              return accumulator +
                                                     cell.vertices.size();
                                            }));

      std::vector<std::size_t> cell_vertices_offsets;
      cell_vertices_offsets.reserve(cells.size() + 1);
      cell_vertices_offsets.push_back(0);

      std::vector<ReferenceCell> cell_types;
      cell_types.reserve(cells.size());

      // loop over cells and create CRS
      for (const auto &cell : cells)
        {
          if constexpr (running_in_debug_mode())
            {
              auto vertices_unique = cell.vertices;
              std::sort(vertices_unique.begin(), vertices_unique.end());
              vertices_unique.erase(std::unique(vertices_unique.begin(),
                                                vertices_unique.end()),
                                    vertices_unique.end());

              Assert(
                vertices_unique.size() == cell.vertices.size(),
                ExcMessage(
                  "The definition of a cell refers to the same vertex several "
                  "times. This is not possible. A common reason is that "
                  "CellData::vertices has a size that does not match the "
                  "size expected from the reference cell. Please resize "
                  "CellData::vertices or use the appropriate constructor of "
                  "CellData."));
            }

          cell_types.push_back(
            ReferenceCell::n_vertices_to_type(dim, cell.vertices.size()));

          // create CRS of vertices (to remove template argument dim)
          cell_vertices.insert(cell_vertices.end(),
                               cell.vertices.begin(),
                               cell.vertices.end());
          cell_vertices_offsets.push_back(cell_vertices.size());
        }

      // do the actual work
      return build_connectivity(dim,
                                cell_types,
                                {std::move(cell_vertices_offsets),
                                 std::move(cell_vertices)});
    }
  } // namespace TriangulationImplementation
} // namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif

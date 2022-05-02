// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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

#ifndef dealii_tria_connectivity_h
#define dealii_tria_connectivity_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/ndarray.h>

#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria_description.h>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace TriangulationImplementation
  {
    /**
     * Interface of geometric cell entities with the focus on creating a
     * reduced connectivity table.
     */
    struct CellTypeBase
    {
      /**
       * Default destructor.
       */
      virtual ~CellTypeBase() = default;

      /**
       * Number of sub-entities of dimension @p d.
       */
      virtual unsigned int
      n_entities(const unsigned int d) const
      {
        Assert(false, ExcNotImplemented());
        (void)d;

        return 0;
      }

      /**
       * Number of vertices of the @p e-th sub-entity of dimension @p d.
       */
      virtual dealii::ArrayView<const unsigned int>
      vertices_of_entity(const unsigned int d, const unsigned int e) const
      {
        Assert(false, ExcNotImplemented());
        (void)d;
        (void)e;

        return {};
      }

      /**
       * Geometric entity type of the @p e-th sub-entity of dimension @p d.
       */
      virtual dealii::ReferenceCell
      type_of_entity(const unsigned int d, const unsigned int e) const
      {
        Assert(false, ExcNotImplemented());
        (void)d;
        (void)e;

        return dealii::ReferenceCells::Vertex;
      }

      /**
       * Number of lines of @p face-th surface.
       */
      virtual unsigned int
      n_lines_of_surface(const unsigned int face) const
      {
        Assert(false, ExcNotImplemented());
        (void)face;

        return 0;
      }

      /**
       * Index of the @p line-th lines of @p face-th surface.
       */
      virtual unsigned int
      nth_line_of_surface(const unsigned int line,
                          const unsigned int face) const
      {
        Assert(false, ExcNotImplemented());
        (void)line;
        (void)face;

        return 0;
      }

      /**
       * Vertex indices of the @p line-th lines of @p face-th surface.
       */
      virtual const std::array<unsigned int, 2> &
      vertices_of_nth_line_of_surface(const unsigned int line,
                                      const unsigned int face) const
      {
        Assert(false, ExcNotImplemented());
        (void)line;
        (void)face;

        const static std::array<unsigned int, 2> table = {};

        return table;
      }
    };



    /**
     * Implementation for lines.
     */
    struct CellTypeLine : public CellTypeBase
    {
      dealii::ArrayView<const unsigned int>
      vertices_of_entity(const unsigned int d,
                         const unsigned int e) const override
      {
        (void)e;

        if (d == 1)
          {
            static const std::array<unsigned int, 2> table = {{0, 1}};

            AssertDimension(e, 0);

            return {table};
          }

        Assert(false, ExcNotImplemented());

        return {};
      }

      dealii::ReferenceCell
      type_of_entity(const unsigned int d, const unsigned int e) const override
      {
        (void)e;

        if (d == 1)
          return dealii::ReferenceCells::Line;

        Assert(false, ExcNotImplemented());

        return dealii::ReferenceCells::Vertex;
      }

      unsigned int
      n_entities(const unsigned int d) const override
      {
        static std::array<unsigned int, 2> table = {{2, 1}};
        return table[d];
      }
    };



    /**
     * Implementation for triangles.
     */
    struct CellTypeTri : public CellTypeBase
    {
      dealii::ArrayView<const unsigned int>
      vertices_of_entity(const unsigned int d,
                         const unsigned int e) const override
      {
        if (d == 2)
          {
            static const std::array<unsigned int, 3> table = {{0, 1, 2}};

            AssertDimension(e, 0);

            return {table};
          }

        if (d == 1)
          {
            static const dealii::ndarray<unsigned int, 3, 2> table = {
              {{{0, 1}}, {{1, 2}}, {{2, 0}}}};

            return {table[e]};
          }

        Assert(false, ExcNotImplemented());

        return {};
      }

      virtual dealii::ReferenceCell
      type_of_entity(const unsigned int d, const unsigned int e) const override
      {
        (void)e;

        if (d == 2)
          return dealii::ReferenceCells::Triangle;

        if (d == 1)
          return dealii::ReferenceCells::Line;

        Assert(false, ExcNotImplemented());

        return dealii::ReferenceCells::Vertex;
      }

      unsigned int
      n_entities(const unsigned int d) const override
      {
        static std::array<unsigned int, 3> table = {{3, 3, 1}};
        return table[d];
      }
    };



    /**
     * Implementation for quadrilaterals.
     */
    struct CellTypeQuad : public CellTypeBase
    {
      dealii::ArrayView<const unsigned int>
      vertices_of_entity(const unsigned int d,
                         const unsigned int e) const override
      {
        if (d == 2)
          {
            static const std::array<unsigned int, 4> table = {{0, 1, 2, 3}};

            AssertDimension(e, 0);

            return {table};
          }

        if (d == 1)
          {
            static const dealii::ndarray<unsigned int, 4, 2> table = {
              {{{0, 2}}, {{1, 3}}, {{0, 1}}, {{2, 3}}}};

            return {table[e]};
          }

        Assert(false, ExcNotImplemented());

        return {};
      }

      virtual dealii::ReferenceCell
      type_of_entity(const unsigned int d, const unsigned int e) const override
      {
        (void)e;

        if (d == 2)
          return dealii::ReferenceCells::Quadrilateral;

        if (d == 1)
          return dealii::ReferenceCells::Line;

        Assert(false, ExcNotImplemented());

        return dealii::ReferenceCells::Vertex;
      }

      unsigned int
      n_entities(const unsigned int d) const override
      {
        static std::array<unsigned int, 3> table = {{4, 4, 1}};
        return table[d];
      }
    };



    /**
     * Implementation for tetrahedrons.
     */
    struct CellTypeTet : public CellTypeBase
    {
      dealii::ArrayView<const unsigned int>
      vertices_of_entity(const unsigned int d,
                         const unsigned int e) const override
      {
        if (d == 3)
          {
            static const std::array<unsigned int, 4> table = {{0, 1, 2, 3}};

            AssertDimension(e, 0);

            return {table};
          }

        if (d == 2)
          {
            static const dealii::ndarray<unsigned int, 4, 3> table = {
              {{{0, 1, 2}}, {{1, 0, 3}}, {{0, 2, 3}}, {{2, 1, 3}}}};

            return {table[e]};
          }

        if (d == 1)
          {
            static const dealii::ndarray<unsigned int, 6, 2> table = {
              {{{0, 1}}, {{1, 2}}, {{2, 0}}, {{0, 3}}, {{1, 3}}, {{2, 3}}}};

            return {table[e]};
          }

        Assert(false, ExcNotImplemented());

        return {};
      }

      virtual dealii::ReferenceCell
      type_of_entity(const unsigned int d, const unsigned int e) const override
      {
        (void)e;

        if (d == 3)
          return dealii::ReferenceCells::Tetrahedron;

        if (d == 2)
          return dealii::ReferenceCells::Triangle;

        if (d == 1)
          return dealii::ReferenceCells::Line;

        Assert(false, ExcNotImplemented());

        return dealii::ReferenceCells::Vertex;
      }

      unsigned int
      n_entities(const unsigned int d) const override
      {
        static std::array<unsigned int, 4> table = {{4, 6, 4, 1}};
        return table[d];
      }

      unsigned int
      n_lines_of_surface(const unsigned int line) const override
      {
        (void)line;
        return 3;
      }

      unsigned int
      nth_line_of_surface(const unsigned int line,
                          const unsigned int face) const override
      {
        const static dealii::ndarray<unsigned int, 4, 3> table = {
          {{{0, 1, 2}}, {{0, 3, 4}}, {{2, 5, 3}}, {{1, 4, 5}}}};

        return table[face][line];
      }

      const std::array<unsigned int, 2> &
      vertices_of_nth_line_of_surface(const unsigned int line,
                                      const unsigned int face) const override
      {
        const static dealii::ndarray<unsigned int, 4, 3, 2> table = {
          {{{{{0, 1}}, {{1, 2}}, {{2, 0}}}},
           {{{{1, 0}}, {{0, 3}}, {{3, 1}}}},
           {{{{0, 2}}, {{2, 3}}, {{3, 0}}}},
           {{{{2, 1}}, {{1, 3}}, {{3, 2}}}}}};

        return table[face][line];
      }
    };


    /**
     * Implementation for pyramids.
     */

    struct CellTypePyramid : public CellTypeBase
    {
      dealii::ArrayView<const unsigned int>
      vertices_of_entity(const unsigned int d,
                         const unsigned int e) const override
      {
        if (d == 3)
          {
            static const std::array<unsigned int, 5> table = {{0, 1, 2, 3, 4}};

            AssertDimension(e, 0);

            return {table};
          }

        if (d == 2)
          {
            if (e == 0)
              {
                static const std::array<unsigned int, 4> table = {{0, 1, 2, 3}};
                return {table};
              }

            static const dealii::ndarray<unsigned int, 4, 3> table = {
              {{{0, 2, 4}}, {{3, 1, 4}}, {{1, 0, 4}}, {{2, 3, 4}}}};

            return {table[e - 1]};
          }

        if (d == 1)
          {
            static const dealii::ndarray<unsigned int, 8, 2> table = {
              {{{0, 2}},
               {{1, 3}},
               {{0, 1}},
               {{2, 3}},
               {{0, 4}},
               {{1, 4}},
               {{2, 4}},
               {{3, 4}}}};

            return {table[e]};
          }

        Assert(false, ExcNotImplemented());

        return {};
      }

      virtual dealii::ReferenceCell
      type_of_entity(const unsigned int d, const unsigned int e) const override
      {
        (void)e;

        if (d == 3)
          return dealii::ReferenceCells::Pyramid;

        if (d == 2 && e == 0)
          return dealii::ReferenceCells::Quadrilateral;
        else if (d == 2)
          return dealii::ReferenceCells::Triangle;

        if (d == 1)
          return dealii::ReferenceCells::Line;

        Assert(false, ExcNotImplemented());

        return dealii::ReferenceCells::Vertex;
      }

      unsigned int
      n_entities(const unsigned int d) const override
      {
        static std::array<unsigned int, 4> table = {{5, 8, 5, 1}};
        return table[d];
      }

      unsigned int
      n_lines_of_surface(const unsigned int surface) const override
      {
        if (surface == 0)
          return 4;

        return 3;
      }

      unsigned int
      nth_line_of_surface(const unsigned int line,
                          const unsigned int face) const override
      {
        const static dealii::ndarray<unsigned int, 5, 4> table = {
          {{{0, 1, 2, 3}},
           {{0, 6, 4, numbers::invalid_unsigned_int}},
           {{1, 5, 7, numbers::invalid_unsigned_int}},
           {{2, 4, 5, numbers::invalid_unsigned_int}},
           {{3, 7, 6, numbers::invalid_unsigned_int}}}};

        return table[face][line];
      }

      const std::array<unsigned int, 2> &
      vertices_of_nth_line_of_surface(const unsigned int line,
                                      const unsigned int face) const override
      {
        static const unsigned int X = static_cast<unsigned int>(-1);

        const static dealii::ndarray<unsigned int, 5, 4, 2> table = {
          {{{{{0, 2}}, {{1, 3}}, {{0, 1}}, {{2, 3}}}},
           {{{{0, 2}}, {{2, 4}}, {{4, 0}}, {{X, X}}}},
           {{{{3, 1}}, {{1, 4}}, {{4, 3}}, {{X, X}}}},
           {{{{1, 0}}, {{0, 4}}, {{4, 1}}, {{X, X}}}},
           {{{{2, 3}}, {{3, 4}}, {{4, 2}}, {{X, X}}}}}};

        return table[face][line];
      }
    };



    /**
     * Implementation for wedges.
     */
    struct CellTypeWedge : public CellTypeBase
    {
      dealii::ArrayView<const unsigned int>
      vertices_of_entity(const unsigned int d,
                         const unsigned int e) const override
      {
        if (d == 3)
          {
            static const std::array<unsigned int, 6> table = {
              {0, 1, 2, 3, 4, 5}};

            AssertDimension(e, 0);

            return {table};
          }

        if (d == 2)
          {
            if (e == 0 || e == 1)
              {
                static const dealii::ndarray<unsigned int, 2, 3> table = {
                  {{{1, 0, 2}}, {{3, 4, 5}}}};

                return {table[e]};
              }

            static const dealii::ndarray<unsigned int, 3, 4> table = {
              {{{0, 1, 3, 4}}, {{1, 2, 4, 5}}, {{2, 0, 5, 3}}}};

            return {table[e - 2]};
          }

        if (d == 1)
          {
            static const dealii::ndarray<unsigned int, 9, 2> table = {
              {{{0, 1}},
               {{1, 2}},
               {{2, 0}},
               {{3, 4}},
               {{4, 5}},
               {{5, 3}},
               {{0, 3}},
               {{1, 4}},
               {{2, 5}}}};

            return {table[e]};
          }

        Assert(false, ExcNotImplemented());

        return {};
      }

      virtual dealii::ReferenceCell
      type_of_entity(const unsigned int d, const unsigned int e) const override
      {
        (void)e;

        if (d == 3)
          return dealii::ReferenceCells::Wedge;

        if (d == 2 && e > 1)
          return dealii::ReferenceCells::Quadrilateral;
        else if (d == 2)
          return dealii::ReferenceCells::Triangle;

        if (d == 1)
          return dealii::ReferenceCells::Line;

        Assert(false, ExcNotImplemented());

        return dealii::ReferenceCells::Vertex;
      }

      unsigned int
      n_entities(const unsigned int d) const override
      {
        static std::array<unsigned int, 4> table = {{6, 9, 5, 1}};
        return table[d];
      }

      unsigned int
      n_lines_of_surface(const unsigned int surface) const override
      {
        if (surface > 1)
          return 4;

        return 3;
      }

      unsigned int
      nth_line_of_surface(const unsigned int line,
                          const unsigned int face) const override
      {
        static const unsigned int X = static_cast<unsigned int>(-1);

        const static dealii::ndarray<unsigned int, 5, 4> table = {
          {{{0, 2, 1, X}},
           {{3, 4, 5, X}},
           {{6, 7, 0, 3}},
           {{7, 8, 1, 4}},
           {{8, 6, 5, 2}}}};

        return table[face][line];
      }

      const std::array<unsigned int, 2> &
      vertices_of_nth_line_of_surface(const unsigned int line,
                                      const unsigned int face) const override
      {
        static const unsigned int X = static_cast<unsigned int>(-1);

        const static dealii::ndarray<unsigned int, 5, 4, 2> table = {
          {{{{{1, 0}}, {{0, 2}}, {{2, 1}}, {{X, X}}}},
           {{{{3, 4}}, {{4, 5}}, {{5, 3}}, {{X, X}}}},
           {{{{0, 3}}, {{1, 4}}, {{0, 1}}, {{3, 4}}}},
           {{{{1, 4}}, {{2, 5}}, {{1, 2}}, {{4, 5}}}},
           {{{{2, 5}}, {{0, 3}}, {{2, 0}}, {{5, 3}}}}}};

        return table[face][line];
      }
    };



    /**
     * Implementation for hexahedra.
     */
    struct CellTypeHex : public CellTypeBase
    {
      dealii::ArrayView<const unsigned int>
      vertices_of_entity(const unsigned int d,
                         const unsigned int e) const override
      {
        if (d == 3)
          {
            static const std::array<unsigned int, 8> table = {
              {0, 1, 2, 3, 4, 5, 6, 7}};

            AssertDimension(e, 0);

            return {table};
          }

        if (d == 2)
          {
            static const dealii::ndarray<unsigned int, 6, 4> table = {
              {{{0, 2, 4, 6}},
               {{1, 3, 5, 7}},
               {{0, 4, 1, 5}},
               {{2, 6, 3, 7}},
               {{0, 1, 2, 3}},
               {{4, 5, 6, 7}}}};

            return {table[e]};
          }

        if (d == 1)
          {
            static const dealii::ndarray<unsigned int, 12, 2> table = {
              {{{0, 2}},
               {{1, 3}},
               {{0, 1}},
               {{2, 3}},
               {{4, 6}},
               {{5, 7}},
               {{4, 5}},
               {{6, 7}},
               {{0, 4}},
               {{1, 5}},
               {{2, 6}},
               {{3, 7}}}};

            return {table[e]};
          }

        Assert(false, ExcNotImplemented());

        return {};
      }

      virtual dealii::ReferenceCell
      type_of_entity(const unsigned int d, const unsigned int e) const override
      {
        (void)e;

        if (d == 3)
          return dealii::ReferenceCells::Hexahedron;

        if (d == 2)
          return dealii::ReferenceCells::Quadrilateral;

        if (d == 1)
          return dealii::ReferenceCells::Line;

        Assert(false, ExcNotImplemented());

        return dealii::ReferenceCells::Vertex;
      }

      unsigned int
      n_entities(const unsigned int d) const override
      {
        static std::array<unsigned int, 4> table = {{8, 12, 6, 1}};
        return table[d];
      }

      unsigned int
      n_lines_of_surface(const unsigned int surface) const override
      {
        (void)surface;
        return 4;
      }

      unsigned int
      nth_line_of_surface(const unsigned int line,
                          const unsigned int face) const override
      {
        const static dealii::ndarray<unsigned int, 6, 4> table = {
          {{{8, 10, 0, 4}},
           {{9, 11, 1, 5}},
           {{2, 6, 8, 9}},
           {{3, 7, 10, 11}},
           {{0, 1, 2, 3}},
           {{4, 5, 6, 7}}}};

        return table[face][line];
      }

      const std::array<unsigned int, 2> &
      vertices_of_nth_line_of_surface(const unsigned int line,
                                      const unsigned int face) const override
      {
        const static dealii::ndarray<unsigned int, 6, 4, 2> table = {
          {{{{{0, 4}}, {{2, 6}}, {{0, 2}}, {{4, 6}}}},
           {{{{1, 5}}, {{3, 7}}, {{1, 3}}, {{5, 7}}}},
           {{{{0, 1}}, {{4, 5}}, {{0, 4}}, {{1, 5}}}},
           {{{{2, 3}}, {{6, 7}}, {{2, 6}}, {{3, 7}}}},
           {{{{0, 2}}, {{1, 3}}, {{0, 1}}, {{2, 3}}}},
           {{{{4, 6}}, {{5, 7}}, {{4, 5}}, {{6, 7}}}}}};

        return table[face][line];
      }
    };



    /**
     * Compressed row storage sparse matrix. This class is similar to
     * SparsityPattern but reduced to the bare minimum as needed here - in the
     * context of setting up the connectivity - and allowing direct simplified
     * access to the entries.
     */
    template <typename T = unsigned int>
    struct CRS
    {
      /**
       * Default constructor.
       */
      CRS()
        : ptr{0} {};

      /**
       * Constructor which allows to set the internal fields directly.
       */
      CRS(const std::vector<std::size_t> &ptr, const std::vector<T> &col)
        : ptr(ptr)
        , col(col)
      {}

      // row index
      std::vector<std::size_t> ptr;

      // column index
      std::vector<T> col;
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
     *  - quad (2 - 3D), line (1 - 2D/3D) to vertices (0) to be able to process
     *    the user provided SubCellData during
     *    Triangulation::create_triangulation().
     * We call a table, which computes the corresponding entries of a full
     * connectivity table a reduced table.
     *
     * The entries of the reduced table are as follows for 1D-3D:
     *
     * 1D :    | 0 1    2D:    | 0 1 2    3D:    | 0 1 2 3
     *      ---+-----       ---+-------       ---+--------
     *       0 |             0 |               0 |
     *       1 | x n         1 | x             1 | x
     *                       2 | s x n         2 | s x
     *                                         3 |     x n
     *
     * with markers highlighting the reason for the entry x:=bounding entities;
     * n:= neighboring entities; s:=sub-cell data
     */
    template <typename T = unsigned int>
    struct Connectivity
    {
      Connectivity(const unsigned int                        dim,
                   const std::vector<dealii::ReferenceCell> &cell_types)
        : dim(dim)
        , cell_types(cell_types)
      {}

      inline std::vector<unsigned char> &
      entity_orientations(const unsigned int structdim)
      {
        if (structdim == 1)
          return line_orientation;

        AssertDimension(structdim, 2);

        return quad_orientation;
      }

      inline const std::vector<unsigned char> &
      entity_orientations(const unsigned int structdim) const
      {
        if (structdim == 1)
          return line_orientation;

        AssertDimension(structdim, 2);

        return quad_orientation;
      }

      inline std::vector<dealii::ReferenceCell> &
      entity_types(const unsigned int structdim)
      {
        if (structdim == dim)
          return cell_types;

        // for vertices/lines the entity types are clear (0/1)
        AssertDimension(structdim, 2);
        AssertDimension(dim, 3);

        return quad_types;
      }

      inline const std::vector<dealii::ReferenceCell> &
      entity_types(const unsigned int structdim) const
      {
        if (structdim == dim)
          return cell_types;

        // for vertices/lines the entity types are clear (0/1)
        AssertDimension(structdim, 2);
        AssertDimension(dim, 3);

        return quad_types;
      }

      inline CRS<T> &
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

        Assert(false, ExcNotImplemented());

        return cell_entities;
      }

      inline const CRS<T> &
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

        Assert(false, ExcNotImplemented());

        return cell_entities;
      }

    private:
      const unsigned int                 dim;
      std::vector<dealii::ReferenceCell> cell_types;

      CRS<T> line_vertices;

      std::vector<unsigned char> line_orientation;

      CRS<T> quad_vertices;
      CRS<T> quad_lines;

      std::vector<unsigned char> quad_orientation;

      CRS<T> cell_entities;
      CRS<T> neighbors;

      std::vector<dealii::ReferenceCell> quad_types;
    };



    /**
     * Determine the neighbors of all cells.
     *
     * @p con_cf connectivity cell-face
     * @p con_cc connectivity cell-cell (for each cell-face it contains the
     *   index of the neighboring cell or -1 for boundary face)
     */
    template <typename T>
    void
    determine_neighbors(const CRS<T> &con_cf, CRS<T> &con_cc)
    {
      const auto &col_cf = con_cf.col;
      const auto &ptr_cf = con_cf.ptr;

      auto &col_cc = con_cc.col;
      auto &ptr_cc = con_cc.ptr;

      const unsigned int n_faces =
        *std::max_element(col_cf.begin(), col_cf.end()) + 1;

      // clear and initialize with -1 (assume that all faces are at the
      // boundary)
      col_cc = std::vector<T>(col_cf.size(), -1);
      ptr_cc = ptr_cf;

      std::vector<std::pair<T, unsigned int>> neighbors(n_faces, {-1, -1});

      // loop over all cells
      for (unsigned int i_0 = 0; i_0 < ptr_cf.size() - 1; ++i_0)
        {
          // ... and all its faces
          for (std::size_t j_0 = ptr_cf[i_0]; j_0 < ptr_cf[i_0 + 1]; ++j_0)
            {
              if (neighbors[col_cf[j_0]].first == static_cast<unsigned int>(-1))
                {
                  // face is visited the first time -> save the visiting cell
                  // and the face pointer
                  neighbors[col_cf[j_0]] = std::pair<T, unsigned int>(i_0, j_0);
                }
              else
                {
                  // face is visited the second time -> now we know the cells
                  // on both sides of the face and we can determine for both
                  // cells the neighbor
                  col_cc[j_0] = neighbors[col_cf[j_0]].first;
                  col_cc[neighbors[col_cf[j_0]].second] = i_0;
                }
            }
        }
    }



    /**
     * Build entities of dimension d (with 0<d<dim). Entities are described by
     * a set of vertices.
     *
     * Furthermore, the function determines for each cell of which d-dimensional
     * entity it consists of and its orientation relative to the cell.
     */
    template <int key_length, typename FU>
    void
    build_entity_templated(
      const unsigned int                                d,
      const std::vector<std::shared_ptr<CellTypeBase>> &cell_types,
      const std::vector<dealii::ReferenceCell> &        cell_types_index,
      const CRS<unsigned int> &                         crs,
      CRS<unsigned int> &                               crs_d,        // result
      CRS<unsigned int> &                               crs_0,        // result
      std::vector<unsigned char> &                      orientations, // result
      const FU &                                        second_key_function)
    {
      const bool compatibility_mode = true;

      const std::vector<std::size_t> & cell_ptr      = crs.ptr;
      const std::vector<unsigned int> &cell_vertices = crs.col;
      std::vector<std::size_t> &       ptr_d         = crs_d.ptr;
      std::vector<unsigned int> &      col_d         = crs_d.col;

      // note: we do not pre-allocate memory for these arrays because it turned
      // out that counting unique entities is more expensive than push_back().
      std::vector<std::size_t> & ptr_0 = crs_0.ptr;
      std::vector<unsigned int> &col_0 = crs_0.col;

      // clear
      ptr_0 = {};
      col_0 = {};

      unsigned int n_entities = 0;

      for (const auto &c : cell_types_index)
        n_entities +=
          cell_types[static_cast<types::geometric_entity_type>(c)]->n_entities(
            d);

      // step 1: store each d-dimensional entity of a cell (described by their
      // vertices) into a vector and create a key for them
      //
      // note: it turned out to be more efficient to have a vector of tuples
      // than to have two vectors (sorting becomes inefficient)
      std::vector<
        std::tuple<std::array<unsigned int, key_length>, unsigned int>>
        keys; // key (sorted vertices), cell-entity index

      std::vector<std::array<unsigned int, key_length>> ad_entity_vertices;
      std::vector<dealii::ReferenceCell>                ad_entity_types;
      std::vector<std::array<unsigned int, key_length>> ad_compatibility;

      keys.reserve(n_entities);
      ad_entity_vertices.reserve(n_entities);
      ad_entity_types.reserve(n_entities);
      ad_compatibility.reserve(n_entities);

      ptr_d.resize(cell_types_index.size() + 1);
      ptr_d[0] = 0;

      static const unsigned int offset = 1;

      // loop over all cells
      for (unsigned int c = 0, counter = 0; c < cell_types_index.size(); ++c)
        {
          const auto &cell_type =
            cell_types[static_cast<types::geometric_entity_type>(
              cell_types_index[c])];
          ptr_d[c + 1] = ptr_d[c] + cell_type->n_entities(d);

          // ... collect vertices of cell
          const dealii::ArrayView<const unsigned int> cell_vertice(
            cell_vertices.data() + cell_ptr[c], cell_ptr[c + 1] - cell_ptr[c]);

          // ... loop over all its entities
          for (unsigned int e = 0; e < cell_type->n_entities(d); ++e)
            {
              // ... determine global entity vertices
              const auto &local_entity_vertices =
                cell_type->vertices_of_entity(d, e);

              std::array<unsigned int, key_length> entity_vertices;
              std::fill(entity_vertices.begin(), entity_vertices.end(), 0);

              for (unsigned int i = 0; i < local_entity_vertices.size(); ++i)
                entity_vertices[i] =
                  cell_vertice[local_entity_vertices[i]] + offset;

              // ... create key
              std::array<unsigned int, key_length> key = entity_vertices;
              std::sort(key.begin(), key.end());
              keys.emplace_back(key, counter++);

              ad_entity_vertices.emplace_back(entity_vertices);

              ad_entity_types.emplace_back(cell_type->type_of_entity(d, e));

              if (compatibility_mode)
                ad_compatibility.emplace_back(
                  second_key_function(entity_vertices, cell_type, c, e));
            }
        }

      col_d.resize(keys.size());
      orientations.resize(keys.size());

      // step 2: sort according to key so that entities with same key can be
      // merged
      std::sort(keys.begin(), keys.end());


      if (compatibility_mode)
        {
          unsigned int n_unique_entities        = 0;
          unsigned int n_unique_entity_vertices = 0;

          std::array<unsigned int, key_length> ref_key, new_key;
          std::fill(ref_key.begin(), ref_key.end(), 0);
          for (unsigned int i = 0; i < keys.size(); ++i)
            {
              const auto offset_i = std::get<1>(keys[i]);

              if (ref_key != std::get<0>(keys[i]))
                {
                  ref_key = std::get<0>(keys[i]);

                  n_unique_entities++;
                  n_unique_entity_vertices +=
                    cell_types[static_cast<types::geometric_entity_type>(
                                 ad_entity_types[offset_i])]
                      ->n_entities(0);

                  new_key = ad_compatibility[offset_i];
                }

              std::get<0>(keys[i]) = new_key;
            }

          std::sort(keys.begin(), keys.end());

          ptr_0.reserve(n_unique_entities);
          col_0.reserve(n_unique_entity_vertices);
        }


      std::array<unsigned int, key_length> ref_key;
      std::array<unsigned int, key_length> ref_indices;
      std::fill(ref_key.begin(), ref_key.end(), 0);

      for (unsigned int i = 0, counter = dealii::numbers::invalid_unsigned_int;
           i < keys.size();
           i++)
        {
          const auto offset_i = std::get<1>(keys[i]);

          if (ref_key != std::get<0>(keys[i]))
            {
              // new key
              counter++;
              ref_key     = std::get<0>(keys[i]);
              ref_indices = ad_entity_vertices[offset_i];

              ptr_0.push_back(col_0.size());
              for (const auto j : ad_entity_vertices[offset_i])
                if (j != 0)
                  col_0.push_back(j - offset);

              // take its orientation as default
              col_d[offset_i]        = counter;
              orientations[offset_i] = 1;
            }
          else
            {
              col_d[offset_i] = counter;
              orientations[offset_i] =
                ad_entity_types[offset_i].compute_orientation(
                  ad_entity_vertices[offset_i], ref_indices);
            }
        }
      ptr_0.push_back(col_0.size());
    }



    /**
     * Call the right templated function to be able to use std::array instead
     * of std::vector.
     */
    template <typename FU>
    void
    build_entity(const unsigned int                                d,
                 const std::vector<std::shared_ptr<CellTypeBase>> &cell_types,
                 const std::vector<dealii::ReferenceCell> &cell_types_index,
                 const CRS<unsigned int> &                 crs,
                 CRS<unsigned int> &                       crs_d,
                 CRS<unsigned int> &                       crs_0,
                 std::vector<unsigned char> &              orientations,
                 const FU &                                second_key_function)
    {
      std::size_t key_length = 0;

      for (const auto &c : cell_types_index)
        {
          const auto &cell_type =
            cell_types[static_cast<types::geometric_entity_type>(c)];
          for (unsigned int e = 0; e < cell_type->n_entities(d); ++e)
            key_length =
              std::max(key_length, cell_type->vertices_of_entity(d, e).size());
        }

      if (key_length == 2)
        build_entity_templated<2>(d,
                                  cell_types,
                                  cell_types_index,
                                  crs,
                                  crs_d,
                                  crs_0,
                                  orientations,
                                  second_key_function);
      else if (key_length == 3)
        build_entity_templated<3>(d,
                                  cell_types,
                                  cell_types_index,
                                  crs,
                                  crs_d,
                                  crs_0,
                                  orientations,
                                  second_key_function);
      else if (key_length == 4)
        build_entity_templated<4>(d,
                                  cell_types,
                                  cell_types_index,
                                  crs,
                                  crs_d,
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
    build_intersection(
      const std::vector<std::shared_ptr<CellTypeBase>> &cell_types,
      const std::vector<dealii::ReferenceCell> &        cell_types_index,
      const CRS<unsigned int> &                         con_cv,
      const CRS<unsigned int> &                         con_cl,
      const CRS<unsigned int> &                         con_lv,
      const CRS<unsigned int> &                         con_cq,
      const CRS<unsigned int> &                         con_qv,
      const std::vector<unsigned char> &                ori_cq,
      CRS<unsigned int> &                               con_ql,   // result
      std::vector<unsigned char> &                      ori_ql,   // result
      std::vector<dealii::ReferenceCell> &              quad_t_id // result
    )
    {
      // reset output
      ori_ql     = {};
      con_ql.ptr = {};
      con_ql.col = {};

      con_ql.ptr.resize(con_qv.ptr.size());
      con_ql.ptr[0] = 0;

      quad_t_id.resize(con_qv.ptr.size() - 1);

      // count the number of lines of each face
      for (unsigned int c = 0; c < con_cq.ptr.size() - 1; ++c)
        {
          const auto &cell_type =
            cell_types[static_cast<types::geometric_entity_type>(
              cell_types_index[c])];

          // loop over faces
          for (unsigned int f_ = con_cq.ptr[c], f_index = 0;
               f_ < con_cq.ptr[c + 1];
               ++f_, ++f_index)
            {
              const unsigned int f = con_cq.col[f_];

              con_ql.ptr[f + 1] = cell_type->n_lines_of_surface(f_index);
            }
        }

      // use the counts to determine the offsets -> prefix sum
      for (unsigned int i = 0; i < con_ql.ptr.size() - 1; ++i)
        con_ql.ptr[i + 1] += con_ql.ptr[i];

      // allocate memory
      con_ql.col.resize(con_ql.ptr.back());
      ori_ql.resize(con_ql.ptr.back());

      // loop over cells
      for (unsigned int c = 0; c < con_cq.ptr.size() - 1; ++c)
        {
          const auto &cell_type =
            cell_types[static_cast<types::geometric_entity_type>(
              cell_types_index[c])];

          // loop over faces
          for (unsigned int f_ = con_cq.ptr[c], f_index = 0;
               f_ < con_cq.ptr[c + 1];
               ++f_, ++f_index)
            {
              const unsigned int f = con_cq.col[f_];

              // only faces with default orientation have to do something
              if (ori_cq[f_] != 1)
                continue;

              // determine entity type of face
              quad_t_id[f] = cell_type->type_of_entity(2, f_index);

              // loop over lines
              for (unsigned int l = 0;
                   l < cell_type->n_lines_of_surface(f_index);
                   ++l)
                {
                  // determine global index of line
                  const unsigned int local_line_index =
                    cell_type->nth_line_of_surface(l, f_index);
                  const unsigned int global_line_index =
                    con_cl.col[con_cl.ptr[c] + local_line_index];
                  con_ql.col[con_ql.ptr[f] + l] = global_line_index;

                  // determine orientation of line
                  const auto line_vertices_1_ref =
                    cell_type->vertices_of_nth_line_of_surface(l, f_index);

                  bool same = true;
                  for (unsigned int v = 0; v < line_vertices_1_ref.size(); ++v)
                    if (con_cv.col[con_cv.ptr[c] + line_vertices_1_ref[v]] !=
                        con_lv.col[con_lv.ptr[global_line_index] + v])
                      {
                        same = false;
                        break;
                      }

                  // ... comparison gives orientation
                  ori_ql[con_ql.ptr[f] + l] = (same ? 1 : 0);
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
    template <typename T>
    Connectivity<T>
    build_connectivity(const unsigned int                                dim,
                       const std::vector<std::shared_ptr<CellTypeBase>> &cell_t,
                       const std::vector<dealii::ReferenceCell> &cell_t_id,
                       const CRS<T> &                            con_cv)
    {
      Connectivity<T> connectivity(dim, cell_t_id);

      CRS<T> temp1; // needed for 3D

      if (dim == 1)
        connectivity.entity_to_entities(1, 0) = con_cv;

      if (dim == 2 || dim == 3) // build lines
        {
          std::vector<unsigned char> dummy;

          build_entity(1,
                       cell_t,
                       connectivity.entity_types(dim),
                       con_cv,
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
          build_entity(
            2,
            cell_t,
            connectivity.entity_types(3),
            con_cv,
            connectivity.entity_to_entities(3, 2),
            connectivity.entity_to_entities(2, 0),
            connectivity.entity_orientations(2),
            [&](auto key, const auto &cell_type, const auto &c, const auto &f) {
              //  to ensure same enumeration as in deal.II
              AssertIndexRange(cell_type->n_lines_of_surface(f),
                               key.size() + 1);

              unsigned int l = 0;

              for (; l < cell_type->n_lines_of_surface(f); ++l)
                key[l] =
                  temp1
                    .col[temp1.ptr[c] + cell_type->nth_line_of_surface(l, f)] +
                  1 /*offset!*/;

              for (; l < key.size(); ++l)
                key[l] = 0;

              return key;
            });

          // create connectivity: quad -> line
          build_intersection(cell_t,
                             connectivity.entity_types(3),
                             con_cv,
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
      determine_neighbors(connectivity.entity_to_entities(dim, dim - 1),
                          connectivity.entity_to_entities(dim, dim));

      return connectivity;
    }



    /**
     * Preprocessing step to remove the template argument dim.
     */
    template <typename T, int dim>
    Connectivity<T>
    build_connectivity(const std::vector<CellData<dim>> &cells)
    {
      AssertThrow(cells.size() > 0, ExcMessage("No cells have been provided!"));

      // vector of possible cell entity types
      std::vector<std::shared_ptr<CellTypeBase>> cell_types_impl(8);

      cell_types_impl[static_cast<types::geometric_entity_type>(
                        dealii::ReferenceCells::Line)]
        .reset(new CellTypeLine());
      cell_types_impl[static_cast<types::geometric_entity_type>(
                        dealii::ReferenceCells::Triangle)]
        .reset(new CellTypeTri());
      cell_types_impl[static_cast<types::geometric_entity_type>(
                        dealii::ReferenceCells::Quadrilateral)]
        .reset(new CellTypeQuad());
      cell_types_impl[static_cast<types::geometric_entity_type>(
                        dealii::ReferenceCells::Tetrahedron)]
        .reset(new CellTypeTet());
      cell_types_impl[static_cast<types::geometric_entity_type>(
                        dealii::ReferenceCells::Pyramid)]
        .reset(new CellTypePyramid());
      cell_types_impl[static_cast<types::geometric_entity_type>(
                        dealii::ReferenceCells::Wedge)]
        .reset(new CellTypeWedge());
      cell_types_impl[static_cast<types::geometric_entity_type>(
                        dealii::ReferenceCells::Hexahedron)]
        .reset(new CellTypeHex());

      // determine cell types and process vertices
      std::vector<T> cell_vertices;
      cell_vertices.reserve(
        std::accumulate(cells.begin(),
                        cells.end(),
                        0,
                        [](const auto &result, const auto &cell) {
                          return result + cell.vertices.size();
                        }));

      std::vector<std::size_t> cell_vertices_ptr;
      cell_vertices_ptr.reserve(cells.size() + 1);
      cell_vertices_ptr.push_back(0);

      std::vector<dealii::ReferenceCell> cell_types_indices;
      cell_types_indices.reserve(cells.size());

      // loop over cells and create CRS
      for (const auto &cell : cells)
        {
#ifdef DEBUG
          auto vertices_unique = cell.vertices;
          std::sort(vertices_unique.begin(), vertices_unique.end());
          vertices_unique.erase(std::unique(vertices_unique.begin(),
                                            vertices_unique.end()),
                                vertices_unique.end());

          Assert(vertices_unique.size() == cell.vertices.size(),
                 ExcMessage(
                   "The definition of a cell refers to the same vertex several "
                   "times. This is not possible. A common reason is that "
                   "CellData::vertices has a size that does not match the "
                   "size expected from the reference cell. Please resize "
                   "CellData::vertices or use the appropriate constructor of "
                   "CellData."));
#endif

          const dealii::ReferenceCell reference_cell =
            dealii::ReferenceCell::n_vertices_to_type(dim,
                                                      cell.vertices.size());

          Assert(reference_cell != dealii::ReferenceCells::Invalid,
                 ExcNotImplemented());
          AssertIndexRange(static_cast<types::geometric_entity_type>(
                             reference_cell),
                           cell_types_impl.size());
          Assert(cell_types_impl[static_cast<types::geometric_entity_type>(
                                   reference_cell)]
                     .get() != nullptr,
                 ExcNotImplemented());

          cell_types_indices.push_back(reference_cell);

          // create CRS of vertices (to remove template argument dim)
          for (const auto &vertex : cell.vertices)
            cell_vertices.push_back(vertex);

          cell_vertices_ptr.push_back(cell_vertices.size());
        }

      // do the actual work
      return build_connectivity<T>(dim,
                                   cell_types_impl,
                                   cell_types_indices,
                                   {cell_vertices_ptr, cell_vertices});
    }
  } // namespace TriangulationImplementation
} // namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif

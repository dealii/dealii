// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

#ifndef dealii_tria_reference_cell_h
#define dealii_tria_reference_cell_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>


DEAL_II_NAMESPACE_OPEN


/**
 * A namespace for reference cells.
 */
namespace ReferenceCell
{
  /**
   * Supported reference cell types.
   */
  enum class Type : std::uint8_t
  {
    Vertex  = 0,
    Line    = 1,
    Tri     = 2,
    Quad    = 3,
    Tet     = 4,
    Pyramid = 5,
    Wedge   = 6,
    Hex     = 7,
    Invalid = static_cast<std::uint8_t>(-1)
  };

  /**
   * Return the correct simplex reference cell type for the given dimension
   * @p dim.
   */
  inline Type
  get_simplex(const unsigned int dim)
  {
    switch (dim)
      {
        case 0:
          return Type::Vertex;
        case 1:
          return Type::Line;
        case 2:
          return Type::Tri;
        case 3:
          return Type::Tet;
        default:
          Assert(false, ExcNotImplemented());
          return Type::Invalid;
      }
  }

  /**
   * Return the correct hypercube reference cell type for the given dimension
   * @p dim.
   */
  inline Type
  get_hypercube(const unsigned int dim)
  {
    switch (dim)
      {
        case 0:
          return Type::Vertex;
        case 1:
          return Type::Line;
        case 2:
          return Type::Quad;
        case 3:
          return Type::Hex;
        default:
          Assert(false, ExcNotImplemented());
          return Type::Invalid;
      }
  }

  namespace internal
  {
    /**
     * Check if the bit at position @p n in @p number is set.
     */
    inline static bool
    get_bit(const unsigned char number, const unsigned int n)
    {
      AssertIndexRange(n, 8);

      // source:
      // https://stackoverflow.com/questions/47981/how-do-you-set-clear-and-toggle-a-single-bit
      // "Checking a bit"
      return (number >> n) & 1U;
    }



    /**
     * Set the bit at position @p n in @p number to value @p x.
     */
    inline static void
    set_bit(unsigned char &number, const unsigned int n, const bool x)
    {
      AssertIndexRange(n, 8);

      // source:
      // https://stackoverflow.com/questions/47981/how-do-you-set-clear-and-toggle-a-single-bit
      // "Changing the nth bit to x"
      number ^= (-static_cast<unsigned char>(x) ^ number) & (1U << n);
    }

    /**
     * A namespace for geometric information on reference cells.
     */
    namespace Info
    {
      /**
       * Interface to be used in TriaAccessor/TriaCellAccessor to access
       * sub-entities of dimension d' of geometric entities of dimension d, with
       * 0<=d'<d<=3.
       */
      struct Base
      {
        /**
         * Destructor.
         */
        virtual ~Base() = default;

        /**
         * Number of vertices.
         */
        virtual unsigned int
        n_vertices() const
        {
          Assert(false, ExcNotImplemented());
          return 0;
        }

        /**
         * Number of lines.
         */
        virtual unsigned int
        n_lines() const
        {
          Assert(false, ExcNotImplemented());
          return 0;
        }


        /**
         * Number of faces.
         */
        virtual unsigned int
        n_faces() const
        {
          Assert(false, ExcNotImplemented());
          return 0;
        }

        /**
         * Return an object that can be thought of as an array containing all
         * indices from zero to n_vertices().
         */
        inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
        vertex_indices() const
        {
          return {0U, n_vertices()};
        }

        /**
         * Return an object that can be thought of as an array containing all
         * indices from zero to n_lines().
         */
        inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
        line_indices() const
        {
          return {0U, n_lines()};
        }

        /**
         * Return an object that can be thought of as an array containing all
         * indices from zero to n_faces().
         */
        inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
        face_indices() const
        {
          return {0U, n_faces()};
        }

        /**
         * Standard decomposition of vertex index into face and face-vertex
         * index.
         */
        virtual std::array<unsigned int, 2>
        standard_vertex_to_face_and_vertex_index(
          const unsigned int vertex) const
        {
          Assert(false, ExcNotImplemented());

          (void)vertex;

          return {{0u, 0u}};
        }

        /**
         * Standard decomposition of line index into face and face-line index.
         */
        virtual std::array<unsigned int, 2>
        standard_line_to_face_and_line_index(const unsigned int line) const
        {
          Assert(false, ExcNotImplemented());

          (void)line;

          return {{0, 0}};
        }

        /**
         * Correct vertex index depending on face orientation.
         */
        virtual unsigned int
        standard_to_real_face_vertex(const unsigned int  vertex,
                                     const unsigned int  face,
                                     const unsigned char face_orientation) const
        {
          Assert(false, ExcNotImplemented());

          (void)vertex;
          (void)face;
          (void)face_orientation;

          return 0;
        }

        /**
         * Correct line index depending on face orientation.
         */
        virtual unsigned int
        standard_to_real_face_line(const unsigned int  line,
                                   const unsigned int  face,
                                   const unsigned char face_orientation) const
        {
          Assert(false, ExcNotImplemented());

          (void)line;
          (void)face;
          (void)face_orientation;

          return 0;
        }

        /**
         * Combine face and line orientation.
         */
        virtual bool
        combine_face_and_line_orientation(
          const unsigned int  line,
          const unsigned char face_orientation,
          const unsigned char line_orientation) const
        {
          Assert(false, ExcNotImplemented());

          (void)line;
          (void)face_orientation;
          (void)line_orientation;

          return true;
        }

        /**
         * Return reference-cell type of face @p face_no.
         */
        virtual ReferenceCell::Type
        face_reference_cell_type(const unsigned int face_no) const
        {
          Assert(false, ExcNotImplemented());
          (void)face_no;

          return ReferenceCell::Type::Invalid;
        }

        /**
         * Map face line number to cell line number.
         */
        virtual unsigned int
        face_to_cell_lines(const unsigned int  face,
                           const unsigned int  line,
                           const unsigned char face_orientation) const
        {
          Assert(false, ExcNotImplemented());
          (void)face;
          (void)line;
          (void)face_orientation;

          return 0;
        }

        /**
         * Map face vertex number to cell vertex number.
         */
        virtual unsigned int
        face_to_cell_vertices(const unsigned int  face,
                              const unsigned int  vertex,
                              const unsigned char face_orientation) const
        {
          Assert(false, ExcNotImplemented());
          (void)face;
          (void)vertex;
          (void)face_orientation;

          return 0;
        }
      };


      /**
       * Base class for tensor-product geometric entities.
       */
      template <int dim>
      struct TensorProductBase : Base
      {
        unsigned int
        n_vertices() const override
        {
          return GeometryInfo<dim>::vertices_per_cell;
        }

        unsigned int
        n_lines() const override
        {
          return GeometryInfo<dim>::lines_per_cell;
        }

        unsigned int
        n_faces() const override
        {
          return GeometryInfo<dim>::faces_per_cell;
        }

        unsigned int
        face_to_cell_lines(const unsigned int  face,
                           const unsigned int  line,
                           const unsigned char face_orientation) const override
        {
          return GeometryInfo<dim>::face_to_cell_lines(
            face,
            line,
            get_bit(face_orientation, 0),
            get_bit(face_orientation, 2),
            get_bit(face_orientation, 1));
        }

        unsigned int
        face_to_cell_vertices(
          const unsigned int  face,
          const unsigned int  vertex,
          const unsigned char face_orientation) const override
        {
          return GeometryInfo<dim>::face_to_cell_vertices(
            face,
            vertex,
            get_bit(face_orientation, 0),
            get_bit(face_orientation, 2),
            get_bit(face_orientation, 1));
        }
      };



      /*
       * Vertex.
       */
      struct Vertex : public TensorProductBase<0>
      {
        ReferenceCell::Type
        face_reference_cell_type(const unsigned int face_no) const override
        {
          (void)face_no;
          return ReferenceCell::Type::Invalid;
        }
      };



      /*
       * Line.
       */
      struct Line : public TensorProductBase<1>
      {
        ReferenceCell::Type
        face_reference_cell_type(const unsigned int face_no) const override
        {
          (void)face_no;
          return ReferenceCell::Type::Vertex;
        }
      };



      /**
       * Triangle.
       */
      struct Tri : public Base
      {
        unsigned int
        n_vertices() const override
        {
          return 3;
        }

        unsigned int
        n_lines() const override
        {
          return 3;
        }

        unsigned int
        n_faces() const override
        {
          return this->n_lines();
        }

        std::array<unsigned int, 2>
        standard_vertex_to_face_and_vertex_index(
          const unsigned int vertex) const override
        {
          AssertIndexRange(vertex, 3);

          static const std::array<std::array<unsigned int, 2>, 3> table = {
            {{{0, 0}}, {{0, 1}}, {{1, 1}}}};

          return table[vertex];
        }

        unsigned int
        standard_to_real_face_vertex(
          const unsigned int  vertex,
          const unsigned int  face,
          const unsigned char line_orientation) const override
        {
          (void)face;

          static const std::array<std::array<unsigned int, 2>, 2> table = {
            {{{1, 0}}, {{0, 1}}}};

          return table[line_orientation][vertex];
        }

        ReferenceCell::Type
        face_reference_cell_type(const unsigned int face_no) const override
        {
          (void)face_no;

          AssertIndexRange(face_no, n_faces());

          return ReferenceCell::Type::Line;
        }

        unsigned int
        face_to_cell_lines(const unsigned int  face,
                           const unsigned int  line,
                           const unsigned char face_orientation) const override
        {
          AssertIndexRange(face, n_faces());
          AssertDimension(line, 0);

          (void)line;
          (void)face_orientation;

          return face;
        }

        unsigned int
        face_to_cell_vertices(
          const unsigned int  face,
          const unsigned int  vertex,
          const unsigned char face_orientation) const override
        {
          static const std::array<std::array<unsigned int, 2>, 3> table = {
            {{0, 1}, {1, 2}, {2, 0}}};

          return table[face][face_orientation ? vertex : (1 - vertex)];
        }
      };



      /**
       * Quad.
       */
      struct Quad : public TensorProductBase<2>
      {
        std::array<unsigned int, 2>
        standard_vertex_to_face_and_vertex_index(
          const unsigned int vertex) const override
        {
          return GeometryInfo<2>::standard_quad_vertex_to_line_vertex_index(
            vertex);
        }

        unsigned int
        standard_to_real_face_vertex(
          const unsigned int  vertex,
          const unsigned int  face,
          const unsigned char line_orientation) const override
        {
          (void)face;

          return GeometryInfo<2>::standard_to_real_line_vertex(
            vertex, line_orientation);
        }

        ReferenceCell::Type
        face_reference_cell_type(const unsigned int face_no) const override
        {
          (void)face_no;
          return ReferenceCell::Type::Line;
        }
      };



      /**
       * Tet.
       */
      struct Tet : public TensorProductBase<3>
      {
        unsigned int
        n_vertices() const override
        {
          return 4;
        }

        unsigned int
        n_lines() const override
        {
          return 6;
        }

        unsigned int
        n_faces() const override
        {
          return 4;
        }

        std::array<unsigned int, 2>
        standard_line_to_face_and_line_index(
          const unsigned int line) const override
        {
          static const std::array<unsigned int, 2> table[6] = {
            {{0, 0}}, {{0, 1}}, {{0, 2}}, {{1, 1}}, {{1, 2}}, {{2, 1}}};

          return table[line];
        }

        unsigned int
        standard_to_real_face_line(
          const unsigned int  line,
          const unsigned int  face,
          const unsigned char face_orientation) const override
        {
          (void)face;

          static const std::array<std::array<unsigned int, 3>, 6> table = {
            {{{2, 1, 0}},
             {{0, 1, 2}},
             {{1, 2, 0}},
             {{0, 2, 1}},
             {{1, 0, 2}},
             {{2, 0, 1}}}};

          return table[face_orientation][line];
        }

        bool
        combine_face_and_line_orientation(
          const unsigned int  line,
          const unsigned char face_orientation_raw,
          const unsigned char line_orientation) const override
        {
          (void)line;
          (void)face_orientation_raw;

          return line_orientation;
        }

        std::array<unsigned int, 2>
        standard_vertex_to_face_and_vertex_index(
          const unsigned int vertex) const override
        {
          AssertIndexRange(vertex, 4);

          static const std::array<unsigned int, 2> table[4] = {{{0, 0}},
                                                               {{0, 1}},
                                                               {{0, 2}},
                                                               {{1, 2}}};

          return table[vertex];
        }

        unsigned int
        standard_to_real_face_vertex(
          const unsigned int  vertex,
          const unsigned int  face,
          const unsigned char face_orientation) const override
        {
          AssertIndexRange(face_orientation, 6);
          (void)face;

          static const std::array<std::array<unsigned int, 3>, 6> table = {
            {{{0, 2, 1}},
             {{0, 1, 2}},
             {{1, 2, 0}},
             {{1, 0, 2}},
             {{2, 1, 0}},
             {{2, 0, 1}}}};

          return table[face_orientation][vertex];
        }

        ReferenceCell::Type
        face_reference_cell_type(const unsigned int face_no) const override
        {
          (void)face_no;

          AssertIndexRange(face_no, n_faces());

          return ReferenceCell::Type::Tri;
        }

        unsigned int
        face_to_cell_lines(const unsigned int  face,
                           const unsigned int  line,
                           const unsigned char face_orientation) const override
        {
          AssertIndexRange(face, n_faces());

          const static std::array<std::array<unsigned int, 3>, 4> table = {
            {{0, 1, 2}, {0, 3, 4}, {2, 5, 3}, {1, 4, 5}}};

          return table[face][standard_to_real_face_line(
            line, face, face_orientation)];
        }

        unsigned int
        face_to_cell_vertices(
          const unsigned int  face,
          const unsigned int  vertex,
          const unsigned char face_orientation) const override
        {
          static const std::array<std::array<unsigned int, 3>, 4> table = {
            {{0, 1, 2}, {1, 0, 3}, {0, 2, 3}, {2, 1, 3}}};

          return table[face][standard_to_real_face_vertex(
            vertex, face, face_orientation)];
        }
      };



      /**
       * Pyramid.
       */
      struct Pyramid : public TensorProductBase<3>
      {
        unsigned int
        n_vertices() const override
        {
          return 5;
        }

        unsigned int
        n_lines() const override
        {
          return 8;
        }

        unsigned int
        n_faces() const override
        {
          return 5;
        }

        std::array<unsigned int, 2>
        standard_line_to_face_and_line_index(
          const unsigned int line) const override
        {
          Assert(false, ExcNotImplemented());

          static const std::array<unsigned int, 2> table[6] = {
            {{0, 0}}, {{0, 1}}, {{0, 2}}, {{1, 1}}, {{1, 2}}, {{2, 1}}};

          return table[line];
        }

        unsigned int
        standard_to_real_face_line(
          const unsigned int  line,
          const unsigned int  face,
          const unsigned char face_orientation) const override
        {
          Assert(false, ExcNotImplemented());

          (void)face;

          static const std::array<std::array<unsigned int, 3>, 6> table = {
            {{{2, 1, 0}},
             {{0, 1, 2}},
             {{1, 2, 0}},
             {{0, 2, 1}},
             {{1, 0, 2}},
             {{2, 0, 1}}}};

          return table[face_orientation][line];
        }

        bool
        combine_face_and_line_orientation(
          const unsigned int  line,
          const unsigned char face_orientation_raw,
          const unsigned char line_orientation) const override
        {
          (void)line;
          (void)face_orientation_raw;

          return line_orientation;
        }

        std::array<unsigned int, 2>
        standard_vertex_to_face_and_vertex_index(
          const unsigned int vertex) const override
        {
          static const std::array<unsigned int, 2> table[5] = {
            {{0, 0}}, {{0, 1}}, {{0, 2}}, {{0, 3}}, {{1, 2}}};

          return table[vertex];
        }

        unsigned int
        standard_to_real_face_vertex(
          const unsigned int  vertex,
          const unsigned int  face,
          const unsigned char face_orientation) const override
        {
          if (face == 0) // Quad
            {
              return GeometryInfo<3>::standard_to_real_face_vertex(
                vertex,
                get_bit(face_orientation, 0),
                get_bit(face_orientation, 2),
                get_bit(face_orientation, 1));
            }
          else // Tri
            {
              static const std::array<std::array<unsigned int, 3>, 6> table = {
                {{{0, 2, 1}},
                 {{0, 1, 2}},
                 {{1, 2, 0}},
                 {{1, 0, 2}},
                 {{2, 1, 0}},
                 {{2, 0, 1}}}};

              return table[face_orientation][vertex];
            }
        }

        ReferenceCell::Type
        face_reference_cell_type(const unsigned int face_no) const override
        {
          AssertIndexRange(face_no, n_faces());

          if (face_no == 1)
            return ReferenceCell::Type::Quad;
          else
            return ReferenceCell::Type::Tri;
        }
      };



      /**
       * Wedge.
       */
      struct Wedge : public TensorProductBase<3>
      {
        unsigned int
        n_vertices() const override
        {
          return 6;
        }

        unsigned int
        n_lines() const override
        {
          return 9;
        }

        unsigned int
        n_faces() const override
        {
          return 6;
        }

        std::array<unsigned int, 2>
        standard_line_to_face_and_line_index(
          const unsigned int line) const override
        {
          Assert(false, ExcNotImplemented());

          static const std::array<unsigned int, 2> table[6] = {
            {{0, 0}}, {{0, 1}}, {{0, 2}}, {{1, 1}}, {{1, 2}}, {{2, 1}}};

          return table[line];
        }

        unsigned int
        standard_to_real_face_line(
          const unsigned int  line,
          const unsigned int  face,
          const unsigned char face_orientation) const override
        {
          Assert(false, ExcNotImplemented());

          (void)face;

          static const std::array<std::array<unsigned int, 3>, 6> table = {
            {{{2, 1, 0}},
             {{0, 1, 2}},
             {{1, 2, 0}},
             {{0, 2, 1}},
             {{1, 0, 2}},
             {{2, 0, 1}}}};

          return table[face_orientation][line];
        }

        bool
        combine_face_and_line_orientation(
          const unsigned int  line,
          const unsigned char face_orientation_raw,
          const unsigned char line_orientation) const override
        {
          (void)line;
          (void)face_orientation_raw;

          return line_orientation;
        }

        std::array<unsigned int, 2>
        standard_vertex_to_face_and_vertex_index(
          const unsigned int vertex) const override
        {
          static const std::array<std::array<unsigned int, 2>, 6> table = {
            {{{0, 1}}, {{0, 0}}, {{0, 2}}, {{1, 0}}, {{1, 1}}, {{1, 2}}}};

          return table[vertex];
        }

        unsigned int
        standard_to_real_face_vertex(
          const unsigned int  vertex,
          const unsigned int  face,
          const unsigned char face_orientation) const override
        {
          if (face > 1) // QUAD
            {
              return GeometryInfo<3>::standard_to_real_face_vertex(
                vertex,
                get_bit(face_orientation, 0),
                get_bit(face_orientation, 2),
                get_bit(face_orientation, 1));
            }
          else // TRI
            {
              static const std::array<std::array<unsigned int, 3>, 6> table = {
                {{{0, 2, 1}},
                 {{0, 1, 2}},
                 {{1, 2, 0}},
                 {{1, 0, 2}},
                 {{2, 1, 0}},
                 {{2, 0, 1}}}};

              return table[face_orientation][vertex];
            }
        }

        ReferenceCell::Type
        face_reference_cell_type(const unsigned int face_no) const override
        {
          AssertIndexRange(face_no, n_faces());

          if (face_no > 1)
            return ReferenceCell::Type::Quad;
          else
            return ReferenceCell::Type::Tri;
        }
      };



      /**
       * Hex.
       */
      struct Hex : public TensorProductBase<3>
      {
        std::array<unsigned int, 2>
        standard_line_to_face_and_line_index(
          const unsigned int line) const override
        {
          return GeometryInfo<3>::standard_hex_line_to_quad_line_index(line);
        }

        unsigned int
        standard_to_real_face_line(
          const unsigned int  line,
          const unsigned int  face,
          const unsigned char face_orientation) const override
        {
          (void)face;

          return GeometryInfo<3>::standard_to_real_face_line(
            line,
            get_bit(face_orientation, 0),
            get_bit(face_orientation, 2),
            get_bit(face_orientation, 1));
        }

        bool
        combine_face_and_line_orientation(
          const unsigned int  line,
          const unsigned char face_orientation_raw,
          const unsigned char line_orientation) const override
        {
          static const bool bool_table[2][2][2][2] = {
            {{{true, false},    // lines 0/1, face_orientation=false,
                                // face_flip=false, face_rotation=false and true
              {false, true}},   // lines 0/1, face_orientation=false,
                                // face_flip=true, face_rotation=false and true
             {{true, true},     // lines 0/1, face_orientation=true,
                                // face_flip=false, face_rotation=false and true
              {false, false}}}, // lines 0/1, face_orientation=true,
                                // face_flip=true, face_rotation=false and true

            {{{true, true}, // lines 2/3 ...
              {false, false}},
             {{true, false}, {false, true}}}};

          const bool face_orientation = get_bit(face_orientation_raw, 0);
          const bool face_flip        = get_bit(face_orientation_raw, 2);
          const bool face_rotation    = get_bit(face_orientation_raw, 1);

          return (
            static_cast<bool>(line_orientation) ==
            bool_table[line / 2][face_orientation][face_flip][face_rotation]);
        }

        std::array<unsigned int, 2>
        standard_vertex_to_face_and_vertex_index(
          const unsigned int vertex) const override
        {
          return GeometryInfo<3>::standard_hex_vertex_to_quad_vertex_index(
            vertex);
        }

        unsigned int
        standard_to_real_face_vertex(
          const unsigned int  vertex,
          const unsigned int  face,
          const unsigned char face_orientation) const override
        {
          (void)face;

          return GeometryInfo<3>::standard_to_real_face_vertex(
            vertex,
            get_bit(face_orientation, 0),
            get_bit(face_orientation, 2),
            get_bit(face_orientation, 1));
        }

        ReferenceCell::Type
        face_reference_cell_type(const unsigned int face_no) const override
        {
          (void)face_no;
          return ReferenceCell::Type::Quad;
        }
      };

      /**
       * Return for a given reference-cell type @p the right Info.
       */
      inline const ReferenceCell::internal::Info::Base &
      get_cell(const ReferenceCell::Type &type)
      {
        static ReferenceCell::internal::Info::Base    gei_invalid;
        static ReferenceCell::internal::Info::Vertex  gei_vertex;
        static ReferenceCell::internal::Info::Line    gei_line;
        static ReferenceCell::internal::Info::Tri     gei_tri;
        static ReferenceCell::internal::Info::Quad    gei_quad;
        static ReferenceCell::internal::Info::Tet     gei_tet;
        static ReferenceCell::internal::Info::Pyramid gei_pyramid;
        static ReferenceCell::internal::Info::Wedge   gei_wedge;
        static ReferenceCell::internal::Info::Hex     gei_hex;

        switch (type)
          {
            case ReferenceCell::Type::Vertex:
              return gei_vertex;
            case ReferenceCell::Type::Line:
              return gei_line;
            case ReferenceCell::Type::Tri:
              return gei_tri;
            case ReferenceCell::Type::Quad:
              return gei_quad;
            case ReferenceCell::Type::Tet:
              return gei_tet;
            case ReferenceCell::Type::Pyramid:
              return gei_pyramid;
            case ReferenceCell::Type::Wedge:
              return gei_wedge;
            case ReferenceCell::Type::Hex:
              return gei_hex;
            default:
              Assert(false, StandardExceptions::ExcNotImplemented());
              return gei_invalid;
          }
      }

      /**
       * Return for a given reference-cell type @p and face number @p face_no the
       * right Info of the @p face_no-th face.
       */
      inline const ReferenceCell::internal::Info::Base &
      get_face(const ReferenceCell::Type &type, const unsigned int face_no)
      {
        return get_cell(get_cell(type).face_reference_cell_type(face_no));
      }

    } // namespace Info
  }   // namespace internal
} // namespace ReferenceCell


DEAL_II_NAMESPACE_CLOSE

#endif

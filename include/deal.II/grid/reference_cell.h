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

#ifndef dealii_tria_referce_cell_h
#define dealii_tria_referce_cell_h

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
      };



      /*
       * Vertex.
       */
      struct Vertex : public TensorProductBase<0>
      {};



      /*
       * Line.
       */
      struct Line : public TensorProductBase<1>
      {};



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
      };

    } // namespace Info
  }   // namespace internal
} // namespace ReferenceCell


DEAL_II_NAMESPACE_CLOSE

#endif

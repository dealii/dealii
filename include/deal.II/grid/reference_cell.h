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

#include <deal.II/base/array_view.h>
#include <deal.II/base/geometry_info.h>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
class Triangulation;
template <int dim, int spacedim>
class Mapping;
template <int dim>
class Quadrature;
template <int dim>
class ScalarPolynomialsBase;
#endif

/**
 * A namespace for reference cells.
 */
namespace ReferenceCell
{
  /**
   * A type that describes the kinds of reference cells that can be used.
   * This includes quadrilaterals and hexahedra (i.e., "hypercubes"),
   * triangles and tetrahedra (simplices), and the pyramids and wedges
   * necessary when using mixed 3d meshes.
   */
  class Type
  {
  public:
    enum CellKinds : std::uint8_t
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
     * Default constructor. Initialize this object as an invalid object.
     */
    Type();

    /**
     * Constructor.
     */
    Type(const CellKinds kind);

    /**
     * Conversion operator to an integer.
     */
    operator std::uint8_t() const;

    /**
     * Operator for equality comparison.
     */
    bool
    operator==(const Type &type) const;

    /**
     * Operator for inequality comparison.
     */
    bool
    operator!=(const Type &type) const;

    /**
     * Operator for equality comparison.
     */
    bool
    operator==(const CellKinds &type) const;

    /**
     * Operator for inequality comparison.
     */
    bool
    operator!=(const CellKinds &type) const;

    /**
     * Write and read the data of this object from a stream for the purpose
     * of serialization using the [BOOST serialization
     * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
     */
    template <class Archive>
    void
    serialize(Archive &archive, const unsigned int /*version*/);

    /**
     * Return a vector of faces a @p vertex belongs to.
     */
    ArrayView<const unsigned int>
    faces_for_given_vertex(const unsigned int vertex) const;

  private:
    /**
     * The variable that stores what this object actually corresponds to.
     */
    CellKinds kind;
  };



  inline Type::Type()
    : Type(Invalid)
  {}



  inline Type::Type(const CellKinds kind)
    : kind(kind)
  {}



  inline Type::operator std::uint8_t() const
  {
    return kind;
  }



  inline bool
  Type::operator==(const Type &type) const
  {
    return kind == type.kind;
  }



  inline bool
  Type::operator!=(const Type &type) const
  {
    return kind != type.kind;
  }



  inline bool
  Type::operator==(const CellKinds &type) const
  {
    return kind == type;
  }



  inline bool
  Type::operator!=(const CellKinds &type) const
  {
    return kind != type;
  }



  template <class Archive>
  inline void
  Type::serialize(Archive &archive, const unsigned int /*version*/)
  {
    // Serialize the state as an 8-bit int. When saving the state, the
    // last of the following 3 lines is a no-op. When loading, the first
    // of these lines is a no-op.
    std::uint8_t kind_as_int = static_cast<std::uint8_t>(kind);
    archive &    kind_as_int;
    kind = static_cast<CellKinds>(kind_as_int);
  }



  inline ArrayView<const unsigned int>
  Type::faces_for_given_vertex(const unsigned int vertex) const
  {
    if (kind == Type::Line)
      {
        AssertIndexRange(vertex, GeometryInfo<1>::vertices_per_cell);
        return {&GeometryInfo<2>::vertex_to_face[vertex][0], 1};
      }
    else if (kind == Type::Quad)
      {
        AssertIndexRange(vertex, GeometryInfo<2>::vertices_per_cell);
        return {&GeometryInfo<2>::vertex_to_face[vertex][0], 2};
      }
    else if (kind == Type::Hex)
      {
        AssertIndexRange(vertex, GeometryInfo<3>::vertices_per_cell);
        return {&GeometryInfo<3>::vertex_to_face[vertex][0], 3};
      }
    else if (kind == Type::Tri)
      {
        AssertIndexRange(vertex, 3);
        static const std::array<std::array<unsigned int, 2>, 3> table = {
          {{{0, 2}}, {{0, 1}}, {{1, 2}}}};

        return table[vertex];
      }
    else if (kind == Type::Tet)
      {
        AssertIndexRange(vertex, 4);
        static const std::array<std::array<unsigned int, 3>, 4> table = {
          {{{0, 1, 2}}, {{0, 1, 3}}, {{0, 2, 3}}, {{1, 2, 3}}}};

        return table[vertex];
      }
    else if (kind == Type::Wedge)
      {
        AssertIndexRange(vertex, 6);
        static const std::array<std::array<unsigned int, 3>, 6> table = {
          {{{0, 2, 4}},
           {{0, 2, 3}},
           {{0, 3, 4}},
           {{1, 2, 4}},
           {{1, 2, 3}},
           {{1, 3, 4}}}};

        return table[vertex];
      }
    else if (kind == Type::Pyramid)
      {
        AssertIndexRange(vertex, 5);
        static const unsigned int X = numbers::invalid_unsigned_int;
        static const std::array<std::array<unsigned int, 4>, 5> table = {
          {{{0, 1, 3, X}},
           {{0, 2, 3, X}},
           {{0, 1, 4, X}},
           {{0, 2, 4, X}},
           {{1, 2, 3, 4}}}};

        return {&table[vertex][0], vertex == 4 ? 4u : 3u};
      }

    Assert(false, ExcNotImplemented());

    return {};
  }



  /**
   * Return the dimension of the given reference-cell type @p type.
   */
  inline unsigned int
  get_dimension(const Type &type)
  {
    switch (type)
      {
        case Type::Vertex:
          return 0;
        case Type::Line:
          return 1;
        case Type::Tri:
        case Type::Quad:
          return 2;
        case Type::Tet:
        case Type::Pyramid:
        case Type::Wedge:
        case Type::Hex:
          return 3;
        default:
          return numbers::invalid_unsigned_int;
      }
  }

  /**
   * Convert the given reference cell type to a string.
   */
  inline std::string
  to_string(const Type &type)
  {
    switch (type)
      {
        case Type::Vertex:
          return "Vertex";
        case Type::Line:
          return "Line";
        case Type::Tri:
          return "Tri";
        case Type::Quad:
          return "Quad";
        case Type::Tet:
          return "Tet";
        case Type::Pyramid:
          return "Pyramid";
        case Type::Wedge:
          return "Wedge";
        case Type::Hex:
          return "Hex";
        case Type::Invalid:
          return "Invalid";
        default:
          Assert(false, ExcNotImplemented());
      }

    return "Invalid";
  }

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

  /**
   * Retrieve the correct ReferenceCell::Type for a given structural dimension
   * and number of vertices.
   */
  inline Type
  n_vertices_to_type(const int dim, const unsigned int n_vertices)
  {
    AssertIndexRange(dim, 4);
    AssertIndexRange(n_vertices, 9);
    const auto X = Type::Invalid;

    static constexpr std::array<std::array<ReferenceCell::Type::CellKinds, 9>,
                                4>
      table = {
        {// dim 0
         {{X, Type::Vertex, X, X, X, X, X, X, X}},
         // dim 1
         {{X, X, Type::Line, X, X, X, X, X, X}},
         // dim 2
         {{X, X, X, Type::Tri, Type::Quad, X, X, X, X}},
         // dim 3
         {{X, X, X, X, Type::Tet, Type::Pyramid, Type::Wedge, X, Type::Hex}}}};
    Assert(table[dim][n_vertices] != Type::Invalid,
           ExcMessage("The combination of dim = " + std::to_string(dim) +
                      " and n_vertices = " + std::to_string(n_vertices) +
                      " does not correspond to a known reference cell type."));
    return table[dim][n_vertices];
  }


  /**
   * Compute the value of the $i$-th linear shape function at location $\xi$ for
   * a given reference-cell type.
   */
  template <int dim>
  inline double
  d_linear_shape_function(const Type &       reference_cell,
                          const Point<dim> & xi,
                          const unsigned int i)
  {
    if (reference_cell == get_hypercube(dim))
      return GeometryInfo<dim>::d_linear_shape_function(xi, i);

    if (reference_cell ==
        Type::Tri) // see also Simplex::ScalarPolynomial::compute_value
      {
        switch (i)
          {
            case 0:
              return 1.0 - xi[std::min(0, dim - 1)] - xi[std::min(1, dim - 1)];
            case 1:
              return xi[std::min(0, dim - 1)];
            case 2:
              return xi[std::min(1, dim - 1)];
          }
      }

    if (reference_cell ==
        Type::Tet) // see also Simplex::ScalarPolynomial::compute_value
      {
        switch (i)
          {
            case 0:
              return 1.0 - xi[std::min(0, dim - 1)] - xi[std::min(1, dim - 1)] -
                     xi[std::min(2, dim - 1)];
            case 1:
              return xi[std::min(0, dim - 1)];
            case 2:
              return xi[std::min(1, dim - 1)];
            case 3:
              return xi[std::min(2, dim - 1)];
          }
      }

    if (reference_cell ==
        Type::Wedge) // see also Simplex::ScalarWedgePolynomial::compute_value
      {
        return d_linear_shape_function<2>(Type::Tri,
                                          Point<2>(xi[std::min(0, dim - 1)],
                                                   xi[std::min(1, dim - 1)]),
                                          i % 3) *
               d_linear_shape_function<1>(Type::Line,
                                          Point<1>(xi[std::min(2, dim - 1)]),
                                          i / 3);
      }

    if (reference_cell ==
        Type::Pyramid) // see also
                       // Simplex::ScalarPyramidPolynomial::compute_value
      {
        const double Q14 = 0.25;
        double       ration;

        const double r = xi[std::min(0, dim - 1)];
        const double s = xi[std::min(1, dim - 1)];
        const double t = xi[std::min(2, dim - 1)];

        if (fabs(t - 1.0) > 1.0e-14)
          {
            ration = (r * s * t) / (1.0 - t);
          }
        else
          {
            ration = 0.0;
          }

        if (i == 0)
          return Q14 * ((1.0 - r) * (1.0 - s) - t + ration);
        if (i == 1)
          return Q14 * ((1.0 + r) * (1.0 - s) - t - ration);
        if (i == 2)
          return Q14 * ((1.0 - r) * (1.0 + s) - t - ration);
        if (i == 3)
          return Q14 * ((1.0 + r) * (1.0 + s) - t + ration);
        else
          return t;
      }

    Assert(false, ExcNotImplemented());

    return 0.0;
  }

  /**
   * Compute the gradient of the $i$-th linear shape function at location $\xi$
   * for a given reference-cell type.
   */
  template <int dim>
  inline Tensor<1, dim>
  d_linear_shape_function_gradient(const Type &       reference_cell,
                                   const Point<dim> & xi,
                                   const unsigned int i)
  {
    if (reference_cell == get_hypercube(dim))
      return GeometryInfo<dim>::d_linear_shape_function_gradient(xi, i);

    if (reference_cell ==
        Type::Tri) // see also Simplex::ScalarPolynomial::compute_grad
      {
        switch (i)
          {
            case 0:
              return Point<dim>(-1.0, -1.0);
            case 1:
              return Point<dim>(+1.0, +0.0);
            case 2:
              return Point<dim>(+0.0, +1.0);
          }
      }

    Assert(false, ExcNotImplemented());

    return Point<dim>(+0.0, +0.0, +0.0);
  }

  /**
   * Return i-th unit tangential vector of a face of the reference cell.
   * The vectors are arranged such that the
   * cross product between the two vectors returns the unit normal vector.
   */
  template <int dim>
  inline Tensor<1, dim>
  unit_tangential_vectors(const Type &       reference_cell,
                          const unsigned int face_no,
                          const unsigned int i)
  {
    AssertDimension(dim, get_dimension(reference_cell));
    AssertIndexRange(i, dim - 1);

    if (reference_cell == get_hypercube(dim))
      {
        AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);
        return GeometryInfo<dim>::unit_tangential_vectors[face_no][i];
      }
    else if (reference_cell == Type::Tri)
      {
        AssertIndexRange(face_no, 3);
        static const std::array<Tensor<1, dim>, 3> table = {
          {Point<dim>(1, 0),
           Point<dim>(-std::sqrt(0.5), +std::sqrt(0.5)),
           Point<dim>(0, -1)}};

        return table[face_no];
      }
    else if (reference_cell == Type::Tet)
      {
        AssertIndexRange(face_no, 4);
        static const std::array<std::array<Tensor<1, dim>, 2>, 4> table = {
          {{{Point<dim>(0, 1, 0), Point<dim>(1, 0, 0)}},
           {{Point<dim>(1, 0, 0), Point<dim>(0, 0, 1)}},
           {{Point<dim>(0, 0, 1), Point<dim>(0, 1, 0)}},
           {{Point<dim>(-std::pow(1.0 / 3.0, 1.0 / 4.0),
                        +std::pow(1.0 / 3.0, 1.0 / 4.0),
                        0),
             Point<dim>(-std::pow(1.0 / 3.0, 1.0 / 4.0),
                        0,
                        +std::pow(1.0 / 3.0, 1.0 / 4.0))}}}};

        return table[face_no][i];
      }
    else if (reference_cell == Type::Wedge)
      {
        AssertIndexRange(face_no, 5);
        static const std::array<std::array<Tensor<1, dim>, 2>, 5> table = {
          {{{Point<dim>(0, 1, 0), Point<dim>(1, 0, 0)}},
           {{Point<dim>(1, 0, 0), Point<dim>(0, 0, 1)}},
           {{Point<dim>(-1 / std::sqrt(2.0), +1 / std::sqrt(2.0), 0),
             Point<dim>(0, 0, 1)}},
           {{Point<dim>(0, 0, 1), Point<dim>(0, 1, 0)}},
           {{Point<dim>(1, 0, 0), Point<dim>(0, 0, 1)}}}};

        return table[face_no][i];
      }
    else if (reference_cell == Type::Pyramid)
      {
        AssertIndexRange(face_no, 5);
        static const std::array<std::array<Tensor<1, dim>, 2>, 5> table = {
          {{{Point<dim>(0, 1, 0), Point<dim>(1, 0, 0)}},
           {{Point<dim>(+1.0 / sqrt(2.0), 0, +1.0 / sqrt(2.0)),
             Point<dim>(0, 1, 0)}},
           {{Point<dim>(+1.0 / sqrt(2.0), 0, -1.0 / sqrt(2.0)),
             Point<dim>(0, 1, 0)}},
           {{Point<dim>(1, 0, 0),
             Point<dim>(0, +1.0 / sqrt(2.0), +1.0 / sqrt(2.0))}},
           {{Point<dim>(1, 0, 0),
             Point<dim>(0, +1.0 / sqrt(2.0), -1.0 / sqrt(2.0))}}}};

        return table[face_no][i];
      }

    Assert(false, ExcNotImplemented());

    return {};
  }

  /**
   * Return the unit normal vector of a face of the reference cell.
   */
  template <int dim>
  inline Tensor<1, dim>
  unit_normal_vectors(const Type &reference_cell, const unsigned int face_no)
  {
    AssertDimension(dim, get_dimension(reference_cell));

    if (reference_cell == get_hypercube(dim))
      {
        AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);
        return GeometryInfo<dim>::unit_normal_vector[face_no];
      }
    else if (dim == 2)
      {
        const auto tangential =
          unit_tangential_vectors<dim>(reference_cell, face_no, 0);

        Tensor<1, dim> result;

        result[0] = tangential[1];
        result[1] = -tangential[0];

        return result;
      }
    else if (dim == 3)
      {
        return cross_product_3d(
          unit_tangential_vectors<dim>(reference_cell, face_no, 0),
          unit_tangential_vectors<dim>(reference_cell, face_no, 1));
      }

    Assert(false, ExcNotImplemented());

    return {};
  }

  /*
   * Create a (coarse) grid with a single cell of the shape of the provided
   * reference cell.
   */
  template <int dim, int spacedim>
  void
  make_triangulation(const Type &                  reference_cell,
                     Triangulation<dim, spacedim> &tria);

  /**
   * Return a default linear mapping matching the given reference cell.
   * If this reference cell is a hypercube, then the returned mapping
   * is a MappingQ1; otherwise, it is an object of type MappingFE
   * initialized with Simplex::FE_P (if the reference cell is a triangle and
   * tetrahedron), with Simplex::FE_PyramidP (if the reference cell is a
   * pyramid), or with Simplex::FE_WedgeP (if the reference cell is a wedge). In
   * other words, the term "linear" in the name of the function has to be
   * understood as $d$-linear (i.e., bilinear or trilinear) for some of the
   * coordinate directions.
   */
  template <int dim, int spacedim>
  const Mapping<dim, spacedim> &
  get_default_linear_mapping(const Type &reference_cell);

  /**
   * Return a default linear mapping that works for the given triangulation.
   * Internally, this function calls the function above for the reference
   * cell used by the given triangulation, assuming that the triangulation
   * uses only a single cell type. If the triangulation uses mixed cell
   * types, then this function will trigger an exception.
   */
  template <int dim, int spacedim>
  const Mapping<dim, spacedim> &
  get_default_linear_mapping(const Triangulation<dim, spacedim> &triangulation);

  /**
   * Return a Gauss-type quadrature matching the given reference cell(QGauss,
   * Simplex::QGauss, Simplex::QGaussPyramid, Simplex::QGaussWedge) and
   * @p n_points_1D the number of quadrature points in each direction (QGuass)
   * or the indication of what polynomial degree to be integrated exactly.
   */
  template <int dim>
  Quadrature<dim>
  get_gauss_type_quadrature(const Type &   reference_cell,
                            const unsigned n_points_1D);

  /**
   * Return a quadrature rule with the support points of the given reference
   * cell.
   *
   * @note The weights are not filled.
   */
  template <int dim>
  Quadrature<dim> &
  get_nodal_type_quadrature(const Type &reference_cell);

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

        /**
         * Map an ExodusII vertex number to a deal.II vertex number.
         */
        virtual unsigned int
        exodusii_vertex_to_deal_vertex(const unsigned int vertex_n) const
        {
          Assert(false, ExcNotImplemented());
          (void)vertex_n;

          return 0;
        }

        /**
         * Map an ExodusII face number to a deal.II face number.
         */
        virtual unsigned int
        exodusii_face_to_deal_face(const unsigned int face_n) const
        {
          Assert(false, ExcNotImplemented());
          (void)face_n;

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

        virtual unsigned int
        exodusii_face_to_deal_face(const unsigned int face_n) const override
        {
          (void)face_n;
          AssertIndexRange(face_n, n_faces());

          return 0;
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

        virtual unsigned int
        exodusii_vertex_to_deal_vertex(
          const unsigned int vertex_n) const override
        {
          AssertIndexRange(vertex_n, n_vertices());
          return vertex_n;
        }

        virtual unsigned int
        exodusii_face_to_deal_face(const unsigned int face_n) const override
        {
          AssertIndexRange(face_n, n_faces());
          return face_n;
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
            {{{0, 1}}, {{1, 2}}, {{2, 0}}}};

          return table[face][face_orientation ? vertex : (1 - vertex)];
        }

        virtual unsigned int
        exodusii_vertex_to_deal_vertex(
          const unsigned int vertex_n) const override
        {
          AssertIndexRange(vertex_n, n_vertices());
          return vertex_n;
        }

        virtual unsigned int
        exodusii_face_to_deal_face(const unsigned int face_n) const override
        {
          AssertIndexRange(face_n, n_faces());
          return face_n;
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

        virtual unsigned int
        exodusii_vertex_to_deal_vertex(
          const unsigned int vertex_n) const override
        {
          AssertIndexRange(vertex_n, n_vertices());
          constexpr std::array<unsigned int, 4> exodus_to_deal{{0, 1, 3, 2}};
          return exodus_to_deal[vertex_n];
        }

        virtual unsigned int
        exodusii_face_to_deal_face(const unsigned int face_n) const override
        {
          AssertIndexRange(face_n, n_faces());
          constexpr std::array<unsigned int, 4> exodus_to_deal{{2, 1, 3, 0}};
          return exodus_to_deal[face_n];
        }
      };



      /**
       * Tet.
       */
      struct Tet : public Base
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
             {{1, 0, 2}},
             {{1, 2, 0}},
             {{0, 2, 1}},
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
             {{2, 1, 0}},
             {{1, 2, 0}},
             {{1, 0, 2}},
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
            {{{0, 1, 2}}, {{0, 3, 4}}, {{2, 5, 3}}, {{1, 4, 5}}}};

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
            {{{0, 1, 2}}, {{1, 0, 3}}, {{0, 2, 3}}, {{2, 1, 3}}}};

          return table[face][standard_to_real_face_vertex(
            vertex, face, face_orientation)];
        }

        virtual unsigned int
        exodusii_vertex_to_deal_vertex(
          const unsigned int vertex_n) const override
        {
          AssertIndexRange(vertex_n, n_vertices());
          return vertex_n;
        }

        virtual unsigned int
        exodusii_face_to_deal_face(const unsigned int face_n) const override
        {
          AssertIndexRange(face_n, n_faces());
          constexpr std::array<unsigned int, 4> exodus_to_deal{{1, 3, 2, 0}};
          return exodus_to_deal[face_n];
        }
      };



      /**
       * Pyramid.
       */
      struct Pyramid : public Base
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
          static const std::array<unsigned int, 2> table[8] = {{{0, 0}},
                                                               {{0, 1}},
                                                               {{0, 2}},
                                                               {{0, 3}},
                                                               {{1, 2}},
                                                               {{2, 1}},
                                                               {{1, 1}},
                                                               {{2, 2}}};

          return table[line];
        }

        unsigned int
        standard_to_real_face_line(
          const unsigned int  line,
          const unsigned int  face,
          const unsigned char face_orientation) const override
        {
          if (face == 0) // QUAD
            {
              return GeometryInfo<3>::standard_to_real_face_line(
                line,
                get_bit(face_orientation, 0),
                get_bit(face_orientation, 2),
                get_bit(face_orientation, 1));
            }
          else // TRI
            {
              static const std::array<std::array<unsigned int, 3>, 6> table = {
                {{{2, 1, 0}},
                 {{0, 1, 2}},
                 {{1, 0, 2}},
                 {{1, 2, 0}},
                 {{0, 2, 1}},
                 {{2, 0, 1}}}};

              return table[face_orientation][line];
            }
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
                 {{2, 1, 0}},
                 {{1, 2, 0}},
                 {{1, 0, 2}},
                 {{2, 0, 1}}}};

              return table[face_orientation][vertex];
            }
        }

        ReferenceCell::Type
        face_reference_cell_type(const unsigned int face_no) const override
        {
          AssertIndexRange(face_no, n_faces());

          if (face_no == 0)
            return ReferenceCell::Type::Quad;
          else
            return ReferenceCell::Type::Tri;
        }

        unsigned int
        face_to_cell_vertices(
          const unsigned int  face,
          const unsigned int  vertex,
          const unsigned char face_orientation) const override
        {
          AssertIndexRange(face, n_faces());
          if (face == 0)
            {
              AssertIndexRange(vertex, 4);
            }
          else
            {
              AssertIndexRange(vertex, 3);
            }
          constexpr auto X = numbers::invalid_unsigned_int;
          static const std::array<std::array<unsigned int, 4>, 5> table = {
            {{{0, 1, 2, 3}},
             {{0, 2, 4, X}},
             {{3, 1, 4, X}},
             {{1, 0, 4, X}},
             {{2, 3, 4, X}}}};

          return table[face][standard_to_real_face_vertex(
            vertex, face, face_orientation)];
        }

        virtual unsigned int
        exodusii_vertex_to_deal_vertex(
          const unsigned int vertex_n) const override
        {
          AssertIndexRange(vertex_n, n_vertices());
          constexpr std::array<unsigned int, 5> exodus_to_deal{{0, 1, 3, 2, 4}};
          return exodus_to_deal[vertex_n];
        }

        virtual unsigned int
        exodusii_face_to_deal_face(const unsigned int face_n) const override
        {
          AssertIndexRange(face_n, n_faces());
          constexpr std::array<unsigned int, 5> exodus_to_deal{{3, 2, 4, 1, 0}};
          return exodus_to_deal[face_n];
        }
      };



      /**
       * Wedge.
       */
      struct Wedge : public Base
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
          return 5;
        }

        std::array<unsigned int, 2>
        standard_line_to_face_and_line_index(
          const unsigned int line) const override
        {
          static const std::array<unsigned int, 2> table[9] = {{{0, 0}},
                                                               {{0, 2}},
                                                               {{0, 1}},
                                                               {{1, 0}},
                                                               {{1, 1}},
                                                               {{1, 2}},
                                                               {{2, 0}},
                                                               {{2, 1}},
                                                               {{3, 1}}};

          return table[line];
        }

        unsigned int
        standard_to_real_face_line(
          const unsigned int  line,
          const unsigned int  face,
          const unsigned char face_orientation) const override
        {
          if (face > 1) // QUAD
            {
              return GeometryInfo<3>::standard_to_real_face_line(
                line,
                get_bit(face_orientation, 0),
                get_bit(face_orientation, 2),
                get_bit(face_orientation, 1));
            }
          else // TRI
            {
              static const std::array<std::array<unsigned int, 3>, 6> table = {
                {{{2, 1, 0}},
                 {{0, 1, 2}},
                 {{1, 0, 2}},
                 {{1, 2, 0}},
                 {{0, 2, 1}},
                 {{2, 0, 1}}}};

              return table[face_orientation][line];
            }
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
                 {{2, 1, 0}},
                 {{1, 2, 0}},
                 {{1, 0, 2}},
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

        unsigned int
        face_to_cell_vertices(
          const unsigned int  face,
          const unsigned int  vertex,
          const unsigned char face_orientation) const override
        {
          AssertIndexRange(face, n_faces());
          if (face < 2)
            {
              AssertIndexRange(vertex, 3);
            }
          else
            {
              AssertIndexRange(vertex, 4);
            }
          constexpr auto X = numbers::invalid_unsigned_int;
          static const std::array<std::array<unsigned int, 4>, 6> table = {
            {{{1, 0, 2, X}},
             {{3, 4, 5, X}},
             {{0, 1, 3, 4}},
             {{1, 2, 4, 5}},
             {{2, 0, 5, 3}}}};

          return table[face][standard_to_real_face_vertex(
            vertex, face, face_orientation)];
        }

        virtual unsigned int
        exodusii_vertex_to_deal_vertex(
          const unsigned int vertex_n) const override
        {
          AssertIndexRange(vertex_n, n_vertices());
          constexpr std::array<unsigned int, 6> exodus_to_deal{
            {2, 1, 0, 5, 4, 3}};
          return exodus_to_deal[vertex_n];
        }

        virtual unsigned int
        exodusii_face_to_deal_face(const unsigned int face_n) const override
        {
          AssertIndexRange(face_n, n_faces());
          constexpr std::array<unsigned int, 6> exodus_to_deal{{3, 4, 2, 0, 1}};
          return exodus_to_deal[face_n];
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

        virtual unsigned int
        exodusii_vertex_to_deal_vertex(
          const unsigned int vertex_n) const override
        {
          AssertIndexRange(vertex_n, n_vertices());
          constexpr std::array<unsigned int, 8> exodus_to_deal{
            {0, 1, 3, 2, 4, 5, 7, 6}};
          return exodus_to_deal[vertex_n];
        }

        virtual unsigned int
        exodusii_face_to_deal_face(const unsigned int face_n) const override
        {
          AssertIndexRange(face_n, n_faces());
          constexpr std::array<unsigned int, 6> exodus_to_deal{
            {2, 1, 3, 0, 4, 5}};
          return exodus_to_deal[face_n];
        }
      };

      /**
       * Return for a given reference-cell type the right Info.
       */
      inline const ReferenceCell::internal::Info::Base &
      get_cell(const ReferenceCell::Type &type)
      {
        static const std::
          array<std::unique_ptr<ReferenceCell::internal::Info::Base>, 8>
            gei{{std::make_unique<ReferenceCell::internal::Info::Vertex>(),
                 std::make_unique<ReferenceCell::internal::Info::Line>(),
                 std::make_unique<ReferenceCell::internal::Info::Tri>(),
                 std::make_unique<ReferenceCell::internal::Info::Quad>(),
                 std::make_unique<ReferenceCell::internal::Info::Tet>(),
                 std::make_unique<ReferenceCell::internal::Info::Pyramid>(),
                 std::make_unique<ReferenceCell::internal::Info::Wedge>(),
                 std::make_unique<ReferenceCell::internal::Info::Hex>()}};
        AssertIndexRange(static_cast<std::uint8_t>(type), 8);
        return *gei[static_cast<std::uint8_t>(type)];
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



  namespace internal
  {
    template <typename T, std::size_t N>
    class NoPermutation : public dealii::ExceptionBase
    {
    public:
      /**
       * Constructor.
       */
      NoPermutation(const ReferenceCell::Type &entity_type,
                    const std::array<T, N> &   vertices_0,
                    const std::array<T, N> &   vertices_1)
        : entity_type(entity_type)
        , vertices_0(vertices_0)
        , vertices_1(vertices_1)
      {}

      /**
       * Destructor.
       */
      virtual ~NoPermutation() noexcept override = default;

      /**
       * Print error message to @p out.
       */
      virtual void
      print_info(std::ostream &out) const override
      {
        out << "[";

        const unsigned int n_vertices =
          ReferenceCell::internal::Info::get_cell(entity_type).n_vertices();

        for (unsigned int i = 0; i < n_vertices; ++i)
          {
            out << vertices_0[i];
            if (i + 1 != n_vertices)
              out << ",";
          }

        out << "] is no permutation of [";

        for (unsigned int i = 0; i < n_vertices; ++i)
          {
            out << vertices_1[i];
            if (i + 1 != n_vertices)
              out << ",";
          }

        out << "]." << std::endl;
      }

      /**
       * Entity type.
       */
      const ReferenceCell::Type entity_type;

      /**
       * First set of values.
       */
      const std::array<T, N> vertices_0;

      /**
       * Second set of values.
       */
      const std::array<T, N> vertices_1;
    };
  } // namespace internal



  /**
   * Determine the orientation of an entity of @p type described by its
   * vertices @p var_1 relative to an entity described by @p var_0.
   */
  template <typename T, std::size_t N>
  inline unsigned char
  compute_orientation(const ReferenceCell::Type entity_type,
                      const std::array<T, N> &  vertices_0,
                      const std::array<T, N> &  vertices_1)
  {
    AssertIndexRange(
      ReferenceCell::internal::Info::get_cell(entity_type).n_vertices(), N + 1);
    if (entity_type == ReferenceCell::Type::Line)
      {
        const std::array<T, 2> i{{vertices_0[0], vertices_0[1]}};
        const std::array<T, 2> j{{vertices_1[0], vertices_1[1]}};

        // line_orientation=true
        if (i == std::array<T, 2>{{j[0], j[1]}})
          return 1;

        // line_orientation=false
        if (i == std::array<T, 2>{{j[1], j[0]}})
          return 0;
      }
    else if (entity_type == ReferenceCell::Type::Tri)
      {
        const std::array<T, 3> i{{vertices_0[0], vertices_0[1], vertices_0[2]}};
        const std::array<T, 3> j{{vertices_1[0], vertices_1[1], vertices_1[2]}};

        // face_orientation=true, face_rotation=false, face_flip=false
        if (i == std::array<T, 3>{{j[0], j[1], j[2]}})
          return 1;

        // face_orientation=true, face_rotation=true, face_flip=false
        if (i == std::array<T, 3>{{j[1], j[2], j[0]}})
          return 3;

        // face_orientation=true, face_rotation=false, face_flip=true
        if (i == std::array<T, 3>{{j[2], j[0], j[1]}})
          return 5;

        // face_orientation=false, face_rotation=false, face_flip=false
        if (i == std::array<T, 3>{{j[0], j[2], j[1]}})
          return 0;

        // face_orientation=false, face_rotation=true, face_flip=false
        if (i == std::array<T, 3>{{j[2], j[1], j[0]}})
          return 2;

        // face_orientation=false, face_rotation=false, face_flip=true
        if (i == std::array<T, 3>{{j[1], j[0], j[2]}})
          return 4;
      }
    else if (entity_type == ReferenceCell::Type::Quad)
      {
        const std::array<T, 4> i{
          {vertices_0[0], vertices_0[1], vertices_0[2], vertices_0[3]}};
        const std::array<T, 4> j{
          {vertices_1[0], vertices_1[1], vertices_1[2], vertices_1[3]}};

        // face_orientation=true, face_rotation=false, face_flip=false
        if (i == std::array<T, 4>{{j[0], j[1], j[2], j[3]}})
          return 1;

        // face_orientation=true, face_rotation=true, face_flip=false
        if (i == std::array<T, 4>{{j[2], j[0], j[3], j[1]}})
          return 3;

        // face_orientation=true, face_rotation=false, face_flip=true
        if (i == std::array<T, 4>{{j[3], j[2], j[1], j[0]}})
          return 5;

        // face_orientation=true, face_rotation=true, face_flip=true
        if (i == std::array<T, 4>{{j[1], j[3], j[0], j[2]}})
          return 7;

        // face_orientation=false, face_rotation=false, face_flip=false
        if (i == std::array<T, 4>{{j[0], j[2], j[1], j[3]}})
          return 0;

        // face_orientation=false, face_rotation=true, face_flip=false
        if (i == std::array<T, 4>{{j[2], j[3], j[0], j[1]}})
          return 2;

        // face_orientation=false, face_rotation=false, face_flip=true
        if (i == std::array<T, 4>{{j[3], j[1], j[2], j[0]}})
          return 4;

        // face_orientation=false, face_rotation=true, face_flip=true
        if (i == std::array<T, 4>{{j[1], j[0], j[3], j[2]}})
          return 6;
      }

    Assert(
      false,
      (internal::NoPermutation<T, N>(entity_type, vertices_0, vertices_1)));

    return -1;
  }

  /**
   * Inverse function of compute_orientation().
   */
  template <typename T, std::size_t N>
  inline std::array<T, N>
  permute_according_orientation(const ReferenceCell::Type entity_type,
                                const std::array<T, N> &  vertices,
                                const unsigned int        orientation)
  {
    std::array<T, 4> temp;

    if (entity_type == ReferenceCell::Type::Line)
      {
        switch (orientation)
          {
            case 1:
              temp = {{vertices[0], vertices[1]}};
              break;
            case 0:
              temp = {{vertices[1], vertices[0]}};
              break;
            default:
              Assert(false, ExcNotImplemented());
          }
      }
    else if (entity_type == ReferenceCell::Type::Tri)
      {
        switch (orientation)
          {
            case 1:
              temp = {{vertices[0], vertices[1], vertices[2]}};
              break;
            case 3:
              temp = {{vertices[1], vertices[2], vertices[0]}};
              break;
            case 5:
              temp = {{vertices[2], vertices[0], vertices[1]}};
              break;
            case 0:
              temp = {{vertices[0], vertices[2], vertices[1]}};
              break;
            case 2:
              temp = {{vertices[2], vertices[1], vertices[0]}};
              break;
            case 4:
              temp = {{vertices[1], vertices[0], vertices[2]}};
              break;
            default:
              Assert(false, ExcNotImplemented());
          }
      }
    else if (entity_type == ReferenceCell::Type::Quad)
      {
        switch (orientation)
          {
            case 1:
              temp = {{vertices[0], vertices[1], vertices[2], vertices[3]}};
              break;
            case 3:
              temp = {{vertices[2], vertices[0], vertices[3], vertices[1]}};
              break;
            case 5:
              temp = {{vertices[3], vertices[2], vertices[1], vertices[0]}};
              break;
            case 7:
              temp = {{vertices[1], vertices[3], vertices[0], vertices[2]}};
              break;
            case 0:
              temp = {{vertices[0], vertices[2], vertices[1], vertices[3]}};
              break;
            case 2:
              temp = {{vertices[2], vertices[3], vertices[0], vertices[1]}};
              break;
            case 4:
              temp = {{vertices[3], vertices[1], vertices[2], vertices[0]}};
              break;
            case 6:
              temp = {{vertices[1], vertices[0], vertices[3], vertices[2]}};
              break;
            default:
              Assert(false, ExcNotImplemented());
          }
      }
    else
      {
        AssertThrow(false, ExcNotImplemented());
      }

    std::array<T, N> temp_;
    std::copy_n(temp.begin(), N, temp_.begin());

    return temp_;
  }

} // namespace ReferenceCell


DEAL_II_NAMESPACE_CLOSE

#endif

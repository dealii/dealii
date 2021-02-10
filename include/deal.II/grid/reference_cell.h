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
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
class Mapping;

template <int dim>
class Quadrature;

class ReferenceCell;
#endif


namespace internal
{
  namespace ReferenceCell
  {
    /**
     * A helper function to create a ReferenceCell object from an
     * integer. ReferenceCell objects are "singletons" (actually,
     * "multitons" -- there are multiple, but they are only a handful and
     * these are all that can be used). What is then necessary is to
     * have a way to create these with their internal id to distinguish
     * the few possible ones in existence. We could do this via a public
     * constructor of ReferenceCell, but that would allow users
     * to create ones outside the range we envision, and we don't want to do
     * that. Rather, the constructor that takes an integer is made `private`
     * but we have this one function in an internal namespace that is a friend
     * of the class and can be used to create the objects.
     */
    DEAL_II_CONSTEXPR dealii::ReferenceCell
                      make_reference_cell_from_int(const std::uint8_t kind);
  } // namespace ReferenceCell
} // namespace internal



/**
 * A type that describes the kinds of reference cells that can be used.
 * This includes quadrilaterals and hexahedra (i.e., "hypercubes"),
 * triangles and tetrahedra (simplices), and the pyramids and wedges
 * necessary when using mixed 3d meshes.
 *
 * Objects of this type should not be created in user code, and as a
 * consequence the class does not have a user-accessible constructor
 * other than the default constructor (which creates an invalid object).
 * Rather, there is a finite number of specific reference cell objects
 * defined in the ReferenceCells namespace that completely enumerate
 * all of the possible values. User codes should therefore rely
 * exclusively on assigning ReferenceCell objects from these special
 * objects, and comparing against those special objects.
 */
class ReferenceCell
{
public:
  /**
   * Return the correct ReferenceCell for a given structural
   * dimension and number of vertices. For example, if `dim==2` and
   * `n_vertices==4`, this function will return ReferenceCells::Quadrilateral.
   * But if `dim==3` and `n_vertices==4`, it will return
   * ReferenceCells::Tetrahedron.
   */
  static ReferenceCell
  n_vertices_to_type(const int dim, const unsigned int n_vertices);

  /**
   * Default constructor. Initialize this object as an invalid object. The
   * end result is that the current object equals ReferenceCells::Invalid.
   *
   * Generally, ReferenceCell objects are created by assignment from
   * the special objects in namespace ReferenceCells, which is the only
   * way to obtain a valid object.
   */
  DEAL_II_CONSTEXPR
  ReferenceCell();

  /**
   * @name Querying information about the kind of reference cells
   * @{
   */

  /**
   * Return `true` if the object is a ReferenceCells::Vertex,
   * ReferenceCells::Line, ReferenceCells::Quadrilateral, or
   * ReferenceCells::Hexahedron.
   */
  bool
  is_hyper_cube() const;

  /**
   * Return true if the object is a Vertex, Line, Triangle, or Tetrahedron.
   */
  bool
  is_simplex() const;

  /**
   * Return the dimension of the reference cell represented by the current
   * object.
   */
  unsigned int
  get_dimension() const;

  /**
   * @}
   */

  /**
   * @name Shape functions, mappings, quadratures defined on a reference cell
   * @{
   */

  /**
   * Compute the value of the $i$-th linear shape function at location $\xi$
   * for the current reference-cell type.
   */
  template <int dim>
  double
  d_linear_shape_function(const Point<dim> &xi, const unsigned int i) const;

  /**
   * Compute the gradient of the $i$-th linear shape function at location
   * $\xi$ for the current reference-cell type.
   */
  template <int dim>
  Tensor<1, dim>
  d_linear_shape_function_gradient(const Point<dim> & xi,
                                   const unsigned int i) const;

  /**
   * Return a default mapping of degree @p degree matching the current
   * reference cell. If this reference cell is a hypercube, then the returned
   * mapping is a MappingQGeneric; otherwise, it is an object of type
   * MappingFE initialized with Simplex::FE_P (if the reference cell is a
   * triangle or tetrahedron), with Simplex::FE_PyramidP (if the reference
   * cell is a pyramid), or with Simplex::FE_WedgeP (if the reference cell is
   * a wedge).
   */
  template <int dim, int spacedim>
  std::unique_ptr<Mapping<dim, spacedim>>
  get_default_mapping(const unsigned int degree) const;

  /**
   * Return a default linear mapping matching the current reference cell.
   * If this reference cell is a hypercube, then the returned mapping
   * is a MappingQ1; otherwise, it is an object of type MappingFE
   * initialized with Simplex::FE_P (if the reference cell is a triangle or
   * tetrahedron), with Simplex::FE_PyramidP (if the reference cell is a
   * pyramid), or with Simplex::FE_WedgeP (if the reference cell is a wedge).
   * In other words, the term "linear" in the name of the function has to be
   * understood as $d$-linear (i.e., bilinear or trilinear) for some of the
   * coordinate directions.
   */
  template <int dim, int spacedim>
  const Mapping<dim, spacedim> &
  get_default_linear_mapping() const;

  /**
   * Return a Gauss-type quadrature matching the given reference cell (QGauss,
   * Simplex::QGauss, Simplex::QGaussPyramid, Simplex::QGaussWedge).
   *
   * @param[in] n_points_1D The number of quadrature points in each direction
   * (QGauss) or an indication of what polynomial degree needs to be
   * integrated exactly for the other types.
   */
  template <int dim>
  Quadrature<dim>
  get_gauss_type_quadrature(const unsigned n_points_1D) const;

  /**
   * Return a quadrature rule with the support points of the given reference
   * cell.
   *
   * @note The weights of the quadrature object are left unfilled.
   */
  template <int dim>
  const Quadrature<dim> &
  get_nodal_type_quadrature() const;

  /**
   * @}
   */

  /**
   * @name Querying the number of building blocks of a reference cell
   * @{
   */

  /**
   * Return the number of vertices that make up the reference
   * cell in question. A vertex is a "corner" (a zero-dimensional
   * object) of the reference cell.
   */
  unsigned int
  n_vertices() const;

  /**
   * Return an object that can be thought of as an array containing all
   * indices from zero to n_vertices().
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  vertex_indices() const;

  /**
   * Return the number of lines that make up the reference
   * cell in question. A line is an "edge" (a one-dimensional
   * object) of the reference cell.
   */
  unsigned int
  n_lines() const;

  /**
   * Return an object that can be thought of as an array containing all
   * indices from zero to n_lines().
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  line_indices() const;

  /**
   * Return the number of faces that make up the reference
   * cell in question. A face is a `(dim-1)`-dimensional
   * object bounding the reference cell.
   */
  unsigned int
  n_faces() const;

  /**
   * Return an object that can be thought of as an array containing all
   * indices from zero to n_faces().
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  face_indices() const;

  /**
   * Return the reference-cell type of face @p face_no of the current
   * object. For example, if the current object is
   * ReferenceCells::Tetrahedron, then `face_no` must be between
   * in the interval $[0,4)$ and the function will always return
   * ReferenceCells::Triangle. If the current object is
   * ReferenceCells::Hexahedron, then `face_no` must be between
   * in the interval $[0,6)$ and the function will always return
   * ReferenceCells::Quadrilateral. For wedges and pyramids, the
   * returned object may be either ReferenceCells::Triangle or
   * ReferenceCells::Quadrilateral, depending on the given index.
   */
  ReferenceCell
  face_reference_cell(const unsigned int face_no) const;

  /**
   * @}
   */

  /**
   * @name Relationships between objects in the cell and on faces
   * @{
   */

  /**
   * Return which child cells are adjacent to a certain face of the
   * mother cell.
   *
   * For example, in 2D the layout of a quadrilateral cell is as follows:
   * @verbatim
   * .      3
   * .   2-->--3
   * .   |     |
   * . 0 ^     ^ 1
   * .   |     |
   * .   0-->--1
   * .      2
   * @endverbatim
   * Vertices and faces are indicated with their numbers, faces also with
   * their directions.
   *
   * Now, when refined, the layout is like this:
   * @verbatim
   * *--*--*
   * | 2|3 |
   * *--*--*
   * | 0|1 |
   * *--*--*
   * @endverbatim
   *
   * Thus, the child cells on face 0 are (ordered in the direction of the
   * face) 0 and 2, on face 3 they are 2 and 3, etc.
   *
   * For three spatial dimensions, the exact order of the children is laid
   * down in the general documentation of this class.
   */
  unsigned int
  child_cell_on_face(const unsigned int face_n,
                     const unsigned int subface_n) const;

  /**
   * For a given vertex in a cell, return a pair of a face index and a
   * vertex index within this face.
   *
   * @note In practice, a vertex is of course generally part of more than one
   *   face, and one could return different faces and the corresponding
   *   index within. Which face this function chooses is often not of
   *   importance (and not exposed by this function on purpose).
   */
  std::array<unsigned int, 2>
  standard_vertex_to_face_and_vertex_index(const unsigned int vertex) const;

  /**
   * For a given line in a cell, return a pair of a face index and a
   * line index within this face.
   *
   * @note In practice, a line is of course generally part of more than one
   *   face, and one could return different faces and the corresponding
   *   index within. Which face this function chooses is often not of
   *   importance (and not exposed by this function on purpose).
   */
  std::array<unsigned int, 2>
  standard_line_to_face_and_line_index(const unsigned int line) const;

  /**
   * Map face line number to cell line number.
   */
  unsigned int
  face_to_cell_lines(const unsigned int  face,
                     const unsigned int  line,
                     const unsigned char face_orientation) const;

  /**
   * Map face vertex number to cell vertex number.
   */
  unsigned int
  face_to_cell_vertices(const unsigned int  face,
                        const unsigned int  vertex,
                        const unsigned char face_orientation) const;

  /**
   * Correct vertex index depending on face orientation.
   */
  unsigned int
  standard_to_real_face_vertex(const unsigned int  vertex,
                               const unsigned int  face,
                               const unsigned char face_orientation) const;

  /**
   * Correct line index depending on face orientation.
   */
  unsigned int
  standard_to_real_face_line(const unsigned int  line,
                             const unsigned int  face,
                             const unsigned char face_orientation) const;

  /**
   * @}
   */

  /**
   * @name Geometric properties of reference cells
   * @name Querying the number of building blocks of a reference cell
   * @{
   */

  /*
   * Return $i$-th unit tangential vector of a face of the reference cell.
   * The vectors are arranged such that the
   * cross product between the two vectors returns the unit normal vector.
   *
   * @pre $i$ must be between zero and `dim-1`.
   */
  template <int dim>
  Tensor<1, dim>
  unit_tangential_vectors(const unsigned int face_no,
                          const unsigned int i) const;

  /**
   * Return the unit normal vector of a face of the reference cell.
   */
  template <int dim>
  Tensor<1, dim>
  unit_normal_vectors(const unsigned int face_no) const;

  /**
   * Determine the orientation of the current entity described by its
   * vertices @p var_1 relative to an entity described by @p var_0.
   */
  template <typename T, std::size_t N>
  unsigned char
  compute_orientation(const std::array<T, N> &vertices_0,
                      const std::array<T, N> &vertices_1) const;

  /**
   * Inverse function of compute_orientation().
   */
  template <typename T, std::size_t N>
  std::array<T, N>
  permute_according_orientation(const std::array<T, N> &vertices,
                                const unsigned int      orientation) const;

  /**
   * Return a vector of faces a given @p vertex_index belongs to.
   */
  ArrayView<const unsigned int>
  faces_for_given_vertex(const unsigned int vertex_index) const;

  /**
   * @}
   */

  /**
   * @name Other functions
   * @{
   */

  /**
   * Return a text representation of the reference cell represented by the
   * current object.
   */
  std::string
  to_string() const;

  /**
   * Conversion operator to an integer.
   */
  constexpr operator std::uint8_t() const;

  /**
   * Operator for equality comparison.
   */
  constexpr bool
  operator==(const ReferenceCell &type) const;

  /**
   * Operator for inequality comparison.
   */
  constexpr bool
  operator!=(const ReferenceCell &type) const;

  /**
   * Write and read the data of this object from a stream for the purpose
   * of serialization using the [BOOST serialization
   * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
   */
  template <class Archive>
  void
  serialize(Archive &archive, const unsigned int /*version*/);

  /**
   * @}
   */

private:
  /**
   * The variable that stores what this object actually corresponds to.
   */
  std::uint8_t kind;

  /**
   * Constructor. This is the constructor used to create the different
   * `static` member variables of this class. It is `private` but can
   * be called by a function in an internal namespace that is a `friend`
   * of this class.
   */
  constexpr ReferenceCell(const std::uint8_t kind);

  /**
   * A kind of constructor -- not quite private because it can be
   * called by anyone, but at least hidden in an internal namespace.
   */
  friend DEAL_II_CONSTEXPR ReferenceCell
                           internal::ReferenceCell::make_reference_cell_from_int(const std::uint8_t);
};



inline constexpr ReferenceCell::ReferenceCell(const std::uint8_t kind)
  : kind(kind)
{}



inline constexpr ReferenceCell::operator std::uint8_t() const
{
  return kind;
}



inline constexpr bool
ReferenceCell::operator==(const ReferenceCell &type) const
{
  return kind == type.kind;
}



inline constexpr bool
ReferenceCell::operator!=(const ReferenceCell &type) const
{
  return kind != type.kind;
}



namespace internal
{
  namespace ReferenceCell
  {
    inline DEAL_II_CONSTEXPR dealii::ReferenceCell
                             make_reference_cell_from_int(const std::uint8_t kind)
    {
      // Make sure these are the only indices from which objects can be
      // created.
      Assert((kind == static_cast<std::uint8_t>(-1)) || (kind < 8),
             ExcInternalError());

      // Call the private constructor, which we can from here because this
      // function is a 'friend'.
      return {kind};
    }
  } // namespace ReferenceCell
} // namespace internal



/**
 * A namespace in which we define objects that correspond to specific
 * reference cells. The objects defined here are a complete enumeration
 * of all possible reference cells that can be used in deal.II.
 *
 * @relates ReferenceCell
 */
namespace ReferenceCells
{
  DEAL_II_CONSTEXPR const ReferenceCell Vertex =
    internal::ReferenceCell::make_reference_cell_from_int(0);
  DEAL_II_CONSTEXPR const ReferenceCell Line =
    internal::ReferenceCell::make_reference_cell_from_int(1);
  DEAL_II_CONSTEXPR const ReferenceCell Triangle =
    internal::ReferenceCell::make_reference_cell_from_int(2);
  DEAL_II_CONSTEXPR const ReferenceCell Quadrilateral =
    internal::ReferenceCell::make_reference_cell_from_int(3);
  DEAL_II_CONSTEXPR const ReferenceCell Tetrahedron =
    internal::ReferenceCell::make_reference_cell_from_int(4);
  DEAL_II_CONSTEXPR const ReferenceCell Pyramid =
    internal::ReferenceCell::make_reference_cell_from_int(5);
  DEAL_II_CONSTEXPR const ReferenceCell Wedge =
    internal::ReferenceCell::make_reference_cell_from_int(6);
  DEAL_II_CONSTEXPR const ReferenceCell Hexahedron =
    internal::ReferenceCell::make_reference_cell_from_int(7);
  DEAL_II_CONSTEXPR const ReferenceCell Invalid =
    internal::ReferenceCell::make_reference_cell_from_int(
      static_cast<std::uint8_t>(-1));

  /**
   * Return the correct simplex reference cell type for the given dimension
   * `dim`. Depending on the template argument `dim`, this function returns a
   * reference to either Vertex, Triangle, or Tetrahedron.
   */
  template <int dim>
  constexpr const ReferenceCell &
  get_simplex();

  /**
   * Return the correct hypercube reference cell type for the given dimension
   * `dim`. Depending on the template argument `dim`, this function returns a
   * reference to either Vertex, Quadrilateral, or Hexahedron.
   */
  template <int dim>
  constexpr const ReferenceCell &
  get_hypercube();
} // namespace ReferenceCells



inline DEAL_II_CONSTEXPR
ReferenceCell::ReferenceCell()
  : ReferenceCell(ReferenceCells::Invalid)
{}



template <class Archive>
inline void
ReferenceCell::serialize(Archive &archive, const unsigned int /*version*/)
{
  archive &kind;
}



inline ArrayView<const unsigned int>
ReferenceCell::faces_for_given_vertex(const unsigned int vertex) const
{
  if (*this == ReferenceCells::Line)
    {
      AssertIndexRange(vertex, GeometryInfo<1>::vertices_per_cell);
      return {&GeometryInfo<2>::vertex_to_face[vertex][0], 1};
    }
  else if (*this == ReferenceCells::Quadrilateral)
    {
      AssertIndexRange(vertex, GeometryInfo<2>::vertices_per_cell);
      return {&GeometryInfo<2>::vertex_to_face[vertex][0], 2};
    }
  else if (*this == ReferenceCells::Hexahedron)
    {
      AssertIndexRange(vertex, GeometryInfo<3>::vertices_per_cell);
      return {&GeometryInfo<3>::vertex_to_face[vertex][0], 3};
    }
  else if (*this == ReferenceCells::Triangle)
    {
      AssertIndexRange(vertex, 3);
      static const std::array<std::array<unsigned int, 2>, 3> table = {
        {{{0, 2}}, {{0, 1}}, {{1, 2}}}};

      return table[vertex];
    }
  else if (*this == ReferenceCells::Tetrahedron)
    {
      AssertIndexRange(vertex, 4);
      static const std::array<std::array<unsigned int, 3>, 4> table = {
        {{{0, 1, 2}}, {{0, 1, 3}}, {{0, 2, 3}}, {{1, 2, 3}}}};

      return table[vertex];
    }
  else if (*this == ReferenceCells::Wedge)
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
  else if (*this == ReferenceCells::Pyramid)
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



inline bool
ReferenceCell::is_hyper_cube() const
{
  return (*this == ReferenceCells::Vertex || *this == ReferenceCells::Line ||
          *this == ReferenceCells::Quadrilateral ||
          *this == ReferenceCells::Hexahedron);
}



inline bool
ReferenceCell::is_simplex() const
{
  return (*this == ReferenceCells::Vertex || *this == ReferenceCells::Line ||
          *this == ReferenceCells::Triangle ||
          *this == ReferenceCells::Tetrahedron);
}



inline unsigned int
ReferenceCell::get_dimension() const
{
  if (*this == ReferenceCells::Vertex)
    return 0;
  else if (*this == ReferenceCells::Line)
    return 1;
  else if ((*this == ReferenceCells::Triangle) ||
           (*this == ReferenceCells::Quadrilateral))
    return 2;
  else if ((*this == ReferenceCells::Tetrahedron) ||
           (*this == ReferenceCells::Pyramid) ||
           (*this == ReferenceCells::Wedge) ||
           (*this == ReferenceCells::Hexahedron))
    return 3;

  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}



inline unsigned int
ReferenceCell::n_vertices() const
{
  if (*this == ReferenceCells::Vertex)
    return 1;
  else if (*this == ReferenceCells::Line)
    return 2;
  else if (*this == ReferenceCells::Triangle)
    return 3;
  else if (*this == ReferenceCells::Quadrilateral)
    return 4;
  else if (*this == ReferenceCells::Tetrahedron)
    return 4;
  else if (*this == ReferenceCells::Pyramid)
    return 5;
  else if (*this == ReferenceCells::Wedge)
    return 6;
  else if (*this == ReferenceCells::Hexahedron)
    return 8;

  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}



inline unsigned int
ReferenceCell::n_lines() const
{
  if (*this == ReferenceCells::Vertex)
    return 0;
  else if (*this == ReferenceCells::Line)
    return 1;
  else if (*this == ReferenceCells::Triangle)
    return 3;
  else if (*this == ReferenceCells::Quadrilateral)
    return 4;
  else if (*this == ReferenceCells::Tetrahedron)
    return 6;
  else if (*this == ReferenceCells::Pyramid)
    return 7;
  else if (*this == ReferenceCells::Wedge)
    return 9;
  else if (*this == ReferenceCells::Hexahedron)
    return 12;

  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}



inline unsigned int
ReferenceCell::n_faces() const
{
  if (*this == ReferenceCells::Vertex)
    return 0;
  else if (*this == ReferenceCells::Line)
    return 2;
  else if (*this == ReferenceCells::Triangle)
    return 3;
  else if (*this == ReferenceCells::Quadrilateral)
    return 4;
  else if (*this == ReferenceCells::Tetrahedron)
    return 4;
  else if (*this == ReferenceCells::Pyramid)
    return 5;
  else if (*this == ReferenceCells::Wedge)
    return 5;
  else if (*this == ReferenceCells::Hexahedron)
    return 6;

  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}



inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
ReferenceCell::vertex_indices() const
{
  return {0U, n_vertices()};
}



inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
ReferenceCell::line_indices() const
{
  return {0U, n_lines()};
}



inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
ReferenceCell::face_indices() const
{
  return {0U, n_faces()};
}



inline ReferenceCell
ReferenceCell::face_reference_cell(const unsigned int face_no) const
{
  AssertIndexRange(face_no, n_faces());

  if (*this == ReferenceCells::Vertex)
    return ReferenceCells::Invalid;
  else if (*this == ReferenceCells::Line)
    return ReferenceCells::Vertex;
  else if (*this == ReferenceCells::Triangle)
    return ReferenceCells::Line;
  else if (*this == ReferenceCells::Quadrilateral)
    return ReferenceCells::Line;
  else if (*this == ReferenceCells::Tetrahedron)
    return ReferenceCells::Triangle;
  else if (*this == ReferenceCells::Pyramid)
    {
      if (face_no == 0)
        return ReferenceCells::Quadrilateral;
      else
        return ReferenceCells::Triangle;
    }
  else if (*this == ReferenceCells::Wedge)
    {
      if (face_no > 1)
        return ReferenceCells::Quadrilateral;
      else
        return ReferenceCells::Triangle;
    }
  else if (*this == ReferenceCells::Hexahedron)
    return ReferenceCells::Quadrilateral;

  Assert(false, ExcNotImplemented());
  return ReferenceCells::Invalid;
}



inline unsigned int
ReferenceCell::child_cell_on_face(const unsigned int face,
                                  const unsigned int subface) const
{
  AssertIndexRange(face, n_faces());

  if (*this == ReferenceCells::Vertex)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Line)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Triangle)
    {
      static constexpr unsigned int subcells[3][2] = {{0, 1}, {1, 2}, {2, 0}};

      return subcells[face][subface];
    }
  else if (*this == ReferenceCells::Quadrilateral)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Tetrahedron)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Pyramid)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Wedge)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Hexahedron)
    {
      Assert(false, ExcNotImplemented());
    }

  Assert(false, ExcNotImplemented());
  return {};
}



inline std::array<unsigned int, 2>
ReferenceCell::standard_vertex_to_face_and_vertex_index(
  const unsigned int vertex) const
{
  AssertIndexRange(vertex, n_vertices());

  if (*this == ReferenceCells::Vertex)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Line)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Triangle)
    {
      static const std::array<std::array<unsigned int, 2>, 3> table = {
        {{{0, 0}}, {{0, 1}}, {{1, 1}}}};

      return table[vertex];
    }
  else if (*this == ReferenceCells::Quadrilateral)
    {
      return GeometryInfo<2>::standard_quad_vertex_to_line_vertex_index(vertex);
    }
  else if (*this == ReferenceCells::Tetrahedron)
    {
      static const std::array<unsigned int, 2> table[4] = {{{0, 0}},
                                                           {{0, 1}},
                                                           {{0, 2}},
                                                           {{1, 2}}};

      return table[vertex];
    }
  else if (*this == ReferenceCells::Pyramid)
    {
      static const std::array<unsigned int, 2> table[5] = {
        {{0, 0}}, {{0, 1}}, {{0, 2}}, {{0, 3}}, {{1, 2}}};

      return table[vertex];
    }
  else if (*this == ReferenceCells::Wedge)
    {
      static const std::array<std::array<unsigned int, 2>, 6> table = {
        {{{0, 1}}, {{0, 0}}, {{0, 2}}, {{1, 0}}, {{1, 1}}, {{1, 2}}}};

      return table[vertex];
    }
  else if (*this == ReferenceCells::Hexahedron)
    {
      return GeometryInfo<3>::standard_hex_vertex_to_quad_vertex_index(vertex);
    }

  Assert(false, ExcNotImplemented());
  return {};
}



inline std::array<unsigned int, 2>
ReferenceCell::standard_line_to_face_and_line_index(
  const unsigned int line) const
{
  AssertIndexRange(line, n_lines());

  if (*this == ReferenceCells::Vertex)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Line)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Triangle)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Quadrilateral)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Tetrahedron)
    {
      static const std::array<unsigned int, 2> table[6] = {
        {{0, 0}}, {{0, 1}}, {{0, 2}}, {{1, 1}}, {{1, 2}}, {{2, 1}}};

      return table[line];
    }
  else if (*this == ReferenceCells::Pyramid)
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
  else if (*this == ReferenceCells::Wedge)
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
  else if (*this == ReferenceCells::Hexahedron)
    {
      return GeometryInfo<3>::standard_hex_line_to_quad_line_index(line);
    }

  Assert(false, ExcNotImplemented());
  return {};
}



inline unsigned int
ReferenceCell::face_to_cell_lines(const unsigned int  face,
                                  const unsigned int  line,
                                  const unsigned char face_orientation) const
{
  AssertIndexRange(face, n_faces());
  AssertIndexRange(line, face_reference_cell(face).n_lines());

  if (*this == ReferenceCells::Vertex)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Line)
    {
      return GeometryInfo<1>::face_to_cell_lines(
        face,
        line,
        Utilities::get_bit(face_orientation, 0),
        Utilities::get_bit(face_orientation, 2),
        Utilities::get_bit(face_orientation, 1));
    }
  else if (*this == ReferenceCells::Triangle)
    {
      return face;
    }
  else if (*this == ReferenceCells::Quadrilateral)
    {
      return GeometryInfo<2>::face_to_cell_lines(
        face,
        line,
        Utilities::get_bit(face_orientation, 0),
        Utilities::get_bit(face_orientation, 2),
        Utilities::get_bit(face_orientation, 1));
    }
  else if (*this == ReferenceCells::Tetrahedron)
    {
      const static unsigned int table[4][3] = {{0, 1, 2},
                                               {0, 3, 4},
                                               {2, 5, 3},
                                               {1, 4, 5}};

      return table[face]
                  [standard_to_real_face_line(line, face, face_orientation)];
    }
  else if (*this == ReferenceCells::Pyramid)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Wedge)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Hexahedron)
    {
      return GeometryInfo<3>::face_to_cell_lines(
        face,
        line,
        Utilities::get_bit(face_orientation, 0),
        Utilities::get_bit(face_orientation, 2),
        Utilities::get_bit(face_orientation, 1));
    }

  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}



inline unsigned int
ReferenceCell::face_to_cell_vertices(const unsigned int  face,
                                     const unsigned int  vertex,
                                     const unsigned char face_orientation) const
{
  AssertIndexRange(face, n_faces());
  AssertIndexRange(vertex, face_reference_cell(face).n_vertices());

  if (*this == ReferenceCells::Vertex)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Line)
    {
      return GeometryInfo<1>::face_to_cell_vertices(
        face,
        vertex,
        Utilities::get_bit(face_orientation, 0),
        Utilities::get_bit(face_orientation, 2),
        Utilities::get_bit(face_orientation, 1));
    }
  else if (*this == ReferenceCells::Triangle)
    {
      static const unsigned int table[3][2] = {{0, 1}, {1, 2}, {2, 0}};

      return table[face][face_orientation ? vertex : (1 - vertex)];
    }
  else if (*this == ReferenceCells::Quadrilateral)
    {
      return GeometryInfo<2>::face_to_cell_vertices(
        face,
        vertex,
        Utilities::get_bit(face_orientation, 0),
        Utilities::get_bit(face_orientation, 2),
        Utilities::get_bit(face_orientation, 1));
    }
  else if (*this == ReferenceCells::Tetrahedron)
    {
      static const unsigned int table[4][3] = {{0, 1, 2},
                                               {1, 0, 3},
                                               {0, 2, 3},
                                               {2, 1, 3}};

      return table[face][standard_to_real_face_vertex(
        vertex, face, face_orientation)];
    }
  else if (*this == ReferenceCells::Pyramid)
    {
      constexpr auto            X           = numbers::invalid_unsigned_int;
      static const unsigned int table[5][4] = {
        {0, 1, 2, 3}, {0, 2, 4, X}, {3, 1, 4, X}, {1, 0, 4, X}, {2, 3, 4, X}};

      return table[face][standard_to_real_face_vertex(
        vertex, face, face_orientation)];
    }
  else if (*this == ReferenceCells::Wedge)
    {
      constexpr auto            X           = numbers::invalid_unsigned_int;
      static const unsigned int table[6][4] = {
        {1, 0, 2, X}, {3, 4, 5, X}, {0, 1, 3, 4}, {1, 2, 4, 5}, {2, 0, 5, 3}};

      return table[face][standard_to_real_face_vertex(
        vertex, face, face_orientation)];
    }
  else if (*this == ReferenceCells::Hexahedron)
    {
      return GeometryInfo<3>::face_to_cell_vertices(
        face,
        vertex,
        Utilities::get_bit(face_orientation, 0),
        Utilities::get_bit(face_orientation, 2),
        Utilities::get_bit(face_orientation, 1));
    }

  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}



inline unsigned int
ReferenceCell::standard_to_real_face_vertex(
  const unsigned int  vertex,
  const unsigned int  face,
  const unsigned char face_orientation) const
{
  AssertIndexRange(face, n_faces());
  AssertIndexRange(vertex, face_reference_cell(face).n_vertices());

  if (*this == ReferenceCells::Vertex)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Line)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Triangle)
    {
      static const unsigned int table[2][2] = {{1, 0}, {0, 1}};

      return table[face_orientation][vertex];
    }
  else if (*this == ReferenceCells::Quadrilateral)
    {
      return GeometryInfo<2>::standard_to_real_line_vertex(vertex,
                                                           face_orientation);
    }
  else if (*this == ReferenceCells::Tetrahedron)
    {
      static const unsigned int table[6][3] = {
        {0, 2, 1}, {0, 1, 2}, {2, 1, 0}, {1, 2, 0}, {1, 0, 2}, {2, 0, 1}};

      return table[face_orientation][vertex];
    }
  else if (*this == ReferenceCells::Pyramid)
    {
      if (face == 0) // Quad
        {
          return GeometryInfo<3>::standard_to_real_face_vertex(
            vertex,
            Utilities::get_bit(face_orientation, 0),
            Utilities::get_bit(face_orientation, 2),
            Utilities::get_bit(face_orientation, 1));
        }
      else // Tri
        {
          static const unsigned int table[6][3] = {
            {0, 2, 1}, {0, 1, 2}, {2, 1, 0}, {1, 2, 0}, {1, 0, 2}, {2, 0, 1}};

          return table[face_orientation][vertex];
        }
    }
  else if (*this == ReferenceCells::Wedge)
    {
      if (face > 1) // QUAD
        {
          return GeometryInfo<3>::standard_to_real_face_vertex(
            vertex,
            Utilities::get_bit(face_orientation, 0),
            Utilities::get_bit(face_orientation, 2),
            Utilities::get_bit(face_orientation, 1));
        }
      else // TRI
        {
          static const unsigned int table[6][3] = {
            {0, 2, 1}, {0, 1, 2}, {2, 1, 0}, {1, 2, 0}, {1, 0, 2}, {2, 0, 1}};

          return table[face_orientation][vertex];
        }
    }
  else if (*this == ReferenceCells::Hexahedron)
    {
      return GeometryInfo<3>::standard_to_real_face_vertex(
        vertex,
        Utilities::get_bit(face_orientation, 0),
        Utilities::get_bit(face_orientation, 2),
        Utilities::get_bit(face_orientation, 1));
    }

  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}



inline unsigned int
ReferenceCell::standard_to_real_face_line(
  const unsigned int  line,
  const unsigned int  face,
  const unsigned char face_orientation) const
{
  AssertIndexRange(face, n_faces());
  AssertIndexRange(line, face_reference_cell(face).n_lines());

  if (*this == ReferenceCells::Vertex)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Line)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Triangle)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Quadrilateral)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Tetrahedron)
    {
      static const unsigned int table[6][3] = {
        {2, 1, 0}, {0, 1, 2}, {1, 0, 2}, {1, 2, 0}, {0, 2, 1}, {2, 0, 1}};

      return table[face_orientation][line];
    }
  else if (*this == ReferenceCells::Pyramid)
    {
      if (face == 0) // QUAD
        {
          return GeometryInfo<3>::standard_to_real_face_line(
            line,
            Utilities::get_bit(face_orientation, 0),
            Utilities::get_bit(face_orientation, 2),
            Utilities::get_bit(face_orientation, 1));
        }
      else // TRI
        {
          static const unsigned int table[6][3] = {
            {2, 1, 0}, {0, 1, 2}, {1, 0, 2}, {1, 2, 0}, {0, 2, 1}, {2, 0, 1}};

          return table[face_orientation][line];
        }
    }
  else if (*this == ReferenceCells::Wedge)
    {
      if (face > 1) // QUAD
        {
          return GeometryInfo<3>::standard_to_real_face_line(
            line,
            Utilities::get_bit(face_orientation, 0),
            Utilities::get_bit(face_orientation, 2),
            Utilities::get_bit(face_orientation, 1));
        }
      else // TRI
        {
          static const unsigned int table[6][3] = {
            {2, 1, 0}, {0, 1, 2}, {1, 0, 2}, {1, 2, 0}, {0, 2, 1}, {2, 0, 1}};

          return table[face_orientation][line];
        }
    }
  else if (*this == ReferenceCells::Hexahedron)
    {
      return GeometryInfo<3>::standard_to_real_face_line(
        line,
        Utilities::get_bit(face_orientation, 0),
        Utilities::get_bit(face_orientation, 2),
        Utilities::get_bit(face_orientation, 1));
    }

  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}



namespace ReferenceCells
{
  template <int dim>
  inline constexpr const ReferenceCell &
  get_simplex()
  {
    switch (dim)
      {
        case 0:
          return ReferenceCells::Vertex;
        case 1:
          return ReferenceCells::Line;
        case 2:
          return ReferenceCells::Triangle;
        case 3:
          return ReferenceCells::Tetrahedron;
        default:
          Assert(false, ExcNotImplemented());
          return ReferenceCells::Invalid;
      }
  }



  template <int dim>
  inline constexpr const ReferenceCell &
  get_hypercube()
  {
    switch (dim)
      {
        case 0:
          return ReferenceCells::Vertex;
        case 1:
          return ReferenceCells::Line;
        case 2:
          return ReferenceCells::Quadrilateral;
        case 3:
          return ReferenceCells::Hexahedron;
        default:
          Assert(false, ExcNotImplemented());
          return ReferenceCells::Invalid;
      }
  }
} // namespace ReferenceCells


inline ReferenceCell
ReferenceCell::n_vertices_to_type(const int dim, const unsigned int n_vertices)
{
  AssertIndexRange(dim, 4);
  AssertIndexRange(n_vertices, 9);

  const auto X = ReferenceCells::Invalid;
  static const std::array<std::array<ReferenceCell, 9>,
                          4>
    table = {{// dim 0
              {{X, ReferenceCells::Vertex, X, X, X, X, X, X, X}},
              // dim 1
              {{X, X, ReferenceCells::Line, X, X, X, X, X, X}},
              // dim 2
              {{X,
                X,
                X,
                ReferenceCells::Triangle,
                ReferenceCells::Quadrilateral,
                X,
                X,
                X,
                X}},
              // dim 3
              {{X,
                X,
                X,
                X,
                ReferenceCells::Tetrahedron,
                ReferenceCells::Pyramid,
                ReferenceCells::Wedge,
                X,
                ReferenceCells::Hexahedron}}}};
  Assert(table[dim][n_vertices] != ReferenceCells::Invalid,
         ExcMessage("The combination of dim = " + std::to_string(dim) +
                    " and n_vertices = " + std::to_string(n_vertices) +
                    " does not correspond to a known reference cell type."));
  return table[dim][n_vertices];
}



template <int dim>
inline double
ReferenceCell::d_linear_shape_function(const Point<dim> & xi,
                                       const unsigned int i) const
{
  AssertDimension(dim, get_dimension());
  if (*this == ReferenceCells::get_hypercube<dim>())
    return GeometryInfo<dim>::d_linear_shape_function(xi, i);

  if (*this ==
      ReferenceCells::Triangle) // see also
                                // Simplex::ScalarPolynomial::compute_value
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

  if (*this ==
      ReferenceCells::Tetrahedron) // see also
                                   // Simplex::ScalarPolynomial::compute_value
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

  if (*this ==
      ReferenceCells::Wedge) // see also
                             // Simplex::ScalarWedgePolynomial::compute_value
    {
      return ReferenceCell(ReferenceCells::Triangle)
               .d_linear_shape_function<2>(Point<2>(xi[std::min(0, dim - 1)],
                                                    xi[std::min(1, dim - 1)]),
                                           i % 3) *
             ReferenceCell(ReferenceCells::Line)
               .d_linear_shape_function<1>(Point<1>(xi[std::min(2, dim - 1)]),
                                           i / 3);
    }

  if (*this == ReferenceCells::
                 Pyramid) // see also
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



template <int dim>
inline Tensor<1, dim>
ReferenceCell::d_linear_shape_function_gradient(const Point<dim> & xi,
                                                const unsigned int i) const
{
  AssertDimension(dim, get_dimension());
  if (*this == ReferenceCells::get_hypercube<dim>())
    return GeometryInfo<dim>::d_linear_shape_function_gradient(xi, i);

  if (*this ==
      ReferenceCells::Triangle) // see also
                                // Simplex::ScalarPolynomial::compute_grad
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


template <int dim>
inline Tensor<1, dim>
ReferenceCell::unit_tangential_vectors(const unsigned int face_no,
                                       const unsigned int i) const
{
  AssertDimension(dim, get_dimension());
  AssertIndexRange(i, dim - 1);

  if (*this == ReferenceCells::get_hypercube<dim>())
    {
      AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);
      return GeometryInfo<dim>::unit_tangential_vectors[face_no][i];
    }
  else if (*this == ReferenceCells::Triangle)
    {
      AssertIndexRange(face_no, 3);
      static const std::array<Tensor<1, dim>, 3> table = {
        {Point<dim>(1, 0),
         Point<dim>(-std::sqrt(0.5), +std::sqrt(0.5)),
         Point<dim>(0, -1)}};

      return table[face_no];
    }
  else if (*this == ReferenceCells::Tetrahedron)
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
  else if (*this == ReferenceCells::Wedge)
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
  else if (*this == ReferenceCells::Pyramid)
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



template <int dim>
inline Tensor<1, dim>
ReferenceCell::unit_normal_vectors(const unsigned int face_no) const
{
  AssertDimension(dim, this->get_dimension());

  if (is_hyper_cube())
    {
      AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);
      return GeometryInfo<dim>::unit_normal_vector[face_no];
    }
  else if (dim == 2)
    {
      Assert(*this == ReferenceCells::Triangle, ExcInternalError());

      // Return the rotated vector
      return cross_product_2d(unit_tangential_vectors<dim>(face_no, 0));
    }
  else if (dim == 3)
    {
      return cross_product_3d(unit_tangential_vectors<dim>(face_no, 0),
                              unit_tangential_vectors<dim>(face_no, 1));
    }

  Assert(false, ExcNotImplemented());

  return {};
}


namespace internal
{
  /**
   * A namespace for geometric information on reference cells.
   */
  namespace ReferenceCell
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

    private:
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
      standard_vertex_to_face_and_vertex_index(const unsigned int vertex) const
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

    public:
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

    private:
      /**
       * Return reference-cell type of face @p face_no.
       */
      virtual dealii::ReferenceCell
      face_reference_cell(const unsigned int face_no) const
      {
        Assert(false, ExcNotImplemented());
        (void)face_no;

        return ReferenceCells::Invalid;
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

    public:
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

    private:
      /**
       * Indices of child cells that are adjacent to a certain face of the
       * mother cell.
       */
      virtual unsigned int
      child_cell_on_face(const unsigned int face_n,
                         const unsigned int subface_n) const
      {
        Assert(false, ExcNotImplemented());
        (void)face_n;
        (void)subface_n;

        return numbers::invalid_unsigned_int;
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
          Utilities::get_bit(face_orientation, 0),
          Utilities::get_bit(face_orientation, 2),
          Utilities::get_bit(face_orientation, 1));
      }

      unsigned int
      face_to_cell_vertices(const unsigned int  face,
                            const unsigned int  vertex,
                            const unsigned char face_orientation) const override
      {
        return GeometryInfo<dim>::face_to_cell_vertices(
          face,
          vertex,
          Utilities::get_bit(face_orientation, 0),
          Utilities::get_bit(face_orientation, 2),
          Utilities::get_bit(face_orientation, 1));
      }
    };



    /*
     * Vertex.
     */
    struct Vertex : public TensorProductBase<0>
    {
      dealii::ReferenceCell
      face_reference_cell(const unsigned int face_no) const override
      {
        (void)face_no;
        return ReferenceCells::Invalid;
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
      dealii::ReferenceCell
      face_reference_cell(const unsigned int face_no) const override
      {
        (void)face_no;
        return ReferenceCells::Vertex;
      }

      virtual unsigned int
      exodusii_vertex_to_deal_vertex(const unsigned int vertex_n) const override
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
    struct Triangle : public Base
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

      dealii::ReferenceCell
      face_reference_cell(const unsigned int face_no) const override
      {
        (void)face_no;

        AssertIndexRange(face_no, n_faces());

        return dealii::ReferenceCells::Line;
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
      face_to_cell_vertices(const unsigned int  face,
                            const unsigned int  vertex,
                            const unsigned char face_orientation) const override
      {
        static const std::array<std::array<unsigned int, 2>, 3> table = {
          {{{0, 1}}, {{1, 2}}, {{2, 0}}}};

        return table[face][face_orientation ? vertex : (1 - vertex)];
      }

      virtual unsigned int
      exodusii_vertex_to_deal_vertex(const unsigned int vertex_n) const override
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

      virtual unsigned int
      child_cell_on_face(const unsigned int face_n,
                         const unsigned int subface_n) const override
      {
        static constexpr unsigned int subcells[3][2] = {{0, 1}, {1, 2}, {2, 0}};

        return subcells[face_n][subface_n];
      }
    };



    /**
     * Quadrilateral
     */
    struct Quadrilateral : public TensorProductBase<2>
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

        return GeometryInfo<2>::standard_to_real_line_vertex(vertex,
                                                             line_orientation);
      }

      dealii::ReferenceCell
      face_reference_cell(const unsigned int face_no) const override
      {
        (void)face_no;
        return dealii::ReferenceCells::Line;
      }

      virtual unsigned int
      exodusii_vertex_to_deal_vertex(const unsigned int vertex_n) const override
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
     * Tetrahedron
     */
    struct Tetrahedron : public Base
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

      dealii::ReferenceCell
      face_reference_cell(const unsigned int face_no) const override
      {
        (void)face_no;

        AssertIndexRange(face_no, n_faces());

        return dealii::ReferenceCells::Triangle;
      }

      unsigned int
      face_to_cell_lines(const unsigned int  face,
                         const unsigned int  line,
                         const unsigned char face_orientation) const override
      {
        AssertIndexRange(face, n_faces());

        const static std::array<std::array<unsigned int, 3>, 4> table = {
          {{{0, 1, 2}}, {{0, 3, 4}}, {{2, 5, 3}}, {{1, 4, 5}}}};

        return table[face]
                    [standard_to_real_face_line(line, face, face_orientation)];
      }

      unsigned int
      face_to_cell_vertices(const unsigned int  face,
                            const unsigned int  vertex,
                            const unsigned char face_orientation) const override
      {
        static const std::array<std::array<unsigned int, 3>, 4> table = {
          {{{0, 1, 2}}, {{1, 0, 3}}, {{0, 2, 3}}, {{2, 1, 3}}}};

        return table[face][standard_to_real_face_vertex(
          vertex, face, face_orientation)];
      }

      virtual unsigned int
      exodusii_vertex_to_deal_vertex(const unsigned int vertex_n) const override
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
              Utilities::get_bit(face_orientation, 0),
              Utilities::get_bit(face_orientation, 2),
              Utilities::get_bit(face_orientation, 1));
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
              Utilities::get_bit(face_orientation, 0),
              Utilities::get_bit(face_orientation, 2),
              Utilities::get_bit(face_orientation, 1));
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

      dealii::ReferenceCell
      face_reference_cell(const unsigned int face_no) const override
      {
        AssertIndexRange(face_no, n_faces());

        if (face_no == 0)
          return dealii::ReferenceCells::Quadrilateral;
        else
          return dealii::ReferenceCells::Triangle;
      }

      unsigned int
      face_to_cell_vertices(const unsigned int  face,
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
      exodusii_vertex_to_deal_vertex(const unsigned int vertex_n) const override
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
              Utilities::get_bit(face_orientation, 0),
              Utilities::get_bit(face_orientation, 2),
              Utilities::get_bit(face_orientation, 1));
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
              Utilities::get_bit(face_orientation, 0),
              Utilities::get_bit(face_orientation, 2),
              Utilities::get_bit(face_orientation, 1));
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

      dealii::ReferenceCell
      face_reference_cell(const unsigned int face_no) const override
      {
        AssertIndexRange(face_no, n_faces());

        if (face_no > 1)
          return dealii::ReferenceCells::Quadrilateral;
        else
          return dealii::ReferenceCells::Triangle;
      }

      unsigned int
      face_to_cell_vertices(const unsigned int  face,
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
      exodusii_vertex_to_deal_vertex(const unsigned int vertex_n) const override
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
    struct Hexahedron : public TensorProductBase<3>
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
          Utilities::get_bit(face_orientation, 0),
          Utilities::get_bit(face_orientation, 2),
          Utilities::get_bit(face_orientation, 1));
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

        const bool face_orientation =
          Utilities::get_bit(face_orientation_raw, 0);
        const bool face_flip     = Utilities::get_bit(face_orientation_raw, 2);
        const bool face_rotation = Utilities::get_bit(face_orientation_raw, 1);

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
          Utilities::get_bit(face_orientation, 0),
          Utilities::get_bit(face_orientation, 2),
          Utilities::get_bit(face_orientation, 1));
      }

      dealii::ReferenceCell
      face_reference_cell(const unsigned int face_no) const override
      {
        (void)face_no;
        return dealii::ReferenceCells::Quadrilateral;
      }

      virtual unsigned int
      exodusii_vertex_to_deal_vertex(const unsigned int vertex_n) const override
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
    inline const internal::ReferenceCell::Base &
    get_cell(const dealii::ReferenceCell &type)
    {
      static const std::array<std::unique_ptr<internal::ReferenceCell::Base>, 8>
        gei{{std::make_unique<internal::ReferenceCell::Vertex>(),
             std::make_unique<internal::ReferenceCell::Line>(),
             std::make_unique<internal::ReferenceCell::Triangle>(),
             std::make_unique<internal::ReferenceCell::Quadrilateral>(),
             std::make_unique<internal::ReferenceCell::Tetrahedron>(),
             std::make_unique<internal::ReferenceCell::Pyramid>(),
             std::make_unique<internal::ReferenceCell::Wedge>(),
             std::make_unique<internal::ReferenceCell::Hexahedron>()}};
      AssertIndexRange(static_cast<std::uint8_t>(type), 8);
      return *gei[static_cast<std::uint8_t>(type)];
    }

    /**
     * Return for a given reference-cell type @p and face number @p face_no the
     * right Info of the @p face_no-th face.
     */
    inline const internal::ReferenceCell::Base &
    get_face(const dealii::ReferenceCell &type, const unsigned int face_no)
    {
      return get_cell(type.face_reference_cell(face_no));
    }

  } // namespace ReferenceCell
} // namespace internal



namespace internal
{
  template <typename T, std::size_t N>
  class NoPermutation : public dealii::ExceptionBase
  {
  public:
    /**
     * Constructor.
     */
    NoPermutation(const dealii::ReferenceCell &entity_type,
                  const std::array<T, N> &     vertices_0,
                  const std::array<T, N> &     vertices_1)
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

      const unsigned int n_vertices = entity_type.n_vertices();

      for (unsigned int i = 0; i < n_vertices; ++i)
        {
          out << vertices_0[i];
          if (i + 1 != n_vertices)
            out << ",";
        }

      out << "] is not a permutation of [";

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
    const dealii::ReferenceCell entity_type;

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



template <typename T, std::size_t N>
inline unsigned char
ReferenceCell::compute_orientation(const std::array<T, N> &vertices_0,
                                   const std::array<T, N> &vertices_1) const
{
  AssertIndexRange(n_vertices(), N + 1);
  if (*this == ReferenceCells::Line)
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
  else if (*this == ReferenceCells::Triangle)
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
  else if (*this == ReferenceCells::Quadrilateral)
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

  Assert(false, (internal::NoPermutation<T, N>(*this, vertices_0, vertices_1)));

  return -1;
}



template <typename T, std::size_t N>
inline std::array<T, N>
ReferenceCell::permute_according_orientation(
  const std::array<T, N> &vertices,
  const unsigned int      orientation) const
{
  std::array<T, 4> temp;

  if (*this == ReferenceCells::Line)
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
  else if (*this == ReferenceCells::Triangle)
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
  else if (*this == ReferenceCells::Quadrilateral)
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


DEAL_II_NAMESPACE_CLOSE

#endif

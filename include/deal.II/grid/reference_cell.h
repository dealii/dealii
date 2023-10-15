// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2023 by the deal.II authors
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
#include <deal.II/base/ndarray.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/tria_orientation.h>

#include <boost/container/small_vector.hpp>

#include <iosfwd>
#include <string>

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
  /**
   * A helper function to create a ReferenceCell object from an integer.
   * ReferenceCell objects are "singletons" (actually, "multitons" -- there are
   * multiple, but they are only a handful and these are all that can be used).
   * What is then necessary is to have a way to create these with their internal
   * id to distinguish the few possible ones in existence. We could do this via
   * a public constructor of ReferenceCell, but that would allow users to create
   * ones outside the range we envision, and we don't want to do that. Rather,
   * the constructor that takes an integer is made `private` but we have this
   * one function in an internal namespace that is a friend of the class and can
   * be used to create the objects.
   */
  constexpr ReferenceCell
  make_reference_cell_from_int(const std::uint8_t kind);
} // namespace internal



/**
 * A type that describes the kinds of reference cells that can be
 * used.  This includes quadrilaterals and hexahedra (i.e.,
 * "hypercubes"), triangles and tetrahedra (simplices), and the
 * pyramids and wedges necessary when using mixed 3d meshes. This
 * class then describes geometric, topological, and other kinds of
 * information about these kinds of reference cells. This includes how
 * many vertices or faces a certain kind of reference cell has
 * (topological information), where these vertices lie, what the
 * cell's volume or center of mass is (geometric information), and how
 * to output these cells in various output formats or what appropriate
 * quadrature rules are. The documentation of this class is separated
 * into a number of sections to group the many member functions into
 * different categories such as those mentioned above.
 *
 * Objects of this type should not be created in user code, and as a
 * consequence the class does not have a user-accessible constructor
 * other than the default constructor (which creates an invalid object).
 * Rather, there is a finite number of specific reference cell objects
 * defined in the ReferenceCells namespace that completely enumerate
 * all of the possible values. User codes should therefore rely
 * exclusively on assigning ReferenceCell objects from these special
 * objects, and comparing against those special objects.
 *
 * The purposes and intents of this class are described in the
 * @ref GlossReferenceCell "reference cell"
 * glossary entry.
 *
 * @ingroup grid geomprimitives aniso reordering
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
  constexpr ReferenceCell();

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
  d_linear_shape_function_gradient(const Point<dim>  &xi,
                                   const unsigned int i) const;

  /**
   * Return a default mapping of degree @p degree matching the current
   * reference cell. If this reference cell is a hypercube, then the returned
   * mapping is a MappingQ; otherwise, it is an object of type
   * MappingFE initialized with FE_SimplexP (if the reference cell is a
   * triangle or tetrahedron), with FE_PyramidP (if the reference
   * cell is a pyramid), or with FE_WedgeP (if the reference cell is
   * a wedge).
   */
  template <int dim, int spacedim = dim>
  std::unique_ptr<Mapping<dim, spacedim>>
  get_default_mapping(const unsigned int degree) const;

  /**
   * Return a default linear mapping matching the current reference cell.
   * If this reference cell is a hypercube, then the returned mapping
   * is a MappingQ1; otherwise, it is an object of type MappingFE
   * initialized with FE_SimplexP (if the reference cell is a triangle or
   * tetrahedron), with FE_PyramidP (if the reference cell is a
   * pyramid), or with FE_WedgeP (if the reference cell is a wedge).
   * In other words, the term "linear" in the name of the function has to be
   * understood as $d$-linear (i.e., bilinear or trilinear) for some of the
   * coordinate directions.
   */
  template <int dim, int spacedim = dim>
  const Mapping<dim, spacedim> &
  get_default_linear_mapping() const;

  /**
   * Return a Gauss-type quadrature matching the given reference cell (QGauss,
   * QGaussSimplex, QGaussPyramid, QGaussWedge).
   *
   * @param[in] n_points_1d The number of quadrature points in each direction
   * (QGauss) or an indication of what polynomial degree needs to be
   * integrated exactly for the other types.
   */
  template <int dim>
  Quadrature<dim>
  get_gauss_type_quadrature(const unsigned n_points_1d) const;

  /**
   * Return a quadrature object that has a single quadrature point at the
   * barycenter of the cell with quadrature weight equal to the volume of the
   * reference cell. This quadrature formula is exact for integrals of constant
   * and linear integrands.
   *
   * The object returned by this function generalizes what the QMidpoint class
   * represents to other reference cells.
   */
  template <int dim>
  Quadrature<dim>
  get_midpoint_quadrature() const;

  /**
   * Return a quadrature rule whose quadrature points are the vertices of the
   * given reference cell. For 1d line segments, this corresponds to the
   * quadrature points of the trapezoidal rule, which by taking tensor products
   * easily generalizes also to other hypercube elements (see also QTrapezoid).
   * For all reference cell shapes, the quadrature points are ordered
   * in the same order as the vertices of the reference cell.
   *
   * @note The weights of the quadrature object are left unfilled and
   *   consequently the object cannot usefully be used for actually
   *   computing integrals. This is in contrast to, for example, the QTrapezoid
   *   class that correctly sets quadrature weights.
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
   * Return the location of the `v`th vertex of the reference
   * cell that corresponds to the current object.
   *
   * Because the ReferenceCell class does not have a `dim` argument,
   * it has to be explicitly specified in the call to this function.
   */
  template <int dim>
  Point<dim>
  vertex(const unsigned int v) const;

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
   * Return the number of cells one would get by isotropically
   * refining the current cell. Here, "isotropic refinement"
   * means that we subdivide in each "direction" of a cell.
   * For example, a square would be refined into four children
   * by introducing new vertices along each edge and a new
   * vertex in the cell center. For triangles, one would introduce
   * new vertices at the center of each edge, and connect them to
   * obtain four children. Similar constructions can be done for
   * the other reference cell types.
   */
  unsigned int
  n_isotropic_children() const;

  /**
   * Return an object that can be thought of as an array containing all
   * indices from zero to n_isotropic_children().
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  isotropic_child_indices() const;

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
   * Return the default combined face orientation flag (i.e., the default set of
   * orientations, defined by orientation, rotate, and flip for a face in 3d).
   */
  static constexpr unsigned char
  default_combined_face_orientation();

  /**
   * Return the reversed (non-default orientation) line orientation flag. As
   * lines only have two possible orientations, this function and
   * ReferenceCell::default_combined_face_orientation() encode all of its
   * possible orientation states.
   *
   * @note Line orientations are typically stored as booleans, but to better
   * enable dimension-independent programming relevant functions typically
   * present these values as unsigned chars.
   */
  static constexpr unsigned char
  reversed_combined_line_orientation();

  /**
   * Return which child cells are adjacent to a certain face of the
   * mother cell.
   *
   * For example, in 2d the layout of a quadrilateral cell is as follows:
   * @verbatim
   *      3
   *   2-->--3
   *   |     |
   * 0 ^     ^ 1
   *   |     |
   *   0-->--1
   *      2
   * @endverbatim
   * Vertices and faces are indicated with their numbers, faces also with
   * their directions.
   *
   * Now, when refined, the layout is like this:
   * @verbatim
   * *---*---*
   * | 2 | 3 |
   * *---*---*
   * | 0 | 1 |
   * *---*---*
   * @endverbatim
   *
   * Thus, the child cells on face 0 are (ordered in the direction of the
   * face) 0 and 2, on face 3 they are 2 and 3, etc.
   *
   * For three spatial dimensions, the exact order of the children is laid
   * down in the general documentation of this class.
   *
   * The <tt>face_orientation</tt> argument is meant exclusively for
   * quadrilaterals and hexahedra at the moment. It determines how this function
   * handles faces oriented in the standard and non-standard orientation. It
   * represents a bit-code for the overall <tt>face_orientation</tt>,
   * <tt>face_flip</tt> and <tt>face_rotation</tt> and defaults to the standard
   * orientation. The concept of face orientations is explained in this
   * @ref GlossFaceOrientation "glossary"
   * entry.
   */
  unsigned int
  child_cell_on_face(const unsigned int  face,
                     const unsigned int  subface,
                     const unsigned char face_orientation =
                       default_combined_face_orientation()) const;

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
   * Map line vertex number to cell vertex number, i.e., return the cell vertex
   * number of the <tt>vertex</tt>th vertex of line <tt>line</tt>.
   *
   * The order of the lines, as well as their direction (which in turn
   * determines which vertices are first and second on a line) is the canonical
   * one in deal.II, as described in the general documentation of this class.
   */
  unsigned int
  line_to_cell_vertices(const unsigned int line,
                        const unsigned int vertex) const;

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
   * For a given face, in standard orientation, return the location of one
   * of its vertices in the ambient space of the cell. For example, for a
   * square or triangular 2d cell, the zeroth vertex of its zeroth face
   * is located at $(0,0)$ -- a location in 2d space.
   *
   * @param[in] face The number of face. This number must be between zero
   *   and `n_faces()`.
   * @param[in] vertex The number of the vertex within the face. This number
   *   must be between zero and `face_reference_cell(face).n_vertices()`.
   * @return The location of the vertex so identified in `dim`-dimensional
   *   space.
   *
   * @post The output of calling `reference_cell.face_vertex_location<dim>(f,v)`
   *   is identical to calling
   *   `reference_cell.vertex<dim>(
   *       reference_cell.face_to_cell_vertices(
                   f, v, ReferenceCell::default_combined_face_orientation()))`.
   */
  template <int dim>
  Point<dim>
  face_vertex_location(const unsigned int face,
                       const unsigned int vertex) const;

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
   * Return whether the line with index @p line is oriented in
   * standard direction within a cell, given the @p face_orientation of
   * the face within the current cell, and @p line_orientation flag
   * for the line within that face. @p true indicates that the line is
   * oriented from vertex 0 to vertex 1, whereas it is the other way
   * around otherwise. In 1d and 2d, this is always @p true, but in 3d
   * it may be different, see the respective discussion in the
   * documentation of the GeometryInfo class.
   */
  bool
  standard_vs_true_line_orientation(const unsigned int  line,
                                    const unsigned int  face,
                                    const unsigned char face_orientation,
                                    const bool          line_orientation) const;

  /**
   * @}
   */

  /**
   * @name Geometric properties of reference cells
   * @name Querying the number of building blocks of a reference cell
   * @{
   */

  /**
   * Return the $d$-dimensional volume of the reference cell that corresponds
   * to the current object, where $d$ is the dimension of the space it lives
   * in. For example, since the quadrilateral reference cell is $[0,1]^2$,
   * its volume is one, whereas the volume of the reference triangle is
   * 0.5 because it occupies the area $\{0 \le x,y \le 1, x+y\le 1\}$.
   *
   * For ReferenceCells::Vertex, the reference cell is a zero-dimensional
   * point in a zero-dimensional space. As a consequence, one cannot
   * meaningfully define a volume for it. The function returns one for
   * this case, because this makes it possible to define useful quadrature
   * rules based on the center of a reference cell and its volume.
   */
  double
  volume() const;

  /**
   * Return the barycenter (i.e., the center of mass) of the reference
   * cell that corresponds to the current object. The function is not
   * called `center()` because one can define the center of an object
   * in a number of different ways whereas the barycenter of a
   * reference cell $K$ is unambiguously defined as
   * @f[
   *   \mathbf x_K = \frac{1}{V} \int_K \mathbf x \; dx
   * @f]
   * where $V$ is the volume of the reference cell (see also the volume()
   * function).
   */
  template <int dim>
  Point<dim>
  barycenter() const;

  /**
   * Return true if the given point is inside the reference cell of the present
   * space dimension up to some tolerance. This function accepts an additional
   * parameter (which defaults to zero) which specifies by how much the point
   * position may actually be outside the true reference cell. This is useful
   * because in practice we may often not be able to compute the coordinates of
   * a point in reference coordinates exactly, but only up to numerical
   * roundoff. For example, strictly speaking one would expect that for points
   * on the boundary of the reference cell, the function would return `true` if
   * the tolerance was zero. But in practice, this may or may not actually be
   * true; for example, the point $(1/3, 2/3)$ is on the boundary of the
   * reference triangle because $1/3+2/3 \le 1$, but since neither of its
   * coordinates are exactly representable in floating point arithmetic, the
   * floating point representations of $1/3$ and $2/3$ may or may not add up to
   * anything that is less than or equal to one.
   *
   * The tolerance parameter may be less than zero, indicating that the point
   * should be safely inside the cell.
   */
  template <int dim>
  bool
  contains_point(const Point<dim> &p, const double tolerance = 0) const;

  /**
   * Return the point on the surface of the reference cell closest (in the
   * Euclidean norm) to @p p.
   */
  template <int dim>
  Point<dim>
  closest_point(const Point<dim> &p) const;

  /**
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
   * Return the number of orientations for a face in the ReferenceCell. For
   * example, for hexahedra this is 8 for every face since quadrilaterals have
   * 8 possible orientations.
   */
  unsigned int
  n_face_orientations(const unsigned int face_no) const;

  /**
   * Determine the orientation of the current entity described by its
   * vertices @p vertices_1 relative to an entity described by @p vertices_0.
   * The two arrays given as arguments can be arrays of global vertex
   * indices or local vertex indices, arrays of vertex locations, or
   * arrays of any other objects identifying the vertices and the order
   * in which they are encountered in a cell.
   *
   * The size of the arrays, i.e., the template argument `N`,
   * must be equal to or larger than the number of vertices of the current
   * entity. If it is larger, only those elements of the input and output
   * arrays are read from or written to that correspond to valid vertex
   * indices.
   *
   * @deprecated Use get_combined_orientation() instead.
   */
  template <typename T, std::size_t N>
  DEAL_II_DEPRECATED unsigned char
  compute_orientation(const std::array<T, N> &vertices_0,
                      const std::array<T, N> &vertices_1) const;

  /**
   * Determine the relative orientation of the current entity described by its
   * vertices @p vertices_1 relative to an entity described by @p vertices_0.
   * Relative orientations are special cases of permutations since every
   * vertex has to appear in the list of vertices of a reoriented cell
   * as well; however, not every permutation can denote the same cell:
   * For example, a square's vertices can be rotated by 90, 180, or
   * 270 degrees, and the cell can be inverted (in essence looking
   * at it from the other side), but one can't just exchange the order of two
   * adjacent vertices because then the resulting cell is no longer a square
   * but an object with two edges that cross each other.
   *
   * The two arrays given as arguments can be arrays of global vertex
   * indices or local vertex indices, arrays of vertex locations, or
   * arrays of any other objects identifying the vertices and the order
   * in which they are encountered in a cell.
   *
   * The size of the input arrays must be equal to the number of vertices of
   * the current entity.
   *
   * @returns A number that describes a relative orientation. How exactly
   *   this index is defined is not important, but it is consistent with the
   *   understanding the permute_by_combined_orientation() has of
   *   these orientation indices.
   */
  template <typename T>
  unsigned char
  get_combined_orientation(const ArrayView<const T> &vertices_0,
                           const ArrayView<const T> &vertices_1) const;


  /**
   * Inverse function of compute_orientation(): Given a set of
   * vertex-associated objects (such as vertex indices, locations, etc.) and
   * a desired orientation permutation, return the permuted vertex information.
   *
   * The size of the input and output arrays, i.e., the template argument `N`,
   * must be equal to or larger than the number of vertices of the current
   * entity. If it is larger, only those elements of the input and output
   * arrays are read from or written to that correspond to valid vertex
   * indices.
   *
   * @deprecated Use permute_by_combined_orientation() instead.
   */
  template <typename T, std::size_t N>
  DEAL_II_DEPRECATED std::array<T, N>
  permute_according_orientation(const std::array<T, N> &vertices,
                                const unsigned int      orientation) const;

  /**
   * This is the inverse function to get_combined_orientation(): Given a set of
   * vertex-associated objects (such as vertex indices, locations, etc.) and
   * a desired orientation permutation, return the permuted vertex information.
   *
   * The size of the input array must be equal to the number of vertices of
   * the current entity. The output is an array or permuted quantities of
   * the same size. It is a vector that can store up to and including as many
   * elements as cells can have vertices (namely eight, as in the case of
   * hexahedra in 3d).
   */
  template <typename T>
  boost::container::small_vector<T, 8>
  permute_by_combined_orientation(const ArrayView<const T> &vertices,
                                  const unsigned char       orientation) const;

  /**
   * Return a vector of faces a given @p vertex_index belongs to.
   */
  ArrayView<const unsigned int>
  faces_for_given_vertex(const unsigned int vertex_index) const;

  /**
   * @}
   */

  /**
   * @name Translating between deal.II indexing and formats used by other programs
   * @{
   */

  /**
   * Map an ExodusII vertex number to a deal.II vertex number.
   */
  unsigned int
  exodusii_vertex_to_deal_vertex(const unsigned int vertex_n) const;

  /**
   * Map an ExodusII face number to a deal.II face number.
   */
  unsigned int
  exodusii_face_to_deal_face(const unsigned int face_n) const;

  /**
   * Map a UNV vertex number to a deal.II vertex number.
   */
  unsigned int
  unv_vertex_to_deal_vertex(const unsigned int vertex_n) const;

  /**
   * Return a VTK linear shape constant that corresponds to the reference cell.
   */
  unsigned int
  vtk_linear_type() const;

  /**
   * Return a VTK quadratic shape constant that corresponds to the reference
   * cell.
   */
  unsigned int
  vtk_quadratic_type() const;

  /**
   * Return a VTK Lagrange shape constant that corresponds to the reference
   * cell.
   */
  unsigned int
  vtk_lagrange_type() const;

  /**
   * Given a set of node indices of the form $(i)$ or $(i,j)$ or $(i,j,k)$
   * (depending on whether the reference cell is in 1d, 2d, or 3d), return
   * the index the VTK format uses for this node for cells that are
   * subdivided as many times in each of the coordinate directions as
   * described by the second argument. For a uniformly subdivided cell,
   * the second argument is an array whose elements will all be equal.
   *
   * The last argument, @p legacy_format, indicates whether to use the
   * old, VTK legacy format (when `true`) or the new, VTU format (when
   * `false`).
   */
  template <int dim>
  unsigned int
  vtk_lexicographic_to_node_index(
    const std::array<unsigned, dim> &node_indices,
    const std::array<unsigned, dim> &nodes_per_direction,
    const bool                       legacy_format) const;

  /**
   * Map a VTK vertex number to a deal.II vertex number.
   */
  unsigned int
  vtk_vertex_to_deal_vertex(const unsigned int vertex_index) const;

  /**
   * Return the GMSH element type code that corresponds to the reference cell.
   */
  unsigned int
  gmsh_element_type() const;

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
   * Return the number of bytes used by an instance of this class.
   */
  static constexpr std::size_t
  memory_consumption();

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
   * Table containing all vertex permutations for a line.
   */
  static constexpr ndarray<unsigned int, 2, 2> line_vertex_permutations = {
    {{{1, 0}}, {{0, 1}}}};


  /**
   * Table containing all vertex permutations for a triangle.
   */
  static constexpr ndarray<unsigned int, 6, 3> triangle_vertex_permutations = {
    {{{0, 2, 1}},
     {{0, 1, 2}},
     {{2, 1, 0}},
     {{2, 0, 1}},
     {{1, 0, 2}},
     {{1, 2, 0}}}};

  /**
   * Table containing all vertex permutations for a quadrilateral.
   */
  static constexpr ndarray<unsigned int, 8, 4>
    quadrilateral_vertex_permutations = {{{{0, 2, 1, 3}},
                                          {{0, 1, 2, 3}},
                                          {{2, 3, 0, 1}},
                                          {{2, 0, 3, 1}},
                                          {{3, 1, 2, 0}},
                                          {{3, 2, 1, 0}},
                                          {{1, 0, 3, 2}},
                                          {{1, 3, 0, 2}}}};

  /**
   * A kind of constructor -- not quite private because it can be
   * called by anyone, but at least hidden in an internal namespace.
   */
  friend constexpr ReferenceCell
  internal::make_reference_cell_from_int(const std::uint8_t);

  friend std::ostream &
  operator<<(std::ostream &out, const ReferenceCell &reference_cell);

  friend std::istream &
  operator>>(std::istream &in, ReferenceCell &reference_cell);
};


/**
 * Output operator that writes the @p reference_cell object to the stream
 * in a text format in which the object is represented by an integer. The
 * details of which integer value represents each kind of reference cell
 * is unimportant and consequently not specified. If you want a string
 * representation of what a ReferenceCell is, use ReferenceCell::to_string().
 */
std::ostream &
operator<<(std::ostream &out, const ReferenceCell &reference_cell);

/**
 * Input operator that reads the @p reference_cell object from the stream
 * in a text format in which the object is represented by an integer. Which
 * specific integer value represents which reference cell is unspecified,
 * but the function uses the same translation as the corresponding
 * output `operator<<`.
 */
std::istream &
operator>>(std::istream &in, ReferenceCell &reference_cell);


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
  inline constexpr ReferenceCell
  make_reference_cell_from_int(const std::uint8_t kind)
  {
#ifndef DEAL_II_CXX14_CONSTEXPR_BUG
    // Make sure these are the only indices from which objects can be
    // created.
    Assert((kind == std::numeric_limits<std::uint8_t>::max()) || (kind < 8),
           ExcInternalError());
#endif

    // Call the private constructor, which we can from here because this
    // function is a 'friend'.
    return {kind};
  }
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
  constexpr const ReferenceCell Vertex =
    internal::make_reference_cell_from_int(0);
  constexpr const ReferenceCell Line =
    internal::make_reference_cell_from_int(1);
  constexpr const ReferenceCell Triangle =
    internal::make_reference_cell_from_int(2);
  constexpr const ReferenceCell Quadrilateral =
    internal::make_reference_cell_from_int(3);
  constexpr const ReferenceCell Tetrahedron =
    internal::make_reference_cell_from_int(4);
  constexpr const ReferenceCell Pyramid =
    internal::make_reference_cell_from_int(5);
  constexpr const ReferenceCell Wedge =
    internal::make_reference_cell_from_int(6);
  constexpr const ReferenceCell Hexahedron =
    internal::make_reference_cell_from_int(7);
  constexpr const ReferenceCell Invalid =
    internal::make_reference_cell_from_int(
      std::numeric_limits<std::uint8_t>::max());

  /**
   * Return the correct simplex reference cell type for the given dimension
   * `dim`. Depending on the template argument `dim`, this function returns a
   * reference to either Vertex, Line, Triangle, or Tetrahedron.
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



inline constexpr ReferenceCell::ReferenceCell()
  : ReferenceCell(ReferenceCells::Invalid)
{}



template <class Archive>
inline void
ReferenceCell::serialize(Archive &archive, const unsigned int /*version*/)
{
  archive &kind;
}



inline constexpr std::size_t
ReferenceCell::memory_consumption()
{
  return sizeof(ReferenceCells::Invalid);
}



inline ArrayView<const unsigned int>
ReferenceCell::faces_for_given_vertex(const unsigned int vertex) const
{
  AssertIndexRange(vertex, n_vertices());
  switch (this->kind)
    {
      case ReferenceCells::Line:
        return {&GeometryInfo<2>::vertex_to_face[vertex][0], 1};
      case ReferenceCells::Quadrilateral:
        return {&GeometryInfo<2>::vertex_to_face[vertex][0], 2};
      case ReferenceCells::Triangle:
        {
          static constexpr ndarray<unsigned int, 3, 2> table = {
            {{{0, 2}}, {{0, 1}}, {{1, 2}}}};
          return table[vertex];
        }
      case ReferenceCells::Tetrahedron:
        {
          static constexpr ndarray<unsigned int, 4, 3> table = {
            {{{0, 1, 2}}, {{0, 1, 3}}, {{0, 2, 3}}, {{1, 2, 3}}}};

          return table[vertex];
        }
      case ReferenceCells::Pyramid:
        {
          static constexpr unsigned int X = numbers::invalid_unsigned_int;
          static constexpr ndarray<unsigned int, 5, 4> table = {
            {{{0, 1, 3, X}},
             {{0, 2, 3, X}},
             {{0, 1, 4, X}},
             {{0, 2, 4, X}},
             {{1, 2, 3, 4}}}};

          return {&table[vertex][0], vertex == 4 ? 4u : 3u};
        }
      case ReferenceCells::Wedge:
        {
          AssertIndexRange(vertex, 6);
          static constexpr ndarray<unsigned int, 6, 3> table = {{{{0, 2, 4}},
                                                                 {{0, 2, 3}},
                                                                 {{0, 3, 4}},
                                                                 {{1, 2, 4}},
                                                                 {{1, 2, 3}},
                                                                 {{1, 3, 4}}}};

          return table[vertex];
        }
      case ReferenceCells::Hexahedron:
        return {&GeometryInfo<3>::vertex_to_face[vertex][0], 3};
      default:
        Assert(false, ExcNotImplemented());
    }

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
  switch (this->kind)
    {
      case ReferenceCells::Vertex:
        return 0;
      case ReferenceCells::Line:
        return 1;
      case ReferenceCells::Triangle:
      case ReferenceCells::Quadrilateral:
        return 2;
      case ReferenceCells::Tetrahedron:
      case ReferenceCells::Pyramid:
      case ReferenceCells::Wedge:
      case ReferenceCells::Hexahedron:
        return 3;
      default:
        Assert(false, ExcNotImplemented());
    }

  return numbers::invalid_unsigned_int;
}



template <int dim>
Quadrature<dim>
ReferenceCell::get_midpoint_quadrature() const
{
  return Quadrature<dim>(std::vector<Point<dim>>({barycenter<dim>()}),
                         std::vector<double>({volume()}));
}



inline unsigned int
ReferenceCell::n_vertices() const
{
  switch (this->kind)
    {
      case ReferenceCells::Vertex:
        return 1;
      case ReferenceCells::Line:
        return 2;
      case ReferenceCells::Triangle:
        return 3;
      case ReferenceCells::Quadrilateral:
        return 4;
      case ReferenceCells::Tetrahedron:
        return 4;
      case ReferenceCells::Pyramid:
        return 5;
      case ReferenceCells::Wedge:
        return 6;
      case ReferenceCells::Hexahedron:
        return 8;
      default:
        Assert(false, ExcNotImplemented());
    }

  return numbers::invalid_unsigned_int;
}



inline unsigned int
ReferenceCell::n_lines() const
{
  switch (this->kind)
    {
      case ReferenceCells::Vertex:
        return 0;
      case ReferenceCells::Line:
        return 1;
      case ReferenceCells::Triangle:
        return 3;
      case ReferenceCells::Quadrilateral:
        return 4;
      case ReferenceCells::Tetrahedron:
        return 6;
      case ReferenceCells::Pyramid:
        return 8;
      case ReferenceCells::Wedge:
        return 9;
      case ReferenceCells::Hexahedron:
        return 12;
      default:
        Assert(false, ExcNotImplemented());
    }

  return numbers::invalid_unsigned_int;
}



template <int dim>
Point<dim>
ReferenceCell::vertex(const unsigned int v) const
{
  AssertDimension(dim, get_dimension());
  AssertIndexRange(v, n_vertices());

  switch (dim)
    {
      case 0:
        {
          if (*this == ReferenceCells::Vertex)
            return Point<dim>();
          break;
        }
      case 1:
        {
          static const Point<dim> vertices[2] = {
            Point<dim>(),              // the origin
            Point<dim>::unit_vector(0) // unit point along x-axis
          };
          if (*this == ReferenceCells::Line)
            return vertices[v];
          break;
        }
      case 2:
        {
          switch (this->kind)
            {
              case ReferenceCells::Triangle:
                {
                  static const Point<dim> vertices[3] = {
                    Point<dim>(),               // the origin
                    Point<dim>::unit_vector(0), // unit point along x-axis
                    Point<dim>::unit_vector(1)  // unit point along y-axis
                  };
                  return vertices[v];
                }
              case ReferenceCells::Quadrilateral:
                {
                  static const Point<dim> vertices[4] = {
                    // First the two points on the x-axis
                    Point<dim>(),
                    Point<dim>::unit_vector(0),
                    // Then these two points shifted in the y-direction
                    Point<dim>() + Point<dim>::unit_vector(1),
                    Point<dim>::unit_vector(0) + Point<dim>::unit_vector(1)};
                  return vertices[v];
                }
            }
          break;
        }
      case 3:
        {
          switch (this->kind)
            {
              case ReferenceCells::Tetrahedron:
                {
                  static const Point<dim> vertices[4] = {
                    Point<dim>(),               // the origin
                    Point<dim>::unit_vector(0), // unit point along x-axis
                    Point<dim>::unit_vector(1), // unit point along y-axis
                    Point<dim>::unit_vector(2)  // unit point along z-axis
                  };
                  return vertices[v];
                }
              case ReferenceCells::Pyramid:
                {
                  static const Point<dim> vertices[5] = {
                    Point<dim>{-1.0, -1.0, 0.0},
                    Point<dim>{+1.0, -1.0, 0.0},
                    Point<dim>{-1.0, +1.0, 0.0},
                    Point<dim>{+1.0, +1.0, 0.0},
                    Point<dim>{+0.0, +0.0, 1.0}};
                  return vertices[v];
                }
              case ReferenceCells::Wedge:
                {
                  static const Point<dim> vertices[6] = {
                    // First the three points on the triangular base of the
                    // wedge:
                    Point<dim>(),
                    Point<dim>::unit_vector(0),
                    Point<dim>::unit_vector(1),
                    // And now everything shifted in the z-direction again
                    Point<dim>() + Point<dim>::unit_vector(2),
                    Point<dim>::unit_vector(0) + Point<dim>::unit_vector(2),
                    Point<dim>::unit_vector(1) + Point<dim>::unit_vector(2)};
                  return vertices[v];
                }
              case ReferenceCells::Hexahedron:
                {
                  static const Point<dim> vertices[8] = {
                    // First the two points on the x-axis
                    Point<dim>(),
                    Point<dim>::unit_vector(0),
                    // Then these two points shifted in the y-direction
                    Point<dim>() + Point<dim>::unit_vector(1),
                    Point<dim>::unit_vector(0) + Point<dim>::unit_vector(1),
                    // And now all four points shifted in the z-direction
                    Point<dim>() + Point<dim>::unit_vector(2),
                    Point<dim>::unit_vector(0) + Point<dim>::unit_vector(2),
                    Point<dim>() + Point<dim>::unit_vector(1) +
                      Point<dim>::unit_vector(2),
                    Point<dim>::unit_vector(0) + Point<dim>::unit_vector(1) +
                      Point<dim>::unit_vector(2)};
                  return vertices[v];
                }
            }
          break;
        }
      default:
        Assert(false, ExcNotImplemented());
    }

  Assert(false, ExcNotImplemented());
  return Point<dim>();
}


inline unsigned int
ReferenceCell::n_faces() const
{
  switch (this->kind)
    {
      case ReferenceCells::Vertex:
        return 0;
      case ReferenceCells::Line:
        return 2;
      case ReferenceCells::Triangle:
        return 3;
      case ReferenceCells::Quadrilateral:
        return 4;
      case ReferenceCells::Tetrahedron:
        return 4;
      case ReferenceCells::Pyramid:
        return 5;
      case ReferenceCells::Wedge:
        return 5;
      case ReferenceCells::Hexahedron:
        return 6;
      default:
        Assert(false, ExcNotImplemented());
    }

  return numbers::invalid_unsigned_int;
}



inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
ReferenceCell::face_indices() const
{
  return {0U, n_faces()};
}



inline unsigned int
ReferenceCell::n_isotropic_children() const
{
  switch (this->kind)
    {
      case ReferenceCells::Vertex:
        return 0;
      case ReferenceCells::Line:
        return 2;
      case ReferenceCells::Triangle:
        return 4;
      case ReferenceCells::Quadrilateral:
        return 4;
      case ReferenceCells::Tetrahedron:
        return 8;
      case ReferenceCells::Pyramid:
        // We haven't yet decided how to refine pyramids. Update this when we
        // have
        Assert(false, ExcNotImplemented());
        return numbers::invalid_unsigned_int;
      case ReferenceCells::Wedge:
        return 8;
      case ReferenceCells::Hexahedron:
        return 8;
      default:
        Assert(false, ExcNotImplemented());
    }

  return numbers::invalid_unsigned_int;
}



inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
ReferenceCell::isotropic_child_indices() const
{
  return {0U, n_isotropic_children()};
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



inline ReferenceCell
ReferenceCell::face_reference_cell(const unsigned int face_no) const
{
  AssertIndexRange(face_no, n_faces());

  switch (this->kind)
    {
      case ReferenceCells::Vertex:
        return ReferenceCells::Invalid;
      case ReferenceCells::Line:
        return ReferenceCells::Vertex;
      case ReferenceCells::Triangle:
      case ReferenceCells::Quadrilateral:
        return ReferenceCells::Line;
      case ReferenceCells::Tetrahedron:
        return ReferenceCells::Triangle;
      case ReferenceCells::Pyramid:
        if (face_no == 0)
          return ReferenceCells::Quadrilateral;
        else
          return ReferenceCells::Triangle;
      case ReferenceCells::Wedge:
        if (face_no > 1)
          return ReferenceCells::Quadrilateral;
        else
          return ReferenceCells::Triangle;
      case ReferenceCells::Hexahedron:
        return ReferenceCells::Quadrilateral;
      default:
        Assert(false, ExcNotImplemented());
    }

  return ReferenceCells::Invalid;
}



inline constexpr unsigned char
ReferenceCell::default_combined_face_orientation()
{
  // Our convention is that 'orientation' has a default value of true and
  // occupies the least-significant bit while rotate and flip have default
  // values of 'false' and occupy the second and third bits.
  return 0b001;
}



inline constexpr unsigned char
ReferenceCell::reversed_combined_line_orientation()
{
  // For a reversed line 'orientation' is false and neither flip nor rotate are
  // defined.
  return 0b000;
}



inline unsigned int
ReferenceCell::child_cell_on_face(
  const unsigned int  face,
  const unsigned int  subface,
  const unsigned char combined_face_orientation) const
{
  AssertIndexRange(face, n_faces());
  AssertIndexRange(subface, face_reference_cell(face).n_isotropic_children());

  switch (this->kind)
    {
      case ReferenceCells::Vertex:
      case ReferenceCells::Line:
        {
          Assert(false, ExcNotImplemented());
          break;
        }
      case ReferenceCells::Triangle:
        {
          static constexpr ndarray<unsigned int, 3, 2> subcells = {
            {{{0, 1}}, {{1, 2}}, {{2, 0}}}};

          return subcells[face][subface];
        }
      case ReferenceCells::Quadrilateral:
        {
          const auto [face_orientation, face_rotation, face_flip] =
            internal::split_face_orientation(combined_face_orientation);

          return GeometryInfo<2>::child_cell_on_face(
            RefinementCase<2>(RefinementPossibilities<2>::isotropic_refinement),
            face,
            subface,
            face_orientation,
            face_flip,
            face_rotation);
        }
      case ReferenceCells::Tetrahedron:
      case ReferenceCells::Pyramid:
      case ReferenceCells::Wedge:
        {
          Assert(false, ExcNotImplemented());
          break;
        }
      case ReferenceCells::Hexahedron:
        {
          const auto [face_orientation, face_rotation, face_flip] =
            internal::split_face_orientation(combined_face_orientation);

          return GeometryInfo<3>::child_cell_on_face(
            RefinementCase<3>(RefinementPossibilities<3>::isotropic_refinement),
            face,
            subface,
            face_orientation,
            face_flip,
            face_rotation);
        }
      default:
        Assert(false, ExcNotImplemented());
    }

  return numbers::invalid_unsigned_int;
}



inline std::array<unsigned int, 2>
ReferenceCell::standard_vertex_to_face_and_vertex_index(
  const unsigned int vertex) const
{
  AssertIndexRange(vertex, n_vertices());
  // Work around a GCC warning at higher optimization levels by making all of
  // these tables the same size
  constexpr unsigned int X = numbers::invalid_unsigned_int;

  switch (this->kind)
    {
      case ReferenceCells::Vertex:
      case ReferenceCells::Line:
        Assert(false, ExcNotImplemented());
        break;
      case ReferenceCells::Triangle:
        {
          static constexpr ndarray<unsigned int, 6, 2> table = {
            {{{0, 0}}, {{0, 1}}, {{1, 1}}, {{X, X}}, {{X, X}}, {{X, X}}}};

          return table[vertex];
        }
      case ReferenceCells::Quadrilateral:
        {
          return GeometryInfo<2>::standard_quad_vertex_to_line_vertex_index(
            vertex);
        }
      case ReferenceCells::Tetrahedron:
        {
          static constexpr ndarray<unsigned int, 6, 2> table = {
            {{{0, 0}}, {{0, 1}}, {{0, 2}}, {{1, 2}}, {{X, X}}, {{X, X}}}};

          return table[vertex];
        }
      case ReferenceCells::Pyramid:
        {
          static constexpr ndarray<unsigned int, 6, 2> table = {
            {{{0, 0}}, {{0, 1}}, {{0, 2}}, {{0, 3}}, {{1, 2}}, {{X, X}}}};

          return table[vertex];
        }
      case ReferenceCells::Wedge:
        {
          static constexpr ndarray<unsigned int, 6, 2> table = {
            {{{0, 1}}, {{0, 0}}, {{0, 2}}, {{1, 0}}, {{1, 1}}, {{1, 2}}}};

          return table[vertex];
        }
      case ReferenceCells::Hexahedron:
        {
          return GeometryInfo<3>::standard_hex_vertex_to_quad_vertex_index(
            vertex);
        }
      default:
        Assert(false, ExcNotImplemented());
    }

  return {};
}



inline std::array<unsigned int, 2>
ReferenceCell::standard_line_to_face_and_line_index(
  const unsigned int line) const
{
  AssertIndexRange(line, n_lines());

  switch (this->kind)
    {
      case ReferenceCells::Vertex:
      case ReferenceCells::Line:
      case ReferenceCells::Triangle:
      case ReferenceCells::Quadrilateral:
        {
          Assert(false, ExcNotImplemented());
          break;
        }
      case ReferenceCells::Tetrahedron:
        {
          static const std::array<unsigned int, 2> table[6] = {
            {{0, 0}}, {{0, 1}}, {{0, 2}}, {{1, 1}}, {{1, 2}}, {{2, 1}}};

          return table[line];
        }
      case ReferenceCells::Pyramid:
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
      case ReferenceCells::Wedge:
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
      case ReferenceCells::Hexahedron:
        {
          return GeometryInfo<3>::standard_hex_line_to_quad_line_index(line);
        }
      default:
        Assert(false, ExcNotImplemented());
    }

  return {};
}


inline unsigned int
ReferenceCell::line_to_cell_vertices(const unsigned int line,
                                     const unsigned int vertex) const
{
  AssertIndexRange(vertex, 2);
  AssertIndexRange(line, n_lines());

  switch (this->kind)
    {
      case ReferenceCells::Vertex:
      case ReferenceCells::Line:
        return vertex;
      case ReferenceCells::Triangle:
        {
          static constexpr ndarray<unsigned int, 3, 2> table = {
            {{{0, 1}}, {{1, 2}}, {{2, 0}}}};
          return table[line][vertex];
        }
      case ReferenceCells::Quadrilateral:
        {
          static constexpr ndarray<unsigned int, 4, 2> table = {
            {{{0, 2}}, {{1, 3}}, {{0, 1}}, {{2, 3}}}};
          return table[line][vertex];
        }
      case ReferenceCells::Tetrahedron:
        {
          static constexpr ndarray<unsigned int, 6, 2> table = {
            {{{0, 1}}, {{1, 2}}, {{2, 0}}, {{0, 3}}, {{1, 3}}, {{2, 3}}}};
          return table[line][vertex];
        }
      case ReferenceCells::Pyramid:
        {
          static constexpr ndarray<unsigned int, 8, 2> table = {{{{0, 2}},
                                                                 {{1, 3}},
                                                                 {{0, 1}},
                                                                 {{2, 3}},
                                                                 {{4, 0}},
                                                                 {{1, 4}},
                                                                 {{2, 4}},
                                                                 {{4, 3}}}};
          return table[line][vertex];
        }
      case ReferenceCells::Wedge:
        {
          static constexpr ndarray<unsigned int, 9, 2> table = {{{{1, 0}},
                                                                 {{2, 1}},
                                                                 {{0, 2}},
                                                                 {{3, 4}},
                                                                 {{4, 5}},
                                                                 {{5, 3}},
                                                                 {{0, 3}},
                                                                 {{1, 4}},
                                                                 {{2, 5}}}};
          return table[line][vertex];
        }
      case ReferenceCells::Hexahedron:
        {
          // first four lines comprise the bottom face, next four are the top,
          // and the last four are 'bottom to top'
          static constexpr ndarray<unsigned int, 12, 2> table = {{{{0, 2}},
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
          return table[line][vertex];
        }

      default:
        Assert(false, ExcNotImplemented());
    }

  return numbers::invalid_unsigned_int;
}


inline unsigned int
ReferenceCell::face_to_cell_lines(
  const unsigned int  face,
  const unsigned int  line,
  const unsigned char combined_face_orientation) const
{
  AssertIndexRange(face, n_faces());
  AssertIndexRange(line, face_reference_cell(face).n_lines());

  static constexpr unsigned int X = numbers::invalid_unsigned_int;

  switch (this->kind)
    {
      case ReferenceCells::Vertex:
        {
          Assert(false, ExcNotImplemented());
          break;
        }
      case ReferenceCells::Line:
        {
          const auto [face_orientation, face_rotation, face_flip] =
            internal::split_face_orientation(combined_face_orientation);

          return GeometryInfo<1>::face_to_cell_lines(
            face, line, face_orientation, face_flip, face_rotation);
        }
      case ReferenceCells::Triangle:
        {
          return face;
        }
      case ReferenceCells::Quadrilateral:
        {
          const auto [face_orientation, face_rotation, face_flip] =
            internal::split_face_orientation(combined_face_orientation);

          return GeometryInfo<2>::face_to_cell_lines(
            face, line, face_orientation, face_flip, face_rotation);
        }
      case ReferenceCells::Tetrahedron:
        {
          static constexpr ndarray<unsigned int, 4, 3> table = {
            {{{0, 1, 2}}, {{0, 3, 4}}, {{2, 5, 3}}, {{1, 4, 5}}}};

          return table[face][standard_to_real_face_line(
            line, face, combined_face_orientation)];
        }
      case ReferenceCells::Pyramid:
        {
          static constexpr ndarray<unsigned int, 5, 4> table = {
            {{{0, 1, 2, 3}},
             {{0, 6, 4, X}},
             {{1, 5, 7, X}},
             {{2, 4, 5, X}},
             {{3, 7, 6, X}}}};

          return table[face][standard_to_real_face_line(
            line, face, combined_face_orientation)];
        }
      case ReferenceCells::Wedge:
        {
          static constexpr ndarray<unsigned int, 5, 4> table = {
            {{{0, 2, 1, X}},
             {{3, 4, 5, X}},
             {{6, 7, 0, 3}},
             {{7, 8, 1, 4}},
             {{8, 6, 5, 2}}}};

          return table[face][standard_to_real_face_line(
            line, face, combined_face_orientation)];
        }
      case ReferenceCells::Hexahedron:
        {
          const auto [face_orientation, face_rotation, face_flip] =
            internal::split_face_orientation(combined_face_orientation);

          return GeometryInfo<3>::face_to_cell_lines(
            face, line, face_orientation, face_flip, face_rotation);
        }
      default:
        Assert(false, ExcNotImplemented());
    }

  return numbers::invalid_unsigned_int;
}



inline unsigned int
ReferenceCell::face_to_cell_vertices(
  const unsigned int  face,
  const unsigned int  vertex,
  const unsigned char combined_face_orientation) const
{
  AssertIndexRange(face, n_faces());
  AssertIndexRange(vertex, face_reference_cell(face).n_vertices());

  switch (this->kind)
    {
      case ReferenceCells::Vertex:
        {
          Assert(false, ExcNotImplemented());
          break;
        }
      case ReferenceCells::Line:
        {
          const auto [face_orientation, face_rotation, face_flip] =
            internal::split_face_orientation(combined_face_orientation);

          return GeometryInfo<1>::face_to_cell_vertices(
            face, vertex, face_orientation, face_flip, face_rotation);
        }
      case ReferenceCells::Triangle:
        {
          static constexpr ndarray<unsigned int, 3, 2> table = {
            {{{0, 1}}, {{1, 2}}, {{2, 0}}}};

          return table[face][combined_face_orientation !=
                                 reversed_combined_line_orientation() ?
                               vertex :
                               (1 - vertex)];
        }
      case ReferenceCells::Quadrilateral:
        {
          const auto [face_orientation, face_rotation, face_flip] =
            internal::split_face_orientation(combined_face_orientation);

          return GeometryInfo<2>::face_to_cell_vertices(
            face, vertex, face_orientation, face_flip, face_rotation);
        }
      case ReferenceCells::Tetrahedron:
        {
          static constexpr ndarray<unsigned int, 4, 3> table = {
            {{{0, 1, 2}}, {{1, 0, 3}}, {{0, 2, 3}}, {{2, 1, 3}}}};

          return table[face][standard_to_real_face_vertex(
            vertex, face, combined_face_orientation)];
        }
      case ReferenceCells::Pyramid:
        {
          constexpr auto X = numbers::invalid_unsigned_int;
          static constexpr ndarray<unsigned int, 5, 4> table = {
            {{{0, 1, 2, 3}},
             {{0, 2, 4, X}},
             {{3, 1, 4, X}},
             {{1, 0, 4, X}},
             {{2, 3, 4, X}}}};

          return table[face][standard_to_real_face_vertex(
            vertex, face, combined_face_orientation)];
        }
      case ReferenceCells::Wedge:
        {
          constexpr auto X = numbers::invalid_unsigned_int;
          static constexpr ndarray<unsigned int, 6, 4> table = {
            {{{1, 0, 2, X}},
             {{3, 4, 5, X}},
             {{0, 1, 3, 4}},
             {{1, 2, 4, 5}},
             {{2, 0, 5, 3}}}};

          return table[face][standard_to_real_face_vertex(
            vertex, face, combined_face_orientation)];
        }
      case ReferenceCells::Hexahedron:
        {
          const auto [face_orientation, face_rotation, face_flip] =
            internal::split_face_orientation(combined_face_orientation);

          return GeometryInfo<3>::face_to_cell_vertices(
            face, vertex, face_orientation, face_flip, face_rotation);
        }
      default:
        Assert(false, ExcNotImplemented());
    }

  return numbers::invalid_unsigned_int;
}



template <int dim>
Point<dim>
ReferenceCell::face_vertex_location(const unsigned int face,
                                    const unsigned int vertex) const
{
  return this->template vertex<dim>(
    face_to_cell_vertices(face, vertex, default_combined_face_orientation()));
}



inline unsigned int
ReferenceCell::standard_to_real_face_vertex(
  const unsigned int  vertex,
  const unsigned int  face,
  const unsigned char face_orientation) const
{
  AssertIndexRange(face, n_faces());
  AssertIndexRange(vertex, face_reference_cell(face).n_vertices());

  switch (this->kind)
    {
      case ReferenceCells::Vertex:
      case ReferenceCells::Line:
        Assert(false, ExcNotImplemented());
        break;
      case ReferenceCells::Triangle:
      case ReferenceCells::Quadrilateral:
        return line_vertex_permutations[face_orientation][vertex];
      case ReferenceCells::Tetrahedron:
        return triangle_vertex_permutations[face_orientation][vertex];
      case ReferenceCells::Pyramid:
        // face 0 is a quadrilateral
        if (face == 0)
          return quadrilateral_vertex_permutations[face_orientation][vertex];
        else
          return triangle_vertex_permutations[face_orientation][vertex];
      case ReferenceCells::Wedge:
        // faces 0 and 1 are triangles
        if (face > 1)
          return quadrilateral_vertex_permutations[face_orientation][vertex];
        else
          return triangle_vertex_permutations[face_orientation][vertex];
      case ReferenceCells::Hexahedron:
        return quadrilateral_vertex_permutations[face_orientation][vertex];
      default:
        Assert(false, ExcNotImplemented());
    }

  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}



inline unsigned int
ReferenceCell::standard_to_real_face_line(
  const unsigned int  line,
  const unsigned int  face,
  const unsigned char combined_face_orientation) const
{
  AssertIndexRange(face, n_faces());
  AssertIndexRange(line, face_reference_cell(face).n_lines());

  static constexpr ndarray<unsigned int, 6, 3> triangle_table = {{{{2, 1, 0}},
                                                                  {{0, 1, 2}},
                                                                  {{1, 0, 2}},
                                                                  {{2, 0, 1}},
                                                                  {{0, 2, 1}},
                                                                  {{1, 2, 0}}}};


  switch (this->kind)
    {
      case ReferenceCells::Vertex:
      case ReferenceCells::Line:
      case ReferenceCells::Triangle:
      case ReferenceCells::Quadrilateral:
        Assert(false, ExcNotImplemented());
        break;
      case ReferenceCells::Tetrahedron:
        return triangle_table[combined_face_orientation][line];
      case ReferenceCells::Pyramid:
        if (face == 0) // The quadrilateral face
          {
            const auto [face_orientation, face_rotation, face_flip] =
              internal::split_face_orientation(combined_face_orientation);

            return GeometryInfo<3>::standard_to_real_face_line(line,
                                                               face_orientation,
                                                               face_flip,
                                                               face_rotation);
          }
        else // One of the triangular faces
          {
            return triangle_table[combined_face_orientation][line];
          }
      case ReferenceCells::Wedge:
        if (face > 1) // One of the quadrilateral faces
          {
            const auto [face_orientation, face_rotation, face_flip] =
              internal::split_face_orientation(combined_face_orientation);

            return GeometryInfo<3>::standard_to_real_face_line(line,
                                                               face_orientation,
                                                               face_flip,
                                                               face_rotation);
          }
        else // One of the triangular faces
          return triangle_table[combined_face_orientation][line];
      case ReferenceCells::Hexahedron:
        {
          static constexpr ndarray<unsigned int, 8, 4> table = {
            {{{2, 3, 0, 1}},
             {{0, 1, 2, 3}},
             {{0, 1, 3, 2}},
             {{3, 2, 0, 1}},
             {{3, 2, 1, 0}},
             {{1, 0, 3, 2}},
             {{1, 0, 2, 3}},
             {{2, 3, 1, 0}}}};
          return table[combined_face_orientation][line];
        }
      default:
        Assert(false, ExcNotImplemented());
    }

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
      }
    return ReferenceCells::Invalid;
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
      }
    return ReferenceCells::Invalid;
  }
} // namespace ReferenceCells


inline ReferenceCell
ReferenceCell::n_vertices_to_type(const int dim, const unsigned int n_vertices)
{
  AssertIndexRange(dim, 4);
  AssertIndexRange(n_vertices, 9);

  const auto                                X     = ReferenceCells::Invalid;
  static const ndarray<ReferenceCell, 4, 9> table = {
    {// dim 0
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
ReferenceCell::d_linear_shape_function(const Point<dim>  &xi,
                                       const unsigned int i) const
{
  AssertDimension(dim, get_dimension());
  AssertIndexRange(i, n_vertices());
  switch (this->kind)
    {
      case ReferenceCells::Vertex:
      case ReferenceCells::Line:
      case ReferenceCells::Quadrilateral:
      case ReferenceCells::Hexahedron:
        return GeometryInfo<dim>::d_linear_shape_function(xi, i);
      // see also BarycentricPolynomials<2>::compute_value
      case ReferenceCells::Triangle:
        {
          switch (i)
            {
              case 0:
                return 1.0 - xi[std::min(0, dim - 1)] -
                       xi[std::min(1, dim - 1)];
              case 1:
                return xi[std::min(0, dim - 1)];
              case 2:
                return xi[std::min(1, dim - 1)];
              default:
                Assert(false, ExcInternalError());
            }
        }
      // see also BarycentricPolynomials<3>::compute_value
      case ReferenceCells::Tetrahedron:
        {
          switch (i)
            {
              case 0:
                return 1.0 - xi[std::min(0, dim - 1)] -
                       xi[std::min(1, dim - 1)] - xi[std::min(2, dim - 1)];
              case 1:
                return xi[std::min(0, dim - 1)];
              case 2:
                return xi[std::min(1, dim - 1)];
              case 3:
                return xi[std::min(2, dim - 1)];
              default:
                Assert(false, ExcInternalError());
            }
        }
      // see also ScalarLagrangePolynomialPyramid::compute_value()
      case ReferenceCells::Pyramid:
        {
          const double Q14 = 0.25;

          const double r = xi[std::min(0, dim - 1)];
          const double s = xi[std::min(1, dim - 1)];
          const double t = xi[std::min(2, dim - 1)];

          const double ratio =
            (std::fabs(t - 1.0) > 1.0e-14 ? (r * s * t) / (1.0 - t) : 0.0);

          if (i == 0)
            return Q14 * ((1.0 - r) * (1.0 - s) - t + ratio);
          if (i == 1)
            return Q14 * ((1.0 + r) * (1.0 - s) - t - ratio);
          if (i == 2)
            return Q14 * ((1.0 - r) * (1.0 + s) - t - ratio);
          if (i == 3)
            return Q14 * ((1.0 + r) * (1.0 + s) - t + ratio);
          else
            return t;
        }
      // see also ScalarLagrangePolynomialWedge::compute_value()
      case ReferenceCells::Wedge:
        return ReferenceCell(ReferenceCells::Triangle)
                 .d_linear_shape_function<2>(Point<2>(xi[std::min(0, dim - 1)],
                                                      xi[std::min(1, dim - 1)]),
                                             i % 3) *
               ReferenceCell(ReferenceCells::Line)
                 .d_linear_shape_function<1>(Point<1>(xi[std::min(2, dim - 1)]),
                                             i / 3);
      default:
        Assert(false, ExcNotImplemented());
    }

  return 0.0;
}



template <int dim>
inline Tensor<1, dim>
ReferenceCell::d_linear_shape_function_gradient(const Point<dim>  &xi,
                                                const unsigned int i) const
{
  AssertDimension(dim, get_dimension());
  switch (this->kind)
    {
      case ReferenceCells::Vertex:
      case ReferenceCells::Line:
      case ReferenceCells::Quadrilateral:
      case ReferenceCells::Hexahedron:
        return GeometryInfo<dim>::d_linear_shape_function_gradient(xi, i);
        // see also BarycentricPolynomials<2>::compute_grad()
      case ReferenceCells::Triangle:
        switch (i)
          {
            case 0:
              return Point<dim>(-1.0, -1.0);
            case 1:
              return Point<dim>(+1.0, +0.0);
            case 2:
              return Point<dim>(+0.0, +1.0);
            default:
              Assert(false, ExcInternalError());
          }
      default:
        Assert(false, ExcNotImplemented());
    }

  return Point<dim>(+0.0, +0.0, +0.0);
}



inline double
ReferenceCell::volume() const
{
  switch (this->kind)
    {
      case ReferenceCells::Vertex:
        return 0;
      case ReferenceCells::Line:
        return 1;
      case ReferenceCells::Triangle:
        return 1. / 2.;
      case ReferenceCells::Quadrilateral:
        return 1;
      case ReferenceCells::Tetrahedron:
        return 1. / 6.;
      case ReferenceCells::Pyramid:
        return 4. / 3.;
      case ReferenceCells::Wedge:
        return 1. / 2.;
      case ReferenceCells::Hexahedron:
        return 1;
      default:
        Assert(false, ExcNotImplemented());
    }

  return 0.0;
}



template <int dim>
inline Point<dim>
ReferenceCell::barycenter() const
{
  AssertDimension(dim, get_dimension());

  switch (this->kind)
    {
      case ReferenceCells::Vertex:
        return Point<dim>();
      case ReferenceCells::Line:
        return Point<dim>(1. / 2.);
      case ReferenceCells::Triangle:
        return Point<dim>(1. / 3., 1. / 3.);
      case ReferenceCells::Quadrilateral:
        return Point<dim>(1. / 2., 1. / 2.);
      case ReferenceCells::Tetrahedron:
        return Point<dim>(1. / 4., 1. / 4., 1. / 4.);
      case ReferenceCells::Pyramid:
        return Point<dim>(0, 0, 1. / 4.);
      case ReferenceCells::Wedge:
        return Point<dim>(1. / 3, 1. / 3, 1. / 2.);
      case ReferenceCells::Hexahedron:
        return Point<dim>(1. / 2., 1. / 2., 1. / 2.);
      default:
        Assert(false, ExcNotImplemented());
    }

  return Point<dim>();
}



template <int dim>
inline bool
ReferenceCell::contains_point(const Point<dim> &p, const double tolerance) const
{
  AssertDimension(dim, get_dimension());

  // Introduce abbreviations to silence compiler warnings about invalid
  // array accesses (that can't happen, of course, but the compiler
  // doesn't know that).
  constexpr unsigned int x_coordinate = 0;
  constexpr unsigned int y_coordinate = (dim >= 2 ? 1 : 0);
  constexpr unsigned int z_coordinate = (dim >= 3 ? 2 : 0);

  switch (this->kind)
    {
      case ReferenceCells::Vertex:
        {
          // Vertices are special cases in that they do not actually
          // have coordinates. Error out if this function is called
          // with a vertex:
          Assert(false,
                 ExcMessage("Vertices are zero-dimensional objects and "
                            "as a consequence have no coordinates. You "
                            "cannot meaningfully ask whether a point is "
                            "inside a vertex (within a certain tolerance) "
                            "without coordinate values."));
          return false;
        }
      case ReferenceCells::Line:
      case ReferenceCells::Quadrilateral:
      case ReferenceCells::Hexahedron:
        {
          for (unsigned int d = 0; d < dim; ++d)
            if ((p[d] < -tolerance) || (p[d] > 1 + tolerance))
              return false;
          return true;
        }
      case ReferenceCells::Triangle:
      case ReferenceCells::Tetrahedron:
        {
          // First make sure that we are in the first quadrant or octant
          for (unsigned int d = 0; d < dim; ++d)
            if (p[d] < -tolerance)
              return false;

          // Now we also need to make sure that we are below the diagonal line
          // or plane that delineates the simplex. This diagonal is given by
          // sum(p[d])<=1, and a diagonal a distance eps away is given by
          // sum(p[d])<=1+eps*sqrt(d). (For example, the point at (1,1) is a
          // distance of 1/sqrt(2) away from the diagonal. That is, its
          // sum satisfies
          //   sum(p[d]) = 2 <= 1 + (1/sqrt(2)) * sqrt(2)
          // in other words, it satisfies the predicate with eps=1/sqrt(2).)
          double sum = 0;
          for (unsigned int d = 0; d < dim; ++d)
            sum += p[d];
          return (sum <= 1 + tolerance * std::sqrt(1. * dim));
        }
      case ReferenceCells::Pyramid:
        {
          // A pyramid only lives in the upper half-space:
          if (p[z_coordinate] < -tolerance)
            return false;

          // It also only lives in the space below z=1:
          if (p[z_coordinate] > 1 + tolerance)
            return false;

          // Within what's left of the space, a pyramid is a cone that tapers
          // towards the top. First compute the distance of the point to the
          // axis in the max norm (this is the right norm because the vertices
          // of the pyramid are at points +/-1, +/-1):
          const double distance_from_axis =
            std::max(std::fabs(p[x_coordinate]), std::fabs(p[y_coordinate]));

          // We are inside the pyramid if the distance from the axis is less
          // than (1-z)
          return (distance_from_axis <= 1 + tolerance - p[z_coordinate]);
        }
      case ReferenceCells::Wedge:
        {
          // The wedge we use is a triangle extruded into the third
          // dimension by one unit. So we can use the same logic as for
          // triangles above (i.e., for the simplex above, using dim==2)
          // and then check the third dimension separately.

          if ((p[x_coordinate] < -tolerance) || (p[y_coordinate] < -tolerance))
            return false;

          const double sum = p[x_coordinate] + p[y_coordinate];
          if (sum > 1 + tolerance * std::sqrt(2.0))
            return false;

          if (p[z_coordinate] < -tolerance)
            return false;
          if (p[z_coordinate] > 1 + tolerance)
            return false;

          return true;
        }
      default:
        Assert(false, ExcNotImplemented());
    }

  return false;
}



template <int dim>
inline Tensor<1, dim>
ReferenceCell::unit_tangential_vectors(const unsigned int face_no,
                                       const unsigned int i) const
{
  AssertDimension(dim, get_dimension());
  AssertIndexRange(i, dim - 1);

  switch (this->kind)
    {
      case ReferenceCells::Vertex:
      case ReferenceCells::Line:
      case ReferenceCells::Quadrilateral:
      case ReferenceCells::Hexahedron:
        AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);
        return GeometryInfo<dim>::unit_tangential_vectors[face_no][i];
      case ReferenceCells::Triangle:
        {
          AssertIndexRange(face_no, 3);
          static const std::array<Tensor<1, dim>, 3> table = {
            {Point<dim>(1, 0),
             Point<dim>(-std::sqrt(0.5), +std::sqrt(0.5)),
             Point<dim>(0, -1)}};

          return table[face_no];
        }
      case ReferenceCells::Tetrahedron:
        {
          AssertIndexRange(face_no, 4);
          static const ndarray<Tensor<1, dim>, 4, 2> table = {
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
      case ReferenceCells::Pyramid:
        {
          AssertIndexRange(face_no, 5);
          static const ndarray<Tensor<1, dim>, 5, 2> table = {
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
      case ReferenceCells::Wedge:
        {
          AssertIndexRange(face_no, 5);
          static const ndarray<Tensor<1, dim>, 5, 2> table = {
            {{{Point<dim>(0, 1, 0), Point<dim>(1, 0, 0)}},
             {{Point<dim>(1, 0, 0), Point<dim>(0, 1, 0)}},
             {{Point<dim>(1, 0, 0), Point<dim>(0, 0, 1)}},
             {{Point<dim>(-1 / std::sqrt(2.0), +1 / std::sqrt(2.0), 0),
               Point<dim>(0, 0, 1)}},
             {{Point<dim>(0, 0, 1), Point<dim>(0, 1, 0)}}}};

          return table[face_no][i];
        }
      default:
        Assert(false, ExcNotImplemented());
    }


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



inline unsigned int
ReferenceCell::n_face_orientations(const unsigned int face_no) const
{
  AssertIndexRange(face_no, n_faces());
  if (get_dimension() == 1)
    return 1;
  if (get_dimension() == 2)
    return 2;
  else if (face_reference_cell(face_no) == ReferenceCells::Quadrilateral)
    return 8;
  else if (face_reference_cell(face_no) == ReferenceCells::Triangle)
    return 6;

  Assert(false, ExcInternalError());
  return numbers::invalid_unsigned_int;
}



inline bool
ReferenceCell::standard_vs_true_line_orientation(
  const unsigned int  line,
  const unsigned int  face,
  const unsigned char combined_face_orientation,
  const bool          line_orientation) const
{
  if (*this == ReferenceCells::Hexahedron)
    {
      static constexpr dealii::ndarray<bool, 2, 8> bool_table{
        {{{true, true, false, true, false, false, true, false}},
         {{true, true, true, false, false, false, false, true}}}};

      return (line_orientation ==
              bool_table[line / 2][combined_face_orientation]);
    }
  else if (*this == ReferenceCells::Tetrahedron)
    {
      static constexpr unsigned int X = numbers::invalid_unsigned_int;
      static constexpr dealii::ndarray<unsigned int, 4, 3> combined_lines{
        {{{0, 0, 0}}, {{X, 0, 1}}, {{X, 0, X}}, {{X, X, X}}}};

      const auto combined_line = combined_lines[face][line];

      Assert(combined_line != X,
             ExcMessage(
               "This function can only be called for following face-line "
               "combinations: (0,0), (0,1), (0,2), (1,1), (1,2), (2,1),"));

      static constexpr dealii::ndarray<bool, 2, 6> bool_table{
        {{{false, true, false, true, false, true}},
         {{true, false, true, false, true, false}}}};

      return (line_orientation ==
              bool_table[combined_line][combined_face_orientation]);
    }
  else
    // TODO: This might actually be wrong for some of the other
    // kinds of objects. We should check this
    return true;
}



namespace internal
{
  template <typename T>
  class NoPermutation : public dealii::ExceptionBase
  {
  public:
    /**
     * Constructor.
     */
    NoPermutation(const ReferenceCell      &entity_type,
                  const ArrayView<const T> &vertices_0,
                  const ArrayView<const T> &vertices_1)
      : entity_type(entity_type)
      , vertices_0(vertices_0)
      , vertices_1(vertices_1)
    {
      Assert(vertices_0.size() >= entity_type.n_vertices(), ExcInternalError());
      Assert(vertices_1.size() >= entity_type.n_vertices(), ExcInternalError());
    }

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
      out << '[';

      const unsigned int n_vertices = entity_type.n_vertices();

      for (unsigned int i = 0; i < n_vertices; ++i)
        {
          out << vertices_0[i];
          if (i + 1 != n_vertices)
            out << ',';
        }

      out << "] is not a valid permutation of [";

      for (unsigned int i = 0; i < n_vertices; ++i)
        {
          out << vertices_1[i];
          if (i + 1 != n_vertices)
            out << ',';
        }

      out << "]." << std::endl;
    }

    /**
     * Entity type.
     */
    const ReferenceCell entity_type;

    /**
     * First set of values.
     */
    const ArrayView<const T> vertices_0;

    /**
     * Second set of values.
     */
    const ArrayView<const T> vertices_1;
  };

  /**
   * This exception is raised whenever the types of two reference cell objects
   * were assumed to be equal, but were not.
   *
   * Parameters to the constructor are the first and second reference cells,
   * both of type <tt>ReferenceCell</tt>.
   */
  DeclException2(
    ExcNonMatchingReferenceCellTypes,
    ReferenceCell,
    ReferenceCell,
    << "The reference-cell type used on this cell (" << arg1.to_string()
    << ") does not match the reference-cell type of the finite element "
    << "associated with this cell (" << arg2.to_string() << "). "
    << "Did you accidentally use simplex elements on hypercube meshes "
    << "(or the other way around), or are you using a mixed mesh and "
    << "assigned a simplex element to a hypercube cell (or the other "
    << "way around) via the active_fe_index?");
} // namespace internal



template <typename T, std::size_t N>
inline unsigned char
ReferenceCell::compute_orientation(const std::array<T, N> &vertices_0,
                                   const std::array<T, N> &vertices_1) const
{
  Assert(N >= n_vertices(),
         ExcMessage("The number of array elements must be equal to or "
                    "greater than the number of vertices of the cell "
                    "referenced by this object."));

  // Call the non-deprecated function, taking care of calling it only with
  // those array elements that we actually care about (see the note
  // in the documentation about the arguments potentially being
  // larger arrays than necessary).
  return get_combined_orientation(
    make_array_view(vertices_0.begin(), vertices_0.begin() + n_vertices()),
    make_array_view(vertices_1.begin(), vertices_1.begin() + n_vertices()));
}



template <typename T>
unsigned char
ReferenceCell::get_combined_orientation(
  const ArrayView<const T> &vertices_0,
  const ArrayView<const T> &vertices_1) const
{
  Assert(vertices_0.size() == n_vertices(),
         ExcMessage("The number of array elements must be equal to "
                    "the number of vertices of the cell "
                    "referenced by this object."));
  Assert(vertices_1.size() == n_vertices(),
         ExcMessage("The number of array elements must be equal to "
                    "the number of vertices of the cell "
                    "referenced by this object."));

  const auto v0_equals = [&](const std::initializer_list<const T> &list) {
    Assert(list.size() == n_vertices(), ExcInternalError());
    return std::equal(vertices_0.begin(), vertices_0.end(), std::begin(list));
  };

  switch (this->kind)
    {
      case ReferenceCells::Line:
        // line_orientation=true
        if (v0_equals({vertices_1[0], vertices_1[1]}))
          return 1;

        // line_orientation=false
        if (v0_equals({vertices_1[1], vertices_1[0]}))
          return 0;
        break;
      case ReferenceCells::Triangle:
        // face_orientation=true, face_rotation=false, face_flip=false
        if (v0_equals({vertices_1[0], vertices_1[1], vertices_1[2]}))
          return 1;

        // face_orientation=true, face_rotation=true, face_flip=false
        if (v0_equals({vertices_1[1], vertices_1[2], vertices_1[0]}))
          return 5;

        // face_orientation=true, face_rotation=false, face_flip=true
        if (v0_equals({vertices_1[2], vertices_1[0], vertices_1[1]}))
          return 3;

        // face_orientation=false, face_rotation=false, face_flip=false
        if (v0_equals({vertices_1[0], vertices_1[2], vertices_1[1]}))
          return 0;

        // face_orientation=false, face_rotation=true, face_flip=false
        if (v0_equals({vertices_1[2], vertices_1[1], vertices_1[0]}))
          return 2;

        // face_orientation=false, face_rotation=false, face_flip=true
        if (v0_equals({vertices_1[1], vertices_1[0], vertices_1[2]}))
          return 4;
        break;
      case ReferenceCells::Quadrilateral:
        // face_orientation=true, face_rotation=false, face_flip=false
        if (v0_equals(
              {vertices_1[0], vertices_1[1], vertices_1[2], vertices_1[3]}))
          return 1;

        // face_orientation=true, face_rotation=true, face_flip=false
        if (v0_equals(
              {vertices_1[2], vertices_1[0], vertices_1[3], vertices_1[1]}))
          return 3;

        // face_orientation=true, face_rotation=false, face_flip=true
        if (v0_equals(
              {vertices_1[3], vertices_1[2], vertices_1[1], vertices_1[0]}))
          return 5;

        // face_orientation=true, face_rotation=true, face_flip=true
        if (v0_equals(
              {vertices_1[1], vertices_1[3], vertices_1[0], vertices_1[2]}))
          return 7;

        // face_orientation=false, face_rotation=false, face_flip=false
        if (v0_equals(
              {vertices_1[0], vertices_1[2], vertices_1[1], vertices_1[3]}))
          return 0;

        // face_orientation=false, face_rotation=true, face_flip=false
        if (v0_equals(
              {vertices_1[2], vertices_1[3], vertices_1[0], vertices_1[1]}))
          return 2;

        // face_orientation=false, face_rotation=false, face_flip=true
        if (v0_equals(
              {vertices_1[3], vertices_1[1], vertices_1[2], vertices_1[0]}))
          return 4;

        // face_orientation=false, face_rotation=true, face_flip=true
        if (v0_equals(
              {vertices_1[1], vertices_1[0], vertices_1[3], vertices_1[2]}))
          return 6;
        break;
      default:
        Assert(false, ExcNotImplemented());
    }

  Assert(false, (internal::NoPermutation<T>(*this, vertices_0, vertices_1)));
  return std::numeric_limits<unsigned char>::max();
}



template <typename T, std::size_t N>
inline std::array<T, N>
ReferenceCell::permute_according_orientation(
  const std::array<T, N> &vertices,
  const unsigned int      orientation) const
{
  Assert(N >= n_vertices(),
         ExcMessage("The number of array elements must be equal to or "
                    "greater than the number of vertices of the cell "
                    "referenced by this object."));

  // Call the non-deprecated function, taking care of calling it only with
  // those array elements that we actually care about (see the note
  // in the documentation about the arguments potentially being
  // larger arrays than necessary).
  const auto permutation = permute_by_combined_orientation(
    make_array_view(vertices.begin(), vertices.begin() + n_vertices()),
    orientation);

  std::array<T, N> temp;
  std::copy(permutation.begin(), permutation.end(), temp.begin());

  return temp;
}



template <typename T>
boost::container::small_vector<T, 8>
ReferenceCell::permute_by_combined_orientation(
  const ArrayView<const T> &vertices,
  const unsigned char       orientation) const
{
  Assert(vertices.size() == n_vertices(),
         ExcMessage("The number of array elements must be equal to "
                    "the number of vertices of the cell "
                    "referenced by this object."));

  switch (this->kind)
    {
      case ReferenceCells::Line:
        switch (orientation)
          {
            case 1:
              return {vertices[0], vertices[1]};
            case 0:
              return {vertices[1], vertices[0]};
            default:
              Assert(false, ExcNotImplemented());
          }
        break;
      case ReferenceCells::Triangle:
        switch (orientation)
          {
            case 1:
              return {vertices[0], vertices[1], vertices[2]};
            case 3:
              return {vertices[1], vertices[2], vertices[0]};
            case 5:
              return {vertices[2], vertices[0], vertices[1]};
            case 0:
              return {vertices[0], vertices[2], vertices[1]};
            case 2:
              return {vertices[2], vertices[1], vertices[0]};
            case 4:
              return {vertices[1], vertices[0], vertices[2]};
            default:
              Assert(false, ExcNotImplemented());
          }
        break;
      case ReferenceCells::Quadrilateral:
        switch (orientation)
          {
            case 1:
              return {vertices[0], vertices[1], vertices[2], vertices[3]};
            case 3:
              return {vertices[2], vertices[0], vertices[3], vertices[1]};
            case 5:
              return {vertices[3], vertices[2], vertices[1], vertices[0]};
            case 7:
              return {vertices[1], vertices[3], vertices[0], vertices[2]};
            case 0:
              return {vertices[0], vertices[2], vertices[1], vertices[3]};
            case 2:
              return {vertices[2], vertices[3], vertices[0], vertices[1]};
            case 4:
              return {vertices[3], vertices[1], vertices[2], vertices[0]};
            case 6:
              return {vertices[1], vertices[0], vertices[3], vertices[2]};
            default:
              Assert(false, ExcNotImplemented());
          }
        break;
      default:
        AssertThrow(false, ExcNotImplemented());
    }

  return {};
}



template <>
unsigned int
ReferenceCell::vtk_lexicographic_to_node_index<0>(
  const std::array<unsigned, 0> &node_indices,
  const std::array<unsigned, 0> &nodes_per_direction,
  const bool                     legacy_format) const;

template <>
unsigned int
ReferenceCell::vtk_lexicographic_to_node_index<1>(
  const std::array<unsigned, 1> &node_indices,
  const std::array<unsigned, 1> &nodes_per_direction,
  const bool                     legacy_format) const;

template <>
unsigned int
ReferenceCell::vtk_lexicographic_to_node_index<2>(
  const std::array<unsigned, 2> &node_indices,
  const std::array<unsigned, 2> &nodes_per_direction,
  const bool                     legacy_format) const;

template <>
unsigned int
ReferenceCell::vtk_lexicographic_to_node_index<3>(
  const std::array<unsigned, 3> &node_indices,
  const std::array<unsigned, 3> &nodes_per_direction,
  const bool                     legacy_format) const;

DEAL_II_NAMESPACE_CLOSE

#endif

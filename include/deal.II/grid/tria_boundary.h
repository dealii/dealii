// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__tria_boundary_h
#define __deal2__tria_boundary_h


/*----------------------------   boundary-function.h     ---------------------------*/

#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/point.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/manifold.h>

DEAL_II_NAMESPACE_OPEN

template <int dim, int space_dim> class Triangulation;



/**
 * This class is used to represent a boundary to a triangulation. When a
 * triangulation creates a new vertex on the boundary of the domain, it
 * determines the new vertex' coordinates through the following code (here in
 * two dimensions):
 *   @code
 *     ...
 *     Point<2> new_vertex = boundary.get_new_point_on_line (line);
 *     ...
 *   @endcode
 * @p line denotes the line at the boundary that shall be refined and for
 * which we seek the common point of the two child lines.
 *
 * In 3D, a new vertex may be placed on the middle of a line or on the middle
 * of a side. Respectively, the library calls
 *   @code
 *     ...
 *     Point<3> new_line_vertices[4]
 *       = { boundary.get_new_point_on_line (face->line(0)),
 *           boundary.get_new_point_on_line (face->line(1)),
 *           boundary.get_new_point_on_line (face->line(2)),
 *           boundary.get_new_point_on_line (face->line(3))  };
 *     ...
 *   @endcode
 * to get the four midpoints of the lines bounding the quad at the boundary,
 * and after that
 *   @code
 *     ...
 *     Point<3> new_quad_vertex = boundary.get_new_point_on_quad (face);
 *     ...
 *   @endcode
 * to get the midpoint of the face. It is guaranteed that this order (first
 * lines, then faces) holds, so you can use information from the children of
 * the four lines of a face, since these already exist at the time the
 * midpoint of the face is to be computed.
 *
 * Since iterators are passed to the functions, you may use information about
 * boundary indicators and the like, as well as all other information provided
 * by these objects.
 *
 * There are specializations, StraightBoundary, which places the new point
 * right into the middle of the given points, and HyperBallBoundary creating a
 * hyperball with given radius around a given center point.
 *
 * @ingroup boundary
 * @author Wolfgang Bangerth, 1999, 2001, 2009, Ralf Hartmann, 2001, 2008,
 * Luca Heltai, 2014
 */
template <int dim, int spacedim=dim>
class Boundary : public FlatManifold<dim, spacedim>
{
public:

  /**
   * Type keeping information about the normals at the vertices of a face of a
   * cell. Thus, there are <tt>GeometryInfo<dim>::vertices_per_face</tt>
   * normal vectors, that define the tangent spaces of the boundary at the
   * vertices. Note that the vectors stored in this object are not required to
   * be normalized, nor to actually point outward, as one often will only want
   * to check for orthogonality to define the tangent plane; if a function
   * requires the normals to be normalized, then it must do so itself.
   *
   * For obvious reasons, this type is not useful in 1d.
   */
  typedef Tensor<1,spacedim> FaceVertexNormals[GeometryInfo<dim>::vertices_per_face];

  /**
   * Destructor. Does nothing here, but needs to be declared to make it
   * virtual.
   */
  virtual ~Boundary ();


  /**
   * Return intermediate points on a line spaced according to the interior
   * support points of the 1D Gauss-Lobatto quadrature formula.
   *
   * The number of points requested is given by the size of the vector @p
   * points. It is the task of the derived classes to arrange the points in
   * approximately equal distances.
   *
   * This function is called by the @p MappingQ class. This happens on each
   * face line of a cells that has got at least one boundary line.
   *
   * As this function is not needed for @p MappingQ1, it is not made pure
   * virtual, to avoid the need to overload it.  The default implementation
   * throws an error in any case, however.
   */
  virtual
  void
  get_intermediate_points_on_line (const typename Triangulation<dim,spacedim>::line_iterator &line,
                                   std::vector<Point<spacedim> > &points) const;

  /**
   * Return intermediate points on a line spaced according to the tensor
   * product of the interior support points of the 1D Gauss-Lobatto quadrature
   * formula.
   *
   * The number of points requested is given by the size of the vector @p
   * points. It is required that this number is a square of another integer,
   * i.e. <tt>n=points.size()=m*m</tt>. It is the task of the derived classes
   * to arrange the points such they split the quad into <tt>(m+1)(m+1)</tt>
   * approximately equal-sized subquads.
   *
   * This function is called by the <tt>MappingQ<3></tt> class. This happens
   * each face quad of cells in 3d that has got at least one boundary face
   * quad.
   *
   * As this function is not needed for @p MappingQ1, it is not made pure
   * virtual, to avoid the need to overload it.  The default implementation
   * throws an error in any case, however.
   */
  virtual
  void
  get_intermediate_points_on_quad (const typename Triangulation<dim,spacedim>::quad_iterator &quad,
                                   std::vector<Point<spacedim> > &points) const;

  /**
   * Depending on <tt>dim=2</tt> or <tt>dim=3</tt> this function calls the
   * get_intermediate_points_on_line or the get_intermediate_points_on_quad
   * function. It throws an exception for <tt>dim=1</tt>. This wrapper allows
   * dimension independent programming.
   */
  void
  get_intermediate_points_on_face (const typename Triangulation<dim,spacedim>::face_iterator &face,
                                   std::vector<Point<spacedim> > &points) const;

  /**
   * Return the normal vector to the surface at the point p. If p is not in
   * fact on the surface, but only close-by, try to return something
   * reasonable, for example the normal vector at the surface point closest to
   * p.  (The point p will in fact not normally lie on the actual surface, but
   * rather be a quadrature point mapped by some polynomial mapping; the
   * mapped surface, however, will not usually coincide with the actual
   * surface.)
   *
   * The face iterator gives an indication which face this function is
   * supposed to compute the normal vector for.  This is useful if the
   * boundary of the domain is composed of different nondifferential pieces
   * (for example when using the StraightBoundary class to approximate a
   * geometry that is completely described by the coarse mesh, with piecewise
   * (bi-)linear components between the vertices, but where the boundary may
   * have a kink at the vertices itself).
   *
   * @note Implementations of this function should be able to assume that the
   * point p lies within or close to the face described by the first argument.
   * In turn, callers of this function should ensure that this is in fact the
   * case.
   */
  virtual
  Tensor<1,spacedim>
  normal_vector (const typename Triangulation<dim,spacedim>::face_iterator &face,
                 const Point<spacedim> &p) const;

  /**
   * Compute the normal vectors to the boundary at each vertex of the given
   * face. It is not required that the normal vectors be normed somehow.
   * Neither is it required that the normals actually point outward.
   *
   * This function is needed to compute data for C1 mappings. The default
   * implementation is to throw an error, so you need not overload this
   * function in case you do not intend to use C1 mappings.
   *
   * Note that when computing normal vectors at a vertex where the boundary is
   * not differentiable, you have to make sure that you compute the one-sided
   * limits, i.e. limit with respect to points inside the given face.
   */
  virtual
  void
  get_normals_at_vertices (const typename Triangulation<dim,spacedim>::face_iterator &face,
                           FaceVertexNormals &face_vertex_normals) const;

  /**
   * Given a candidate point and a line segment characterized by the iterator,
   * return a point that lies on the surface described by this object. This
   * function is used in some mesh smoothing algorithms that try to move
   * around points in order to improve the mesh quality but need to ensure
   * that points that were on the boundary remain on the boundary.
   *
   * If spacedim==1, then the line represented by the line iterator is the
   * entire space (i.e. it is a cell, not a part of the boundary), and the
   * returned point equals the given input point.
   *
   * Derived classes do not need to implement this function unless mesh
   * smoothing algorithms are used with a particular boundary object. The
   * default implementation of this function throws an exception of type
   * ExcPureFunctionCalled.
   */
  virtual
  Point<spacedim>
  project_to_surface (const typename Triangulation<dim,spacedim>::line_iterator &line,
                      const Point<spacedim> &candidate) const;

  /**
   * Same function as above but for a point that is to be projected onto the
   * area characterized by the given quad.
   *
   * If spacedim<=2, then the surface represented by the quad iterator is the
   * entire space (i.e. it is a cell, not a part of the boundary), and the
   * returned point equals the given input point.
   */
  virtual
  Point<spacedim>
  project_to_surface (const typename Triangulation<dim,spacedim>::quad_iterator &quad,
                      const Point<spacedim> &candidate) const;

  /**
   * Same function as above but for a point that is to be projected onto the
   * area characterized by the given quad.
   *
   * If spacedim<=3, then the manifold represented by the hex iterator is the
   * entire space (i.e. it is a cell, not a part of the boundary), and the
   * returned point equals the given input point.
   */
  virtual
  Point<spacedim>
  project_to_surface (const typename Triangulation<dim,spacedim>::hex_iterator &hex,
                      const Point<spacedim> &candidate) const;

protected:
  /**
   * Returns the support points of the Gauss-Lobatto quadrature formula used
   * for intermediate points.
   *
   * @note Since the boundary description is closely tied to the unit cell
   * support points of MappingQ, new boundary descriptions need to explicitly
   * use these Gauss-Lobatto points and not equidistant points.
   */
  const std::vector<Point<1> > &
  get_line_support_points (const unsigned int n_intermediate_points) const;

private:
  /**
   * Point generator for the intermediate points on a boundary.
   */
  mutable std::vector<std_cxx11::shared_ptr<QGaussLobatto<1> > > points;

  /**
   * Mutex for protecting the points array.
   */
  mutable Threads::Mutex mutex;
};



/**
 * Specialization of Boundary<dim,spacedim>, which places the new point right
 * into the middle of the given points. The middle is defined as the
 * arithmetic mean of the points.
 *
 * This class does not really describe a boundary in the usual sense. By
 * placing new points in the middle of old ones, it rather assumes that the
 * boundary of the domain is given by the polygon/polyhedron defined by the
 * boundary of the initial coarse triangulation.
 *
 * @ingroup boundary
 *
 * @author Wolfgang Bangerth, 1998, 2001, Ralf Hartmann, 2001
 */
template <int dim, int spacedim=dim>
class StraightBoundary : public Boundary<dim,spacedim>
{
public:
  /**
   * Default constructor. Some compilers require this for some reasons.
   */
  StraightBoundary ();

  /**
   * Let the new point be the arithmetic mean of the two vertices of the line.
   *
   * Refer to the general documentation of this class and the documentation of
   * the base class for more information.
   */
  virtual Point<spacedim>
  get_new_point_on_line (const typename Triangulation<dim,spacedim>::line_iterator &line) const;

  /**
   * Let the new point be the arithmetic mean of the four vertices of this
   * quad and the four midpoints of the lines, which are already created at
   * the time of calling this function.
   *
   * Refer to the general documentation of this class and the documentation of
   * the base class for more information.
   */
  virtual
  Point<spacedim>
  get_new_point_on_quad (const typename Triangulation<dim,spacedim>::quad_iterator &quad) const;

  /**
   * Gives <tt>n=points.size()</tt> points that splits the StraightBoundary
   * line into $n+1$ partitions of equal lengths.
   *
   * Refer to the general documentation of this class and the documentation of
   * the base class.
   */
  virtual
  void
  get_intermediate_points_on_line (const typename Triangulation<dim,spacedim>::line_iterator &line,
                                   std::vector<Point<spacedim> > &points) const;

  /**
   * Gives <tt>n=points.size()=m*m</tt> points that splits the
   * StraightBoundary quad into $(m+1)(m+1)$ subquads of equal size.
   *
   * Refer to the general documentation of this class and the documentation of
   * the base class.
   */
  virtual
  void
  get_intermediate_points_on_quad (const typename Triangulation<dim,spacedim>::quad_iterator &quad,
                                   std::vector<Point<spacedim> > &points) const;

  /**
   * Implementation of the function declared in the base class.
   *
   * Refer to the general documentation of this class and the documentation of
   * the base class.
   */
  virtual
  Tensor<1,spacedim>
  normal_vector (const typename Triangulation<dim,spacedim>::face_iterator &face,
                 const Point<spacedim> &p) const;

  /**
   * Compute the normals to the boundary at the vertices of the given face.
   *
   * Refer to the general documentation of this class and the documentation of
   * the base class.
   */
  virtual
  void
  get_normals_at_vertices (const typename Triangulation<dim,spacedim>::face_iterator &face,
                           typename Boundary<dim,spacedim>::FaceVertexNormals &face_vertex_normals) const;

  /**
   * Given a candidate point and a line segment characterized by the iterator,
   * return a point that lies on the surface described by this object. This
   * function is used in some mesh smoothing algorithms that try to move
   * around points in order to improve the mesh quality but need to ensure
   * that points that were on the boundary remain on the boundary.
   *
   * The point returned is the projection of the candidate point onto the line
   * through the two vertices of the given line iterator.
   *
   * If spacedim==1, then the line represented by the line iterator is the
   * entire space (i.e. it is a cell, not a part of the boundary), and the
   * returned point equals the given input point.
   */
  virtual
  Point<spacedim>
  project_to_surface (const typename Triangulation<dim,spacedim>::line_iterator &line,
                      const Point<spacedim> &candidate) const;

  /**
   * Same function as above but for a point that is to be projected onto the
   * area characterized by the given quad.
   *
   * The point returned is the projection of the candidate point onto the
   * bilinear surface spanned by the four vertices of the given quad iterator.
   *
   * If spacedim<=2, then the surface represented by the quad iterator is the
   * entire space (i.e. it is a cell, not a part of the boundary), and the
   * returned point equals the given input point.
   */
  virtual
  Point<spacedim>
  project_to_surface (const typename Triangulation<dim,spacedim>::quad_iterator &quad,
                      const Point<spacedim> &candidate) const;

  /**
   * Same function as above but for a point that is to be projected onto the
   * area characterized by the given quad.
   *
   * The point returned is the projection of the candidate point onto the
   * trilinear manifold spanned by the eight vertices of the given hex
   * iterator.
   *
   * If spacedim<=3, then the manifold represented by the hex iterator is the
   * entire space (i.e. it is a cell, not a part of the boundary), and the
   * returned point equals the given input point.
   */
  virtual
  Point<spacedim>
  project_to_surface (const typename Triangulation<dim,spacedim>::hex_iterator &hex,
                      const Point<spacedim> &candidate) const;
};



/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN

template <>
void
Boundary<1,1>::
get_intermediate_points_on_face (const Triangulation<1,1>::face_iterator &,
                                 std::vector<Point<1> > &) const;

template <>
void
Boundary<1,2>::
get_intermediate_points_on_face (const Triangulation<1,2>::face_iterator &,
                                 std::vector<Point<2> > &) const;

template <>
void
Boundary<1,3>::
get_intermediate_points_on_face (const Triangulation<1,3>::face_iterator &,
                                 std::vector<Point<3> > &) const;
template <>
void
StraightBoundary<1,1>::
get_normals_at_vertices (const Triangulation<1,1>::face_iterator &,
                         Boundary<1,1>::FaceVertexNormals &) const;
template <>
void
StraightBoundary<2,2>::
get_normals_at_vertices (const Triangulation<2,2>::face_iterator &face,
                         Boundary<2,2>::FaceVertexNormals &face_vertex_normals) const;
template <>
void
StraightBoundary<3,3>::
get_normals_at_vertices (const Triangulation<3,3>::face_iterator &face,
                         Boundary<3,3>::FaceVertexNormals &face_vertex_normals) const;

template <>
Point<3>
StraightBoundary<3,3>::
get_new_point_on_quad (const Triangulation<3,3>::quad_iterator &quad) const;

template <>
void
StraightBoundary<3,3>::
get_intermediate_points_on_quad (const Triangulation<3,3>::quad_iterator &quad,
                                 std::vector<Point<3> > &points) const;

template <>
Point<3>
StraightBoundary<1,3>::
project_to_surface (const Triangulation<1, 3>::quad_iterator &quad,
                    const Point<3>  &y) const;


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

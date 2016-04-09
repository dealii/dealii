// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2015 by the deal.II authors
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

#ifndef dealii__tria_manifold_h
#define dealii__tria_manifold_h


/*----------------------------   manifold.h     ---------------------------*/

#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/point.h>
#include <deal.II/base/derivative_form.h>
#include <deal.II/grid/tria.h>

DEAL_II_NAMESPACE_OPEN

/**
 * We collect here some helper functions used in the Manifold<dim,spacedim>
 * classes.
 */
namespace Manifolds
{
  /**
   * Given a general mesh iterator, construct a quadrature object that
   * contains the following points:
   * - If the iterator points to a line, then the quadrature points
   *   are the two vertices of the line. This results in a quadrature
   *   object with two points.
   * - If the iterator points to a quad, then the quadrature points
   *   are the vertices and line mid-points. This results in a quadrature
   *   object with eight (4+4) points.
   * - If the iterator points to a hex, then the quadrature points
   *   are the vertices, the line mid-points, and the face mid-points.
   *   This results in a quadrature object with 26 (8+12+6) points.
   *
   * The quadrature weights for these points are either chosen identically
   * and equal to one over the number of quadrature points (if @p with_laplace
   * is @p false), or in a way that gives points closer to the cell center
   * (measured on the reference cell) a higher weight. These weights correspond
   * to solving a Laplace equation and evaluating the solution at the quadrature
   * points (if @p with_laplace is @p true).
   *
   * The function is primarily used to construct the input argument
   * for the Manifold::get_new_point() function, which computes a new
   * point on a manifold based on a weighted average of "surrounding"
   * points represented by the quadrature points and weights stored in a
   * Quadrature object. This function creates such an object based on
   * the points that "surround" a cell, face, or edge, and weights
   * are chosen in a way appropriate for computing the new "mid-point"
   * of the object pointed to. An example of where this is necessary
   * is for mesh refinement, where (using the 2d situation as an example)
   * we need to first create new edge mid-points, and then a new cell-point.
   *
   * @param[in] iterator A mesh iterator that points to either a line, quad,
   *   or hex.
   * @param[in] with_laplace Whether or not to compute the quadrature weights
   *   by solving a Laplace equation, as discussed above.
   * @tparam MeshIteratorType An iterator type that corresponds to either
   *   Triangulation::cell_iterator (or variants such as
   *   Triangulation::active_cell_iterator or DoFHandler::cell_iterator) or
   *   that is the result of statements such as
   *   <code>cell-@>face(f)</code> or <code>cell-@>line(l)</code>.
   */
  template <typename MeshIteratorType>
  Quadrature<MeshIteratorType::AccessorType::space_dimension>
  get_default_quadrature(const MeshIteratorType &iterator,
                         const bool              with_laplace = false);
}


/**
 * Manifolds are used to describe the geometry of boundaries of domains as
 * well as the geometry of the interior. Manifold objects are therefore
 * associated with cells, faces, and/or edges, either by direct user action
 * or, if a user program does not do this explicitly, a default manifold
 * object is used.
 *
 * Manifolds are best understood by using the language of differential
 * geometry, but their common uses can be easily described simply through
 * examples.
 *
 *
 * <h3>Common use case: Creating a new vertex</h3>
 *
 * In the most essential use of manifolds, manifold descriptions are used
 * to create a "point between other points". For example, when a triangulation
 * creates a new vertex on a cell, face, or edge, it determines the new
 * vertex' coordinates through the following function call:
 *   @code
 *     ...
 *     Point<spacedim> new_vertex = manifold.get_new_point (quadrature);
 *     ...
 *   @endcode
 * Here, @p quadrature is a Quadrature<spacedim> object, which contains a collection
 * of points in @p spacedim dimension, and a collection of weights. The points
 * in this context will then be the vertices of the cell, face, or edge, and
 * the weights are typically one over the number of points when a new midpoint
 * of the cell, face, or edge is needed. Derived classes then will implement the
 * Manifold::get_new_point() function in a way that computes the location of this
 * new point. In the simplest case, for example in the FlatManifold class, the
 * function simply computes the arithmetic average (with given weights) of
 * the given points. However, other classes do something differently; for example,
 * the SphericalManifold class, which is used to describe domains that form (part of) the
 * sphere, will ensure that, given the two vertices of an edge at
 * the boundary, the new returned point will lie on the grand circle that connects
 * the two points, rather than choosing a point that is half-way between the
 * two points in ${\mathbb R}^d$.
 *
 *
 * @note Unlike almost all other cases in the library, we here interpret the points
 * in the quadrature object to be in real space, not on the reference cell.
 *
 * Manifold::get_new_point() has a default implementation that can simplify
 * this process somewhat:
 * Internally, the function calls the Manifold::project_to_manifold()
 * function after computing the weighted average of the quadrature points.
 * This allows derived classes to only overload Manifold::project_to_manifold()
 * for simple situations. This is often useful when describing manifolds that
 * are embedded in higher dimensional space, e.g., the surface of a sphere.
 * In those cases, the desired new point is simply the (weighted) average
 * of the provided point, projected back out onto the sphere.
 *
 *
 * <h3>Common use case: Computing tangent vectors</h3>
 *
 * The second use of this class is in computing directions on domains and
 * boundaries. For example, we may need to compute the normal vector to a
 * face in order to impose the no-flow boundary condition
 * $\mathbf u \cdot \mathbf n = 0$ (see the
 * VectorTools::compute_no_normal_flux_constraints() as an example). Similarly,
 * we may need normal vectors in the computation of the normal component of
 * the gradient of the numerical solution in order to compute the jump in the
 * gradient of the solution in error estimators (see, for example, the
 * KellyErrorEstimator class).
 *
 * To make this possible, the Manifold class provides a member function
 * (to be implemented by derived classes) that computes a "vector tangent
 * to the manifold at one point, in direction of another point" via the
 * Manifold::get_tangent_vector() function. For example, in 2d, one would
 * use this function with the two vertices of an edge at the boundary
 * to compute a "tangential" vector along the edge, and then get the normal
 * vector by rotation by 90 degrees. In 3d, one would compute the two
 * vectors "tangential" to the two edges of a boundary face adjacent to a
 * boundary vertex, and then take the cross product of these two to
 * obtain a vector normal to the boundary.
 *
 * For reasons that are more
 * difficult to understand, these direction vectors are normalized in a very
 * specific way, rather than to have unit norm. See the documentation of
 * Manifold::get_tangent_vector(), as well as below, for more information.
 *
 * In the simplest case (namely, the FlatManifold class), these tangent
 * vectors are just the difference vector between the two given points.
 * However, in more complicated (and more interesting) cases, the direction may
 * be different. For example, for the SphericalManifold case, if the two given
 * points lie on a common grand circle around the origin, then the tangent
 * vector will be tangential to the grand circle, rather than pointing straight
 * from one point to the other.
 *
 *
 * <h3>A unified description</h3>
 *
 * The "real" way to understand what this class does is to see it in the
 * framework of differential geometry. More specifically, differential geometry
 * is fundamentally based on the assumption that two sufficiently close points
 * are connected via a line of "shortest distance". This line is called a
 * "geodesic", and it is selected from all other lines that connect the two
 * points by the property that it is shortest if distances are measured in
 * terms of the "metric" that describes a manifold. To give examples, recall
 * that the geodesics of a flat manifold (implemented in the FlatManifold
 * class) are simply the straight lines connecting two points, whereas for
 * spherical manifolds (see the SphericalManifold class) geodesics between
 * two points of same distance are the grand circles, and are in general
 * curved lines when connecting two lines of different distance from the
 * origin.
 *
 * In the following discussion, and for the purposes of implementing the
 * current class, the concept of "metrics" that is so fundamental to
 * differential geometry is no longer of great importance to us. Rather,
 * everything can simply be described by postulating the existence of
 * geodesics connecting points on a manifold.
 *
 * Given geodesics, the operations discussed in the previous two sections
 * can be described in a more formal way. In essence, they rely on the
 * fact that we can assume that a geodesic is parameterized by a "time"
 * like variable $t$ so that $\mathbf s(t)$ describes the curve and so
 * that $\mathbf s(0)$ is the location of the first and $\mathbf s(1)$
 * the location of the second point. Furthermore, $\mathbf s(t)$ traces
 * out the geodesic at constant speed, covering equal distance in equal
 * time (as measured by the metric). Note that this parameterization
 * uses time, not arc length to denote progress along the geodesic.
 *
 * In this picture, computing a mid-point between points $\mathbf x_1$
 * and $\mathbf x_2$, with weights $w_1$ and $w_2=1-w_1$, simply
 * requires computing the point $\mathbf s(w_1)$. Computing a new
 * point as a weighted average of more than two points can be done
 * by considering pairwise geodetics, finding suitable points on
 * the geodetic between the first two points, then on the geodetic
 * between this new point and the third given point, etc.
 *
 * Likewise, the "tangential" vector described above is simply the
 * velocity vector, $\mathbf s'(t)$, evaluated at one of the end
 * points of a geodesic (i.e., at $t=0$ or $t=1$). In the case of a flat
 * manifold, the geodesic is simply the straight line connecting two points,
 * and the velocity vector is just the connecting vector in that
 * case. On the other hand, for two points on a spherical manifold,
 * the geodesic is a grand circle, and the velocity vector is
 * tangent to the spherical surface.
 *
 * Note that if we wanted to, we could use this to compute the length
 * of the geodesic that connects two points $\mathbf x_1$
 * and $\mathbf x_2$ by computing $\int_0^1 \|\mathbf s'(t)\| dt$
 * along the geodesic that connects them, but this operation will
 * not be of use to us in practice. One could also conceive
 * computing the direction vector using the "new point" operation
 * above, using the formula $\mathbf s'(0)=\lim_{w\rightarrow 0}
 * \frac{\mathbf s(w)-\mathbf s(0)}{w}$ where all we need to do
 * is compute the new point $\mathbf s(w)$ with weights $w$ and
 * $1-w$ along the geodesic connecting $\mathbf x_1$ and $\mathbf x_2$.
 * The default implementation of the function does this, by evaluating
 * the quotient for a small but finite weight $w$.
 * In practice, however, it is almost always possible to explicitly
 * compute the direction vector, i.e., without the need to numerically
 * approximate the limit process, and derived classes should do so.
 *
 *
 * @ingroup manifold
 * @author Luca Heltai, Wolfgang Bangerth, 2014, 2016
 */
template <int dim, int spacedim=dim>
class Manifold : public Subscriptor
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
   * Destructor. Does nothing here, but needs to be declared virtual to make
   * class hierarchies derived from this class possible.
   */
  virtual ~Manifold ();

  /**
   * @name Computing the location of points.
   */
  /// @{

  /**
   * Return the point which shall become the new vertex surrounded by the
   * given points which make up the quadrature. We use a quadrature object,
   * which should be filled with the surrounding points together with
   * appropriate weights.
   *
   * In its default implementation it calls internally the function
   * project_to_manifold. User classes can get away by simply implementing
   * that method.
   */
  virtual
  Point<spacedim>
  get_new_point(const Quadrature<spacedim> &quad) const;

  /**
   * Given a point which lies close to the given manifold, it modifies it and
   * projects it to manifold itself.
   *
   * This class is used by the default implementation of the function
   * get_new_point(). It should be made pure virtual, but for historical
   * reason, derived classes like Boundary<dim, spacedim> do not implement it.
   * The default behavior of this class, however, is to throw an exception
   * when called.
   *
   * If your manifold is simple, you could implement this function only, and
   * the default behavior should work out of the box.
   */
  virtual
  Point<spacedim> project_to_manifold (const std::vector<Point<spacedim> > &surrounding_points,
                                       const Point<spacedim> &candidate) const;

  /**
   * Backward compatibility interface.  Return the point which shall become
   * the new middle vertex of the two children of a regular line. In 2D, this
   * line is a line at the boundary, while in 3d, it is bounding a face at the
   * boundary (the lines therefore is also on the boundary).
   *
   * The default implementation of this function passes its argument to the
   * Manifolds::get_default_quadrature() function, and then calls the
   * Manifold<dim,spacedim>::get_new_point() function. User derived classes
   * can overload Manifold<dim,spacedim>::get_new_point() or
   * Manifold<dim,spacedim>::project_to_manifold(), which is called by the
   * default implementation of Manifold<dim,spacedim>::get_new_point().
   */
  virtual
  Point<spacedim>
  get_new_point_on_line (const typename Triangulation<dim,spacedim>::line_iterator &line) const;

  /**
   * Backward compatibility interface. Return the point which shall become the
   * common point of the four children of a quad at the boundary in three or
   * more spatial dimensions. This function therefore is only useful in at
   * least three dimensions and should not be called for lower dimensions.
   *
   * This function is called after the four lines bounding the given @p quad
   * are refined, so you may want to use the information provided by
   * <tt>quad->line(i)->child(j)</tt>, <tt>i=0...3</tt>, <tt>j=0,1</tt>.
   *
   * The default implementation of this function passes its argument to the
   * Manifolds::get_default_quadrature() function, and then calls the
   * Manifold<dim,spacedim>::get_new_point() function. User derived classes
   * can overload Manifold<dim,spacedim>::get_new_point() or
   * Manifold<dim,spacedim>::project_to_manifold(), which is called by the
   * default implementation of Manifold<dim,spacedim>::get_new_point().
   */
  virtual
  Point<spacedim>
  get_new_point_on_quad (const typename Triangulation<dim,spacedim>::quad_iterator &quad) const;

  /**
   * Backward compatibility interface.  Return the point which shall become
   * the common point of the eight children of a hex in three or spatial
   * dimensions. This function therefore is only useful in at least three
   * dimensions and should not be called for lower dimensions.
   *
   * This function is called after the all the bounding objects of the given
   * @p hex are refined, so you may want to use the information provided by
   * <tt>hex->quad(i)->line(j)->child(k)</tt>, <tt>i=0...5</tt>,
   * <tt>j=0...3</tt>, <tt>k=0,1</tt>.
   *
   * The default implementation of this function passes its argument to the
   * Manifolds::get_default_quadrature() function, and then calls the
   * Manifold<dim,spacedim>::get_new_point() function. User derived classes
   * can overload Manifold<dim,spacedim>::get_new_point() or
   * Manifold<dim,spacedim>::project_to_manifold(), which is called by the
   * default implementation of Manifold<dim,spacedim>::get_new_point().
   */
  virtual
  Point<spacedim>
  get_new_point_on_hex (const typename Triangulation<dim,spacedim>::hex_iterator &hex) const;


  /**
   * Backward compatibility interface. Depending on <tt>dim=2</tt> or
   * <tt>dim=3</tt> this function calls the get_new_point_on_line or the
   * get_new_point_on_quad function. It throws an exception for
   * <tt>dim=1</tt>. This wrapper allows dimension independent programming.
   */
  Point<spacedim>
  get_new_point_on_face (const typename Triangulation<dim,spacedim>::face_iterator &face) const;


  /**
   * Backward compatibility interface.  Depending on <tt>dim=1</tt>,
   * <tt>dim=2</tt> or <tt>dim=3</tt> this function calls the
   * get_new_point_on_line, get_new_point_on_quad or the get_new_point_on_hex
   * function. This wrapper allows dimension independent programming.
   */
  Point<spacedim>
  get_new_point_on_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell) const;

  /// @}

  /**
   * @name Computing tangent vectors
   */
  /// @{

  /**
   * Return a vector that, at $\mathbf x_1$, is tangential to
   * the geodesic that connects two points $\mathbf x_1,\mathbf x_2$. The geodesic
   * is the shortest line between these two points, where "shortest" is defined
   * via a metric specific to a particular implementation of this class in a
   * derived class. For example, in the case of a FlatManifold, the shortest
   * line between two points is just the straight line, and in this case the
   * tangent vector is just the difference $\mathbf d=\mathbf x_2-\mathbf x_1$.
   * On the other hand, for a manifold that describes a surface embedded in
   * a higher dimensional space (e.g., the surface of a sphere), then the
   * tangent vector is tangential to the surface, and consequently may point in
   * a different direction than the straight line that connects the two points.
   *
   * While tangent vectors are often normalized to unit length, the vectors
   * returned by this function are normalized as described in the introduction
   * of this class. Specifically, if $\mathbf s(t)$ traces out the geodesic
   * between the two points where $\mathbf x_1 = \mathbf s(0)$ and
   * $\mathbf x_2 = \mathbf s(1)$, then the returned vector must equal
   * $\mathbf s'(0)$. In other words, the norm of the returned vector also
   * encodes, in some sense, the <i>length</i> of the geodesic because a curve
   * $\mathbf s(t)$ must move "faster" if the two points it connects between
   * arguments $t=0$ and $t=1$ are farther apart.
   *
   * The default implementation of this function approximates
   * $\mathbf s'(0) \approx \frac{$\mathbf s(\epsilon)-\mathbf x_1}{\epsilon}$
   * for a small value of $\epsilon$, and the evaluation of $\mathbf s(\epsilon)$
   * is done by calling get_new_point(). If possible, derived classes should
   * override this function by an implemention of the exact derivative.
   *
   * @param x1 The first point that describes the geodesic, and the one
   *   at which the "direction" is to be evaluated.
   * @param x2 The second point that describes the geodesic.
   * @return A "direction" vector tangential to the geodesic.
   */
  virtual
  Tensor<1,spacedim>
  get_tangent_vector (const Point<spacedim> &x1,
                      const Point<spacedim> &x2) const;

  /// @}

  /**
   * @name Computing normal vectors
   */
  /// @{

  /**
   * Return the normal vector to a face embedded in this manifold, at
   * the point p. If p is not in fact on the surface, but only
   * close-by, try to return something reasonable, for example the
   * normal vector at the surface point closest to p.  (The point p
   * will in fact not normally lie on the actual surface, but rather
   * be a quadrature point mapped by some polynomial mapping; the
   * mapped surface, however, will not usually coincide with the
   * actual surface.)
   *
   * The face iterator gives an indication which face this function is
   * supposed to compute the normal vector for.  This is useful if the
   * boundary of the domain is composed of different nondifferential
   * pieces (for example when using the StraightBoundary class to
   * approximate a geometry that is completely described by the coarse
   * mesh, with piecewise (bi-)linear components between the vertices,
   * but where the boundary may have a kink at the vertices itself).
   *
   * @note The default implementation of this function computes the
   * normal vector by taking the cross product between the tangent
   * vectors from p to the most orthogonal and further non consecutive
   * vertices of the face.
   */
  virtual
  Tensor<1,spacedim>
  normal_vector (const typename Triangulation<dim,spacedim>::face_iterator &face,
                 const Point<spacedim> &p) const;

  /**
   * Compute the normal vectors to the boundary at each vertex of the
   * given face embedded in the Manifold. It is not required that the
   * normal vectors be normed somehow.  Neither is it required that
   * the normals actually point outward.
   *
   * This function is needed to compute data for C1 mappings. The
   * default implementation calls normal_vector() on each vertex.
   *
   * Note that when computing normal vectors at a vertex where the
   * boundary is not differentiable, you have to make sure that you
   * compute the one-sided limits, i.e. limit with respect to points
   * inside the given face.
   */
  virtual
  void
  get_normals_at_vertices (const typename Triangulation<dim,spacedim>::face_iterator &face,
                           FaceVertexNormals &face_vertex_normals) const;

  /// @}
};


/**
 * Specialization of Manifold<dim,spacedim>, which represent a possibly
 * periodic Euclidean space of dimension @p dim embedded in the Euclidean
 * space of @p spacedim dimensions. The main characteristic of this Manifold
 * is the fact that the function
 * FlatManifold<dim,spacedim>::project_to_manifold() is the identity function.
 *
 * @ingroup manifold
 *
 * @author Luca Heltai, 2014
 */
template <int dim, int spacedim=dim>
class FlatManifold : public Manifold<dim, spacedim>
{
public:
  /**
   * Default constructor. The optional argument can be used to specify the
   * periodicity of the spacedim-dimensional manifold (one period per
   * direction). A periodicity value of zero means that along that direction
   * there is no periodicity. By default no periodicity is assumed.
   *
   * Periodicity affects the way a middle point is computed. It is assumed
   * that if two points are more than half period distant, then the distance
   * should be computed by crossing the periodicity boundary, i.e., the
   * average is computed by adding a full period to the sum of the two. For
   * example, if along direction 0 we have 2*pi periodicity, then the average
   * of (2*pi-eps) and (eps) is not pi, but 2*pi (or zero), since, on a
   * periodic manifold, these two points are at distance 2*eps and not (2*pi-
   * eps). Special cases are taken into account, to ensure that the behavior
   * is always as expected. The third argument is used as a relative tolerance
   * when computing distances.
   *
   * Periodicity will be intended in the following way: the domain is
   * considered to be the box contained in [Point<spacedim>(), periodicity)
   * where the right extreme is excluded. If any of the components of this box
   * has zero length, then no periodicity is assumed in that direction.
   * Whenever a function that tries to compute averages is called, an
   * exception will be thrown if one of the points which you are using for the
   * average lies outside the periodicity box. The return points are
   * guaranteed to lie in the periodicity box plus or minus
   * tolerance*periodicity.norm().
   */
  FlatManifold (const Tensor<1,spacedim> &periodicity = Tensor<1,spacedim>(),
                const double tolerance=1e-10);

  /**
   * Let the new point be the average sum of surrounding vertices.
   *
   * This particular implementation constructs the weighted average of the
   * surrounding points, and then calls internally the function
   * project_to_manifold(). The reason why we do it this way, is to allow lazy
   * programmers to implement only the project_to_manifold() function for their
   * own Manifold classes which are small (or trivial) perturbations of a flat
   * manifold. This is the case whenever the coarse mesh is a decent
   * approximation of the manifold geometry. In this case, the middle point of
   * a cell is close to true middle point of the manifold, and a projection
   * may suffice.
   *
   * For most simple geometries, it is possible to get reasonable results by
   * deriving your own Manifold class from FlatManifold, and write a new
   * interface only for the project_to_manifold function. You will have good
   * approximations also with large deformations, as long as in the coarsest
   * mesh size you are trying to refine, the middle point is not too far from
   * the manifold mid point, i.e., as long as the coarse mesh size is small
   * enough.
   */
  virtual
  Point<spacedim>
  get_new_point(const Quadrature<spacedim> &quad) const;


  /**
   * Project to FlatManifold. This is the identity function for flat,
   * Euclidean spaces. Note however that this function can be overloaded by
   * derived classes, which will then benefit from the logic behind the
   * get_new_point() function which are often very similar (if not identical) to
   * the one implemented in this class.
   */
  virtual
  Point<spacedim>
  project_to_manifold (const std::vector<Point<spacedim> > &points,
                       const Point<spacedim> &candidate) const;

  /**
   * Return a vector that, at $\mathbf x_1$, is tangential to
   * the geodesic that connects two points $\mathbf x_1,\mathbf x_2$.
   * For the current class, we assume that the manifold is flat, so
   * the geodesic is the straight line between the two points, and we
   * return $\mathbf x_2-\mathbf x_1$. The normalization of the vector
   * is chosen so that it fits the convention described in
   * Manifold::get_tangent_vector().
   *
   * @note If you use this class as a stepping stone to build a manifold
   *   that only "slightly" deviates from a flat manifold, by overloading
   *   the project_to_manifold() function.
   *
   * @param x1 The first point that describes the geodesic, and the one
   *   at which the "direction" is to be evaluated.
   * @param x2 The second point that describes the geodesic.
   * @return A "direction" vector tangential to the geodesic. Here, this is
   *   $\mathbf x_2-\mathbf x_1$, possibly modified by the periodicity of
   *   the domain as set in the constructor, to use the "shortest" connection
   *   between the points through the periodic boundary as necessary.
   */
  virtual
  Tensor<1,spacedim>
  get_tangent_vector (const Point<spacedim> &x1,
                      const Point<spacedim> &x2) const;

  /**
   * Return the periodicity of this Manifold.
   */
  const Tensor<1,spacedim> &get_periodicity() const;

private:
  /**
   * The periodicity of this Manifold. Periodicity affects the way a middle
   * point is computed. It is assumed that if two points are more than half
   * period distant, then the distance should be computed by crossing the
   * periodicity boundary, i.e., the average is computed by adding a full
   * period to the sum of the two. For example, if along direction 0 we have
   * 2*pi periodicity, then the average of (2*pi-eps) and (eps) is not pi, but
   * 2*pi (or zero), since, on a periodic manifold, these two points are at
   * distance 2*eps and not (2*pi-eps).
   *
   * A periodicity 0 along one direction means no periodicity. This is the
   * default value for all directions.
   */
  const Tensor<1,spacedim> periodicity;

  DeclException3(ExcPeriodicBox, int, Point<spacedim>, double,
                 << "The component number " << arg1 << " of the point [ " << arg2
                 << " ] is not in the interval [ 0, " << arg3 << "), bailing out.");

  /**
   * Relative tolerance. This tolerance is used to compute distances in double
   * precision.
   */
  const double tolerance;
};


/**
 * This class describes mappings that can be expressed in terms of charts.
 * Specifically, this class with its template arguments describes a chart of
 * dimension chartdim, which is part of a Manifold<dim,spacedim> and is used
 * in an object of type Triangulation<dim,spacedim>:  It specializes a
 * Manifold of dimension chartdim embedded in a manifold of dimension
 * spacedim, for which you have explicit pull_back() and push_forward()
 * transformations. Its use is explained in great detail in step-53.
 *
 * This is a helper class which is useful when you have an explicit map from
 * an Euclidean space of dimension chartdim to an Euclidean space of dimension
 * spacedim which represents your manifold, i.e., when your manifold
 * $\mathcal{M}$ can be represented by a map \f[ F: \mathcal{B} \subset
 * R^{\text{chartdim}} \mapsto \mathcal{M} \subset R^{\text{spacedim}} \f]
 * (the push_forward() function) and that admits the inverse transformation
 * \f[ F^{-1}: \mathcal{M} \subset R^{\text{spacedim}} \mapsto \mathcal{B}
 * \subset R^{\text{chartdim}} \f] (the pull_back() function).
 *
 * The get_new_point() function of the ChartManifold class is implemented by
 * calling the pull_back() method for all <tt>surrounding_points</tt>,
 * computing their weighted average in the chartdim Euclidean space, and
 * calling the push_forward() method with the resulting point, i.e., \f[
 * \mathbf x^{\text{new}} = F(\sum_i w_i F^{-1}(\mathbf x_i)).  \f]
 *
 * Derived classes are required to implement the push_forward() and the
 * pull_back() methods. All other functions (with the exception of the
 * push_forward_gradient() function, see below) that are required by mappings
 * will then be provided by this class.
 *
 *
 * <h3>Providing function gradients</h3>
 *
 * In order to compute vectors that are tangent to the manifold (for example,
 * tangent to a surface embedded in higher dimensional space, or simply the
 * three unit vectors of ${\mathbb R}^3$), one needs to also have access
 * to the <i>gradient</i> of the push-forward function $F$. The gradient
 * is the matrix ${\nabla F)_{ij}=\partial_j F_i$, where we take the derivative
 * with regard to the chartdim reference coordinates on the flat Euclidean
 * space in which $\mathcal B$ is located. In other words, at a point
 * $\mathbf x$, $\nabla F(\mathbf x)$ is a matrix of size @p spacedim
 * times @p chartdim.
 *
 * Only the ChartManifold::get_tangent_vector() function uses the gradient
 * of the push-forward, but only a subset of all finite element codes
 * actually require the computation of tangent vectors. Consequently,
 * while derived classes need to implement the abstract virtual push_forward()
 * and pull_back() functions of this class, they do not need to implement
 * the virtual push_forward_gradient() function. Rather, that function has a
 * default implementation (and consequently is not abstract, therefore not
 * forcing derived classes to overload it), but the default implementation
 * clearly can not compute anything useful and therefore simply triggers
 * and exception.
 *
 *
 * <h3>A note on the template arguments</h3>
 *
 * The dimension arguments @p chartdim, @p dim and @p spacedim must satisfy
 * the following relationships:
 *   @code
 *      dim <= spacedim
 *      chartdim <= spacedim
 *   @endcode
 * However, there is no a priori relationship between @p dim and @p chartdim.
 * For example, if you want to describe a mapping for an edge (a 1d object) in
 * a 2d triangulation embedded in 3d space, you could do so by parameterizing
 * it via a line
 *   @f[
 *      F: [0,1] \rightarrow {\mathbb R}^3
 *   @f]
 * in which case @p chartdim is 1. On the other hand, there is no reason why
 * one can't describe this as a mapping
 *   @f[
 *      F: {\mathbb R}^3 \rightarrow {\mathbb R}^3
 *   @f]
 * in such a way that the line $[0,1]\times \{0\}\times \{0\}$ happens to be
 * mapped onto the edge in question. Here, @p chartdim is 3. This may seem
 * cumbersome but satisfies the requirements of an invertible function $F$
 * just fine as long as it is possible to get from the edge to the pull-back
 * space and then back again. Finally, given that we are dealing with a 2d
 * triangulation in 3d, one will often have a mapping from, say, the 2d unit
 * square or unit disk to the domain in 3d space, and the edge in question may
 * simply be the mapped edge of the unit domain in 2d space. In this case, @p
 * chartdim is 2.
 *
 * @ingroup manifold
 *
 * @author Luca Heltai, 2013, 2014
 */
template <int dim, int spacedim=dim, int chartdim=dim>
class ChartManifold : public Manifold<dim,spacedim>
{
public:
  /**
   * Constructor. The optional argument can be used to specify the periodicity
   * of the chartdim-dimensional manifold (one period per direction). A
   * periodicity value of zero means that along that direction there is no
   * periodicity. By default no periodicity is assumed.
   *
   * Periodicity affects the way a middle point is computed. It is assumed
   * that if two points are more than half period distant, then the distance
   * should be computed by crossing the periodicity boundary, i.e., then the
   * average is computed by adding a full period to the sum of the two. For
   * example, if along direction 0 we have 2*pi periodicity, then the average
   * of (2*pi-eps) and (eps) is not pi, but 2*pi (or zero), since, on the
   * manifold, these two points are at distance 2*eps and not (2*pi-eps)
   */
  ChartManifold(const Tensor<1,chartdim> &periodicity = Tensor<1,chartdim>());

  /**
   * Destructor. Does nothing here, but needs to be declared to make it
   * virtual.
   */
  virtual ~ChartManifold ();


  /**
   * Refer to the general documentation of this class and the documentation of
   * the base class for more information.
   */
  virtual
  Point<spacedim>
  get_new_point(const Quadrature<spacedim> &quad) const;

  /**
   * Pull back the given point in spacedim to the Euclidean chartdim
   * dimensional space.
   *
   * Refer to the general documentation of this class for more information.
   */
  virtual
  Point<chartdim>
  pull_back(const Point<spacedim> &space_point) const = 0;

  /**
   * Given a point in the chartdim dimensional Euclidean space, this method
   * returns a point on the manifold embedded in the spacedim Euclidean space.
   *
   * Refer to the general documentation of this class for more information.
   */
  virtual
  Point<spacedim>
  push_forward(const Point<chartdim> &chart_point) const = 0;

  /**
   * Given a point in the chartdim dimensional Euclidean space, this method
   * returns the derivatives of the function $F$ that maps from the
   * chartdim-dimensional to the spacedim-dimensional space. In other
   * words, it is a matrix of size $\text{spacedim}\times\text{chartdim}$.
   *
   * This function is used in the computations required by the
   * get_tangent_vector() function. Since not all users of the Manifold
   * class interface will require calling that function, the current
   * function is implemented but will trigger an exception whenever
   * called. This allows derived classes to avoid implementing the
   * push_forward_gradient function if this functionality is not
   * needed in the user program.
   *
   * Refer to the general documentation of this class for more information.
   */
  virtual
  DerivativeForm<1,chartdim,spacedim>
  push_forward_gradient(const Point<chartdim> &chart_point) const;

  /**
   * Return a vector that, at $\mathbf x_1$, is tangential to
   * the geodesic that connects two points $\mathbf x_1,\mathbf x_2$.
   * See the documentation of the Manifold class and of
   * Manifold::get_tangent_vector() for a more detailed description.
   *
   * For the current class, we assume that this geodesic is the image
   * under the push_forward() operation of a straight line of the
   * pre-images of @p x1 and @p x2 (where pre-images are computed by pulling
   * back the locations @p x1 and @p x2). In other words, if these
   * preimages are $\xi_1=F^{-1}(\mathbf x_1), \xi_2=F^{-1}(\mathbf x_2)$,
   * then the geodesic in preimage (the chartdim-dimensional Euclidean) space
   * is
   * @f{align*}{
   *   \zeta(t) &= \xi_1 +  t (\xi_2-\xi_1)
   *  \\          &= F^{-1}(\mathbf x_1) + t\left[F^{-1}(\mathbf x_2)
   *                                             -F^{-1}(\mathbf x_1)\right]
   * @f}
   * In image space, i.e., in the space in which we operate, this
   * leads to the curve
   * @f{align*}{
   *   \mathbf s(t) &= F(s(t)
   *  \\          &= F(\xi_1 +  t (\xi_2-\xi_1))
   *  \\          &= F\left(F^{-1}(\mathbf x_1) + t\left[F^{-1}(\mathbf x_2)
   *                                     -F^{-1}(\mathbf x_1)\right]\right).
   * @f}
   * What the current function is supposed to return is $\mathbf s'(0)$. By
   * the chain rule, this is equal to
   * @f{align*}{
   *   \mathbf s'(0) &=
   *     \frac{d}{dt}\left. F\left(F^{-1}(\mathbf x_1)
   *                        + t\left[F^{-1}(\mathbf x_2)
   *                                 -F^{-1}(\mathbf x_1)\right]\right)
   *                 \right|_{t=0}
   * \\ &= \nabla_\xi F\left(F^{-1}(\mathbf x_1)\right)
   *                    \left[F^{-1}(\mathbf x_2)
   *                                 -F^{-1}(\mathbf x_1)\right].
   * @f}
   * This formula may then have to be slightly modified by
   * considering any periodicity that was assumed in the call to
   * the constructor.
   *
   * Thus, the computation of tangent vectors also requires the
   * implementation of <i>derivatives</i> $\nabla_\xi F(\xi)$ of
   * the push-forward mapping. Here, $F^{-1}(\mathbf x_2)-F^{-1}(\mathbf x_1)$
   * is a chartdim-dimensional vector, and $\nabla_\xi F\left(F^{-1}(\mathbf x_1)\right)
   * = \nabla_\xi F\left(\xi_1\right)$ is a spacedim-times-chartdim-dimensional
   * matrix. Consequently, and as desired, the operation results in a
   * spacedim-dimensional vector.
   *
   * @param x1 The first point that describes the geodesic, and the one
   *   at which the "direction" is to be evaluated.
   * @param x2 The second point that describes the geodesic.
   * @return A "direction" vector tangential to the geodesic.
   */
  virtual
  Tensor<1,spacedim>
  get_tangent_vector (const Point<spacedim> &x1,
                      const Point<spacedim> &x2) const;

  /**
   * Return the periodicity associated with the submanifold.
   */
  const Tensor<1,chartdim> &get_periodicity() const;

private:
  /**
   * The sub_manifold object is used to compute the average of the points in
   * the chart coordinates system.
   */
  const FlatManifold<dim,chartdim> sub_manifold;
};




/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN

template <>
Point<1>
Manifold<1,1>::
get_new_point_on_face (const Triangulation<1,1>::face_iterator &) const;

template <>
Point<2>
Manifold<1,2>::
get_new_point_on_face (const Triangulation<1,2>::face_iterator &) const;


template <>
Point<3>
Manifold<1,3>::
get_new_point_on_face (const Triangulation<1,3>::face_iterator &) const;


template <>
Point<1>
Manifold<1,1>::
get_new_point_on_quad (const Triangulation<1,1>::quad_iterator &) const;

template <>
Point<2>
Manifold<1,2>::
get_new_point_on_quad (const Triangulation<1,2>::quad_iterator &) const;


template <>
Point<3>
Manifold<1,3>::
get_new_point_on_quad (const Triangulation<1,3>::quad_iterator &) const;


template <>
Point<3>
Manifold<3,3>::
get_new_point_on_hex (const Triangulation<3,3>::hex_iterator &) const;

/*---Templated functions---*/

namespace Manifolds
{
  template <typename MeshIteratorType>
  Quadrature<MeshIteratorType::AccessorType::space_dimension>
  get_default_quadrature(const MeshIteratorType &iterator,
                         const bool              with_laplace)
  {
    const int spacedim = MeshIteratorType::AccessorType::space_dimension;
    const int dim = MeshIteratorType::AccessorType::structure_dimension;

    std::vector<Point<spacedim> > sp;
    std::vector<double> wp;


    // note that the exact weights are chosen such as to minimize the
    // distortion of the four new quads from the optimal shape; their
    // derivation and values is copied over from the
    // @p{MappingQ::set_laplace_on_vector} function
    switch (dim)
      {
      case 1:
        sp.resize(2);
        wp.resize(2);
        sp[0] = iterator->vertex(0);
        wp[0] = .5;
        sp[1] = iterator->vertex(1);
        wp[1] = .5;
        break;
      case 2:
        sp.resize(8);
        wp.resize(8);

        for (unsigned int i=0; i<4; ++i)
          {
            sp[i] = iterator->vertex(i);
            sp[4+i] = ( iterator->line(i)->has_children() ?
                        iterator->line(i)->child(0)->vertex(1) :
                        iterator->line(i)->get_manifold().get_new_point_on_line(iterator->line(i)) );
          }

        if (with_laplace)
          {
            std::fill(wp.begin(), wp.begin()+4, 1.0/16.0);
            std::fill(wp.begin()+4, wp.end(), 3.0/16.0);
          }
        else
          std::fill(wp.begin(), wp.end(), 1.0/8.0);
        break;
      case 3:
      {
        TriaIterator<TriaAccessor<3, 3, 3> > hex
          = static_cast<TriaIterator<TriaAccessor<3, 3, 3> > >(iterator);
        const unsigned int np =
          GeometryInfo<dim>::vertices_per_cell+
          GeometryInfo<dim>::lines_per_cell+
          GeometryInfo<dim>::faces_per_cell;
        sp.resize(np);
        wp.resize(np);
        std::vector<Point<3> > *sp3 = reinterpret_cast<std::vector<Point<3> > *>(&sp);

        unsigned int j=0;

        // note that the exact weights are chosen such as to minimize the
        // distortion of the eight new hexes from the optimal shape; their
        // derivation and values is copied over from the
        // @p{MappingQ::set_laplace_on_vector} function
        for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i, ++j)
          {
            (*sp3)[j] = hex->vertex(i);
            wp[j] = 1.0/128.0;
          }
        for (unsigned int i=0; i<GeometryInfo<dim>::lines_per_cell; ++i, ++j)
          {
            (*sp3)[j] = (hex->line(i)->has_children() ?
                         hex->line(i)->child(0)->vertex(1) :
                         hex->line(i)->get_manifold().get_new_point_on_line(hex->line(i)));
            wp[j] = 7.0/192.0;
          }
        for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i, ++j)
          {
            (*sp3)[j] = (hex->quad(i)->has_children() ?
                         hex->quad(i)->isotropic_child(0)->vertex(3) :
                         hex->quad(i)->get_manifold().get_new_point_on_quad(hex->quad(i)));
            wp[j] = 1.0/12.0;
          }
        // Overwrite the weights with 1/np if we don't want to use
        // laplace vectors.
        if (with_laplace == false)
          std::fill(wp.begin(), wp.end(), 1.0/np);
      }
      break;
      default:
        Assert(false, ExcInternalError());
        break;
      }
    return Quadrature<spacedim>(sp,wp);
  }
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

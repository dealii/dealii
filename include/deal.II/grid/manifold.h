// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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

#ifndef __deal2__tria_manifold_h
#define __deal2__tria_manifold_h


/*----------------------------   manifold.h     ---------------------------*/

#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/point.h>
#include <deal.II/grid/tria.h>

DEAL_II_NAMESPACE_OPEN

template <int dim, int space_dim> class Triangulation;


/** We collect here some helper functions used in the
    Manifold<dim,spacedim> classes.
 */
namespace Manifolds {
  /** Given a hex iterator, construct a quadrature with the Laplace
      weigths, and all relevant points of the hex: vertices, line
      centers and face centers, which can be called when creating
      middle vertices in the manifold routines.*/
  void get_default_quadrature(const TriaIterator<CellAccessor<3, 3> >& hex,
                              Quadrature<3> &quad);
  
  /** Given a general mesh iterator, construct a quadrature with the
      Laplace weigths or with uniform weights according the parameter
      @p with_laplace, and with all relevant points of the iterator:
      vertices, line centers and/or face centers, which can be called
      when creating new vertices in the manifold routines.*/
  template <typename OBJECT, int spacedim>
  void get_default_quadrature(const OBJECT& obj, 
			      Quadrature<spacedim> &quad, 
			      bool with_laplace = false);
}


/**
 *   This class is used to represent a manifold to a triangulation.
 *   When a triangulation creates a new vertex on this manifold, it
 *   determines the new vertex' coordinates through the following
 *   function:
 *
 *   @code
 *     ...
 *     Point<spacedim> new_vertex = manifold.get_new_point (quadrature);
 *     ...
 *   @endcode
 *   @p quadrature is a Quadrature<spacedim> object, which contains a
 *   collection of points in @p spacedim dimension, and a collection of
 *   weights.
 *
 *   Internally, the get_new_point() function calls the
 *   project_to_manifold() function after computing the weighted
 *   average of the quadrature poitns. This allows end users to only
 *   overload project_to_manifold() for simple situations.
 *
 *   Should a finer control be necessary, then get_new_point() can be
 *   overloaded.  For backward compatibility, this function also
 *   offers an interface which is compatible with
 *   Boundary<dim,spacedim>, which are all derived from
 *   FlatManifold<dim,spacedim>, allowing old user codes to keep using
 *   their boundary descriptors as Manifold<dim,spacedim> objects.
 *
 *   The default behavior of these backward compatible interfaces is
 *   to construct a Quadrature<spacedim> object containting the
 *   vertices, midpoints of lines, and midpoints of quads with the
 *   correct weight, and call get_new_point() with this quadrature. If
 *   you need finer tuning for lines, quads or hexes, you can overload
 *   any of the get_new_point_on_* functions. 
 *
 *   FlatManifold is the specialization from which StraigthBoundary is
 *   derived, where the project_to_manifold() function is the identity.
 *
 * @ingroup manifold
 * @author Luca Heltai, 2014
 */
template <int dim, int spacedim=dim>
class Manifold : public Subscriptor
{
public:


  /**
   * Destructor. Does nothing here, but needs to be declared to make it
   * virtual.
   */
  virtual ~Manifold ();

  /**
   * Return the point which shall become the new vertex surrounded by
   * the given points which make up the quadrature. We use a
   * quadrature object, which should be filled with the surrounding
   * points together with appropriate weights.
   *
   * In its default implementation it calls internally the function
   * project to manifold. User classes can get away by simply
   * implementing that method.
   */
  virtual
  Point<spacedim>
  get_new_point(const Quadrature<spacedim> &quad) const;

  /**
   * Given a point which lies close to the given manifold, it modifies
   * it and projects it to manifold itself.
   * 
   * This class is used by the default implementation of the function
   * get_new_point(). It should be made pure virtual, but for
   * historical reason, derived classes like Boundary<dim, spacedim>
   * do not implement it. The default behavior of this class, however,
   * is to throw an exception when called.
   *
   * If your manifold is simple, you could implement this function
   * only, and the default behavior should work out of the box.
   */
  virtual
  Point<spacedim> project_to_manifold (const std::vector<Point<spacedim> > &surrounding_points,
                                       const Point<spacedim> &candidate) const;

  /**
   * Backward compatibility interface.  Return the point which shall
   * become the new middle vertex of the two children of a regular
   * line. In 2D, this line is a line at the boundary, while in 3d, it
   * is bounding a face at the boundary (the lines therefore is also
   * on the boundary).
   *
   * The default implementation of this function
   * passes its argument to the Manifolds::get_default_quadrature()
   * function, and then calls the
   * Manifold<dim,spacedim>::get_new_point() function. User derived
   * classes can overload Manifold<dim,spacedim>::get_new_point() or
   * Manifold<dim,spacedim>::project_to_surface(), which is called by
   * the default implementation of
   * Manifold<dim,spacedim>::get_new_point().
   */
  virtual
  Point<spacedim>
  get_new_point_on_line (const typename Triangulation<dim,spacedim>::line_iterator &line) const;

  /**
   * Backward compatibility interface. Return the point which shall
   * become the common point of the four children of a quad at the
   * boundary in three or more spatial dimensions. This function
   * therefore is only useful in at least three dimensions and should
   * not be called for lower dimensions.
   *
   * This function is called after the four lines bounding the given @p quad
   * are refined, so you may want to use the information provided by
   * <tt>quad->line(i)->child(j)</tt>, <tt>i=0...3</tt>, <tt>j=0,1</tt>.
   *
   * The default implementation of this function passes its argument
   * to the Manifolds::get_default_quadrature() function, and then
   * calls the Manifold<dim,spacedim>::get_new_point() function. User
   * derived classes can overload
   * Manifold<dim,spacedim>::get_new_point() or
   * Manifold<dim,spacedim>::project_to_surface(), which is called by
   * the default implementation of
   * Manifold<dim,spacedim>::get_new_point().
   */
  virtual
  Point<spacedim>
  get_new_point_on_quad (const typename Triangulation<dim,spacedim>::quad_iterator &quad) const;

  /**
   * Backward compatibility interface.  Return the point which shall
   * become the common point of the eight children of a hex in three
   * or spatial dimensions. This function therefore is only useful in
   * at least three dimensions and should not be called for lower
   * dimensions.
   *
   * This function is called after the all the bounding objects of the
   * given @p hex are refined, so you may want to use the information
   * provided by <tt>hex->quad(i)->line(j)->child(k)</tt>,
   * <tt>i=0...5</tt>, <tt>j=0...3</tt>, <tt>k=0,1</tt>.
   *
   * The default implementation of this function passes its argument
   * to the Manifolds::get_default_quadrature() function, and then
   * calls the Manifold<dim,spacedim>::get_new_point() function. User
   * derived classes can overload
   * Manifold<dim,spacedim>::get_new_point() or
   * Manifold<dim,spacedim>::project_to_surface(), which is called by
   * the default implementation of
   * Manifold<dim,spacedim>::get_new_point().
   */
  virtual
  Point<spacedim>
  get_new_point_on_hex (const typename Triangulation<dim,spacedim>::hex_iterator &hex) const;


  /**
   * Backward compatibility interface. Depending on <tt>dim=2</tt> or
   * <tt>dim=3</tt> this function calls the get_new_point_on_line or
   * the get_new_point_on_quad function. It throws an exception for
   * <tt>dim=1</tt>. This wrapper allows dimension independent
   * programming.
   */
  Point<spacedim>
  get_new_point_on_face (const typename Triangulation<dim,spacedim>::face_iterator &face) const;


  /**
   * Backward compatibility interface.  Depending on <tt>dim=1</tt>,
   * <tt>dim=2</tt> or <tt>dim=3</tt> this function calls the
   * get_new_point_on_line, get_new_point_on_quad or the
   * get_new_point_on_hex function. This wrapper allows dimension
   * independent programming.
   */
  Point<spacedim>
  get_new_point_on_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell) const;
};


/**
 *   Specialization of Manifold<dim,spacedim>, which represent a
 *   possibly periodic Euclidean space of dimension @p dim embedded in
 *   the Euclidean space of @p spacedim dimensions. The main
 *   characteristic of this Manifold is the fact that the function
 *   FlatManifold<dim,spacedim>::project_to_manifold() is the identity
 *   function.
 *
 *   @ingroup manifold
 *
 *   @author Luca Heltai, 2014
 */
template <int dim, int spacedim=dim>
class FlatManifold: public Manifold<dim, spacedim>
{
public:
  /**
   * Default constructor.The optional argument can be used to specify
   * the periodicity of the spacedim-dimensional manifold (one period
   * per direction). A peridicity value of zero means that along that
   * direction there is no peridicity. By default no periodicity is
   * assumed.
   *
   * Periodicity affects the way a middle point is computed. It is
   * assumed that if two points are more than half period distant,
   * then the distance should be computed by crossing the periodicity
   * boundary, i.e., the average is computed by adding a full period
   * to the sum of the two. For example, if along direction 0 we have
   * 2*pi periodicity, then the average of (2*pi-eps) and (eps) is not
   * pi, but 2*pi (or zero), since, on a periodic manifold, these two
   * points are at distance 2*eps and not (2*pi-eps). Special cases
   * are taken into account, to ensure that the behavior is always as
   * expected.
   *
   * Periodicity will be intended in the following way: the domain is
   * considered to be the box contained in [Point<spacedim>(),
   * periodicity) where the right extreme is excluded. If any of the
   * components of this box has zero length, then no periodicity is
   * computed in that direction. Whenever a function that tries to
   * compute averages is called, an exception will be thrown if one of
   * the points which you are using for the average lies outside the
   * periodicity box. The return points are garanteed to lie in the
   * perodicity box.
   */
  FlatManifold (const Point<spacedim> periodicity=Point<spacedim>());

  /**
   * Let the new point be the average sum of surrounding vertices.
   *
   * This particular implementation constructs the weighted average of
   * the surrounding points, and then calls internally the function
   * project_to_manifold. The reason why we do it this way, is to
   * allow lazy programmers to implement only the project_to_manifold
   * function for their own Manifold classes which are small (or
   * trivial) perturbations of a flat manifold. For most simple
   * geometries, it is possible to get reasonable results by deriving
   * your own Manifold class from FlatManifold, and write a new
   * interface only for the project_to_manifold function.
   */
    virtual Point<spacedim>
    get_new_point(const Quadrature<spacedim> &quad) const;
  
  
  /**
   *  Project to FlatManifold. This is the identity function for flat,
   *  Euclidean spaces. Note however that this function can be
   *  overloaded by derived classes, which will then benefit from the
   *  logic behind the get_new_point class which are often very
   *  similar (if not identical) to the one implemented in this class.
   */
  virtual
  Point<spacedim> project_to_manifold (const std::vector<Point<spacedim> > &points, 
				       const Point<spacedim> &candidate) const;

private:
  /**
   * The periodicity of this Manifold. Periodicity affects the way a
   * middle point is computed. It is assumed that if two points are
   * more than half period distant, then the distance should be
   * computed by crossing the periodicity boundary, i.e., the average
   * is computed by adding a full period to the sum of the two. For
   * example, if along direction 0 we have 2*pi periodicity, then the
   * average of (2*pi-eps) and (eps) is not pi, but 2*pi (or zero),
   * since, on a periodic manifold, these two points are at distance
   * 2*eps and not (2*pi-eps).
   *
   * A periodicity 0 along one direction means no periodicity. This is
   * the default value for all directions.
   */
  const Point<spacedim> periodicity;
};


/**
 *   A chart of dimension chartdim, which is part of a
 *   Manifold<dim,spacedim>.  This object specializes a Manifold of
 *   dimension chartdim embedded in a manifold of dimension spacedim,
 *   for which you have explicit pull_back and push_forward
 *   transformations. This object only makes sense when chartdim <=
 *   dim, and the constructor throws an exception if this is not the
 *   case. 
 *
 *   This is an helper class which is useful when you have an explicit
 *   map from an Euclidean space of dimension dim to an Euclidean
 *   space of dimension spacedim which represents your manifold, i.e.,
 *   when your manifold \f$\mathcal{M}\f$ can be represented by a map
 *   \f[
 *   F: \mathcal{B} \subset R^{\text{chartdim}} \mapsto \mathcal{M}
 *   \subset R^{\text{spacedim}}
 *   \f]
 *   (the push_forward() function)
 *   which admits the inverse transformation
 *   \f[
 *   F^{-1}: \mathcal{M}
 *   \subset R^{\text{spacedim}} \mapsto
 *   \mathcal{B} \subset R^{\text{chartdim}}
 *   \f]
 *   (the pull_back() function).
 *
 *   The get_new_point() function of the ManifoldChart class is
 *   implemented by calling the pull_back() method for all
 *   #surrounding_points, computing their weighted average in the
 *   chartdim Euclidean space, and calling the push_forward() method
 *   with the resulting point, i.e., \f[ p^{\text{new}} = F(\sum_i w_i
 *   F^{-1}(p_i)).  \f]
 *
 *   Derived classes are required to implement the push_forward() and
 *   the pull_back() methods.
 *
 *   @ingroup manifold
 *
 *   @author Luca Heltai, 2013
 */
template <int dim, int spacedim=dim, int chartdim=dim>
class ManifoldChart: public Manifold<dim,spacedim>
{
public:
  /**
   * Constructor. The optional argument can be used to specify the
   * periodicity of the chartdim-dimensional manifold (one period per
   * direction). A peridicity value of zero means that along that
   * direction there is no peridicity. By default no periodicity is
   * assumed.
   *
   * Periodicity affects the way a middle point is computed. It is
   * assumed that if two points are more than half period distant,
   * then the distance should be computed by crossing the periodicity
   * boundary, i.e., then the average is computed by adding a full
   * period to the sum of the two. For example, if along direction 0
   * we have 2*pi periodicity, then the average of (2*pi-eps) and
   * (eps) is not pi, but 2*pi (or zero), since, on the manifold,
   * these two points are at distance 2*eps and not (2*pi-eps)
   */
  ManifoldChart(const Point<chartdim> periodicity=Point<chartdim>());
  
  /**
   * Destructor. Does nothing here, but needs to be declared to make
   * it virtual.
   */
  virtual ~ManifoldChart ();


  /**
   * Refer to the general documentation of this class and the
   * documentation of the base class for more information.
   */
    virtual Point<spacedim>
    get_new_point(const Quadrature<spacedim> &quad) const;

    /**
     * Pull back the given point in spacedim to the Euclidean dim
     * dimensional space.
     *
     * Refer to the general documentation of this class for more
     * information.
     */
    virtual Point<chartdim>
    pull_back(const Point<spacedim> &space_point) const = 0;

    /**
     * Given a point in the dim dimensianal Euclidean space, this
     * method returns a point on the manifold embedded in the spacedim
     * Euclidean space.
     *
     * Refer to the general documentation of this class for more
     * information.
     */
    virtual Point<spacedim>
    push_forward(const Point<chartdim> &chart_point) const = 0;

  private:
  /**
   * The sub_manifold object is used to compute the average of the
   * points in the chart coordinates system. 
   */
  FlatManifold<dim,chartdim> sub_manifold;
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

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

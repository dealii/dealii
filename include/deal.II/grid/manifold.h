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


/**
 * We collect here some helper functions used in the
 * Manifold<dim,spacedim> classes.
 */
namespace Manifolds
{
  /**
    * Given a general mesh iterator, construct a quadrature with the
    * Laplace weights or with uniform weights according the parameter
    * @p with_laplace, and with all relevant points of the iterator:
    * vertices, line centers and/or face centers, which can be called
    * when creating new vertices in the manifold routines.
    */
  template <typename OBJECT>
  Quadrature<OBJECT::AccessorType::space_dimension>
  get_default_quadrature(const OBJECT &obj, bool with_laplace = false);
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
 *   weights (Note that unlike almost all other cases in the library,
 *   we here interpret the points in the quadrature object to be in
 *   real space, not on the reference cell.)
 *
 *   Internaly, the get_new_point() function calls the
 *   project_to_manifold() function after computing the weighted
 *   average of the quadrature points. This allows end users to only
 *   overload project_to_manifold() for simple situations.
 *
 *   Should a finer control be necessary, then get_new_point() can be
 *   overloaded.
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
   * project_to_manifold. User classes can get away by simply
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
   * Default constructor. The optional argument can be used to specify
   * the periodicity of the spacedim-dimensional manifold (one period
   * per direction). A periodicity value of zero means that along that
   * direction there is no periodicity. By default no periodicity is
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
   * expected. The third argument is used as a relative tolerance when
   * computing distances.
   *
   * Periodicity will be intended in the following way: the domain is
   * considered to be the box contained in [Point<spacedim>(),
   * periodicity) where the right extreme is excluded. If any of the
   * components of this box has zero length, then no periodicity is
   * computed in that direction. Whenever a function that tries to
   * compute averages is called, an exception will be thrown if one of
   * the points which you are using for the average lies outside the
   * periodicity box. The return points are garanteed to lie in the
   * perodicity box plus or minus tolerance*periodicity.norm().
   */
  FlatManifold (const Point<spacedim> periodicity=Point<spacedim>(),
                const double tolerance=1e-10);

  /**
   * Let the new point be the average sum of surrounding vertices.
   *
   * This particular implementation constructs the weighted average of
   * the surrounding points, and then calls internally the function
   * project_to_manifold. The reason why we do it this way, is to
   * allow lazy programmers to implement only the project_to_manifold
   * function for their own Manifold classes which are small (or
   * trivial) perturbations of a flat manifold. This is the case
   * whenever the coarse mesh is a decent approximation of the
   * manifold geometry. In this case, the middle point of a cell is
   * close to true middle point of the manifold, and a projection may
   * suffice.
   *
   * For most simple geometries, it is possible to get reasonable
   * results by deriving your own Manifold class from FlatManifold,
   * and write a new interface only for the project_to_manifold
   * function. You will have good approximations also with large
   * deformations, as long as in the coarsest mesh size you are trying
   * to refine, the middle point is not too far from the manifold mid
   * point, i.e., as long as the coarse mesh size is small enough.
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

  DeclException4(ExcPeriodicBox, int, Point<spacedim>, Point<spacedim>, double,
                 << "The component number " << arg1 << " of the point [ " << arg2
                 << " ]  is not in the interval [ " << -arg4
                 << ", " << arg3[arg4] << "), bailing out.");

  /**
   * Relative tolerance. This tolerance is used to compute distances
   * in double precision.
   */
  const double tolerance;
};


/**
 *   This class describes mappings that can be expressed in terms
 *   of charts. Specifically, this class with its template arguments
 *   describes a chart of dimension chartdim, which is part of a
 *   Manifold<dim,spacedim> and is used in an object of type
 *   Triangulation<dim,spacedim>:  It specializes a Manifold of
 *   dimension chartdim embedded in a manifold of dimension spacedim,
 *   for which you have explicit pull_back() and push_forward()
 *   transformations. Its use is explained in great detail in step-53.
 *
 *   This is a helper class which is useful when you have an explicit
 *   map from an Euclidean space of dimension chartdim to an Euclidean
 *   space of dimension spacedim which represents your manifold, i.e.,
 *   when your manifold $\mathcal{M}$ can be represented by a map
 *   \f[
 *   F: \mathcal{B} \subset R^{\text{chartdim}} \mapsto \mathcal{M}
 *   \subset R^{\text{spacedim}}
 *   \f]
 *   (the push_forward() function)
 *   and that admits the inverse transformation
 *   \f[
 *   F^{-1}: \mathcal{M}
 *   \subset R^{\text{spacedim}} \mapsto
 *   \mathcal{B} \subset R^{\text{chartdim}}
 *   \f]
 *   (the pull_back() function).
 *
 *   The get_new_point() function of the ChartManifold class is
 *   implemented by calling the pull_back() method for all
 *   #surrounding_points, computing their weighted average in the
 *   chartdim Euclidean space, and calling the push_forward() method
 *   with the resulting point, i.e., \f[ p^{\text{new}} = F(\sum_i w_i
 *   F^{-1}(p_i)).  \f]
 *
 *   Derived classes are required to implement the push_forward() and
 *   the pull_back() methods. All other functions required by mappings
 *   will then be provided by this class.
 *
 *   The dimension arguments #chartdim, #dim and #spacedim must
 *   satisfy the following relationships:
 *   @code
 *      dim <= spacedim
 *      chartdim <= spacedim
 *   @endcode
 *   However, there is no a priori relationship between #dim and
 *   #chartdim. For example, if you want to describe a mapping
 *   for an edge (a 1d object) in a 2d triangulation embedded in
 *   3d space, you could do so by parameterizing it via a line
 *   @f[
 *      F: [0,1] \rightarrow {\mathbb R}^3
 *   @f]
 *   in which case #chartdim is 1. On the other hand, there is
 *   no reason why one can't describe this as a mapping
 *   @f[
 *      F: {\mathbb R}^3 \rightarrow {\mathbb R}^3
 *   @f]
 *   in such a way that the line $[0,1]\times \{0\}\times \{0\}$ happens to be
 *   mapped onto the edge in question. Here, #chartdim is 3. This may seem
 *   cumbersome but satisfies the requirements of an invertible function $F$
 *   just fine as long as it is possible to get from the edge to the pull-back
 *   space and then back again. Finally, given that we are dealing with a 2d
 *   triangulation in 3d, one will often have a mapping from, say, the 2d unit
 *   square or unit disk to the domain in 3d space, and the edge in question
 *   may simply be the mapped edge of the unit domain in 2d space. In
 *   this case, #chartdim is 2.
 *
 *   @ingroup manifold
 *
 *   @author Luca Heltai, 2013, 2014
 */
template <int dim, int spacedim=dim, int chartdim=dim>
class ChartManifold: public Manifold<dim,spacedim>
{
public:
  /**
   * Constructor. The optional argument can be used to specify the
   * periodicity of the chartdim-dimensional manifold (one period per
   * direction). A periodicity value of zero means that along that
   * direction there is no periodicity. By default no periodicity is
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
  ChartManifold(const Point<chartdim> periodicity=Point<chartdim>());

  /**
   * Destructor. Does nothing here, but needs to be declared to make
   * it virtual.
   */
  virtual ~ChartManifold ();


  /**
   * Refer to the general documentation of this class and the
   * documentation of the base class for more information.
   */
  virtual Point<spacedim>
  get_new_point(const Quadrature<spacedim> &quad) const;

  /**
   * Pull back the given point in spacedim to the Euclidean chartdim
   * dimensional space.
   *
   * Refer to the general documentation of this class for more
   * information.
   */
  virtual Point<chartdim>
  pull_back(const Point<spacedim> &space_point) const = 0;

  /**
   * Given a point in the chartdim dimensional Euclidean space, this
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

  template <typename OBJECT>
  Quadrature<OBJECT::AccessorType::space_dimension>
  get_default_quadrature(const OBJECT &obj, bool with_laplace)
  {
    const int spacedim = OBJECT::AccessorType::space_dimension;
    const int dim = OBJECT::AccessorType::structure_dimension;

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
        sp[0] = obj->vertex(0);
        wp[0] = .5;
        sp[1] = obj->vertex(1);
        wp[1] = .5;
        break;
      case 2:
        sp.resize(8);
        wp.resize(8);

        for (unsigned int i=0; i<4; ++i)
          {
            sp[i] = obj->vertex(i);
            sp[4+i] = ( obj->line(i)->has_children() ?
                        obj->line(i)->child(0)->vertex(1) :
                        obj->line(i)->get_manifold().get_new_point_on_line(obj->line(i)) );
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
          = static_cast<TriaIterator<TriaAccessor<3, 3, 3> > >(obj);
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
        // Overwrited the weights with 1/np if we don't want to use
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

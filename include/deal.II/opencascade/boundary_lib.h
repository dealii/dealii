// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2015 by the deal.II authors
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


#ifndef dealii__occ_boundary_lib_h
#define dealii__occ_boundary_lib_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_OPENCASCADE

#include <deal.II/opencascade/utilities.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/manifold.h>

// opencascade needs "HAVE_CONFIG_H" to be exported...
#define HAVE_CONFIG_H
#include <BRepAdaptor_Curve.hxx>
#include <Adaptor3d_Curve.hxx>
#undef HAVE_CONFIG_H

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup OpenCASCADE
 * @{
 */

namespace OpenCASCADE
{
  /**
   * A Boundary object based on OpenCASCADE TopoDS_Shape where where new
   * points are first computed by averaging the surrounding points in the same
   * way as FlatManifold does, and are then projected in the normal direction
   * using OpenCASCADE utilities.
   *
   * This class makes no assumptions on the shape you pass to it, and the
   * topological dimension of the Manifold is inferred from the TopoDS_Shape
   * itself. In debug mode there is a sanity check to make sure that the
   * surrounding points (the ones used in project_to_manifold()) actually live
   * on the Manifold, i.e., calling OpenCASCADE::closest_point() on those
   * points leaves them untouched. If this is not the case, an
   * ExcPointNotOnManifold is thrown.
   *
   * This could happen, for example, if you are trying to use a shape of type
   * TopoDS_Edge when projecting on a face. In this case, the vertices of the
   * face would be collapsed to the edge, and your surrounding points would
   * not be lying on the given shape, raising an exception.
   *
   * @author Luca Heltai, Andrea Mola, 2011--2014.
   */
  template <int dim, int spacedim>
  class NormalProjectionBoundary : public Boundary<dim,spacedim>
  {
  public:

    /**
     * The standard constructor takes a generic TopoDS_Shape @p sh, and a
     * tolerance used to compute distances internally.
     *
     * The TopoDS_Shape can be arbitrary, i.e., a collection of shapes, faces,
     * edges or a single face or edge.
     */
    NormalProjectionBoundary(const TopoDS_Shape &sh,
                             const double tolerance=1e-7);

    /**
     * Perform the actual projection onto the manifold. This function, in
     * debug mode, checks that each of the @p surrounding_points is within
     * tolerance from the given TopoDS_Shape. If this is not the case, an
     * exception is thrown.
     *
     * The projected point is computed using OpenCASCADE normal projection
     * algorithms.
     */
    virtual Point<spacedim>
    project_to_manifold (const std::vector<Point<spacedim> > &surrounding_points,
                         const Point<spacedim> &candidate) const;


  private:
    /**
     * The topological shape which is used internally to project points. You
     * can construct such a shape by calling the OpenCASCADE::read_IGES()
     * function, which will create a TopoDS_Shape with the geometry contained
     * in the IGES file.
     */
    const TopoDS_Shape sh;

    /**
     * Relative tolerance used by this class to compute distances.
     */
    const double tolerance;
  };

  /**
   * A Boundary object based on OpenCASCADE TopoDS_Shape where new points are
   * first computed by averaging the surrounding points in the same way as
   * FlatManifold does, and then projecting them onto the manifold along the
   * direction specified at construction time using OpenCASCADE utilities.
   *
   * This class makes no assumptions on the shape you pass to it, and the
   * topological dimension of the Manifold is inferred from the TopoDS_Shape
   * itself. In debug mode there is a sanity check to make sure that the
   * surrounding points (the ones used in project_to_manifold()) actually live
   * on the Manifold, i.e., calling OpenCASCADE::closest_point() on those
   * points leaves them untouched. If this is not the case, an
   * ExcPointNotOnManifold is thrown.
   *
   * Notice that this type of Boundary descriptor may fail to give results if
   * the triangulation to be refined is close to the boundary of the given
   * TopoDS_Shape, or when the direction you use at construction time does not
   * intersect the shape. An exception is thrown when this happens.
   *
   * @author Luca Heltai, Andrea Mola, 2011--2014.
   */
  template <int dim, int spacedim>
  class DirectionalProjectionBoundary : public Boundary<dim,spacedim>
  {
  public:
    /**
     * Construct a Boundary object which will project points on the
     * TopoDS_Shape @p sh, along the given @p direction.
     */
    DirectionalProjectionBoundary(const TopoDS_Shape &sh,
                                  const Tensor<1,spacedim> &direction,
                                  const double tolerance=1e-7);

    /**
     * Perform the actual projection onto the manifold. This function, in
     * debug mode, checks that each of the @p surrounding_points is within
     * tolerance from the given TopoDS_Shape. If this is not the case, an
     * exception is thrown.
     *
     * The projected point is computed using OpenCASCADE directional
     * projection algorithms.
     */
    virtual Point<spacedim>
    project_to_manifold (const std::vector<Point<spacedim> > &surrounding_points,
                         const Point<spacedim> &candidate) const;

  private:
    /**
     * The topological shape which is used internally to project points. You
     * can construct such a shape by calling the OpenCASCADE::read_IGES()
     * function, which will create a TopoDS_Shape with the geometry contained
     * in the IGES file.
     */
    const TopoDS_Shape sh;

    /**
     * Direction used to project new points on the shape.
     */
    const Point<3> direction;

    /**
     * Relative tolerance used by this class to compute distances.
     */
    const double tolerance;
  };


  /**
   * A Boundary object based on OpenCASCADE TopoDS_Shape where new points are
   * first computed by averaging the surrounding points in the same way as
   * FlatManifold does, and then projecting them using OpenCASCADE utilities
   * onto the manifold along a direction which is an estimation of the
   * surrounding points (hence mesh cell) normal.
   *
   * The direction normal to the mesh is particularly useful because it is the
   * direction in which the mesh is missing nodes. For instance, during the
   * refinement of a cell a new node is initially created around the
   * baricenter of the cell. This location somehow ensures a uniform distance
   * from the nodes of the old cell. Projecting such cell baricenter onto the
   * CAD surface in the direction normal to the original cell will then retain
   * uniform distance from the points of the original cell. Of course, at the
   * stage of mesh generation, no dof handler nor finite element are defined,
   * and such direction has to be estimated. For the case in which 8
   * surrounding points are present, 4 different triangles are identified with
   * the points assigned, and the normals of such triangles are averaged to
   * obtain the approximation of the normal to the cell.
   *
   * The case in which 2 surrounding points are present (i.e.:a cell edge is
   * being refined) is of course more tricky. The average of the CAD surface
   * normals at the 2 surrounding points is first computed, and then projected
   * onto the plane normal to the segment linking the surrounding points. This
   * again is an attempt to have the new point with equal distance with
   * respect to the surrounding points
   *
   * This class only operates with CAD faces and makes the assumption that the
   * shape you pass to it contains at least one face. If that is not the case,
   * an Exception is thrown. In debug mode there is a sanity check to make
   * sure that the surrounding points (the ones used in project_to_manifold())
   * actually live on the Manifold, i.e., calling OpenCASCADE::closest_point()
   * on those points leaves them untouched. If this is not the case, an
   * ExcPointNotOnManifold is thrown.
   *
   *
   * Notice that this type of Boundary descriptor may fail to give results if
   * the triangulation to be refined is close to the boundary of the given
   * TopoDS_Shape, or when the normal direction estimated from the surrounding
   * points does not intersect the shape.  An exception is thrown when this
   * happens.
   *
   * @author Luca Heltai, Andrea Mola, 2011--2014.
   */
  template <int dim, int spacedim>
  class NormalToMeshProjectionBoundary : public Boundary<dim,spacedim>
  {
  public:
    /**
     * Construct a Boundary object which will project points on the
     * TopoDS_Shape @p sh, along a direction which is approximately normal to
     * the mesh cell.
     */
    NormalToMeshProjectionBoundary(const TopoDS_Shape &sh,
                                   const double tolerance=1e-7);

    /**
     * Perform the actual projection onto the manifold. This function, in
     * debug mode, checks that each of the @p surrounding_points is within
     * tolerance from the given TopoDS_Shape. If this is not the case, an
     * exception is thrown.
     */
    virtual Point<spacedim>
    project_to_manifold (const std::vector<Point<spacedim> > &surrounding_points,
                         const Point<spacedim> &candidate) const;

  private:
    /**
     * The topological shape which is used internally to project points. You
     * can construct such a shape by calling the OpenCASCADE::read_IGES()
     * function, which will create a TopoDS_Shape with the geometry contained
     * in the IGES file.
     */
    const TopoDS_Shape sh;

    /**
     * Direction used to project new points on the shape.
     */
    const Point<3> direction;

    /**
     * Relative tolerance used by this class to compute distances.
     */
    const double tolerance;
  };

  /**
   * A Boundary object based on OpenCASCADE TopoDS_Shape objects which have
   * topological dimension equal to one (TopoDS_Edge or TopoDS_Wire) where new
   * points are located at the arclength average of the surrounding points. If
   * the given TopoDS_Shape can be casted to a periodic (closed) curve, then
   * this information is used internally to set the periodicity of the base
   * ChartManifold class.
   *
   * This class can only work on TopoDS_Edge or TopoDS_Wire objects, and it
   * only makes sense when spacedim is three. If you use an object of
   * topological dimension different from one, an exception is throw.
   *
   * In debug mode there is an additional sanity check to make sure that the
   * surrounding points actually live on the Manifold, i.e., calling
   * OpenCASCADE::closest_point() on those points leaves them untouched. If
   * this is not the case, an ExcPointNotOnManifold is thrown.
   *
   * @author Luca Heltai, Andrea Mola, 2011--2014.
   */
  template <int dim, int spacedim>
  class ArclengthProjectionLineManifold : public  ChartManifold<dim,spacedim,1>
  {
  public:
    /**
     * Default constructor with a TopoDS_Edge.
     */
    ArclengthProjectionLineManifold(const TopoDS_Shape &sh,
                                    const double tolerance=1e-7);

    /**
     * Given a point on real space, find its arclength parameter. Throws an
     * error in debug mode, if the point is not on the TopoDS_Edge given at
     * construction time.
     */
    virtual Point<1>
    pull_back(const Point<spacedim> &space_point) const;

    /**
     * Given an arclength parameter, find its image in real space.
     */
    virtual Point<spacedim>
    push_forward(const Point<1> &chart_point) const;

  private:
    /**
     * A Curve adaptor. This is the one which is used in the computations, and
     * it points to the right one above.
     */
    Handle_Adaptor3d_HCurve curve;

    /**
     * Relative tolerance used in all internal computations.
     */
    const double tolerance;

    /**
     * The total length of the curve. This is also used as a period if the
     * edge is periodic.
     */
    const double length;
  };

  /**
   * Manifold description for the face of a CAD imported usign OpenCASCADE.
   *
   * @ingroup manifold
   *
   * @author Andrea Mola, Mauro Bardelloni, 2016
   */
  template <int dim, int spacedim>
  class NURBSPatchManifold : public ChartManifold<dim, spacedim, 2>
  {
  public:
    /**
     * The constructor takes an OpenCASCADE TopoDS_Face @p face and an optional
     * @p tolerance. This class uses the interval OpenCASCADE variables @var u,
     * @var v to descrive the manifold.
     */
    NURBSPatchManifold(const TopoDS_Face &face, const double tolerance = 1e-7);

    /**
     * Pull back the given point from the Euclidean space. Will return the uv
     * coordinates associated with the point @p space_point.
     */
    virtual Point<2>
    pull_back(const Point<spacedim> &space_point) const;

    /**
     * Given a @p chart_point in the uv coordinate system, this method returns the
     * Euclidean coordinates associated.
     */
    virtual Point<spacedim>
    push_forward(const Point<2> &chart_point) const;

    /**
     * Given a point in the spacedim dimensional Euclidean space, this
     * method returns the derivatives of the function $F$ that maps from
     * the uv coordinate system to the Euclidean coordinate
     * system. In other words, it is a matrix of size
     * $\text{spacedim}\times\text{chartdim}$.
     *
     * This function is used in the computations required by the
     * get_tangent_vector() function.
     *
     * Refer to the general documentation of this class for more information.
     */
    virtual
    DerivativeForm<1,2,spacedim>
    push_forward_gradient(const Point<2> &chart_point) const;

  private:
    /**
     * Return a tuple representing the minimum and maximum values of u
     * and v.  Precisely, it returns (u_min, u_max, v_min, v_max)
     */
    std_cxx11::tuple<double, double, double, double>
    get_uv_bounds() const;

    /**
     * An OpenCASCADE TopoDS_Face @p face given by the CAD.
     */
    TopoDS_Face face;

    /**
     * Tolerance used by OpenCASCADE to identify points in each
     * operation.
     */
    double tolerance;
  };

}

/*@}*/

DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_WITH_OPENCASCADE

/*------------------------------ occ_boundary_lib.h ------------------------------*/
#endif
/*------------------------------ occ_boundary_lib.h ------------------------------*/

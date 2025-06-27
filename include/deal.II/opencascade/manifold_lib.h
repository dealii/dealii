// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_occ_manifold_lib_h
#define dealii_occ_manifold_lib_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_OPENCASCADE

#  include <deal.II/grid/manifold.h>

#  include <deal.II/opencascade/utilities.h>

// opencascade needs "HAVE_CONFIG_H" to be exported...
#  define HAVE_CONFIG_H
#  include <Adaptor3d_Curve.hxx>
#  if !DEAL_II_OPENCASCADE_VERSION_GTE(7, 6, 0)
#    include <Adaptor3d_HCurve.hxx>
#  endif
#  include <BRepAdaptor_Curve.hxx>
#  undef HAVE_CONFIG_H

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup OpenCASCADE
 * @{
 */

namespace OpenCASCADE
{
  /**
   * A Manifold object based on OpenCASCADE TopoDS_Shape where new points are
   * first computed by averaging the surrounding points in the same way as
   * FlatManifold does, and are then projected in the normal direction using
   * OpenCASCADE utilities.
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
   */
  template <int dim, int spacedim>
  class NormalProjectionManifold : public FlatManifold<dim, spacedim>
  {
  public:
    /**
     * The standard constructor takes a generic TopoDS_Shape @p sh, and a
     * tolerance used to compute distances internally.
     *
     * The TopoDS_Shape can be arbitrary, i.e., a collection of shapes, faces,
     * edges or a single face or edge.
     */
    NormalProjectionManifold(const TopoDS_Shape &sh,
                             const double        tolerance = 1e-7);

    /**
     * Clone the current Manifold.
     */
    virtual std::unique_ptr<Manifold<dim, spacedim>>
    clone() const override;

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
    project_to_manifold(
      const ArrayView<const Point<spacedim>> &surrounding_points,
      const Point<spacedim>                  &candidate) const override;


  protected:
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
   * A Manifold object based on OpenCASCADE TopoDS_Shape where new points are
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
   * Notice that this type of Manifold descriptor may fail to give results if
   * the triangulation to be refined is close to the boundary of the given
   * TopoDS_Shape, or when the direction you use at construction time does not
   * intersect the shape. An exception is thrown when this happens.
   */
  template <int dim, int spacedim>
  class DirectionalProjectionManifold : public FlatManifold<dim, spacedim>
  {
  public:
    /**
     * Construct a Manifold object which will project points on the
     * TopoDS_Shape @p sh, along the given @p direction.
     */
    DirectionalProjectionManifold(const TopoDS_Shape        &sh,
                                  const Tensor<1, spacedim> &direction,
                                  const double               tolerance = 1e-7);

    /**
     * Clone the current Manifold.
     */
    virtual std::unique_ptr<Manifold<dim, spacedim>>
    clone() const override;

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
    project_to_manifold(
      const ArrayView<const Point<spacedim>> &surrounding_points,
      const Point<spacedim>                  &candidate) const override;

  protected:
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
    const Tensor<1, spacedim> direction;

    /**
     * Relative tolerance used by this class to compute distances.
     */
    const double tolerance;
  };


  /**
   * A Manifold object based on OpenCASCADE TopoDS_Shape where new points are
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
   * Notice that this type of Manifold descriptor may fail to give results if
   * the triangulation to be refined is close to the boundary of the given
   * TopoDS_Shape, or when the normal direction estimated from the surrounding
   * points does not intersect the shape.  An exception is thrown when this
   * happens.
   */
  template <int dim, int spacedim>
  class NormalToMeshProjectionManifold : public FlatManifold<dim, spacedim>
  {
  public:
    /**
     * Construct a Manifold object which will project points on the
     * TopoDS_Shape @p sh, along a direction which is approximately normal to
     * the mesh cell.
     */
    NormalToMeshProjectionManifold(const TopoDS_Shape &sh,
                                   const double        tolerance = 1e-7);

    /**
     * Clone the current Manifold.
     */
    virtual std::unique_ptr<Manifold<dim, spacedim>>
    clone() const override;

    /**
     * Perform the actual projection onto the manifold. This function, in
     * debug mode, checks that each of the @p surrounding_points is within
     * tolerance from the given TopoDS_Shape. If this is not the case, an
     * exception is thrown.
     */
    virtual Point<spacedim>
    project_to_manifold(
      const ArrayView<const Point<spacedim>> &surrounding_points,
      const Point<spacedim>                  &candidate) const override;

  protected:
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
   * A Manifold object based on OpenCASCADE TopoDS_Shape objects which have
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
   */
  template <int dim, int spacedim>
  class ArclengthProjectionLineManifold : public ChartManifold<dim, spacedim, 1>
  {
  public:
    /**
     * Default constructor with a TopoDS_Edge.
     */
    ArclengthProjectionLineManifold(const TopoDS_Shape &sh,
                                    const double        tolerance = 1e-7);

    /**
     * Clone the current Manifold.
     */
    virtual std::unique_ptr<Manifold<dim, spacedim>>
    clone() const override;

    /**
     * Given a point on real space, find its arclength parameter. Throws an
     * error in debug mode, if the point is not on the TopoDS_Edge given at
     * construction time.
     */
    virtual Point<1>
    pull_back(const Point<spacedim> &space_point) const override;

    /**
     * Given an arclength parameter, find its image in real space.
     */
    virtual Point<spacedim>
    push_forward(const Point<1> &chart_point) const override;

  protected:
    /**
     * The actual shape used to build this object.
     */
    const TopoDS_Shape sh;

    /**
     * A Curve adaptor. This is the one which is used in the computations, and
     * it points to the right one above.
     */
#  if DEAL_II_OPENCASCADE_VERSION_GTE(7, 6, 0)
    Handle_Adaptor3d_Curve curve;
#  else
    Handle_Adaptor3d_HCurve curve;
#  endif

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
   * Manifold description for the face of a CAD imported using OpenCASCADE.
   *
   * @ingroup manifold
   */
  template <int dim, int spacedim>
  class NURBSPatchManifold : public ChartManifold<dim, spacedim, 2>
  {
  public:
    /**
     * The constructor takes an OpenCASCADE TopoDS_Face @p face and an optional
     * @p tolerance. This class uses the interval OpenCASCADE variables u, v
     * to describe the manifold.
     */
    NURBSPatchManifold(const TopoDS_Face &face, const double tolerance = 1e-7);

    /**
     * Clone the current Manifold.
     */
    virtual std::unique_ptr<Manifold<dim, spacedim>>
    clone() const override;

    /**
     * Pull back the given point from the Euclidean space. Will return the uv
     * coordinates associated with the point @p space_point.
     */
    virtual Point<2>
    pull_back(const Point<spacedim> &space_point) const override;

    /**
     * Given a @p chart_point in the uv coordinate system, this method returns the
     * Euclidean coordinates associated.
     */
    virtual Point<spacedim>
    push_forward(const Point<2> &chart_point) const override;

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
    virtual DerivativeForm<1, 2, spacedim>
    push_forward_gradient(const Point<2> &chart_point) const override;

  protected:
    /**
     * Return a tuple representing the minimum and maximum values of u
     * and v.  Precisely, it returns (u_min, u_max, v_min, v_max)
     */
    std::tuple<double, double, double, double>
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

} // namespace OpenCASCADE

/** @} */

DEAL_II_NAMESPACE_CLOSE


#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_OPENCASCADE
#endif // dealii_occ_manifold_lib_h

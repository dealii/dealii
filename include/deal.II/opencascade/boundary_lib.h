// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2018 by the deal.II authors
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


#ifndef dealii_occ_boundary_lib_h
#  define dealii_occ_boundary_lib_h

#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_OPENCASCADE

#    include <deal.II/opencascade/manifold_lib.h>

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
   *
   * @author Luca Heltai, Andrea Mola, 2011--2014.
   *
   * @deprecated Use NormalProjectionManifold instead, which is identical to
   * this class but satisfies the modern Manifold-based naming convention.
   */
  template <int dim, int spacedim>
  class DEAL_II_DEPRECATED NormalProjectionBoundary
    : public NormalProjectionManifold<dim, spacedim>
  {
  public:
    /**
     * Inherit all constructors.
     */
    using NormalProjectionManifold<dim, spacedim>::NormalProjectionManifold;

    /**
     * Clone the current Manifold.
     */
    virtual std::unique_ptr<Manifold<dim, spacedim>>
    clone() const override;
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
   *
   * @author Luca Heltai, Andrea Mola, 2011--2014.
   *
   * @deprecated Use DirectionalProjectionManifold instead, which is identical
   * to this class but satisfies the modern Manifold-based naming convention.
   */
  template <int dim, int spacedim>
  class DEAL_II_DEPRECATED DirectionalProjectionBoundary
    : public DirectionalProjectionManifold<dim, spacedim>
  {
  public:
    /**
     * Inherit all constructors.
     */
    using DirectionalProjectionManifold<dim, spacedim>::
      DirectionalProjectionManifold;

    /**
     * Clone the current Manifold.
     */
    virtual std::unique_ptr<Manifold<dim, spacedim>>
    clone() const override;
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
   *
   * @author Luca Heltai, Andrea Mola, 2011--2014.
   *
   * @deprecated Use NormalToMeshProjectionManifold instead, which is identical
   * to this class but satisfies the modern Manifold-based naming convention.
   */
  template <int dim, int spacedim>
  class DEAL_II_DEPRECATED NormalToMeshProjectionBoundary
    : public NormalToMeshProjectionManifold<dim, spacedim>
  {
  public:
    /**
     * Inherit all constructors.
     */
    using NormalToMeshProjectionManifold<dim, spacedim>::
      NormalToMeshProjectionManifold;

    /**
     * Clone the current Manifold.
     */
    virtual std::unique_ptr<Manifold<dim, spacedim>>
    clone() const override;
  };
} // namespace OpenCASCADE

/*@}*/

DEAL_II_NAMESPACE_CLOSE


#  endif // DEAL_II_WITH_OPENCASCADE

#endif // dealii_occ_boundary_lib_h
/*---------------------------- occ_boundary_lib.h ---------------------------*/

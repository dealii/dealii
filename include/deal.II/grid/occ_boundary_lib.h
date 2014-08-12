// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 2014 by the deal.II authors
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


#ifndef __deal2__occ_boundary_lib_h
#define __deal2__occ_boundary_lib_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_OPENCASCADE

#include <deal.II/grid/occ_utilities.h>
#include <deal.II/grid/tria_boundary.h>

DEAL_II_NAMESPACE_OPEN

namespace OpenCASCADE 
{
  /**
   * @addtogroup OpenCASCADE
   * @{
   * 
   * A Boundary object based on OpenCASCADE TopoDS_Shape where new
   * points are first computed using the FlatManifold class, and then
   * projected in the normal direction using OpenCASCADE utilities.
   *
   * This class makes no assumptions on the shape you pass to it, and
   * the topological dimension of the Manifold is inferred from the
   * TopoDS_Shape itself. In debug mode there is a sanity check to
   * make sure that the surrounding points (the ones used in
   * project_to_manifold()) actually live on the Manifold, i.e.,
   * calling OpenCASCADE::closest_point() on those points leaves them
   * untouched. If this is not the case, an ExcPointNotOnManifold is
   * thrown.
   *
   * This could happen, for example, if you are trying to use a shape
   * of type TopoDS_Edge when projecting on a face. In this case, the
   * vertices of the face would be collapsed to the edge, and your
   * surrounding points would not be on lying on the given shape,
   * raising an exception.
   * 
   * @author Luca Heltai, Andrea Mola, 2011--2014.
   */
  template <int dim, int spacedim>
  class NormalProjectionBoundary : public Boundary<dim,spacedim> {
    public:
    NormalProjectionBoundary(const TopoDS_Shape sh, 
			     const double tolerance=1e-7);
    
    /**
     * Perform the actual projection onto the manifold. This function,
     * in debug mode, checks that each of the #surrounding_points is
     * within tolerance from the given TopoDS_Shape. If this is not
     * the case, an exception is thrown.
     *
     * The projected point is computed using OpenCASCADE normal
     * projection algorithms.
     */
    virtual Point<spacedim>
    project_to_manifold (const std::vector<Point<spacedim> > &surrounding_points,
			 const Point<spacedim> &candidate) const;


  private:
    /**
     * The topological shape which is used internally to project
     * points. You can construct one such a shape by calling the
     * OpenCASCADE::read_IGES() function, which will create a
     * TopoDS_Shape with the geometry contained in the IGES file.
     */ 
    const TopoDS_Shape sh;

    /**
     * Relative tolerance used by this class to compute distances.
     */
    const double tolerance;
  };

  /**
   * A Boundary object based on OpenCASCADE TopoDS_Shape where new
   * points are first computed using the FlatManifold class, and then
   * projected along the direction given at construction time, using
   * OpenCASCADE utilities.
   *
   * This class makes no assumptions on the shape you pass to it, and
   * the topological dimension of the Manifold is inferred from the
   * TopoDS_Shape itself. In debug mode there is a sanity check to
   * make sure that the surrounding points (the ones used in
   * project_to_manifold()) actually live on the Manifold, i.e.,
   * calling OpenCASCADE::closest_point() on those points leaves them
   * untouched. If this is not the case, an ExcPointNotOnManifold is
   * thrown.
   * 
   * @author Luca Heltai, Andrea Mola, 2011--2014.
   */
  template <int dim, int spacedim>
  class AxisProjectionBoundary : public Boundary<dim,spacedim> {
    public:
    AxisProjectionBoundary(const TopoDS_Shape sh, 
			   const Point<3> direction, 
			   const double tolerance=1e-7);
    
    /**
     * Perform the actual projection onto the manifold. This function,
     * in debug mode, checks that each of the #surrounding_points is
     * within tolerance from the given TopoDS_Shape. If this is not
     * the case, an exception is thrown.
     *
     * The projected point is computed using OpenCASCADE normal
     * projection algorithms.
     */
    virtual Point<spacedim>
    project_to_manifold (const std::vector<Point<spacedim> > &surrounding_points,
			 const Point<spacedim> &candidate) const;

  private:
    /**
     * The topological shape which is used internally to project
     * points. You can construct one such a shape by calling the
     * OpenCASCADE::read_IGES() function, which will create a
     * TopoDS_Shape with the geometry contained in the IGES file.
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

}

/*@}*/

DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_WITH_OPENCASCADE

/*------------------------------ occ_boundary_lib.h ------------------------------*/
#endif
/*------------------------------ occ_boundary_lib.h ------------------------------*/

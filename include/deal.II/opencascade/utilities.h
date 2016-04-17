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


#ifndef dealii__occ_utilities_h
#define dealii__occ_utilities_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_OPENCASCADE

#include <deal.II/grid/tria.h>
#include <deal.II/base/point.h>

#include <string>

// opencascade needs "HAVE_CONFIG_H" to be exported...
#define HAVE_CONFIG_H
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_CompSolid.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Wire.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <gp_Pnt.hxx>
#undef HAVE_CONFIG_H



DEAL_II_NAMESPACE_OPEN

/**
 * We collect in this namespace all utilities which operate on OpenCASCADE
 * entities. OpenCASCADE splits every object into a topological description
 * and a geometrical entity. The basic topological description is a
 * TopoDS_Shape. TopoDS_Shapes are light objects, and can be copied around.
 * The closest deal.II analog is a TriaIterator.
 *
 * The OpenCASCADE topology is designed with reference to the STEP standard
 * ISO-10303-42.  The structure is an oriented one-way graph, where parents
 * refer to their children, and there are no back references. Abstract
 * structure is implemented as C++ classes from the TopoDS package. A
 * TopoDS_Shape is manipulated by value and contains 3 fields: location,
 * orientation and a myTShape handle (of the TopoDS_TShape type). According to
 * OpenCASCADE documentation, myTShape and Location are used to share data
 * between various shapes to save memory. For example, an edge belonging to
 * two faces has equal Locations and myTShape fields but different
 * Orientations (Forward in context of one face and Reversed in one of the
 * other).
 *
 * Valid shapes include collection of other shapes, solids, faces, edges,
 * vertices, etc.
 *
 * Once a topological description is available, if a concrete geometrical
 * object can be created, the BRep classes allow one to extract the actual
 * geometrical information from a shape.
 *
 * This is done by inheriting abstract topology classes from the TopoDS
 * package by those implementing a boundary representation model (from the
 * BRep package). Only 3 types of topological objects have geometric
 * representations – vertex, edge, and face.
 *
 * Every TopoDS_Shape can be queried to figure out what type of shape it is,
 * and actual geometrical objects, like surfaces, curves or points, can be
 * extracted using BRepTools.
 *
 * In this namespace we provide readers and writers that read standard CAD
 * files, and return a TopoDS_Shape, or that write a CAD file, given a
 * TopoDS_Shape. Most of the functions in the OpenCASCADE namespace deal with
 * TopoDS_Shapes of one type or another, and provide interfaces to common
 * deal.II objects, like Triangulation, Manifold, and so on.
 *
 * Notice that these tools are only useful when spacedim is equal to three,
 * since OpenCASCADE only operates in three-dimensional mode.
 *
 * @author Luca Heltai, Andrea Mola, 2011--2014.
 */
namespace OpenCASCADE
{
  /**
   * Count the subobjects of a shape. This function is useful to gather
   * information about the TopoDS_Shape passed as argument. It returns the
   * number of faces, edges and vertices (the only topological entities
   * associated with actual geometries) which are contained in the given
   * shape.
   */
  std_cxx11::tuple<unsigned int, unsigned int, unsigned int>
  count_elements(const TopoDS_Shape &shape);

  /**
   * Read IGES files and translate their content into openCascade topological
   * entities. The option scale_factor is used to compensate for different
   * units being used in the IGES files and in the target application. The
   * standard unit for IGES files is millimiters. The return object is a
   * TopoDS_Shape which contains all objects from the file.
   */
  TopoDS_Shape read_IGES(const std::string &filename,
                         const double scale_factor=1e-3);

  /**
   * Write the given topological shape into an IGES file.
   */
  void write_IGES(const TopoDS_Shape &shape,
                  const std::string &filename);

  /**
   * Read STEP files and translate their content into openCascade topological
   * entities. The option scale_factor is used to compensate for different
   * units being used in the STEP files and in the target application. The
   * standard unit for STEP files is millimiters. The return object is a
   * TopoDS_Shape which contains all objects from the file.
   */
  TopoDS_Shape read_STEP(const std::string &filename,
                         const double scale_factor=1e-3);


  /**
   * Write the given topological shape into an STEP file.
   */
  void write_STEP(const TopoDS_Shape &shape,
                  const std::string &filename);

  /**
   * This function returns the tolerance associated with the shape. Each CAD
   * geometrical object is defined along with a tolerance, which indicates
   * possible inaccuracy of its placement. For instance, the tolerance of a
   * vertex indicates that it can be located in any point contained in a
   * sphere centered in the nominal position and having radius tol. While
   * carrying out an operation such as projecting a point onto a surface
   * (which will in turn have its tolerance) we must keep in mind that the
   * precision of the projection will be limited by the tolerance with which
   * the surface is built.  The tolerance is computed taking the maximum
   * tolerance among the subshapes composing the shape.
   */
  double get_shape_tolerance(const TopoDS_Shape &shape);

  /**
   * Perform the intersection of the given topological shape with the plane
   * $c_x x + c_y y + c_z z +c = 0$. The returned topological shape will
   * contain as few bsplines as possible. An exception is thrown if the
   * intersection produces an empty shape.
   */
  TopoDS_Shape  intersect_plane(const TopoDS_Shape &in_shape,
                                const double c_x,
                                const double c_y,
                                const double c_z,
                                const double c,
                                const double tolerance=1e-7);

  /**
   * Try to join all edges contained in the given TopoDS_Shape into a single
   * TopoDS_Edge, containing as few BSPlines as possible. If the input shape
   * contains faces, they will be ignored by this function. If the contained
   * edges cannot be joined into a single one, i.e., they form disconnected
   * curves, an exception will be thrown.
   */
  TopoDS_Edge join_edges(const TopoDS_Shape &in_shape,
                         const double tolerance=1e-7);

  /**
   * Creates a 3D smooth BSpline curve passing through the points in the
   * assigned vector, and store it in the returned TopoDS_Shape (which is of
   * type TopoDS_Edge). The points are reordered internally according to their
   * scalar product with the direction, if direction is different from zero,
   * otherwise they are used as passed. Notice that this function changes the
   * input points if required by the algorithm.
   *
   * This class is used to interpolate a BsplineCurve passing through an array
   * of points, with a C2 Continuity. If the optional parameter @p closed is
   * set to true, then the curve will be C2 at all points except the first
   * (where only C1 continuity will be given), and it will be a closed curve.
   *
   * The curve is garanteed to be at distance @p tolerance from the input
   * points. If the algorithm fails in generating such a curve, an exception
   * is thrown.
   */
  TopoDS_Edge interpolation_curve(std::vector<Point<3> >  &curve_points,
                                  const Tensor<1,3> &direction=Tensor<1,3>(),
                                  const bool closed=false,
                                  const double tolerance=1e-7);

  /**
   * Extract all subshapes from a TopoDS_Shape, and store the results into
   * standard containers. If the shape does not contain a certain type of
   * shape, the respective container will be empty.
   */
  void extract_geometrical_shapes(const TopoDS_Shape &shape,
                                  std::vector<TopoDS_Face> &faces,
                                  std::vector<TopoDS_Edge> &edges,
                                  std::vector<TopoDS_Vertex> &vertices);

  /**
   * Create a triangulation from a single face. This class extract the first u
   * and v parameter of the parametric surface making up this face, and
   * creates a Triangulation<2,3> containing a single coarse cell reflecting
   * this face. If the surface is not a trimmed surface, the vertices of this
   * cell will coincide with the TopoDS_Vertex vertices of the original
   * TopoDS_Face. This, however, is often not the case, and the user should be
   * careful on how this mesh is used.
   */
  void create_triangulation(const TopoDS_Face &face,
                            Triangulation<2,3> &tria);

  /**
   * Extract all compound shapes from a TopoDS_Shape, and store the results
   * into standard containers. If the shape does not contain a certain type of
   * compound, the respective container will be empty.
   */
  void extract_compound_shapes(const TopoDS_Shape &shape,
                               std::vector<TopoDS_Compound> &compounds,
                               std::vector<TopoDS_CompSolid> &compsolids,
                               std::vector<TopoDS_Solid> &solids,
                               std::vector<TopoDS_Shell> &shells,
                               std::vector<TopoDS_Wire> &wires);

  /**
   * Project the point @p origin on the topological shape given by @p
   * in_shape, and returns the projected point, the subshape which contains
   * the point and the parametric u and v coordinates of the point within the
   * resulting shape. If the shape is not elementary, all its subshapes are
   * iterated, faces first, then edges, and the returned shape is the closest
   * one to the point @p origin. If the returned shape is an edge, then only
   * the u coordinate is filled with sensible information, and the v
   * coordinate is set to zero.
   *
   * This function returns a tuple containing the projected point, the shape,
   * the u coordinate and the v coordinate (which is different from zero only
   * if the resulting shape is a face).
   */
  std_cxx11::tuple<Point<3>, TopoDS_Shape, double, double>
  project_point_and_pull_back(const TopoDS_Shape &in_shape,
                              const Point<3> &origin,
                              const double tolerance=1e-7);

  /**
   * Return the projection of the point @p origin on the topological shape
   * given by @p in_shape. If the shape is not elementary, all its subshapes
   * are iterated, faces first, then edges, and the returned point is the
   * closest one to the @p in_shape, regardless of its type.
   */
  Point<3> closest_point(const TopoDS_Shape &in_shape,
                         const Point<3> &origin,
                         const double tolerance=1e-7);

  /**
   * Given an elementary shape @p in_shape and the reference coordinates
   * within the shape, returns the corresponding point in real space. If the
   * shape is a TopoDS_Edge, the @p v coordinate is ignored. Only edges or
   * faces, as returned by the function project_point_and_pull_back(), can be
   * used as input to this function. If this is not the case, an Exception is
   * thrown.
   */
  Point<3> push_forward(const TopoDS_Shape &in_shape,
                        const double u,
                        const double v);


  /**
   * Given a TopoDS_Face @p face and the reference coordinates within this
   * face, returns the corresponding point in real space, the normal to the
   * surface at that point and the min and max curvatures as a tuple.
   */
  std_cxx11::tuple<Point<3>,  Tensor<1,3>, double, double>
  push_forward_and_differential_forms(const TopoDS_Face &face,
                                      const double u,
                                      const double v,
                                      const double tolerance=1e-7);


  /**
   * Get the closest point to the given topological shape, together with the
   * normal and the min and max curvatures at that point. If the shape is not
   * elementary, all its sub-faces (only the faces) are iterated, faces first,
   * and only the closest point is returned. This function will throw an
   * exception if the @p in_shape does not contain at least one face.
   */
  std_cxx11::tuple<Point<3>,  Tensor<1,3>, double, double>
  closest_point_and_differential_forms(const TopoDS_Shape &in_shape,
                                       const Point<3> &origin,
                                       const double tolerance=1e-7);


  /**
   * Intersect a line passing through the given @p origin point along @p
   * direction and the given topological shape. If there is more than one
   * intersection, it will return the closest one.
   *
   * The optional @p tolerance parameter is used to compute distances.
   */
  Point<3> line_intersection(const TopoDS_Shape &in_shape,
                             const Point<3> &origin,
                             const Tensor<1,3> &direction,
                             const double tolerance=1e-7);


  /**
   * Convert OpenCASCADE point into a Point<3>.
   */
  Point<3> point(const gp_Pnt &p);


  /**
   * Convert Point<3> into OpenCASCADE point.
   */
  gp_Pnt point(const Point<3> &p);


  /**
   * Sort two points according to their scalar product with direction. If the
   * norm of the direction is zero, then use lexicographical ordering. The
   * optional parameter is used as a relative tolerance when comparing
   * objects.
   */
  bool point_compare(const Point<3> &p1, const Point<3> &p2,
                     const Tensor<1,3> &direction=Tensor<1,3>(),
                     const double tolerance=1e-10);


  /**
   * Exception thrown when the point specified as argument does not lie
   * between @p tolerance from the given TopoDS_Shape.
   */
  DeclException1 (ExcPointNotOnManifold,
                  Point<3>,
                  <<"The point [ "<<arg1<<" ] is not on the manifold.");

  /**
   * Exception thrown when the point specified as argument cannot be projected
   * to the manifold.
   */
  DeclException1 (ExcProjectionFailed,
                  Point<3>,
                  <<"Projection of point [ "<< arg1
                  << " ] failed.");

  /**
   * Thrown when internal OpenCASCADE utilities fail to return the OK status.
   */
  DeclException1 (ExcOCCError,
                  IFSelect_ReturnStatus,
                  <<"An OpenCASCADE routine failed with return status "
                  <<arg1);

  /**
   * Trying to make curve operations on a degenerate edge.
   */
  DeclException0(ExcEdgeIsDegenerate);

  /**
   * Trying to make operations on the wrong type of shapes.
   */
  DeclException0(ExcUnsupportedShape);
}


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_OPENCASCADE

/*------------------------------ occ_utilities.h ------------------------------*/
#endif
/*------------------------------ occ_utilities.h ------------------------------*/

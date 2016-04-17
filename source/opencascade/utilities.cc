// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2016 by the deal.II authors
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

#include <deal.II/opencascade/utilities.h>

#ifdef DEAL_II_WITH_OPENCASCADE

#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/exceptions.h>

#include <boost/bind.hpp>

#include <cstdio>
#include <iostream>
#include <set>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS

#include <IGESControl_Controller.hxx>
#include <IGESControl_Reader.hxx>
#include <IGESControl_Writer.hxx>

#include <STEPControl_Controller.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>

#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopExp_Explorer.hxx>

#include <Handle_Standard_Transient.hxx>

#include <TColStd_SequenceOfTransient.hxx>
#include <TColStd_HSequenceOfTransient.hxx>
#include <TColgp_HArray1OfPnt.hxx>

#include <gp_Pnt.hxx>
#include <gp_Lin.hxx>
#include <gp_Vec.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <IntCurvesFace_ShapeIntersector.hxx>

#include <BRepTools.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_HCurve.hxx>
#include <BRepAdaptor_HCompCurve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepAlgo_Section.hxx>

#include <Geom_Plane.hxx>
#include <Geom_BoundedCurve.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <GeomConvert_CompCurveToBSplineCurve.hxx>
#include <GeomLProp_SLProps.hxx>

#include <GCPnts_AbscissaPoint.hxx>
#include <ShapeAnalysis_Surface.hxx>

DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <vector>
#include <algorithm>

DEAL_II_NAMESPACE_OPEN

namespace OpenCASCADE
{
  std_cxx11::tuple<unsigned int, unsigned int, unsigned int>
  count_elements(const TopoDS_Shape &shape)
  {
    TopExp_Explorer exp;
    unsigned int n_faces=0, n_edges=0, n_vertices=0;
    for (exp.Init(shape, TopAbs_FACE);
         exp.More(); exp.Next(), ++n_faces)
      {}
    for (exp.Init(shape, TopAbs_EDGE);
         exp.More(); exp.Next(), ++n_edges)
      {}
    for (exp.Init(shape, TopAbs_VERTEX);
         exp.More(); exp.Next(), ++n_vertices)
      {}
    return std_cxx11::tuple<unsigned int, unsigned int, unsigned int>(n_faces, n_edges, n_vertices);
  }

  void extract_geometrical_shapes(const TopoDS_Shape &shape,
                                  std::vector<TopoDS_Face> &faces,
                                  std::vector<TopoDS_Edge> &edges,
                                  std::vector<TopoDS_Vertex> &vertices)
  {
    faces.resize(0);
    edges.resize(0);
    vertices.resize(0);

    TopExp_Explorer exp;
    for (exp.Init(shape, TopAbs_FACE); exp.More(); exp.Next())
      {
        faces.push_back(TopoDS::Face(exp.Current()));
      }
    for (exp.Init(shape, TopAbs_EDGE); exp.More(); exp.Next())
      {
        edges.push_back(TopoDS::Edge(exp.Current()));
      }
    for (exp.Init(shape, TopAbs_VERTEX); exp.More(); exp.Next())
      {
        vertices.push_back(TopoDS::Vertex(exp.Current()));
      }
  }


  void extract_compound_shapes(const TopoDS_Shape &shape,
                               std::vector<TopoDS_Compound> &compounds,
                               std::vector<TopoDS_CompSolid> &compsolids,
                               std::vector<TopoDS_Solid> &solids,
                               std::vector<TopoDS_Shell> &shells,
                               std::vector<TopoDS_Wire> &wires)
  {
    compounds.resize(0);
    compsolids.resize(0);
    solids.resize(0);
    shells.resize(0);
    wires.resize(0);

    TopExp_Explorer exp;
    for (exp.Init(shape, TopAbs_COMPOUND); exp.More(); exp.Next())
      {
        compounds.push_back(TopoDS::Compound(exp.Current()));
      }
    for (exp.Init(shape, TopAbs_COMPSOLID); exp.More(); exp.Next())
      {
        compsolids.push_back(TopoDS::CompSolid(exp.Current()));
      }
    for (exp.Init(shape, TopAbs_SOLID); exp.More(); exp.Next())
      {
        solids.push_back(TopoDS::Solid(exp.Current()));
      }
    for (exp.Init(shape, TopAbs_SHELL); exp.More(); exp.Next())
      {
        shells.push_back(TopoDS::Shell(exp.Current()));
      }
    for (exp.Init(shape, TopAbs_WIRE); exp.More(); exp.Next())
      {
        wires.push_back(TopoDS::Wire(exp.Current()));
      }
  }

  gp_Pnt point(const Point<3> &p)
  {
    return gp_Pnt(p(0), p(1), p(2));
  }


  Point<3> point(const gp_Pnt &p)
  {
    return Point<3>(p.X(), p.Y(), p.Z());
  }

  bool point_compare(const Point<3> &p1, const Point<3> &p2,
                     const Tensor<1,3> &direction,
                     const double tolerance)
  {
    const double rel_tol=std::max(tolerance, std::max(p1.norm(), p2.norm())*tolerance);
    if (direction.norm() > 0.0)
      return (p1*direction < p2*direction-rel_tol);
    else
      for (int d=2; d>=0; --d)
        if (p1[d] < p2[d]-rel_tol)
          return true;
        else if (p2[d] < p1[d]-rel_tol)
          return false;

    // If we got here, for all d, none of the conditions above was
    // satisfied. The two points are equal up to tolerance
    return false;
  }


  TopoDS_Shape read_IGES(const std::string &filename,
                         const double scale_factor)
  {
    IGESControl_Reader reader;
    IFSelect_ReturnStatus stat;
    stat = reader.ReadFile(filename.c_str());
    AssertThrow(stat == IFSelect_RetDone,
                ExcMessage("Error in reading file!"));

    Standard_Boolean failsonly = Standard_False;
    IFSelect_PrintCount mode = IFSelect_ItemsByEntity;
    reader.PrintCheckLoad (failsonly, mode);

    Standard_Integer nRoots = reader.TransferRoots();
    //selects all IGES entities (including non visible ones) in the
    //file and puts them into a list called MyList,

    AssertThrow(nRoots > 0,
                ExcMessage("Read nothing from file."));

    // Handle IGES Scale here.
    gp_Pnt Origin;
    gp_Trsf scale;
    scale.SetScale (Origin, scale_factor);

    TopoDS_Shape sh = reader.OneShape();
    BRepBuilderAPI_Transform trans(sh, scale);

    return trans.Shape();   // this is the actual translation
  }

  void write_IGES(const TopoDS_Shape &shape,
                  const std::string &filename)
  {
    IGESControl_Controller::Init();
    IGESControl_Writer ICW ("MM", 0);
    Standard_Boolean ok = ICW.AddShape (shape);
    AssertThrow(ok, ExcMessage("Failed to add shape to IGES controller."));
    ICW.ComputeModel();
    Standard_Boolean OK = ICW.Write (filename.c_str());
    AssertThrow(OK, ExcMessage("Failed to write IGES file."));
  }

  TopoDS_Shape read_STEP(const std::string &filename,
                         const double scale_factor)
  {
    STEPControl_Reader reader;
    IFSelect_ReturnStatus stat;
    stat = reader.ReadFile(filename.c_str());
    AssertThrow(stat == IFSelect_RetDone,
                ExcMessage("Error in reading file!"));

    Standard_Boolean failsonly = Standard_False;
    IFSelect_PrintCount mode = IFSelect_ItemsByEntity;
    reader.PrintCheckLoad (failsonly, mode);

    Standard_Integer nRoots = reader.TransferRoots();
    //selects all IGES entities (including non visible ones) in the
    //file and puts them into a list called MyList,

    AssertThrow(nRoots > 0,
                ExcMessage("Read nothing from file."));

    // Handle STEP Scale here.
    gp_Pnt Origin;
    gp_Trsf scale;
    scale.SetScale (Origin, scale_factor);

    TopoDS_Shape sh = reader.OneShape();
    BRepBuilderAPI_Transform trans(sh, scale);

    return trans.Shape();   // this is the actual translation
  }

  void write_STEP(const TopoDS_Shape &shape,
                  const std::string &filename)
  {
    STEPControl_Controller::Init();
    STEPControl_Writer SCW;
    IFSelect_ReturnStatus status;
    status = SCW.Transfer(shape, STEPControl_AsIs);
    AssertThrow(status == IFSelect_RetDone, ExcMessage("Failed to add shape to STEP controller."));

    status = SCW.Write(filename.c_str());

    AssertThrow(status == IFSelect_RetDone, ExcMessage("Failed to write translated shape to STEP file."));
  }

  double get_shape_tolerance(const TopoDS_Shape &shape)
  {
    double tolerance = 0.0;

    std::vector<TopoDS_Face> faces;
    std::vector<TopoDS_Edge> edges;
    std::vector<TopoDS_Vertex> vertices;

    extract_geometrical_shapes(shape,
                               faces,
                               edges,
                               vertices);

    for (unsigned int i=0; i<vertices.size(); ++i)
      tolerance = fmax(tolerance,BRep_Tool::Tolerance(vertices[i]));

    for (unsigned int i=0; i<edges.size(); ++i)
      tolerance = fmax(tolerance,BRep_Tool::Tolerance(edges[i]));

    for (unsigned int i=0; i<faces.size(); ++i)
      tolerance = fmax(tolerance,BRep_Tool::Tolerance(faces[i]));


    return tolerance;
  }

  TopoDS_Shape intersect_plane(const TopoDS_Shape &in_shape,
                               const double c_x,
                               const double c_y,
                               const double c_z,
                               const double c,
                               const double /*tolerance*/)
  {
    Handle(Geom_Plane) plane = new Geom_Plane(c_x,c_y,c_z,c);
    BRepAlgo_Section section(in_shape, plane);
    TopoDS_Shape edges = section.Shape();
    return edges;
  }

  TopoDS_Edge join_edges(const TopoDS_Shape &in_shape,
                         const double tolerance)
  {
    TopoDS_Edge out_shape;
    TopoDS_Shape edges = in_shape;
    std::vector<Handle_Geom_BoundedCurve> intersections;
    TopLoc_Location L;
    Standard_Real First;
    Standard_Real Last;
    gp_Pnt PIn(0.0,0.0,0.0);
    gp_Pnt PFin(0.0,0.0,0.0);
    gp_Pnt PMid(0.0,0.0,0.0);
    TopExp_Explorer edgeExplorer(edges , TopAbs_EDGE);
    TopoDS_Edge edge;
    while (edgeExplorer.More())
      {
        edge = TopoDS::Edge(edgeExplorer.Current());
        Handle(Geom_Curve) curve = BRep_Tool::Curve(edge,L,First,Last);
        intersections.push_back(Handle(Geom_BoundedCurve)::DownCast(curve));
        edgeExplorer.Next();
      }

    // Now we build a single bspline out of all the geometrical
    // curves, in Lexycographical order
    unsigned int numIntersEdges = intersections.size();
    Assert(numIntersEdges>0, ExcMessage("No curves to process!"));

    GeomConvert_CompCurveToBSplineCurve convert_bspline(intersections[0]);

    bool check = false, one_added = true, one_failed=true;
    std::vector<bool> added(numIntersEdges, false);
    added[0] = true;
    while (one_added == true)
      {
        one_added = false;
        one_failed = false;
        for (unsigned int i=1; i<numIntersEdges; ++i)
          if (added[i] == false)
            {
              Handle(Geom_Curve) curve = intersections[i];
              Handle(Geom_BoundedCurve) bcurve = Handle(Geom_BoundedCurve)::DownCast(curve);
              check = convert_bspline.Add(bcurve,tolerance,0,1,0);
              if (check == false)  // If we failed, try again with the reversed curve
                {
                  curve->Reverse();
                  Handle(Geom_BoundedCurve) bcurve = Handle(Geom_BoundedCurve)::DownCast(curve);
                  check = convert_bspline.Add(bcurve,tolerance,0,1,0);
                }
              one_failed = one_failed || (check == false);
              one_added = one_added || (check == true);
              added[i] = check;
            }
      }

    Assert(one_failed == false,
           ExcMessage("Joining some of the Edges failed."));

    Handle(Geom_Curve) bspline = convert_bspline.BSplineCurve();

    out_shape = BRepBuilderAPI_MakeEdge(bspline);
    return out_shape;
  }


  Point<3> line_intersection(const TopoDS_Shape &in_shape,
                             const Point<3> &origin,
                             const Tensor<1,3> &direction,
                             const double tolerance)
  {
    // translating original Point<dim> to gp point

    gp_Pnt P0 = point(origin);
    gp_Ax1 gpaxis(P0, gp_Dir(direction[0], direction[1], direction[2]));
    gp_Lin line(gpaxis);

    // destination point
    gp_Pnt Pproj(0.0,0.0,0.0);

    // we prepare now the surface for the projection we get the whole
    // shape from the iges model
    IntCurvesFace_ShapeIntersector Inters;
    Inters.Load(in_shape,tolerance);

    // Keep in mind: PerformNearest sounds pretty but DOESN'T WORK!!!
    // The closest point must be found by hand
    Inters.Perform(line,-RealLast(),+RealLast());
    Assert(Inters.IsDone(), ExcMessage("Could not project point."));

    double minDistance = 1e7;
    double distance;
    Point<3> result;
    for (int i=0; i<Inters.NbPnt(); ++i)
      {
        distance = point(origin).Distance(Inters.Pnt(i+1));
        //cout<<"Point "<<i<<": "<<point(Inters.Pnt(i+1))<<"  distance: "<<distance<<endl;
        if (distance < minDistance)
          {
            minDistance = distance;
            result = point(Inters.Pnt(i+1));
          }
      }

    return result;
  }

  TopoDS_Edge interpolation_curve(std::vector<Point<3> > &curve_points,
                                  const Tensor<1,3> &direction,
                                  const bool closed,
                                  const double tolerance)
  {

    unsigned int n_vertices = curve_points.size();

    if (direction*direction > 0)
      {
        std::sort(curve_points.begin(), curve_points.end(),
                  boost::bind(&OpenCASCADE::point_compare, _1, _2, direction, tolerance));
      }

    // set up array of vertices
    Handle(TColgp_HArray1OfPnt) vertices = new TColgp_HArray1OfPnt(1,n_vertices);
    for (unsigned int vertex=0; vertex<n_vertices; ++vertex)
      {
        vertices->SetValue(vertex+1,point(curve_points[vertex]));
      }


    GeomAPI_Interpolate bspline_generator(vertices, closed, tolerance);
    bspline_generator.Perform();
    Assert( (bspline_generator.IsDone()), ExcMessage("Interpolated bspline generation failed"));

    Handle(Geom_BSplineCurve) bspline = bspline_generator.Curve();
    TopoDS_Edge out_shape = BRepBuilderAPI_MakeEdge(bspline);
    return out_shape;
  }

  std_cxx11::tuple<Point<3>, TopoDS_Shape, double, double>
  project_point_and_pull_back(const TopoDS_Shape &in_shape,
                              const Point<3> &origin,
                              const double tolerance)
  {
    TopExp_Explorer exp;
    gp_Pnt Pproj = point(origin);

    double minDistance = 1e7;
    gp_Pnt tmp_proj(0.0,0.0,0.0);

    unsigned int counter = 0;
    unsigned int face_counter = 0;

    TopoDS_Shape out_shape;
    double u=0;
    double v=0;

    for (exp.Init(in_shape, TopAbs_FACE); exp.More(); exp.Next())
      {
        TopoDS_Face face = TopoDS::Face(exp.Current());

        // the projection function needs a surface, so we obtain the
        // surface upon which the face is defined
        Handle(Geom_Surface) SurfToProj = BRep_Tool::Surface(face);

        ShapeAnalysis_Surface projector(SurfToProj);
        gp_Pnt2d proj_params = projector.ValueOfUV(point(origin), tolerance);

        SurfToProj->D0(proj_params.X(),proj_params.Y(),tmp_proj);

        double distance = point(tmp_proj).distance(origin);
        if (distance < minDistance)
          {
            minDistance = distance;
            Pproj = tmp_proj;
            out_shape = face;
            u=proj_params.X();
            v=proj_params.Y();
            ++counter;
          }
        ++face_counter;
      }

    // face counter tells us if the shape contained faces: if it does, there is no need
    // to loop on edges. Even if the closest point lies on the boundary of a parametric surface,
    // we need in fact to retain the face and both u and v, if we want to use this method to
    // retrieve the surface normal
    if (face_counter==0)
      for (exp.Init(in_shape, TopAbs_EDGE); exp.More(); exp.Next())
        {
          TopoDS_Edge edge = TopoDS::Edge(exp.Current());
          if (!BRep_Tool::Degenerated(edge))
            {
              TopLoc_Location L;
              Standard_Real First;
              Standard_Real Last;

              // the projection function needs a Curve, so we obtain the
              // curve upon which the edge is defined
              Handle(Geom_Curve) CurveToProj = BRep_Tool::Curve(edge,L,First,Last);

              GeomAPI_ProjectPointOnCurve Proj(point(origin),CurveToProj);
              unsigned int num_proj_points = Proj.NbPoints();
              if ((num_proj_points > 0) && (Proj.LowerDistance() < minDistance))
                {
                  minDistance = Proj.LowerDistance();
                  Pproj = Proj.NearestPoint();
                  out_shape = edge;
                  u=Proj.LowerDistanceParameter();
                  ++counter;
                }
            }
        }

    Assert(counter > 0, ExcMessage("Could not find projection points."));
    return std_cxx11::tuple<Point<3>, TopoDS_Shape, double, double>
           (point(Pproj),out_shape, u, v);
  }


  Point<3> closest_point(const TopoDS_Shape &in_shape,
                         const Point<3> &origin,
                         const double tolerance)
  {
    std_cxx11::tuple<Point<3>, TopoDS_Shape, double, double>
    ref = project_point_and_pull_back(in_shape, origin, tolerance);
    return std_cxx11::get<0>(ref);
  }

  std_cxx11::tuple<Point<3>,  Tensor<1,3>, double, double>
  closest_point_and_differential_forms(const TopoDS_Shape &in_shape,
                                       const Point<3> &origin,
                                       const double tolerance)

  {
    std_cxx11::tuple<Point<3>, TopoDS_Shape, double, double>
    shape_and_params = project_point_and_pull_back(in_shape,
                                                   origin,
                                                   tolerance);

    TopoDS_Shape &out_shape = std_cxx11::get<1>(shape_and_params);
    double &u = std_cxx11::get<2>(shape_and_params);
    double &v = std_cxx11::get<3>(shape_and_params);

    // just a check here: the number of faces in out_shape must be 1, otherwise
    // something is wrong
    std_cxx11::tuple<unsigned int, unsigned int, unsigned int> numbers =
      count_elements(out_shape);
    (void)numbers;

    Assert(std_cxx11::get<0>(numbers) > 0,
           ExcMessage("Could not find normal: the shape containing the closest point has 0 faces."));
    Assert(std_cxx11::get<0>(numbers) < 2,
           ExcMessage("Could not find normal: the shape containing the closest point has more than 1 face."));


    TopExp_Explorer exp;
    exp.Init(out_shape, TopAbs_FACE);
    TopoDS_Face face = TopoDS::Face(exp.Current());
    return push_forward_and_differential_forms(face, u, v, tolerance);
  }

  Point<3> push_forward(const TopoDS_Shape &in_shape,
                        const double u,
                        const double v)
  {
    switch (in_shape.ShapeType())
      {
      case TopAbs_FACE:
      {
        BRepAdaptor_Surface surf(TopoDS::Face(in_shape));
        return point(surf.Value(u,v));
      }
      case TopAbs_EDGE:
      {
        BRepAdaptor_Curve curve(TopoDS::Edge(in_shape));
        return point(curve.Value(u));
      }
      default:
        Assert(false, ExcUnsupportedShape());
      }
    return Point<3>();
  }

  std_cxx11::tuple<Point<3>,  Tensor<1,3>, double, double>
  push_forward_and_differential_forms(const TopoDS_Face &face,
                                      const double u,
                                      const double v,
                                      const double /*tolerance*/)
  {
    Handle(Geom_Surface) SurfToProj = BRep_Tool::Surface(face);
    GeomLProp_SLProps props(SurfToProj, u, v, 1, 1e-7);
    gp_Pnt Value = props.Value();
    Assert(props.IsNormalDefined(), ExcMessage("Normal is not well defined!"));
    gp_Dir Normal = props.Normal();
    Assert(props.IsCurvatureDefined(), ExcMessage("Curvature is not well defined!"));
    Standard_Real Min_Curvature = props.MinCurvature();
    Standard_Real Max_Curvature = props.MaxCurvature();
    Tensor<1,3> normal = Point<3>(Normal.X(),Normal.Y(),Normal.Z());

    // In the case your manifold changes from convex to concave or viceversa
    // the normal could jump from "inner" to "outer" normal.
    // However, you should be able to change the normal sense preserving
    // the manifold orientation:
    if (face.Orientation()==TopAbs_REVERSED)
      {
        normal *= -1;
        Min_Curvature *= -1;
        Max_Curvature *= -1;
      }

    return std_cxx11::tuple<Point<3>, Tensor<1,3>, double, double>(point(Value), normal, Min_Curvature, Max_Curvature);
  }



  void create_triangulation(const TopoDS_Face &face,
                            Triangulation<2,3> &tria)
  {
    BRepAdaptor_Surface surf(face);
    const double u0 = surf.FirstUParameter();
    const double u1 = surf.LastUParameter();
    const double v0 = surf.FirstVParameter();
    const double v1 = surf.LastVParameter();

    std::vector<CellData<2> > cells;
    std::vector<Point<3> > vertices;
    SubCellData t;

    vertices.push_back(point(surf.Value(u0,v0)));
    vertices.push_back(point(surf.Value(u1,v0)));
    vertices.push_back(point(surf.Value(u0,v1)));
    vertices.push_back(point(surf.Value(u1,v1)));

    CellData<2> cell;
    for (unsigned int i=0; i<4; ++i)
      cell.vertices[i] = i;

    cells.push_back(cell);
    tria.create_triangulation(vertices, cells, t);
  }

} // end namespace

DEAL_II_NAMESPACE_CLOSE

#endif

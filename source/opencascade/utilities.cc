// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/opencascade/utilities.h>

#ifdef DEAL_II_WITH_OPENCASCADE

#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/point.h>
#  include <deal.II/base/utilities.h>

#  include <IGESControl_Controller.hxx>
#  include <IGESControl_Reader.hxx>
#  include <IGESControl_Writer.hxx>
#  include <STEPControl_Controller.hxx>
#  include <STEPControl_Reader.hxx>
#  include <STEPControl_Writer.hxx>
#  include <Standard_Version.hxx>
#  include <TopExp_Explorer.hxx>
#  include <TopoDS.hxx>
#  include <TopoDS_Edge.hxx>
#  include <TopoDS_Face.hxx>
#  include <TopoDS_Shape.hxx>

#  include <cstdio>
#  include <iostream>
#  include <set>
#  if (OCC_VERSION_MAJOR < 7)
#    include <Handle_Standard_Transient.hxx>
#  else
#    include <Standard_Transient.hxx>
#  endif

#  include <BRepAdaptor_Curve.hxx>
#  include <BRepAdaptor_HCompCurve.hxx>
#  include <BRepAdaptor_HCurve.hxx>
#  include <BRepAdaptor_Surface.hxx>
#  include <BRepAlgo_Section.hxx>
#  include <BRepBndLib.hxx>
#  include <BRepBuilderAPI_MakeEdge.hxx>
#  include <BRepBuilderAPI_Sewing.hxx>
#  include <BRepBuilderAPI_Transform.hxx>
#  include <BRepMesh_IncrementalMesh.hxx>
#  include <BRepTools.hxx>
#  include <BRep_Builder.hxx>
#  include <GCPnts_AbscissaPoint.hxx>
#  include <GeomAPI_Interpolate.hxx>
#  include <GeomAPI_ProjectPointOnCurve.hxx>
#  include <GeomAPI_ProjectPointOnSurf.hxx>
#  include <GeomConvert_CompCurveToBSplineCurve.hxx>
#  include <GeomLProp_SLProps.hxx>
#  include <Geom_BoundedCurve.hxx>
#  include <Geom_Plane.hxx>
#  include <IntCurvesFace_ShapeIntersector.hxx>
#  include <Poly_Triangulation.hxx>
#  include <ShapeAnalysis_Surface.hxx>
#  include <StlAPI_Reader.hxx>
#  include <StlAPI_Writer.hxx>
#  include <TColStd_HSequenceOfTransient.hxx>
#  include <TColStd_SequenceOfTransient.hxx>
#  include <TColgp_HArray1OfPnt.hxx>
#  include <TopLoc_Location.hxx>
#  include <gp_Lin.hxx>
#  include <gp_Pnt.hxx>
#  include <gp_Vec.hxx>

#  include <algorithm>
#  include <vector>


DEAL_II_NAMESPACE_OPEN

namespace OpenCASCADE
{
  std::tuple<unsigned int, unsigned int, unsigned int>
  count_elements(const TopoDS_Shape &shape)
  {
    TopExp_Explorer exp;
    unsigned int    n_faces = 0, n_edges = 0, n_vertices = 0;
    for (exp.Init(shape, TopAbs_FACE); exp.More(); exp.Next(), ++n_faces)
      {
      }
    for (exp.Init(shape, TopAbs_EDGE); exp.More(); exp.Next(), ++n_edges)
      {
      }
    for (exp.Init(shape, TopAbs_VERTEX); exp.More(); exp.Next(), ++n_vertices)
      {
      }
    return std::tuple<unsigned int, unsigned int, unsigned int>(n_faces,
                                                                n_edges,
                                                                n_vertices);
  }

  void
  extract_geometrical_shapes(const TopoDS_Shape &        shape,
                             std::vector<TopoDS_Face> &  faces,
                             std::vector<TopoDS_Edge> &  edges,
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


  void
  extract_compound_shapes(const TopoDS_Shape &           shape,
                          std::vector<TopoDS_Compound> & compounds,
                          std::vector<TopoDS_CompSolid> &compsolids,
                          std::vector<TopoDS_Solid> &    solids,
                          std::vector<TopoDS_Shell> &    shells,
                          std::vector<TopoDS_Wire> &     wires)
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

  template <int spacedim>
  gp_Pnt
  point(const Point<spacedim> &p)
  {
    switch (spacedim)
      {
        case 1:
          return gp_Pnt(p[0], 0, 0);
        case 2:
          return gp_Pnt(p[0], p[1], 0);
        case 3:
          return gp_Pnt(p[0], p[1], p[2]);
      }
    AssertThrow(false, ExcNotImplemented());
    return {};
  }

  template <int spacedim>
  Point<spacedim>
  point(const gp_Pnt &p, const double tolerance)
  {
    (void)tolerance;
    switch (spacedim)
      {
        case 1:
          Assert(std::abs(p.Y()) <= tolerance,
                 ExcMessage(
                   "Cannot convert OpenCASCADE point to 1d if p.Y() != 0."));
          Assert(std::abs(p.Z()) <= tolerance,
                 ExcMessage(
                   "Cannot convert OpenCASCADE point to 1d if p.Z() != 0."));
          return Point<spacedim>(p.X());
        case 2:
          Assert(std::abs(p.Z()) <= tolerance,
                 ExcMessage(
                   "Cannot convert OpenCASCADE point to 2d if p.Z() != 0."));
          return Point<spacedim>(p.X(), p.Y());
        case 3:
          return Point<spacedim>(p.X(), p.Y(), p.Z());
      }
    AssertThrow(false, ExcNotImplemented());
    return {};
  }

  template <int dim>
  bool
  point_compare(const Point<dim> &    p1,
                const Point<dim> &    p2,
                const Tensor<1, dim> &direction,
                const double          tolerance)
  {
    const double rel_tol =
      std::max(tolerance, std::max(p1.norm(), p2.norm()) * tolerance);
    if (direction.norm() > 0.0)
      return (p1 * direction < p2 * direction - rel_tol);
    else
      for (int d = dim; d >= 0; --d)
        if (p1[d] < p2[d] - rel_tol)
          return true;
        else if (p2[d] < p1[d] - rel_tol)
          return false;

    // If we got here, for all d, none of the conditions above was
    // satisfied. The two points are equal up to tolerance
    return false;
  }


  TopoDS_Shape
  read_IGES(const std::string &filename, const double scale_factor)
  {
    IGESControl_Reader    reader;
    IFSelect_ReturnStatus stat;
    stat = reader.ReadFile(filename.c_str());
    AssertThrow(stat == IFSelect_RetDone, ExcMessage("Error in reading file!"));

    Standard_Boolean    failsonly = Standard_False;
    IFSelect_PrintCount mode      = IFSelect_ItemsByEntity;
    reader.PrintCheckLoad(failsonly, mode);

    Standard_Integer nRoots = reader.TransferRoots();
    // selects all IGES entities (including non visible ones) in the
    // file and puts them into a list called MyList,

    AssertThrow(nRoots > 0, ExcMessage("Read nothing from file."));

    // Handle IGES Scale here.
    gp_Pnt  Origin;
    gp_Trsf scale;
    scale.SetScale(Origin, scale_factor);

    TopoDS_Shape             sh = reader.OneShape();
    BRepBuilderAPI_Transform trans(sh, scale);

    return trans.Shape(); // this is the actual translation
  }

  void
  write_IGES(const TopoDS_Shape &shape, const std::string &filename)
  {
    IGESControl_Controller::Init();
    IGESControl_Writer ICW("MM", 0);
    Standard_Boolean   ok = ICW.AddShape(shape);
    AssertThrow(ok, ExcMessage("Failed to add shape to IGES controller."));
    ICW.ComputeModel();
    Standard_Boolean OK = ICW.Write(filename.c_str());
    AssertThrow(OK, ExcMessage("Failed to write IGES file."));
  }


  TopoDS_Shape
  read_STL(const std::string &filename)
  {
    StlAPI_Reader reader;
    TopoDS_Shape  shape;
    reader.Read(shape, filename.c_str());
    return shape;
  }


  void
  write_STL(const TopoDS_Shape &shape,
            const std::string & filename,
            const double        deflection,
            const bool          sew_different_faces,
            const double        sewer_tolerance,
            const bool          is_relative,
            const double        angular_deflection,
            const bool          in_parallel)
  {
    TopLoc_Location            Loc;
    std::vector<TopoDS_Vertex> vertices;
    std::vector<TopoDS_Edge>   edges;
    std::vector<TopoDS_Face>   faces;
    OpenCASCADE::extract_geometrical_shapes(shape, faces, edges, vertices);
    const bool mesh_is_present =
      std::none_of(faces.begin(), faces.end(), [&Loc](const TopoDS_Face &face) {
        Handle(Poly_Triangulation) theTriangulation =
          BRep_Tool::Triangulation(face, Loc);
        return theTriangulation.IsNull();
      });
    TopoDS_Shape shape_to_be_written = shape;
    if (!mesh_is_present)
      {
        if (sew_different_faces)
          {
            BRepBuilderAPI_Sewing sewer(sewer_tolerance);
            sewer.Add(shape_to_be_written);
            sewer.Perform();
            shape_to_be_written = sewer.SewedShape();
          }
        else
          shape_to_be_written = shape;
        // BRepMesh_IncrementalMesh automatically calls the perform method to
        // create the triangulation which is stored in the argument
        // `shape_to_be_written`.
        BRepMesh_IncrementalMesh mesh_im(shape_to_be_written,
                                         deflection,
                                         is_relative,
                                         angular_deflection,
                                         in_parallel);
      }

    StlAPI_Writer writer;

#  if ((OCC_VERSION_MAJOR * 100 + OCC_VERSION_MINOR * 10) >= 690)
    // opencascade versions 6.9.0 onwards return an error status
    const auto error = writer.Write(shape_to_be_written, filename.c_str());

    // which is a custom type between 6.9.0 and 7.1.0
#    if ((OCC_VERSION_MAJOR * 100 + OCC_VERSION_MINOR * 10) < 720)
    AssertThrow(error == StlAPI_StatusOK,
                ExcMessage("Error writing STL from shape."));
#    else
    // and a boolean from version 7.2.0 onwards
    AssertThrow(error == true, ExcMessage("Error writing STL from shape."));
#    endif

#  else

    // for opencascade versions 6.8.0 and older the return value is void
    writer.Write(shape_to_be_written, filename.c_str());
#  endif
  }

  TopoDS_Shape
  read_STEP(const std::string &filename, const double scale_factor)
  {
    STEPControl_Reader    reader;
    IFSelect_ReturnStatus stat;
    stat = reader.ReadFile(filename.c_str());
    AssertThrow(stat == IFSelect_RetDone, ExcMessage("Error in reading file!"));

    Standard_Boolean    failsonly = Standard_False;
    IFSelect_PrintCount mode      = IFSelect_ItemsByEntity;
    reader.PrintCheckLoad(failsonly, mode);

    Standard_Integer nRoots = reader.TransferRoots();
    // selects all IGES entities (including non visible ones) in the
    // file and puts them into a list called MyList,

    AssertThrow(nRoots > 0, ExcMessage("Read nothing from file."));

    // Handle STEP Scale here.
    gp_Pnt  Origin;
    gp_Trsf scale;
    scale.SetScale(Origin, scale_factor);

    TopoDS_Shape             sh = reader.OneShape();
    BRepBuilderAPI_Transform trans(sh, scale);

    return trans.Shape(); // this is the actual translation
  }

  void
  write_STEP(const TopoDS_Shape &shape, const std::string &filename)
  {
    STEPControl_Controller::Init();
    STEPControl_Writer    SCW;
    IFSelect_ReturnStatus status;
    status = SCW.Transfer(shape, STEPControl_AsIs);
    AssertThrow(status == IFSelect_RetDone,
                ExcMessage("Failed to add shape to STEP controller."));

    status = SCW.Write(filename.c_str());

    AssertThrow(status == IFSelect_RetDone,
                ExcMessage("Failed to write translated shape to STEP file."));
  }

  double
  get_shape_tolerance(const TopoDS_Shape &shape)
  {
    double tolerance = 0.0;

    std::vector<TopoDS_Face>   faces;
    std::vector<TopoDS_Edge>   edges;
    std::vector<TopoDS_Vertex> vertices;

    extract_geometrical_shapes(shape, faces, edges, vertices);

    for (const auto &vertex : vertices)
      tolerance = std::fmax(tolerance, BRep_Tool::Tolerance(vertex));

    for (const auto &edge : edges)
      tolerance = std::fmax(tolerance, BRep_Tool::Tolerance(edge));

    for (const auto &face : faces)
      tolerance = std::fmax(tolerance, BRep_Tool::Tolerance(face));


    return tolerance;
  }

  TopoDS_Shape
  intersect_plane(const TopoDS_Shape &in_shape,
                  const double        c_x,
                  const double        c_y,
                  const double        c_z,
                  const double        c,
                  const double /*tolerance*/)
  {
    Handle(Geom_Plane) plane = new Geom_Plane(c_x, c_y, c_z, c);
    BRepAlgo_Section section(in_shape, plane);
    TopoDS_Shape     edges = section.Shape();
    return edges;
  }

  TopoDS_Edge
  join_edges(const TopoDS_Shape &in_shape, const double tolerance)
  {
    TopoDS_Edge                           out_shape;
    const TopoDS_Shape &                  edges = in_shape;
    std::vector<Handle_Geom_BoundedCurve> intersections;
    TopLoc_Location                       L;
    Standard_Real                         First;
    Standard_Real                         Last;
    gp_Pnt                                PIn(0.0, 0.0, 0.0);
    gp_Pnt                                PFin(0.0, 0.0, 0.0);
    gp_Pnt                                PMid(0.0, 0.0, 0.0);
    TopExp_Explorer                       edgeExplorer(edges, TopAbs_EDGE);
    TopoDS_Edge                           edge;
    while (edgeExplorer.More())
      {
        edge                     = TopoDS::Edge(edgeExplorer.Current());
        Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, L, First, Last);
        intersections.push_back(Handle(Geom_BoundedCurve)::DownCast(curve));
        edgeExplorer.Next();
      }

    // Now we build a single bspline out of all the geometrical
    // curves, in Lexycographical order
    unsigned int numIntersEdges = intersections.size();
    Assert(numIntersEdges > 0, ExcMessage("No curves to process!"));

    GeomConvert_CompCurveToBSplineCurve convert_bspline(intersections[0]);

    bool              check = false, one_added = true, one_failed = true;
    std::vector<bool> added(numIntersEdges, false);
    added[0] = true;
    while (one_added == true)
      {
        one_added  = false;
        one_failed = false;
        for (unsigned int i = 1; i < numIntersEdges; ++i)
          if (added[i] == false)
            {
              Handle(Geom_Curve) curve = intersections[i];
              Handle(Geom_BoundedCurve) bcurve =
                Handle(Geom_BoundedCurve)::DownCast(curve);
              check = convert_bspline.Add(bcurve, tolerance, false, true, 0);
              if (check ==
                  false) // If we failed, try again with the reversed curve
                {
                  curve->Reverse();
                  Handle(Geom_BoundedCurve) bcurve =
                    Handle(Geom_BoundedCurve)::DownCast(curve);
                  check =
                    convert_bspline.Add(bcurve, tolerance, false, true, 0);
                }
              one_failed = one_failed || (check == false);
              one_added  = one_added || (check == true);
              added[i]   = check;
            }
      }

    Assert(one_failed == false,
           ExcMessage("Joining some of the Edges failed."));

    Handle(Geom_Curve) bspline = convert_bspline.BSplineCurve();

    out_shape = BRepBuilderAPI_MakeEdge(bspline);
    return out_shape;
  }

  template <int dim>
  Point<dim>
  line_intersection(const TopoDS_Shape &  in_shape,
                    const Point<dim> &    origin,
                    const Tensor<1, dim> &direction,
                    const double          tolerance)
  {
    // translating original Point<dim> to gp point

    gp_Pnt P0 = point(origin);
    gp_Ax1 gpaxis(P0,
                  gp_Dir(direction[0],
                         dim > 1 ? direction[1] : 0,
                         dim > 2 ? direction[2] : 0));
    gp_Lin line(gpaxis);

    // destination point
    gp_Pnt Pproj(0.0, 0.0, 0.0);

    // we prepare now the surface for the projection we get the whole
    // shape from the iges model
    IntCurvesFace_ShapeIntersector Inters;
    Inters.Load(in_shape, tolerance);

    // Keep in mind: PerformNearest sounds pretty but DOESN'T WORK!!!
    // The closest point must be found by hand
    Inters.Perform(line, -RealLast(), +RealLast());
    Assert(Inters.IsDone(), ExcMessage("Could not project point."));

    double     minDistance = 1e7;
    Point<dim> result;
    for (int i = 0; i < Inters.NbPnt(); ++i)
      {
        const double distance = point(origin).Distance(Inters.Pnt(i + 1));
        // cout<<"Point "<<i<<": "<<point(Inters.Pnt(i+1))<<"  distance:
        // "<<distance<<endl;
        if (distance < minDistance)
          {
            minDistance = distance;
            result      = point<dim>(Inters.Pnt(i + 1));
          }
      }

    return result;
  }

  template <int dim>
  TopoDS_Edge
  interpolation_curve(std::vector<Point<dim>> &curve_points,
                      const Tensor<1, dim> &   direction,
                      const bool               closed,
                      const double             tolerance)
  {
    unsigned int n_vertices = curve_points.size();

    if (direction * direction > 0)
      {
        std::sort(curve_points.begin(),
                  curve_points.end(),
                  [&](const Point<dim> &p1, const Point<dim> &p2) {
                    return OpenCASCADE::point_compare(p1,
                                                      p2,
                                                      direction,
                                                      tolerance);
                  });
      }

    // set up array of vertices
    Handle(TColgp_HArray1OfPnt) vertices =
      new TColgp_HArray1OfPnt(1, n_vertices);
    for (unsigned int vertex = 0; vertex < n_vertices; ++vertex)
      {
        vertices->SetValue(vertex + 1, point(curve_points[vertex]));
      }


    GeomAPI_Interpolate bspline_generator(vertices, closed, tolerance);
    bspline_generator.Perform();
    Assert((bspline_generator.IsDone()),
           ExcMessage("Interpolated bspline generation failed"));

    Handle(Geom_BSplineCurve) bspline = bspline_generator.Curve();
    TopoDS_Edge out_shape             = BRepBuilderAPI_MakeEdge(bspline);
    out_shape.Closed(closed);
    return out_shape;
  }



  template <int spacedim>
  std::vector<TopoDS_Edge>
  create_curves_from_triangulation_boundary(
    const Triangulation<2, spacedim> &triangulation,
    const Mapping<2, spacedim> &      mapping)

  {
    // store maps from global vertex index to pairs of global face   indices
    //        and from global face   index to pairs of global vertex indices
    std::map<unsigned int, std::pair<unsigned int, unsigned int>> vert_to_faces;
    std::map<unsigned int, std::pair<unsigned int, unsigned int>> face_to_verts;
    std::map<unsigned int, bool>                                  visited_faces;
    std::map<unsigned int, Point<spacedim>>                       vert_to_point;

    unsigned int face_index;

    for (const auto &cell : triangulation.active_cell_iterators())
      for (unsigned int f : GeometryInfo<2>::face_indices())
        if (cell->face(f)->at_boundary())
          {
            // get global face and vertex indices
            face_index                       = cell->face(f)->index();
            const unsigned int v0            = cell->face(f)->vertex_index(0);
            const unsigned int v1            = cell->face(f)->vertex_index(1);
            face_to_verts[face_index].first  = v0;
            face_to_verts[face_index].second = v1;
            visited_faces[face_index]        = false;

            // extract mapped vertex locations
            std::array<Point<spacedim>, GeometryInfo<2>::vertices_per_cell>
              verts           = mapping.get_vertices(cell);
            vert_to_point[v0] = verts[GeometryInfo<2>::face_to_cell_vertices(
              f, 0, true, false, false)];
            vert_to_point[v1] = verts[GeometryInfo<2>::face_to_cell_vertices(
              f, 1, true, false, false)];

            // distribute indices into maps
            if (vert_to_faces.find(v0) == vert_to_faces.end())
              {
                vert_to_faces[v0].first = face_index;
              }
            else
              {
                vert_to_faces[v0].second = face_index;
              }
            if (vert_to_faces.find(v1) == vert_to_faces.end())
              {
                vert_to_faces[v1].first = face_index;
              }
            else
              {
                vert_to_faces[v1].second = face_index;
              }
          }

    // run through maps in an orderly fashion, i.e., through the
    // boundary in one cycle and add points to pointlist.
    std::vector<TopoDS_Edge> interpolation_curves;
    bool                     finished = (face_to_verts.size() == 0);
    face_index = finished ? 0 : face_to_verts.begin()->first;

    while (finished == false)
      {
        const unsigned int start_point_index = face_to_verts[face_index].first;
        unsigned int       point_index       = start_point_index;

        // point_index and face_index always run together
        std::vector<Point<spacedim>> pointlist;
        do
          {
            visited_faces[face_index] = true;
            auto current_point        = vert_to_point[point_index];
            pointlist.push_back(current_point);

            // Get next point
            if (face_to_verts[face_index].first != point_index)
              point_index = face_to_verts[face_index].first;
            else
              point_index = face_to_verts[face_index].second;

            // Get next face
            if (vert_to_faces[point_index].first != face_index)
              face_index = vert_to_faces[point_index].first;
            else
              face_index = vert_to_faces[point_index].second;
          }
        while (point_index != start_point_index);

        interpolation_curves.push_back(
          interpolation_curve(pointlist, Tensor<1, spacedim>(), true));

        finished = true;
        for (const auto &f : visited_faces)
          if (f.second == false)
            {
              face_index = f.first;
              finished   = false;
              break;
            }
      }
    return interpolation_curves;
  }


  template <int dim>
  std::tuple<Point<dim>, TopoDS_Shape, double, double>
  project_point_and_pull_back(const TopoDS_Shape &in_shape,
                              const Point<dim> &  origin,
                              const double        tolerance)
  {
    TopExp_Explorer exp;
    gp_Pnt          Pproj = point(origin);

    double minDistance = 1e7;
    gp_Pnt tmp_proj(0.0, 0.0, 0.0);

    unsigned int counter      = 0;
    unsigned int face_counter = 0;

    TopoDS_Shape out_shape;
    double       u = 0;
    double       v = 0;

    for (exp.Init(in_shape, TopAbs_FACE); exp.More(); exp.Next())
      {
        TopoDS_Face face = TopoDS::Face(exp.Current());

        // the projection function needs a surface, so we obtain the
        // surface upon which the face is defined
        Handle(Geom_Surface) SurfToProj = BRep_Tool::Surface(face);

        ShapeAnalysis_Surface projector(SurfToProj);
        gp_Pnt2d proj_params = projector.ValueOfUV(point(origin), tolerance);

        SurfToProj->D0(proj_params.X(), proj_params.Y(), tmp_proj);

        double distance = point<dim>(tmp_proj).distance(origin);
        if (distance < minDistance)
          {
            minDistance = distance;
            Pproj       = tmp_proj;
            out_shape   = face;
            u           = proj_params.X();
            v           = proj_params.Y();
            ++counter;
          }
        ++face_counter;
      }

    // face counter tells us if the shape contained faces: if it does, there is
    // no need to loop on edges. Even if the closest point lies on the boundary
    // of a parametric surface, we need in fact to retain the face and both u
    // and v, if we want to use this method to retrieve the surface normal
    if (face_counter == 0)
      for (exp.Init(in_shape, TopAbs_EDGE); exp.More(); exp.Next())
        {
          TopoDS_Edge edge = TopoDS::Edge(exp.Current());
          if (!BRep_Tool::Degenerated(edge))
            {
              TopLoc_Location L;
              Standard_Real   First;
              Standard_Real   Last;

              // the projection function needs a Curve, so we obtain the
              // curve upon which the edge is defined
              Handle(Geom_Curve) CurveToProj =
                BRep_Tool::Curve(edge, L, First, Last);

              GeomAPI_ProjectPointOnCurve Proj(point(origin), CurveToProj);
              unsigned int                num_proj_points = Proj.NbPoints();
              if ((num_proj_points > 0) && (Proj.LowerDistance() < minDistance))
                {
                  minDistance = Proj.LowerDistance();
                  Pproj       = Proj.NearestPoint();
                  out_shape   = edge;
                  u           = Proj.LowerDistanceParameter();
                  ++counter;
                }
            }
        }

    Assert(counter > 0, ExcMessage("Could not find projection points."));
    return std::tuple<Point<dim>, TopoDS_Shape, double, double>(
      point<dim>(Pproj), out_shape, u, v);
  }


  template <int dim>
  Point<dim>
  closest_point(const TopoDS_Shape &in_shape,
                const Point<dim> &  origin,
                const double        tolerance)
  {
    std::tuple<Point<dim>, TopoDS_Shape, double, double> ref =
      project_point_and_pull_back(in_shape, origin, tolerance);
    return std::get<0>(ref);
  }

  std::tuple<Point<3>, Tensor<1, 3>, double, double>
  closest_point_and_differential_forms(const TopoDS_Shape &in_shape,
                                       const Point<3> &    origin,
                                       const double        tolerance)

  {
    std::tuple<Point<3>, TopoDS_Shape, double, double> shape_and_params =
      project_point_and_pull_back(in_shape, origin, tolerance);

    TopoDS_Shape &out_shape = std::get<1>(shape_and_params);
    double &      u         = std::get<2>(shape_and_params);
    double &      v         = std::get<3>(shape_and_params);

    // just a check here: the number of faces in out_shape must be 1, otherwise
    // something is wrong
    std::tuple<unsigned int, unsigned int, unsigned int> numbers =
      count_elements(out_shape);
    (void)numbers;

    Assert(
      std::get<0>(numbers) > 0,
      ExcMessage(
        "Could not find normal: the shape containing the closest point has 0 faces."));
    Assert(
      std::get<0>(numbers) < 2,
      ExcMessage(
        "Could not find normal: the shape containing the closest point has more than 1 face."));


    TopExp_Explorer exp;
    exp.Init(out_shape, TopAbs_FACE);
    TopoDS_Face face = TopoDS::Face(exp.Current());
    return push_forward_and_differential_forms(face, u, v, tolerance);
  }

  template <int dim>
  Point<dim>
  push_forward(const TopoDS_Shape &in_shape, const double u, const double v)
  {
    switch (in_shape.ShapeType())
      {
        case TopAbs_FACE:
          {
            BRepAdaptor_Surface surf(TopoDS::Face(in_shape));
            return point<dim>(surf.Value(u, v));
          }
        case TopAbs_EDGE:
          {
            BRepAdaptor_Curve curve(TopoDS::Edge(in_shape));
            return point<dim>(curve.Value(u));
          }
        default:
          Assert(false, ExcUnsupportedShape());
      }
    return Point<dim>();
  }

  std::tuple<Point<3>, Tensor<1, 3>, double, double>
  push_forward_and_differential_forms(const TopoDS_Face &face,
                                      const double       u,
                                      const double       v,
                                      const double /*tolerance*/)
  {
    Handle(Geom_Surface) SurfToProj = BRep_Tool::Surface(face);
    GeomLProp_SLProps props(SurfToProj, u, v, 1, 1e-7);
    gp_Pnt            Value = props.Value();
    Assert(props.IsNormalDefined(), ExcMessage("Normal is not well defined!"));
    gp_Dir Normal = props.Normal();
    Assert(props.IsCurvatureDefined(),
           ExcMessage("Curvature is not well defined!"));
    Standard_Real Min_Curvature = props.MinCurvature();
    Standard_Real Max_Curvature = props.MaxCurvature();
    Tensor<1, 3>  normal        = Point<3>(Normal.X(), Normal.Y(), Normal.Z());

    // In the case your manifold changes from convex to concave or viceversa
    // the normal could jump from "inner" to "outer" normal.
    // However, you should be able to change the normal sense preserving
    // the manifold orientation:
    if (face.Orientation() == TopAbs_REVERSED)
      {
        normal *= -1;
        Min_Curvature *= -1;
        Max_Curvature *= -1;
      }

    return std::tuple<Point<3>, Tensor<1, 3>, double, double>(point<3>(Value),
                                                              normal,
                                                              Min_Curvature,
                                                              Max_Curvature);
  }



  template <int spacedim>
  void
  create_triangulation(const TopoDS_Face &         face,
                       Triangulation<2, spacedim> &tria)
  {
    BRepAdaptor_Surface surf(face);
    const double        u0 = surf.FirstUParameter();
    const double        u1 = surf.LastUParameter();
    const double        v0 = surf.FirstVParameter();
    const double        v1 = surf.LastVParameter();

    std::vector<CellData<2>>     cells;
    std::vector<Point<spacedim>> vertices;
    SubCellData                  t;

    vertices.push_back(point<spacedim>(surf.Value(u0, v0)));
    vertices.push_back(point<spacedim>(surf.Value(u1, v0)));
    vertices.push_back(point<spacedim>(surf.Value(u0, v1)));
    vertices.push_back(point<spacedim>(surf.Value(u1, v1)));

    CellData<2> cell;
    for (unsigned int i = 0; i < 4; ++i)
      cell.vertices[i] = i;

    cells.push_back(cell);
    tria.create_triangulation(vertices, cells, t);
  }

#  include "utilities.inst"

} // namespace OpenCASCADE

DEAL_II_NAMESPACE_CLOSE

#endif

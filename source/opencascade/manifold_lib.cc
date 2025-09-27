// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/config.h>

#include <deal.II/opencascade/manifold_lib.h>

#ifdef DEAL_II_WITH_OPENCASCADE


#  include <BRepAdaptor_CompCurve.hxx>
#  include <BRepAdaptor_Curve.hxx>
#  if !DEAL_II_OPENCASCADE_VERSION_GTE(7, 6, 0)
#    include <BRepAdaptor_HCompCurve.hxx>
#    include <BRepAdaptor_HCurve.hxx>
#  endif
#  include <BRepTools.hxx>
#  include <BRep_Tool.hxx>
#  include <GCPnts_AbscissaPoint.hxx>
#  include <ShapeAnalysis_Curve.hxx>
#  include <ShapeAnalysis_Surface.hxx>
#  include <TopoDS.hxx>
#  if !DEAL_II_OPENCASCADE_VERSION_GTE(7, 0, 0)
#    include <Handle_Adaptor3d_HCurve.hxx>
#  endif


DEAL_II_NAMESPACE_OPEN


namespace OpenCASCADE
{
  namespace
  {
    /**
     * Return a Geometrical curve representation for the given
     * TopoDS_Shape. This function will fail when the given shape is
     * not of topological dimension one.
     */
#  if DEAL_II_OPENCASCADE_VERSION_GTE(7, 6, 0)
    Handle_Adaptor3d_Curve
    curve_adaptor(const TopoDS_Shape &shape)
    {
      Assert((shape.ShapeType() == TopAbs_WIRE) ||
               (shape.ShapeType() == TopAbs_EDGE),
             ExcUnsupportedShape());
      if (shape.ShapeType() == TopAbs_WIRE)
        return Handle(BRepAdaptor_CompCurve)(
          new BRepAdaptor_CompCurve(TopoDS::Wire(shape)));
      else if (shape.ShapeType() == TopAbs_EDGE)
        return Handle(BRepAdaptor_Curve)(
          new BRepAdaptor_Curve(TopoDS::Edge(shape)));

      DEAL_II_ASSERT_UNREACHABLE();
      return Handle(BRepAdaptor_Curve)(new BRepAdaptor_Curve());
    }
#  else
    Handle_Adaptor3d_HCurve
    curve_adaptor(const TopoDS_Shape &shape)
    {
      Assert((shape.ShapeType() == TopAbs_WIRE) ||
               (shape.ShapeType() == TopAbs_EDGE),
             ExcUnsupportedShape());
      if (shape.ShapeType() == TopAbs_WIRE)
        return Handle(BRepAdaptor_HCompCurve)(
          new BRepAdaptor_HCompCurve(TopoDS::Wire(shape)));
      else if (shape.ShapeType() == TopAbs_EDGE)
        return Handle(BRepAdaptor_HCurve)(
          new BRepAdaptor_HCurve(TopoDS::Edge(shape)));

      DEAL_II_ASSERT_UNREACHABLE();
      return Handle(BRepAdaptor_HCurve)(new BRepAdaptor_HCurve());
    }
#  endif



    // Helper internal functions.
    double
    shape_length(const TopoDS_Shape &sh)
    {
#  if DEAL_II_OPENCASCADE_VERSION_GTE(7, 6, 0)
      Handle_Adaptor3d_Curve adapt = curve_adaptor(sh);
      return GCPnts_AbscissaPoint::Length(*adapt);
#  else
      Handle_Adaptor3d_HCurve adapt = curve_adaptor(sh);
      return GCPnts_AbscissaPoint::Length(adapt->GetCurve());
#  endif
    }
  } // namespace

  /*======================= NormalProjectionManifold =========================*/
  template <int dim, int spacedim>
  NormalProjectionManifold<dim, spacedim>::NormalProjectionManifold(
    const TopoDS_Shape &sh,
    const double        tolerance)
    : sh(sh)
    , tolerance(tolerance)
  {
    Assert(spacedim == 3, ExcNotImplemented());
  }



  template <int dim, int spacedim>
  std::unique_ptr<Manifold<dim, spacedim>>
  NormalProjectionManifold<dim, spacedim>::clone() const
  {
    return std::unique_ptr<Manifold<dim, spacedim>>(
      new NormalProjectionManifold(sh, tolerance));
  }



  template <int dim, int spacedim>
  Point<spacedim>
  NormalProjectionManifold<dim, spacedim>::project_to_manifold(
    const ArrayView<const Point<spacedim>> &surrounding_points,
    const Point<spacedim>                  &candidate) const
  {
    (void)surrounding_points;
    if constexpr (running_in_debug_mode())
      {
        for (unsigned int i = 0; i < surrounding_points.size(); ++i)
          Assert(closest_point(sh, surrounding_points[i], tolerance)
                     .distance(surrounding_points[i]) <
                   std::max(tolerance * surrounding_points[i].norm(),
                            tolerance),
                 ExcPointNotOnManifold<spacedim>(surrounding_points[i]));
      }
    return closest_point(sh, candidate, tolerance);
  }


  /*===================== DirectionalProjectionManifold ======================*/
  template <int dim, int spacedim>
  DirectionalProjectionManifold<dim, spacedim>::DirectionalProjectionManifold(
    const TopoDS_Shape        &sh,
    const Tensor<1, spacedim> &direction,
    const double               tolerance)
    : sh(sh)
    , direction(direction)
    , tolerance(tolerance)
  {
    Assert(spacedim == 3, ExcNotImplemented());
  }



  template <int dim, int spacedim>
  std::unique_ptr<Manifold<dim, spacedim>>
  DirectionalProjectionManifold<dim, spacedim>::clone() const
  {
    return std::unique_ptr<Manifold<dim, spacedim>>(
      new DirectionalProjectionManifold(sh, direction, tolerance));
  }



  template <int dim, int spacedim>
  Point<spacedim>
  DirectionalProjectionManifold<dim, spacedim>::project_to_manifold(
    const ArrayView<const Point<spacedim>> &surrounding_points,
    const Point<spacedim>                  &candidate) const
  {
    (void)surrounding_points;
    if constexpr (running_in_debug_mode())
      {
        for (unsigned int i = 0; i < surrounding_points.size(); ++i)
          Assert(closest_point(sh, surrounding_points[i], tolerance)
                     .distance(surrounding_points[i]) <
                   std::max(tolerance * surrounding_points[i].norm(),
                            tolerance),
                 ExcPointNotOnManifold<spacedim>(surrounding_points[i]));
      }
    return line_intersection(sh, candidate, direction, tolerance);
  }



  /*===================== NormalToMeshProjectionManifold =====================*/
  template <int dim, int spacedim>
  NormalToMeshProjectionManifold<dim, spacedim>::NormalToMeshProjectionManifold(
    const TopoDS_Shape &sh,
    const double        tolerance)
    : sh(sh)
    , tolerance(tolerance)
  {
    Assert(spacedim == 3, ExcNotImplemented());
    Assert(
      std::get<0>(count_elements(sh)) > 0,
      ExcMessage(
        "NormalToMeshProjectionManifold needs a shape containing faces to operate."));
  }

  template <int dim, int spacedim>
  std::unique_ptr<Manifold<dim, spacedim>>
  NormalToMeshProjectionManifold<dim, spacedim>::clone() const
  {
    return std::unique_ptr<Manifold<dim, spacedim>>(
      new NormalToMeshProjectionManifold<dim, spacedim>(sh, tolerance));
  }


  namespace
  {
    template <int spacedim>
    Point<spacedim>
    internal_project_to_manifold(const TopoDS_Shape &,
                                 const double,
                                 const ArrayView<const Point<spacedim>> &,
                                 const Point<spacedim> &)
    {
      DEAL_II_NOT_IMPLEMENTED();
      return {};
    }

    template <>
    Point<3>
    internal_project_to_manifold(
      const TopoDS_Shape              &sh,
      const double                     tolerance,
      const ArrayView<const Point<3>> &surrounding_points,
      const Point<3>                  &candidate)
    {
      constexpr int       spacedim = 3;
      TopoDS_Shape        out_shape;
      Tensor<1, spacedim> average_normal;
      if constexpr (running_in_debug_mode())
        {
          for (const auto &point : surrounding_points)
            {
              Assert(closest_point(sh, point, tolerance).distance(point) <
                       std::max(tolerance * point.norm(), tolerance),
                     ExcPointNotOnManifold<spacedim>(point));
            }
        }

      switch (surrounding_points.size())
        {
          case 2:
            {
              for (const auto &point : surrounding_points)
                {
                  std::tuple<Point<3>, Tensor<1, 3>, double, double>
                    p_and_diff_forms =
                      closest_point_and_differential_forms(sh,
                                                           point,
                                                           tolerance);
                  average_normal += std::get<1>(p_and_diff_forms);
                }

              average_normal /= 2.0;

              Assert(
                average_normal.norm() > 1e-4,
                ExcMessage(
                  "Failed to refine cell: the average of the surface normals at the surrounding edge turns out to be a null vector, making the projection direction undetermined."));

              Tensor<1, 3> T = surrounding_points[0] - surrounding_points[1];
              T /= T.norm();
              average_normal = average_normal - (average_normal * T) * T;
              average_normal /= average_normal.norm();
              break;
            }
          case 4:
            {
              Tensor<1, 3> u = surrounding_points[1] - surrounding_points[0];
              Tensor<1, 3> v = surrounding_points[2] - surrounding_points[0];
              const double n1_coords[3] = {u[1] * v[2] - u[2] * v[1],
                                           u[2] * v[0] - u[0] * v[2],
                                           u[0] * v[1] - u[1] * v[0]};
              Tensor<1, 3> n1(n1_coords);
              n1 = n1 / n1.norm();
              u  = surrounding_points[2] - surrounding_points[3];
              v  = surrounding_points[1] - surrounding_points[3];
              const double n2_coords[3] = {u[1] * v[2] - u[2] * v[1],
                                           u[2] * v[0] - u[0] * v[2],
                                           u[0] * v[1] - u[1] * v[0]};
              Tensor<1, 3> n2(n2_coords);
              n2 = n2 / n2.norm();

              average_normal = (n1 + n2) / 2.0;

              Assert(
                average_normal.norm() > tolerance,
                ExcMessage(
                  "Failed to refine cell: the normal estimated via the surrounding points turns out to be a null vector, making the projection direction undetermined."));

              average_normal /= average_normal.norm();
              break;
            }
          case 8:
            {
              Tensor<1, 3> u = surrounding_points[1] - surrounding_points[0];
              Tensor<1, 3> v = surrounding_points[2] - surrounding_points[0];
              const double n1_coords[3] = {u[1] * v[2] - u[2] * v[1],
                                           u[2] * v[0] - u[0] * v[2],
                                           u[0] * v[1] - u[1] * v[0]};
              Tensor<1, 3> n1(n1_coords);
              n1 = n1 / n1.norm();
              u  = surrounding_points[2] - surrounding_points[3];
              v  = surrounding_points[1] - surrounding_points[3];
              const double n2_coords[3] = {u[1] * v[2] - u[2] * v[1],
                                           u[2] * v[0] - u[0] * v[2],
                                           u[0] * v[1] - u[1] * v[0]};
              Tensor<1, 3> n2(n2_coords);
              n2 = n2 / n2.norm();
              u  = surrounding_points[4] - surrounding_points[7];
              v  = surrounding_points[6] - surrounding_points[7];
              const double n3_coords[3] = {u[1] * v[2] - u[2] * v[1],
                                           u[2] * v[0] - u[0] * v[2],
                                           u[0] * v[1] - u[1] * v[0]};
              Tensor<1, 3> n3(n3_coords);
              n3 = n3 / n3.norm();
              u  = surrounding_points[6] - surrounding_points[7];
              v  = surrounding_points[5] - surrounding_points[7];
              const double n4_coords[3] = {u[1] * v[2] - u[2] * v[1],
                                           u[2] * v[0] - u[0] * v[2],
                                           u[0] * v[1] - u[1] * v[0]};
              Tensor<1, 3> n4(n4_coords);
              n4 = n4 / n4.norm();

              average_normal = (n1 + n2 + n3 + n4) / 4.0;

              Assert(
                average_normal.norm() > tolerance,
                ExcMessage(
                  "Failed to refine cell: the normal estimated via the surrounding points turns out to be a null vector, making the projection direction undetermined."));

              average_normal /= average_normal.norm();
              break;
            }
          default:
            {
              // Given an arbitrary number of points we compute all the possible
              // normal vectors
              for (unsigned int i = 0; i < surrounding_points.size(); ++i)
                for (unsigned int j = 0; j < surrounding_points.size(); ++j)
                  if (j != i)
                    for (unsigned int k = 0; k < surrounding_points.size(); ++k)
                      if (k != j && k != i)
                        {
                          Tensor<1, 3> u =
                            surrounding_points[i] - surrounding_points[j];
                          Tensor<1, 3> v =
                            surrounding_points[i] - surrounding_points[k];
                          const double n_coords[3] = {u[1] * v[2] - u[2] * v[1],
                                                      u[2] * v[0] - u[0] * v[2],
                                                      u[0] * v[1] -
                                                        u[1] * v[0]};
                          Tensor<1, 3> n1(n_coords);
                          if (n1.norm() > tolerance)
                            {
                              n1 = n1 / n1.norm();
                              if (average_normal.norm() < tolerance)
                                average_normal = n1;
                              else
                                {
                                  auto dot_prod = n1 * average_normal;
                                  // We check that the direction of the normal
                                  // vector w.r.t the current average, and make
                                  // sure we flip it if it is opposite
                                  if (dot_prod > 0)
                                    average_normal += n1;
                                  else
                                    average_normal -= n1;
                                }
                            }
                        }
              Assert(
                average_normal.norm() > tolerance,
                ExcMessage(
                  "Failed to compute a normal: the normal estimated via the surrounding points turns out to be a null vector, making the projection direction undetermined."));
              average_normal = average_normal / average_normal.norm();
              break;
            }
        }

      return line_intersection(sh, candidate, average_normal, tolerance);
    }
  } // namespace


  template <int dim, int spacedim>
  Point<spacedim>
  NormalToMeshProjectionManifold<dim, spacedim>::project_to_manifold(
    const ArrayView<const Point<spacedim>> &surrounding_points,
    const Point<spacedim>                  &candidate) const
  {
    return internal_project_to_manifold(sh,
                                        tolerance,
                                        surrounding_points,
                                        candidate);
  }


  /*==================== ArclengthProjectionLineManifold =====================*/
  template <int dim, int spacedim>
  ArclengthProjectionLineManifold<dim, spacedim>::
    ArclengthProjectionLineManifold(const TopoDS_Shape &sh,
                                    const double        tolerance)
    :

    ChartManifold<dim, spacedim, 1>(sh.Closed() ? Point<1>(shape_length(sh)) :
                                                  Point<1>())
    , sh(sh)
    , curve(curve_adaptor(sh))
    , tolerance(tolerance)
    , length(shape_length(sh))
  {
    Assert(spacedim >= 2, ExcImpossibleInDimSpacedim(dim, spacedim));
  }



  template <int dim, int spacedim>
  std::unique_ptr<Manifold<dim, spacedim>>
  ArclengthProjectionLineManifold<dim, spacedim>::clone() const
  {
    return std::unique_ptr<Manifold<dim, spacedim>>(
      new ArclengthProjectionLineManifold(sh, tolerance));
  }



  template <int dim, int spacedim>
  Point<1>
  ArclengthProjectionLineManifold<dim, spacedim>::pull_back(
    const Point<spacedim> &space_point) const
  {
    double              t(0.0);
    ShapeAnalysis_Curve curve_analysis;
    gp_Pnt              proj;

    const double dist = curve_analysis.Project(
#  if DEAL_II_OPENCASCADE_VERSION_GTE(7, 6, 0)
      *curve, point(space_point), tolerance, proj, t, true);
#  else
      curve->GetCurve(), point(space_point), tolerance, proj, t, true);
#  endif

    (void)dist;
    Assert(dist < tolerance * length,
           ExcPointNotOnManifold<spacedim>(space_point));

    return Point<1>(GCPnts_AbscissaPoint::Length(
#  if DEAL_II_OPENCASCADE_VERSION_GTE(7, 6, 0)
      *curve, curve->FirstParameter(), t));
#  else
      curve->GetCurve(), curve->GetCurve().FirstParameter(), t));
#  endif
  }



  template <int dim, int spacedim>
  Point<spacedim>
  ArclengthProjectionLineManifold<dim, spacedim>::push_forward(
    const Point<1> &chart_point) const
  {
#  if DEAL_II_OPENCASCADE_VERSION_GTE(7, 6, 0)
    GCPnts_AbscissaPoint AP(*curve, chart_point[0], curve->FirstParameter());
    gp_Pnt               P = curve->Value(AP.Parameter());
#  else
    GCPnts_AbscissaPoint AP(curve->GetCurve(),
                            chart_point[0],
                            curve->GetCurve().FirstParameter());
    gp_Pnt               P = curve->GetCurve().Value(AP.Parameter());
#  endif

    return point<spacedim>(P);
  }

  template <int dim, int spacedim>
  NURBSPatchManifold<dim, spacedim>::NURBSPatchManifold(const TopoDS_Face &face,
                                                        const double tolerance)
    : face(face)
    , tolerance(tolerance)
  {}



  template <int dim, int spacedim>
  std::unique_ptr<Manifold<dim, spacedim>>
  NURBSPatchManifold<dim, spacedim>::clone() const
  {
    return std::unique_ptr<Manifold<dim, spacedim>>(
      new NURBSPatchManifold<dim, spacedim>(face, tolerance));
  }



  template <int dim, int spacedim>
  Point<2>
  NURBSPatchManifold<dim, spacedim>::pull_back(
    const Point<spacedim> &space_point) const
  {
    Handle(Geom_Surface) SurfToProj = BRep_Tool::Surface(face);

    ShapeAnalysis_Surface projector(SurfToProj);
    gp_Pnt2d proj_params = projector.ValueOfUV(point(space_point), tolerance);

    double u = proj_params.X();
    double v = proj_params.Y();

    return {u, v};
  }

  template <int dim, int spacedim>
  Point<spacedim>
  NURBSPatchManifold<dim, spacedim>::push_forward(
    const Point<2> &chart_point) const
  {
    return ::dealii::OpenCASCADE::push_forward<spacedim>(face,
                                                         chart_point[0],
                                                         chart_point[1]);
  }

  template <int dim, int spacedim>
  DerivativeForm<1, 2, spacedim>
  NURBSPatchManifold<dim, spacedim>::push_forward_gradient(
    const Point<2> &chart_point) const
  {
    DerivativeForm<1, 2, spacedim> DX;
    Handle(Geom_Surface) surf = BRep_Tool::Surface(face);

    gp_Pnt q;
    gp_Vec Du, Dv;
    surf->D1(chart_point[0], chart_point[1], q, Du, Dv);

    DX[0][0] = Du.X();
    DX[1][0] = Du.Y();
    if (spacedim > 2)
      DX[2][0] = Du.Z();
    else
      Assert(std::abs(Du.Z()) < tolerance,
             ExcMessage(
               "Expecting derivative along Z to be zero! Bailing out."));
    DX[0][1] = Dv.X();
    DX[1][1] = Dv.Y();
    if (spacedim > 2)
      DX[2][1] = Dv.Z();
    else
      Assert(std::abs(Dv.Z()) < tolerance,
             ExcMessage(
               "Expecting derivative along Z to be zero! Bailing out."));
    return DX;
  }

  template <int dim, int spacedim>
  std::tuple<double, double, double, double>
  NURBSPatchManifold<dim, spacedim>::get_uv_bounds() const
  {
    Standard_Real umin, umax, vmin, vmax;
    BRepTools::UVBounds(face, umin, umax, vmin, vmax);
    return std::make_tuple(umin, umax, vmin, vmax);
  }

// We don't build the .inst file if deal.II isn't configured
// with GMSH, but doxygen doesn't know that and tries to find that
// file anyway for parsing -- which then of course it fails on. So
// exclude the following from doxygen consideration.
#  ifndef DOXYGEN
#    include "opencascade/manifold_lib.inst"
#  endif
} // end namespace OpenCASCADE

DEAL_II_NAMESPACE_CLOSE

#endif

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
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

#include <deal.II/cgal/intersections.h>

#include <algorithm>

#ifdef DEAL_II_WITH_CGAL

#  include <deal.II/base/quadrature_lib.h>
#  include <deal.II/base/utilities.h>

#  include <deal.II/fe/mapping.h>

#  include <deal.II/grid/tria.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <CGAL/Boolean_set_operations_2.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#  include <deal.II/cgal/utilities.h>

#  include <CGAL/Cartesian.h>
#  include <CGAL/Circular_kernel_intersections.h>
#  include <CGAL/Constrained_Delaunay_triangulation_2.h>
#  include <CGAL/Delaunay_mesh_face_base_2.h>
#  include <CGAL/Delaunay_mesh_size_criteria_2.h>
#  include <CGAL/Delaunay_mesher_2.h>
#  include <CGAL/Delaunay_triangulation_2.h>
#  include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#  include <CGAL/Kernel_traits.h>
#  include <CGAL/Polygon_2.h>
#  include <CGAL/Polygon_with_holes_2.h>
#  include <CGAL/Projection_traits_xy_3.h>
#  include <CGAL/Segment_3.h>
#  include <CGAL/Simple_cartesian.h>
#  include <CGAL/Surface_mesh/Surface_mesh.h>
#  include <CGAL/Tetrahedron_3.h>
#  include <CGAL/Triangle_2.h>
#  include <CGAL/Triangle_3.h>
#  include <CGAL/Triangulation_2.h>
#  include <CGAL/Triangulation_3.h>
#  include <CGAL/Triangulation_face_base_with_id_2.h>
#  include <CGAL/Triangulation_face_base_with_info_2.h>

#  include <optional>
#  include <variant>

DEAL_II_NAMESPACE_OPEN

namespace CGALWrappers
{
  using K       = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;
  using K_exact = CGAL::Exact_predicates_exact_constructions_kernel;
  using CGALPolygon          = CGAL::Polygon_2<K>;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;
  using CGALTriangle2        = K::Triangle_2;
  using CGALTriangle3        = K::Triangle_3;
  using CGALTriangle3_exact  = K_exact::Triangle_3;
  using CGALPoint2           = K::Point_2;
  using CGALPoint3           = K::Point_3;
  using CGALPoint3_exact     = K_exact::Point_3;
  using CGALSegment2         = K::Segment_2;
  using Surface_mesh         = CGAL::Surface_mesh<K_exact::Point_3>;
  using CGALSegment3         = K::Segment_3;
  using CGALSegment3_exact   = K_exact::Segment_3;
  using CGALTetra            = K::Tetrahedron_3;
  using CGALTetra_exact      = K_exact::Tetrahedron_3;
  using Triangulation2       = CGAL::Triangulation_2<K>;
  using Triangulation3       = CGAL::Triangulation_3<K>;
  using Triangulation3_exact = CGAL::Triangulation_3<K_exact>;

  struct FaceInfo2
  {
    FaceInfo2() = default;
    int nesting_level;
    bool
    in_domain()
    {
      return nesting_level % 2 == 1;
    }
  };

  using Vb       = CGAL::Triangulation_vertex_base_2<K>;
  using Fbb      = CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K>;
  using CFb      = CGAL::Constrained_triangulation_face_base_2<K, Fbb>;
  using Fb       = CGAL::Delaunay_mesh_face_base_2<K, CFb>;
  using Tds      = CGAL::Triangulation_data_structure_2<Vb, Fb>;
  using Itag     = CGAL::Exact_predicates_tag;
  using CDT      = CGAL::Constrained_Delaunay_triangulation_2<K, Tds, Itag>;
  using Criteria = CGAL::Delaunay_mesh_size_criteria_2<CDT>;
  using Vertex_handle = CDT::Vertex_handle;
  using Face_handle   = CDT::Face_handle;

  template <class T, class... Types>
  const T *
  get_if_(const std::variant<Types...> *v)
  {
    return std::get_if<T>(v);
  }

  template <class T, class... Types>
  const T *
  get_if_(const boost::variant<Types...> *v)
  {
    return boost::get<T>(v);
  }

  namespace internal
  {
    namespace
    {
      /**
       * Take a boost::variant object and convert it to a std::variant
       * object by applying the visitor pattern.
       */
      template <typename TargetVariant>
      struct Repackage : boost::static_visitor<TargetVariant>
      {
        template <typename T>
        TargetVariant
        operator()(const T &t) const
        {
          return TargetVariant(t);
        }
      };

      /**
       * Convert a boost::optional<std::variant<...>> to the
       * corresponding C++ type using std::optional and
       * std::variant. The former is what CGAL gives us, the latter is
       * what we want to use because we like to use std data types.
       */
      template <typename... Types>
      std::optional<std::variant<Types...>>
      convert_boost_to_std(const boost::optional<boost::variant<Types...>> &x)
      {
        if (x)
          {
            // The boost::optional object contains an object of type
            // boost::variant. We need to unpack which type the
            // variant contains, and re-package that into a
            // std::variant. This is easily done using a visitor
            // object.
            using std_variant = std::variant<Types...>;
            return boost::apply_visitor(Repackage<std_variant>(), *x);
          }
        else
          {
            // The boost::optional object was empty. Return an empty
            // std::optional object.
            return {};
          }
      }

      template <typename... Types>
      const std::optional<std::variant<Types...>> &
      convert_boost_to_std(const std::optional<std::variant<Types...>> &opt)
      {
        return opt;
      }
    } // namespace



    void
    mark_domains(CDT                  &ct,
                 Face_handle           start,
                 int                   index,
                 std::list<CDT::Edge> &border)
    {
      if (start->info().nesting_level != -1)
        {
          return;
        }
      std::list<Face_handle> queue;
      queue.push_back(start);
      while (!queue.empty())
        {
          Face_handle fh = queue.front();
          queue.pop_front();
          if (fh->info().nesting_level == -1)
            {
              fh->info().nesting_level = index;
              for (int i = 0; i < 3; i++)
                {
                  CDT::Edge   e(fh, i);
                  Face_handle n = fh->neighbor(i);
                  if (n->info().nesting_level == -1)
                    {
                      if (ct.is_constrained(e))
                        border.push_back(e);
                      else
                        queue.push_back(n);
                    }
                }
            }
        }
    }



    void
    mark_domains(CDT &cdt)
    {
      for (CDT::Face_handle f : cdt.all_face_handles())
        {
          f->info().nesting_level = -1;
        }
      std::list<CDT::Edge> border;
      mark_domains(cdt, cdt.infinite_face(), 0, border);
      while (!border.empty())
        {
          CDT::Edge e = border.front();
          border.pop_front();
          Face_handle n = e.first->neighbor(e.second);
          if (n->info().nesting_level == -1)
            {
              mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
            }
        }
    }

    // Collection of utilities that compute intersection between simplices
    // identified by array of points. The return type is the one of
    // CGAL::intersection(), i.e. a std::optional<std::variant<>>.
    // Intersection between 2d and 3d objects and 1d/3d objects are available
    // only with CGAL versions greater or equal than 5.5, hence the
    // corresponding functions are guarded by #ifdef directives. All the
    // signatures follow the convection that the first entity has an intrinsic
    // dimension higher than the second one.

    std::optional<std::variant<CGALPoint2,
                               CGALSegment2,
                               CGALTriangle2,
                               std::vector<CGALPoint2>>>
    compute_intersection_triangle_triangle(
      const ArrayView<const Point<2>> &triangle0,
      const ArrayView<const Point<2>> &triangle1)
    {
      AssertDimension(triangle0.size(), 3);
      AssertDimension(triangle0.size(), triangle1.size());

      std::array<CGALPoint2, 3> pts0, pts1;

      std::transform(triangle0.begin(),
                     triangle0.end(),
                     pts0.begin(),
                     &CGALWrappers::dealii_point_to_cgal_point<CGALPoint2, 2>);

      std::transform(triangle1.begin(),
                     triangle1.end(),
                     pts1.begin(),
                     &CGALWrappers::dealii_point_to_cgal_point<CGALPoint2, 2>);

      CGALTriangle2 cgal_triangle0{pts0[0], pts0[1], pts0[2]};
      CGALTriangle2 cgal_triangle1{pts1[0], pts1[1], pts1[2]};
      return convert_boost_to_std(
        CGAL::intersection(cgal_triangle0, cgal_triangle1));
    }


    std::optional<std::variant<CGALPoint2, CGALSegment2>>
    compute_intersection_triangle_segment(
      const ArrayView<const Point<2>> &triangle,
      const ArrayView<const Point<2>> &segment)
    {
      AssertDimension(triangle.size(), 3);
      AssertDimension(segment.size(), 2);

      std::array<CGALPoint2, 3> pts0;
      std::array<CGALPoint2, 2> pts1;

      std::transform(triangle.begin(),
                     triangle.end(),
                     pts0.begin(),
                     &CGALWrappers::dealii_point_to_cgal_point<CGALPoint2, 2>);

      std::transform(segment.begin(),
                     segment.end(),
                     pts1.begin(),
                     &CGALWrappers::dealii_point_to_cgal_point<CGALPoint2, 2>);

      CGALTriangle2 cgal_triangle{pts0[0], pts0[1], pts0[2]};
      CGALSegment2  cgal_segment{pts1[0], pts1[1]};
      return convert_boost_to_std(
        CGAL::intersection(cgal_segment, cgal_triangle));
    }



    // rectangle-rectangle
    std::vector<Polygon_with_holes_2>
    compute_intersection_rect_rect(const ArrayView<const Point<2>> &rectangle0,
                                   const ArrayView<const Point<2>> &rectangle1)
    {
      AssertDimension(rectangle0.size(), 4);
      AssertDimension(rectangle0.size(), rectangle1.size());

      std::array<CGALPoint2, 4> pts0, pts1;

      std::transform(rectangle0.begin(),
                     rectangle0.end(),
                     pts0.begin(),
                     &CGALWrappers::dealii_point_to_cgal_point<CGALPoint2, 2>);

      std::transform(rectangle1.begin(),
                     rectangle1.end(),
                     pts1.begin(),
                     &CGALWrappers::dealii_point_to_cgal_point<CGALPoint2, 2>);

      const CGALPolygon first_poly{pts0.begin(), pts0.end()};
      const CGALPolygon second_poly{pts1.begin(), pts1.end()};

      std::vector<Polygon_with_holes_2> poly_list;
      CGAL::intersection(first_poly,
                         second_poly,
                         std::back_inserter(poly_list));
      return poly_list;
    }



    std::optional<std::variant<CGALPoint3, CGALSegment3>>
    compute_intersection_tetra_segment(
      const ArrayView<const Point<3>> &tetrahedron,
      const ArrayView<const Point<3>> &segment)
    {
#  if DEAL_II_CGAL_VERSION_GTE(5, 5, 0)

      AssertDimension(tetrahedron.size(), 4);
      AssertDimension(segment.size(), 2);

      std::array<CGALPoint3, 4> pts0;
      std::array<CGALPoint3, 2> pts1;

      std::transform(tetrahedron.begin(),
                     tetrahedron.end(),
                     pts0.begin(),
                     &CGALWrappers::dealii_point_to_cgal_point<CGALPoint3, 3>);

      std::transform(segment.begin(),
                     segment.end(),
                     pts1.begin(),
                     &CGALWrappers::dealii_point_to_cgal_point<CGALPoint3, 3>);

      CGALTetra    cgal_tetrahedron{pts0[0], pts0[1], pts0[2], pts0[3]};
      CGALSegment3 cgal_segment{pts1[0], pts1[1]};
      return convert_boost_to_std(
        CGAL::intersection(cgal_segment, cgal_tetrahedron));
#  else
      Assert(
        false,
        ExcMessage(
          "This function requires a version of CGAL greater or equal than 5.5."));
      (void)tetrahedron;
      (void)segment;
      return {};
#  endif
    }


    // tetra, triangle
    std::optional<std::variant<CGALPoint3,
                               CGALSegment3,
                               CGALTriangle3,
                               std::vector<CGALPoint3>>>
    compute_intersection_tetra_triangle(
      const ArrayView<const Point<3>> &tetrahedron,
      const ArrayView<const Point<3>> &triangle)
    {
#  if DEAL_II_CGAL_VERSION_GTE(5, 5, 0)

      AssertDimension(tetrahedron.size(), 4);
      AssertDimension(triangle.size(), 3);

      std::array<CGALPoint3, 4> pts0;
      std::array<CGALPoint3, 3> pts1;

      std::transform(tetrahedron.begin(),
                     tetrahedron.end(),
                     pts0.begin(),
                     &CGALWrappers::dealii_point_to_cgal_point<CGALPoint3, 3>);

      std::transform(triangle.begin(),
                     triangle.end(),
                     pts1.begin(),
                     &CGALWrappers::dealii_point_to_cgal_point<CGALPoint3, 3>);

      CGALTetra     cgal_tetrahedron{pts0[0], pts0[1], pts0[2], pts0[3]};
      CGALTriangle3 cgal_triangle{pts1[0], pts1[1], pts1[2]};
      return convert_boost_to_std(
        CGAL::intersection(cgal_triangle, cgal_tetrahedron));
#  else

      Assert(
        false,
        ExcMessage(
          "This function requires a version of CGAL greater or equal than 5.5."));
      (void)tetrahedron;
      (void)triangle;
      return {};
#  endif
    }

    // quad-quad
    std::vector<std::array<Point<2>, 3>>
    compute_intersection_quad_quad(const ArrayView<const Point<2>> &quad0,
                                   const ArrayView<const Point<2>> &quad1,
                                   const double                     tol)
    {
      AssertDimension(quad0.size(), 4);
      AssertDimension(quad0.size(), quad1.size());

      const auto intersection_test =
        internal::compute_intersection_rect_rect(quad0, quad1);

      if (!intersection_test.empty())
        {
          const auto        &poly      = intersection_test[0].outer_boundary();
          const unsigned int size_poly = poly.size();
          if (size_poly == 3)
            {
              // intersection is a triangle itself, so directly return its
              // vertices.
              return {
                {{CGALWrappers::cgal_point_to_dealii_point<2>(poly.vertex(0)),
                  CGALWrappers::cgal_point_to_dealii_point<2>(poly.vertex(1)),
                  CGALWrappers::cgal_point_to_dealii_point<2>(
                    poly.vertex(2))}}};
            }
          else if (size_poly >= 4)
            {
              // intersection is a polygon, need to triangulate it.
              std::vector<std::array<Point<2>, 3>> collection;

              CDT cdt;
              cdt.insert_constraint(poly.vertices_begin(),
                                    poly.vertices_end(),
                                    true);

              internal::mark_domains(cdt);

              for (Face_handle f : cdt.finite_face_handles())
                {
                  if (f->info().in_domain() &&
                      CGAL::to_double(cdt.triangle(f).area()) > tol)
                    {
                      collection.push_back(
                        {{CGALWrappers::cgal_point_to_dealii_point<2>(
                            cdt.triangle(f).vertex(0)),
                          CGALWrappers::cgal_point_to_dealii_point<2>(
                            cdt.triangle(f).vertex(1)),
                          CGALWrappers::cgal_point_to_dealii_point<2>(
                            cdt.triangle(f).vertex(2))}});
                    }
                }
              return collection;
            }
          else
            {
              Assert(false, ExcMessage("The polygon is degenerate."));
              return {};
            }
        }
      else
        {
          return {};
        }
    }

    // Specialization for quad \cap line
    std::vector<std::array<Point<2>, 2>>
    compute_intersection_quad_line(const ArrayView<const Point<2>> &quad,
                                   const ArrayView<const Point<2>> &line,
                                   const double                     tol)
    {
      AssertDimension(quad.size(), 4);
      AssertDimension(line.size(), 2);

      std::array<CGALPoint2, 4> pts;

      std::transform(quad.begin(),
                     quad.end(),
                     pts.begin(),
                     &CGALWrappers::dealii_point_to_cgal_point<CGALPoint2, 2>);

      CGALPolygon poly(pts.begin(), pts.end());

      CGALSegment2 segm(
        CGALWrappers::dealii_point_to_cgal_point<CGALPoint2>(line[0]),
        CGALWrappers::dealii_point_to_cgal_point<CGALPoint2>(line[1]));
      CDT cdt;
      cdt.insert_constraint(poly.vertices_begin(), poly.vertices_end(), true);
      std::vector<std::array<Point<2>, 2>> vertices;
      internal::mark_domains(cdt);
      for (Face_handle f : cdt.finite_face_handles())
        {
          if (f->info().in_domain() &&
              CGAL::to_double(cdt.triangle(f).area()) > tol &&
              CGAL::do_intersect(segm, cdt.triangle(f)))
            {
              const auto intersection =
                CGAL::intersection(segm, cdt.triangle(f));
              if (const CGALSegment2 *s = get_if_<CGALSegment2>(&*intersection))
                {
                  vertices.push_back(
                    {{CGALWrappers::cgal_point_to_dealii_point<2>((*s)[0]),
                      CGALWrappers::cgal_point_to_dealii_point<2>((*s)[1])}});
                }
            }
        }

      return vertices;
    }

    // specialization for hex \cap line
    std::vector<std::array<Point<3>, 2>>
    compute_intersection_hexa_line(const ArrayView<const Point<3>> &hexa,
                                   const ArrayView<const Point<3>> &line,
                                   const double                     tol)
    {
#  if DEAL_II_CGAL_VERSION_GTE(5, 5, 0)

      AssertDimension(hexa.size(), 8);
      AssertDimension(line.size(), 2);

      std::array<CGALPoint3_exact, 8> pts;

      std::transform(
        hexa.begin(),
        hexa.end(),
        pts.begin(),
        &CGALWrappers::dealii_point_to_cgal_point<CGALPoint3_exact, 3>);

      CGALSegment3_exact cgal_segment(
        CGALWrappers::dealii_point_to_cgal_point<CGALPoint3_exact>(line[0]),
        CGALWrappers::dealii_point_to_cgal_point<CGALPoint3_exact>(line[1]));

      // Subdivide the hex into tetrahedrons, and intersect each one of them
      // with the line
      std::vector<std::array<Point<3>, 2>> vertices;
      Triangulation3_exact                 cgal_triangulation;
      cgal_triangulation.insert(pts.begin(), pts.end());
      for (const auto &c : cgal_triangulation.finite_cell_handles())
        {
          const auto &cgal_tetrahedron = cgal_triangulation.tetrahedron(c);
          if (CGAL::do_intersect(cgal_segment, cgal_tetrahedron))
            {
              const auto intersection =
                CGAL::intersection(cgal_segment, cgal_tetrahedron);
              if (const CGALSegment3_exact *s =
                    get_if_<CGALSegment3_exact>(&*intersection))
                {
                  if (s->squared_length() > tol * tol)
                    {
                      vertices.push_back(
                        {{CGALWrappers::cgal_point_to_dealii_point<3>(
                            s->vertex(0)),
                          CGALWrappers::cgal_point_to_dealii_point<3>(
                            s->vertex(1))}});
                    }
                }
            }
        }
      return vertices;
#  else
      Assert(
        false,
        ExcMessage(
          "This function requires a version of CGAL greater or equal than 5.5."));
      (void)hexa;
      (void)line;
      (void)tol;
      return {};
#  endif
    }

    std::vector<std::array<Point<3>, 3>>
    compute_intersection_hexa_quad(const ArrayView<const Point<3>> &hexa,
                                   const ArrayView<const Point<3>> &quad,
                                   const double                     tol)
    {
#  if DEAL_II_CGAL_VERSION_GTE(5, 5, 0)

      AssertDimension(hexa.size(), 8);
      AssertDimension(quad.size(), 4);

      std::array<CGALPoint3_exact, 8> pts_hex;
      std::array<CGALPoint3_exact, 4> pts_quad;

      std::transform(
        hexa.begin(),
        hexa.end(),
        pts_hex.begin(),
        &CGALWrappers::dealii_point_to_cgal_point<CGALPoint3_exact, 3>);

      std::transform(
        quad.begin(),
        quad.end(),
        pts_quad.begin(),
        &CGALWrappers::dealii_point_to_cgal_point<CGALPoint3_exact, 3>);

      // Subdivide hex into tetrahedrons
      std::vector<std::array<Point<3>, 3>> vertices;
      Triangulation3_exact                 triangulation_hexa;
      triangulation_hexa.insert(pts_hex.begin(), pts_hex.end());

      // Subdivide quad into triangles
      Triangulation3_exact triangulation_quad;
      triangulation_quad.insert(pts_quad.begin(), pts_quad.end());

      for (const auto &c : triangulation_hexa.finite_cell_handles())
        {
          const auto &tet = triangulation_hexa.tetrahedron(c);

          for (const auto &f : triangulation_quad.finite_facets())
            {
              if (CGAL::do_intersect(tet, triangulation_quad.triangle(f)))
                {
                  const auto intersection =
                    CGAL::intersection(triangulation_quad.triangle(f), tet);

                  if (const CGALTriangle3_exact *t =
                        get_if_<CGALTriangle3_exact>(&*intersection))
                    {
                      if (CGAL::to_double(t->squared_area()) > tol * tol)
                        {
                          vertices.push_back(
                            {{cgal_point_to_dealii_point<3>((*t)[0]),
                              cgal_point_to_dealii_point<3>((*t)[1]),
                              cgal_point_to_dealii_point<3>((*t)[2])}});
                        }
                    }

                  if (const std::vector<CGALPoint3_exact> *vps =
                        get_if_<std::vector<CGALPoint3_exact>>(&*intersection))
                    {
                      Triangulation3_exact tria_inter;
                      tria_inter.insert(vps->begin(), vps->end());

                      for (auto it = tria_inter.finite_facets_begin();
                           it != tria_inter.finite_facets_end();
                           ++it)
                        {
                          const auto triangle = tria_inter.triangle(*it);
                          if (CGAL::to_double(triangle.squared_area()) >
                              tol * tol)
                            {
                              std::array<Point<3>, 3> verts = {
                                {CGALWrappers::cgal_point_to_dealii_point<3>(
                                   triangle[0]),
                                 CGALWrappers::cgal_point_to_dealii_point<3>(
                                   triangle[1]),
                                 CGALWrappers::cgal_point_to_dealii_point<3>(
                                   triangle[2])}};

                              vertices.push_back(verts);
                            }
                        }
                    }
                }
            }
        }

      return vertices;
#  else
      Assert(
        false,
        ExcMessage(
          "This function requires a version of CGAL greater or equal than 5.5."));
      (void)hexa;
      (void)quad;
      (void)tol;
      return {};
#  endif
    }

    std::vector<std::array<Point<3>, 4>>
    compute_intersection_hexa_hexa(const ArrayView<const Point<3>> &hexa0,
                                   const ArrayView<const Point<3>> &hexa1,
                                   const double                     tol)
    {
      AssertDimension(hexa0.size(), 8);
      AssertDimension(hexa0.size(), hexa1.size());

      std::array<CGALPoint3_exact, 8> pts_hex0;
      std::array<CGALPoint3_exact, 8> pts_hex1;

      std::transform(
        hexa0.begin(),
        hexa0.end(),
        pts_hex0.begin(),
        &CGALWrappers::dealii_point_to_cgal_point<CGALPoint3_exact, 3>);

      std::transform(
        hexa1.begin(),
        hexa1.end(),
        pts_hex1.begin(),
        &CGALWrappers::dealii_point_to_cgal_point<CGALPoint3_exact, 3>);

      Surface_mesh surf0, surf1, sm;
      // Subdivide hex into tetrahedrons
      std::vector<std::array<Point<3>, 4>> vertices;
      Triangulation3_exact                 tria0, tria1;

      tria0.insert(pts_hex0.begin(), pts_hex0.end());
      tria1.insert(pts_hex1.begin(), pts_hex1.end());

      for (const auto &c0 : tria0.finite_cell_handles())
        {
          const auto &tet0  = tria1.tetrahedron(c0);
          const auto &tetg0 = CGAL::make_tetrahedron(tet0.vertex(0),
                                                     tet0.vertex(1),
                                                     tet0.vertex(2),
                                                     tet0.vertex(3),
                                                     surf0);
          (void)tetg0; // instead of C++ 17s [[maybe unused]]
          for (const auto &c1 : tria1.finite_cell_handles())
            {
              const auto &tet1  = tria1.tetrahedron(c1);
              const auto &tetg1 = CGAL::make_tetrahedron(tet1.vertex(0),
                                                         tet1.vertex(1),
                                                         tet1.vertex(2),
                                                         tet1.vertex(3),
                                                         surf1);
              (void)tetg1; // instead of C++ 17s [[maybe unused]]
              namespace PMP = CGAL::Polygon_mesh_processing;
              const bool test_intersection =
                PMP::corefine_and_compute_intersection(surf0, surf1, sm);
              if (PMP::volume(sm) > tol && test_intersection)
                {
                  // Collect tetrahedrons
                  Triangulation3_exact triangulation_hexa;
                  triangulation_hexa.insert(sm.points().begin(),
                                            sm.points().end());
                  for (const auto &c : triangulation_hexa.finite_cell_handles())
                    {
                      const auto &tet = triangulation_hexa.tetrahedron(c);
                      vertices.push_back(
                        {{CGALWrappers::cgal_point_to_dealii_point<3>(
                            tet.vertex(0)),
                          CGALWrappers::cgal_point_to_dealii_point<3>(
                            tet.vertex(1)),
                          CGALWrappers::cgal_point_to_dealii_point<3>(
                            tet.vertex(2)),
                          CGALWrappers::cgal_point_to_dealii_point<3>(
                            tet.vertex(3))}});
                    }
                }
              surf1.clear();
              sm.clear();
            }
          surf0.clear();
        }
      return vertices;
    }

  } // namespace internal


  template <int structdim0, int structdim1, int spacedim>
  std::vector<std::array<Point<spacedim>, structdim1 + 1>>
  compute_intersection_of_cells(
    const ArrayView<const Point<spacedim>> &vertices0,
    const ArrayView<const Point<spacedim>> &vertices1,
    const double                            tol)
  {
    const unsigned int n_vertices0 = vertices0.size();
    const unsigned int n_vertices1 = vertices1.size();

    Assert(
      n_vertices0 > 0 || n_vertices1 > 0,
      ExcMessage(
        "The intersection cannot be computed as at least one of the two cells has no vertices."));

    if constexpr (structdim0 == 2 && structdim1 == 2 && spacedim == 2)
      {
        if (n_vertices0 == 4 && n_vertices1 == 4)
          {
            return internal::compute_intersection_quad_quad(vertices0,
                                                            vertices1,
                                                            tol);
          }
      }
    else if constexpr (structdim0 == 2 && structdim1 == 1 && spacedim == 2)
      {
        if (n_vertices0 == 4 && n_vertices1 == 2)
          {
            return internal::compute_intersection_quad_line(vertices0,
                                                            vertices1,
                                                            tol);
          }
      }
    else if constexpr (structdim0 == 3 && structdim1 == 1 && spacedim == 3)
      {
        if (n_vertices0 == 8 && n_vertices1 == 2)
          {
            return internal::compute_intersection_hexa_line(vertices0,
                                                            vertices1,
                                                            tol);
          }
      }
    else if constexpr (structdim0 == 3 && structdim1 == 2 && spacedim == 3)
      {
        if (n_vertices0 == 8 && n_vertices1 == 4)
          {
            return internal::compute_intersection_hexa_quad(vertices0,
                                                            vertices1,
                                                            tol);
          }
      }
    else if constexpr (structdim0 == 3 && structdim1 == 3 && spacedim == 3)
      {
        if (n_vertices0 == 8 && n_vertices1 == 8)
          {
            return internal::compute_intersection_hexa_hexa(vertices0,
                                                            vertices1,
                                                            tol);
          }
      }
    else
      {
        DEAL_II_NOT_IMPLEMENTED();
        return {};
      }
    (void)tol;
    return {};
  }


  template <int structdim0, int structdim1, int spacedim>
  std::vector<std::array<Point<spacedim>, structdim1 + 1>>
  compute_intersection_of_cells(
    const typename Triangulation<structdim0, spacedim>::cell_iterator &cell0,
    const typename Triangulation<structdim1, spacedim>::cell_iterator &cell1,
    const Mapping<structdim0, spacedim>                               &mapping0,
    const Mapping<structdim1, spacedim>                               &mapping1,
    const double                                                       tol)
  {
    Assert(mapping0.get_vertices(cell0).size() ==
             ReferenceCells::get_hypercube<structdim0>().n_vertices(),
           ExcNotImplemented());
    Assert(mapping1.get_vertices(cell1).size() ==
             ReferenceCells::get_hypercube<structdim1>().n_vertices(),
           ExcNotImplemented());

    const auto &vertices0 =
      CGALWrappers::get_vertices_in_cgal_order(cell0, mapping0);
    const auto &vertices1 =
      CGALWrappers::get_vertices_in_cgal_order(cell1, mapping1);

    return compute_intersection_of_cells<structdim0, structdim1, spacedim>(
      vertices0, vertices1, tol);
  }

// Explicit instantiations.
//
// We don't build the instantiations.inst file if deal.II isn't
// configured with CGAL, but doxygen doesn't know that and tries to
// find that file anyway for parsing -- which then of course it fails
// on. So exclude the following from doxygen consideration.
#  ifndef DOXYGEN
#    include "cgal/intersections.inst"
#  endif

} // namespace CGALWrappers

DEAL_II_NAMESPACE_CLOSE

#else

DEAL_II_NAMESPACE_OPEN

template <int structdim0,
          int structdim1,
          int spacedim,
          int n_components0,
          int n_components1>
std::vector<std::array<Point<spacedim>, structdim1 + 1>>
compute_intersection_of_cells(
  const std::array<Point<spacedim>, n_components0> &vertices0,
  const std::array<Point<spacedim>, n_components1> &vertices1,
  const double                                      tol)
{
  (void)vertices0;
  (void)vertices1;
  (void)tol;
  AssertThrow(false, ExcNeedsCGAL());
}

template <int structdim0, int structdim1, int spacedim>
std::vector<std::array<Point<spacedim>, structdim1 + 1>>
compute_intersection_of_cells(
  const typename Triangulation<structdim0, spacedim>::cell_iterator &cell0,
  const typename Triangulation<structdim1, spacedim>::cell_iterator &cell1,
  const Mapping<structdim0, spacedim>                               &mapping0,
  const Mapping<structdim1, spacedim>                               &mapping1,
  const double                                                       tol)
{
  (void)cell0;
  (void)cell1;
  (void)mapping0;
  (void)mapping1;
  (void)tol;
  AssertThrow(false, ExcNeedsCGAL());
}

DEAL_II_NAMESPACE_CLOSE

#endif

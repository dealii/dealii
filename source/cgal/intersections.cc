// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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

#include <deal.II/base/config.h>

#include <deal.II/cgal/intersections.h>

#include <algorithm>

#ifdef DEAL_II_WITH_CGAL

#  include <deal.II/base/quadrature_lib.h>

#  include <deal.II/fe/mapping.h>

#  include <deal.II/grid/tria.h>

#  include <CGAL/Boolean_set_operations_2.h>
#  include <CGAL/Cartesian.h>
#  include <CGAL/Circular_kernel_intersections.h>
#  include <CGAL/Constrained_Delaunay_triangulation_2.h>
#  include <CGAL/Delaunay_mesh_face_base_2.h>
#  include <CGAL/Delaunay_mesh_size_criteria_2.h>
#  include <CGAL/Delaunay_mesher_2.h>
#  include <CGAL/Delaunay_triangulation_2.h>
#  include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#  include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#  include <CGAL/Kernel_traits.h>
#  include <CGAL/Polygon_2.h>
#  include <CGAL/Polygon_with_holes_2.h>
#  include <CGAL/Projection_traits_xy_3.h>
#  include <CGAL/Segment_3.h>
#  include <CGAL/Simple_cartesian.h>
#  include <CGAL/Tetrahedron_3.h>
#  include <CGAL/Triangle_2.h>
#  include <CGAL/Triangle_3.h>
#  include <CGAL/Triangulation_2.h>
#  include <CGAL/Triangulation_3.h>
#  include <CGAL/Triangulation_face_base_with_id_2.h>
#  include <CGAL/Triangulation_face_base_with_info_2.h>
#  include <deal.II/cgal/utilities.h>

#  include <fstream>
#  include <type_traits>

DEAL_II_NAMESPACE_OPEN

namespace CGALWrappers
{
  using K         = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;
  using K_exact   = CGAL::Exact_predicates_exact_constructions_kernel;
  using K_inexact = CGAL::Exact_predicates_inexact_constructions_kernel;
  using CGALPolygon            = CGAL::Polygon_2<K>;
  using Polygon_with_holes_2   = CGAL::Polygon_with_holes_2<K>;
  using CGALTriangle2          = K::Triangle_2;
  using CGALTriangle3          = K::Triangle_3;
  using CGALTriangle3_exact    = K_exact::Triangle_3;
  using CGALPoint2             = K::Point_2;
  using CGALPoint3             = K::Point_3;
  using CGALPoint3_exact       = K_exact::Point_3;
  using CGALPoint3_inexact     = K_inexact::Point_3;
  using CGALSegment2           = K::Segment_2;
  using Surface_mesh           = CGAL::Surface_mesh<K_inexact::Point_3>;
  using CGALSegment3           = K::Segment_3;
  using CGALSegment3_exact     = K_exact::Segment_3;
  using CGALTetra              = K::Tetrahedron_3;
  using CGALTetra_exact        = K_exact::Tetrahedron_3;
  using Triangulation2         = CGAL::Triangulation_2<K>;
  using Triangulation3         = CGAL::Triangulation_3<K>;
  using Triangulation3_exact   = CGAL::Triangulation_3<K_exact>;
  using Triangulation3_inexact = CGAL::Triangulation_3<K_inexact>;

  struct FaceInfo2
  {
    FaceInfo2()
    {}
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

  namespace internal
  {
    void
    mark_domains(CDT &                 ct,
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
    // CGAL::intersection(), i.e. a boost::optional<boost::variant<>>.
    // Intersection between 2D and 3D objects and 1D/3D objects are available
    // only with CGAL versions greater or equal than 5.1.5, hence the
    // corresponding functions are guarded by #ifdef directives. All the
    // signatures follow the convection that the first entity has an intrinsic
    // dimension higher than the second one.

    boost::optional<boost::variant<CGALPoint2,
                                   CGALSegment2,
                                   CGALTriangle2,
                                   std::vector<CGALPoint2>>>
    compute_intersection(const std::array<Point<2>, 3> &first_simplex,
                         const std::array<Point<2>, 3> &second_simplex)
    {
      std::array<CGALPoint2, 3> pts0, pts1;
      std::transform(
        first_simplex.begin(),
        first_simplex.end(),
        pts0.begin(),
        [&](const Point<2> &p) {
          return CGALWrappers::dealii_point_to_cgal_point<CGALPoint2>(p);
        });

      std::transform(
        second_simplex.begin(),
        second_simplex.end(),
        pts1.begin(),
        [&](const Point<2> &p) {
          return CGALWrappers::dealii_point_to_cgal_point<CGALPoint2>(p);
        });

      CGALTriangle2 triangle1{pts0[0], pts0[1], pts0[2]};
      CGALTriangle2 triangle2{pts1[0], pts1[1], pts1[2]};
      return CGAL::intersection(triangle1, triangle2);
    }



    boost::optional<boost::variant<CGALPoint2, CGALSegment2>>
    compute_intersection(const std::array<Point<2>, 3> &first_simplex,
                         const std::array<Point<2>, 2> &second_simplex)
    {
      std::array<CGALPoint2, 3> pts0;
      std::array<CGALPoint2, 2> pts1;
      std::transform(
        first_simplex.begin(),
        first_simplex.end(),
        pts0.begin(),
        [&](const Point<2> &p) {
          return CGALWrappers::dealii_point_to_cgal_point<CGALPoint2>(p);
        });

      std::transform(
        second_simplex.begin(),
        second_simplex.end(),
        pts1.begin(),
        [&](const Point<2> &p) {
          return CGALWrappers::dealii_point_to_cgal_point<CGALPoint2>(p);
        });

      CGALTriangle2 triangle{pts0[0], pts0[1], pts0[2]};
      CGALSegment2  segm{pts1[0], pts1[1]};
      return CGAL::intersection(segm, triangle);
    }



    // rectangle-rectangle
    std::vector<Polygon_with_holes_2>
    compute_intersection(const std::array<Point<2>, 4> &first_simplex,
                         const std::array<Point<2>, 4> &second_simplex)
    {
      std::array<CGALPoint2, 4> pts0, pts1;
      std::transform(
        first_simplex.begin(),
        first_simplex.end(),
        pts0.begin(),
        [&](const Point<2> &p) {
          return CGALWrappers::dealii_point_to_cgal_point<CGALPoint2>(p);
        });
      std::transform(
        second_simplex.begin(),
        second_simplex.end(),
        pts1.begin(),
        [&](const Point<2> &p) {
          return CGALWrappers::dealii_point_to_cgal_point<CGALPoint2>(p);
        });
      const CGALPolygon first{pts0.begin(), pts0.end()};
      const CGALPolygon second{pts1.begin(), pts1.end()};

      std::vector<Polygon_with_holes_2> poly_list;
      CGAL::intersection(first, second, std::back_inserter(poly_list));
      return poly_list;
    }



    boost::optional<boost::variant<CGALPoint3, CGALSegment3>>
    compute_intersection(const std::array<Point<3>, 2> &first_simplex,
                         const std::array<Point<3>, 4> &second_simplex)
    {
#  if DEAL_II_CGAL_VERSION_GTE(5, 1, 5)
      std::array<CGALPoint3, 4> pts0;
      std::array<CGALPoint3, 2> pts1;
      std::transform(
        first_simplex.begin(),
        first_simplex.end(),
        pts0.begin(),
        [&](const Point<3> &p) {
          return CGALWrappers::dealii_point_to_cgal_point<CGALPoint3>(p);
        });


      std::transform(
        second_simplex.begin(),
        second_simplex.end(),
        pts1.begin(),
        [&](const Point<3> &p) {
          return CGALWrappers::dealii_point_to_cgal_point<CGALPoint3>(p);
        });

      CGALTetra    tetra{pts0[0], pts0[1], pts0[2], pts0[3]};
      CGALSegment3 segm{pts1[0], pts1[1]};
      return CGAL::intersection(segm, tetra);
#  else
      Assert(
        false,
        ExcMessage(
          "This function requires a version of CGAL greater or equal than 5.1.5."));
      (void)first_simplex;
      (void)second_simplex;
      return {};
#  endif
    }

    // tetra, triangle
    boost::optional<boost::variant<CGALPoint3,
                                   CGALSegment3,
                                   CGALTriangle3,
                                   std::vector<CGALPoint3>>>
    compute_intersection(const std::array<Point<3>, 3> &first_simplex,
                         const std::array<Point<3>, 4> &second_simplex)
    {
#  if DEAL_II_CGAL_VERSION_GTE(5, 1, 5)
      std::array<CGALPoint3, 4> pts0;
      std::array<CGALPoint3, 3> pts1;
      std::transform(
        first_simplex.begin(),
        first_simplex.end(),
        pts0.begin(),
        [&](const Point<3> &p) {
          return CGALWrappers::dealii_point_to_cgal_point<CGALPoint3>(p);
        });

      std::transform(
        second_simplex.begin(),
        second_simplex.end(),
        pts1.begin(),
        [&](const Point<3> &p) {
          return CGALWrappers::dealii_point_to_cgal_point<CGALPoint3>(p);
        });
      CGALTetra     tetra{pts0[0], pts0[1], pts0[2], pts0[3]};
      CGALTriangle3 triangle{pts1[0], pts1[1], pts1[2]};
      return CGAL::intersection(triangle, tetra);
#  else

      Assert(
        false,
        ExcMessage(
          "This function requires a version of CGAL greater or equal than 5.1.5."));
      (void)first_simplex;
      (void)second_simplex;
      return {};
#  endif
    }
  } // namespace internal



  // Specialization for quads
  template <>
  std::vector<std::array<Point<2>, 3>>
  compute_intersection_of_cells<2, 2, 2, 4, 4>(
    const std::array<Point<2>, 4> &vertices0,
    const std::array<Point<2>, 4> &vertices1,
    const double                   tol)
  {
    const auto intersection_test =
      internal::compute_intersection(vertices0, vertices1);

    if (!intersection_test.empty())
      {
        const auto &       poly      = intersection_test[0].outer_boundary();
        const unsigned int size_poly = poly.size();
        if (size_poly == 3)
          {
            // intersection is a triangle itself, so directly return its
            // vertices.
            return {
              {{CGALWrappers::cgal_point_to_dealii_point<2>(poly.vertex(0)),
                CGALWrappers::cgal_point_to_dealii_point<2>(poly.vertex(1)),
                CGALWrappers::cgal_point_to_dealii_point<2>(poly.vertex(2))}}};
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
            std::array<Point<2>, 3> vertices;

            for (const Face_handle &f : cdt.finite_face_handles())
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
  template <>
  std::vector<std::array<Point<2>, 2>>
  compute_intersection_of_cells<2, 1, 2, 4, 2>(
    const std::array<Point<2>, 4> &vertices0,
    const std::array<Point<2>, 2> &vertices1,
    const double                   tol)
  {
    std::array<CGALPoint2, 4> pts;
    std::transform(
      vertices0.begin(), vertices0.end(), pts.begin(), [&](const Point<2> &p) {
        return CGALWrappers::dealii_point_to_cgal_point<CGALPoint2>(p);
      });

    CGALPolygon poly(pts.begin(), pts.end());

    CGALSegment2 segm(
      CGALWrappers::dealii_point_to_cgal_point<CGALPoint2>(vertices1[0]),
      CGALWrappers::dealii_point_to_cgal_point<CGALPoint2>(vertices1[1]));
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
            const auto intersection = CGAL::intersection(segm, cdt.triangle(f));
            if (const CGALSegment2 *s =
                  boost::get<CGALSegment2>(&*intersection))
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
  template <>
  std::vector<std::array<Point<3>, 2>>
  compute_intersection_of_cells<3, 1, 3, 8, 2>(
    const std::array<Point<3>, 8> &vertices0,
    const std::array<Point<3>, 2> &vertices1,
    const double                   tol)
  {
#  if DEAL_II_CGAL_VERSION_GTE(5, 1, 5)
    std::array<CGALPoint3_exact, 8> pts;
    std::transform(
      vertices0.begin(), vertices0.end(), pts.begin(), [&](const Point<3> &p) {
        return CGALWrappers::dealii_point_to_cgal_point<CGALPoint3_exact>(p);
      });

    CGALSegment3_exact segm(
      CGALWrappers::dealii_point_to_cgal_point<CGALPoint3_exact>(vertices1[0]),
      CGALWrappers::dealii_point_to_cgal_point<CGALPoint3_exact>(vertices1[1]));

    // Subdivide the hex into tetrahedrons, and intersect each one of them with
    // the line
    std::vector<std::array<Point<3>, 2>> vertices;
    Triangulation3_exact                 tria;
    tria.insert(pts.begin(), pts.end());
    for (const auto &c : tria.finite_cell_handles())
      {
        const auto &tet = tria.tetrahedron(c);
        if (CGAL::do_intersect(segm, tet))
          {
            const auto intersection = CGAL::intersection(segm, tet);
            if (const CGALSegment3_exact *s =
                  boost::get<CGALSegment3_exact>(&*intersection))
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
        "This function requires a version of CGAL greater or equal than 5.1.5."));
    (void)vertices0;
    (void)vertices1;
    (void)tol;
    return {};
#  endif
  }

  template <>
  std::vector<std::array<Point<3>, 3>>
  compute_intersection_of_cells<3, 2, 3, 8, 4>(
    const std::array<Point<3>, 8> &vertices0,
    const std::array<Point<3>, 4> &vertices1,
    const double                   tol)
  {
#  if DEAL_II_CGAL_VERSION_GTE(5, 1, 5)
    std::array<CGALPoint3_exact, 8> pts_hex;
    std::array<CGALPoint3_exact, 4> pts_quad;
    std::transform(
      vertices0.begin(),
      vertices0.end(),
      pts_hex.begin(),
      [&](const Point<3> &p) {
        return CGALWrappers::dealii_point_to_cgal_point<CGALPoint3_exact>(p);
      });

    std::transform(
      vertices1.begin(),
      vertices1.end(),
      pts_quad.begin(),
      [&](const Point<3> &p) {
        return CGALWrappers::dealii_point_to_cgal_point<CGALPoint3_exact>(p);
      });

    // Subdivide hex into tetrahedrons
    std::vector<std::array<Point<3>, 3>> vertices;
    Triangulation3_exact                 tria;
    tria.insert(pts_hex.begin(), pts_hex.end());

    // Subdivide quad into triangles
    Triangulation3_exact tria_quad;
    tria_quad.insert(pts_quad.begin(), pts_quad.end());

    for (const auto &c : tria.finite_cell_handles())
      {
        const auto &tet = tria.tetrahedron(c);

        for (const auto &f : tria_quad.finite_facets())
          {
            if (CGAL::do_intersect(tet, tria_quad.triangle(f)))
              {
                const auto intersection =
                  CGAL::intersection(tria_quad.triangle(f), tet);

                if (const CGALTriangle3_exact *t =
                      boost::get<CGALTriangle3_exact>(&*intersection))
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
                      boost::get<std::vector<CGALPoint3_exact>>(&*intersection))
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
        "This function requires a version of CGAL greater or equal than 5.1.5."));
    (void)vertices0;
    (void)vertices1;
    (void)tol;
    return {};
#  endif
  }



  template <>
  std::vector<std::array<Point<3>, 4>>
  compute_intersection_of_cells<3, 3, 3, 8, 8>(
    const std::array<Point<3>, 8> &vertices0,
    const std::array<Point<3>, 8> &vertices1,

    const double tol)
  {
    std::array<CGALPoint3_inexact, 8> pts_hex0;
    std::array<CGALPoint3_inexact, 8> pts_hex1;
    std::transform(
      vertices0.begin(),
      vertices0.end(),
      pts_hex0.begin(),
      [&](const Point<3> &p) {
        return CGALWrappers::dealii_point_to_cgal_point<CGALPoint3_inexact>(p);
      });

    std::transform(
      vertices1.begin(),
      vertices1.end(),
      pts_hex1.begin(),
      [&](const Point<3> &p) {
        return CGALWrappers::dealii_point_to_cgal_point<CGALPoint3_inexact>(p);
      });


    Surface_mesh surf0, surf1, sm;
    // Subdivide hex into tetrahedrons
    std::vector<std::array<Point<3>, 4>> vertices;
    Triangulation3_inexact               tria0, tria1;

    tria0.insert(pts_hex0.begin(), pts_hex0.end());
    tria1.insert(pts_hex1.begin(), pts_hex1.end());

    for (const auto &c0 : tria0.finite_cell_handles())
      {
        const auto &                 tet0 = tria1.tetrahedron(c0);
        [[maybe_unused]] const auto &tetg0 =
          CGAL::make_tetrahedron(tet0.vertex(0),
                                 tet0.vertex(1),
                                 tet0.vertex(2),
                                 tet0.vertex(3),
                                 surf0);
        for (const auto &c1 : tria1.finite_cell_handles())
          {
            const auto &                 tet1 = tria1.tetrahedron(c1);
            [[maybe_unused]] const auto &tetg1 =
              CGAL::make_tetrahedron(tet1.vertex(0),
                                     tet1.vertex(1),
                                     tet1.vertex(2),
                                     tet1.vertex(3),
                                     surf1);
            const bool test_intersection =
              CGAL::PMP::corefine_and_compute_intersection(surf0, surf1, sm);
            if (CGAL::Polygon_mesh_processing::volume(sm) > tol &&
                test_intersection)
              {
                // Collect tetrahedrons
                Triangulation3_inexact tria;
                tria.insert(sm.points().begin(), sm.points().end());
                for (const auto &c : tria.finite_cell_handles())
                  {
                    const auto &tet = tria.tetrahedron(c);
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



  template <int dim0, int dim1, int spacedim>
  std::vector<std::array<Point<spacedim>, dim1 + 1>>
  compute_intersection_of_cells(
    const typename Triangulation<dim0, spacedim>::cell_iterator &cell0,
    const typename Triangulation<dim1, spacedim>::cell_iterator &cell1,
    const Mapping<dim0, spacedim> &                              mapping0,
    const Mapping<dim1, spacedim> &                              mapping1,
    const double                                                 tol)
  {
    Assert(mapping0.get_vertices(cell0).size() == std::pow(2, dim0),
           ExcNotImplemented());
    Assert(mapping1.get_vertices(cell1).size() == std::pow(2, dim1),
           ExcNotImplemented());

    const auto vertices0 =
      CGALWrappers::get_vertices_in_cgal_order<int(std::pow(2, dim0))>(
        cell0, mapping0);
    const auto vertices1 =
      CGALWrappers::get_vertices_in_cgal_order<int(std::pow(2, dim1))>(
        cell1, mapping1);

    return compute_intersection_of_cells<dim0,
                                         dim1,
                                         spacedim,
                                         int(std::pow(2, dim0)),
                                         int(std::pow(2, dim1))>(vertices0,
                                                                 vertices1,
                                                                 tol);
  }

#  include "intersections.inst"

} // namespace CGALWrappers

DEAL_II_NAMESPACE_CLOSE

#else

DEAL_II_NAMESPACE_OPEN

template <int dim0,
          int dim1,
          int spacedim,
          int n_components0,
          int n_components1>
std::vector<std::array<Point<spacedim>, dim1 + 1>>
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

template <int dim0, int dim1, int spacedim>
std::vector<std::array<Point<spacedim>, dim1 + 1>>
compute_intersection_of_cells(
  const typename Triangulation<dim0, spacedim>::cell_iterator &cell0,
  const typename Triangulation<dim1, spacedim>::cell_iterator &cell1,
  const Mapping<dim0, spacedim> &                              mapping0,
  const Mapping<dim1, spacedim> &                              mapping1,
  const double                                                 tol)
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

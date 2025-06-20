// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Convert a deal.II triangulation to a CGAL polygon_2.
#include <deal.II/base/config.h>

#include <deal.II/cgal/polygon.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include "../tests.h"

using namespace CGALWrappers;

using K           = CGAL::Exact_predicates_exact_constructions_kernel;
using CGALPoint2  = CGAL::Point_2<K>;
using CGALPolygon = CGAL::Polygon_2<K>;

void
test_quadrilaterals()
{
  Triangulation<2, 2> tria_in;
  CGALPolygon         poly;

  std::vector<std::pair<std::string, std::string>> names_and_args;

  // only for meshes that have no holes
  // extension to CGAL::Polygon_with_holes possible
  names_and_args = {{"hyper_cube", "0.0 : 1.0 : false"},
                    {"hyper_ball_balanced", "0.,0. : 1. "},
                    {"hyper_L", "0.0 : 1.0 : false"},
                    {"simplex", "0.0, 0.0 ; 1.0 , 0.0 ; 0.0, 1.0"}};


  for (const auto &info_pair : names_and_args)
    {
      auto name = info_pair.first;
      auto args = info_pair.second;
      deallog << "name: " << name << std::endl;
      GridGenerator::generate_from_name_and_arguments(tria_in, name, args);
      tria_in.refine_global(2);
      dealii_tria_to_cgal_polygon(tria_in, poly);

      deallog << "Simple polygon: " << std::boolalpha << poly.is_simple()
              << std::endl;
      deallog << "Counterclockwise oriented polygon: " << std::boolalpha
              << poly.is_counterclockwise_oriented() << std::endl;
      deallog << "deal measure: " << GridTools::volume(tria_in)
              << ", cgal measure: " << poly.area() << std::endl;

      // every vertex of the tria should be inside or on the boundary
      const auto &vertices = tria_in.get_vertices();
      for (unsigned int i = 0; i < vertices.size(); i++)
        {
          auto bounded_side = poly.bounded_side(
            dealii_point_to_cgal_point<CGALPoint2, 2>(vertices[i]));
          Assert(
            bounded_side == CGAL::ON_BOUNDED_SIDE ||
              bounded_side == CGAL::ON_BOUNDARY,
            ExcMessage(
              "The polygon area doesn't contain all triangulation points"));
        }

      tria_in.clear();
      poly.clear();
      deallog << std::endl;
    }
}



void
test_triangles()
{
  Triangulation<2, 2> tria_quad;
  Triangulation<2, 2> tria_simplex;
  CGALPolygon         poly;

  deallog << "name: hyper_cube converted to simplex mesh" << std::endl;
  GridGenerator::hyper_cube(tria_quad, 0., 1., false);
  tria_quad.refine_global(2);
  GridGenerator::convert_hypercube_to_simplex_mesh(tria_quad, tria_simplex);
  dealii_tria_to_cgal_polygon(tria_simplex, poly);

  deallog << "Simple polygon: " << std::boolalpha << poly.is_simple()
          << std::endl;
  deallog << "Counterclockwise oriented polygon: " << std::boolalpha
          << poly.is_counterclockwise_oriented() << std::endl;
  deallog << "deal measure: " << GridTools::volume(tria_simplex)
          << ", cgal measure: " << poly.area() << std::endl;

  // every vertex of the tria should be inside or on the boundary
  const auto &vertices = tria_simplex.get_vertices();
  for (unsigned int i = 0; i < vertices.size(); i++)
    {
      auto bounded_side = poly.bounded_side(
        dealii_point_to_cgal_point<CGALPoint2, 2>(vertices[i]));
      Assert(bounded_side == CGAL::ON_BOUNDED_SIDE ||
               bounded_side == CGAL::ON_BOUNDARY,
             ExcMessage(
               "The polygon area doesn't contain all triangulation points"));
    }

  tria_quad.clear();
  tria_simplex.clear();
  poly.clear();
}



int
main()
{
  initlog();
  test_quadrilaterals();
  test_triangles();
}

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

// Boolean operations on polygons
#include <deal.II/base/config.h>

#include <deal.II/cgal/polygon.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include "../tests.h"

using namespace CGALWrappers;

using K = CGAL::Exact_predicates_exact_constructions_kernel;

using CGALPolygonWithHoles = CGAL::Polygon_with_holes_2<K>;

void
write_volumes(std::vector<CGALPolygonWithHoles> &poly_vec)
{
  deallog << "With volumes: ";
  for (const auto &polygon : poly_vec)
    {
      if (polygon.has_holes())
        {
          for (const auto &hole : polygon.holes())
            {
              deallog << hole.area() << " ";
            }
        }
      deallog << polygon.outer_boundary().area() << " ";
    }
  deallog << std::endl;
}



void
test(unsigned int refinement_1, unsigned int refinement_2)
{
  Triangulation<2, 2> tria_1;
  Triangulation<2, 2> tria_2;
  MappingQ<2, 2>      mapping(1);

  std::vector<std::pair<std::string, std::string>> names_and_args;

  // only for operations on meshes that have no holes
  // extension to CGAL::Polygon_with_holes possible
  names_and_args = {{"hyper_cube", "0.0 : 1.0 : false"},
                    {"hyper_cube", "-0.1 : 0.9 : false"},
                    {"hyper_rectangle", "0.0, -0.1 : 1.0 , 0.5 : false"},
                    {"hyper_rectangle", "-0.1, -0.1 : 0.5 , 1.1 : false"},
                    {"hyper_rectangle", "-0.1, 0.5 : 0.0 , 1.1 : false"},
                    {"simplex", "0.0, 0.0 ; 1.0 , 0.0 ; 0.0, 1.0"},
                    {"simplex", "-0.1, -0.1 ; 0.9 , -0.1 ; -0.1, 0.9"},
                    {"simplex", "-0.5, -0.5 ; 0.5 , 0.5 ; -0.5, 1.5"},
                    {"hyper_ball_balanced", "0.0,0.0 : 1.0"},
                    {"hyper_ball_balanced", "-0.1 ,1.0 : 0.5"}};

  GridGenerator::generate_from_name_and_arguments(tria_1,
                                                  names_and_args[0].first,
                                                  names_and_args[0].second);
  tria_1.refine_global(refinement_1);
  auto poly_1 = dealii_tria_to_cgal_polygon<K>(tria_1, mapping);

  for (const auto &info_pair : names_and_args)
    {
      auto name = info_pair.first;
      auto args = info_pair.second;
      deallog << "name: " << name << std::endl;
      GridGenerator::generate_from_name_and_arguments(tria_2, name, args);
      tria_2.refine_global(refinement_2);
      auto poly_2 = dealii_tria_to_cgal_polygon<K>(tria_2, mapping);

      auto poly_out_vec =
        compute_boolean_operation(poly_1,
                                  poly_2,
                                  BooleanOperation::compute_intersection);
      deallog << "Intersection: " << poly_out_vec.size() << " polygons"
              << std::endl;
      write_volumes(poly_out_vec);

      poly_out_vec =
        compute_boolean_operation(poly_1,
                                  poly_2,
                                  BooleanOperation::compute_difference);
      deallog << "Difference: " << poly_out_vec.size() << " polygons"
              << std::endl;
      write_volumes(poly_out_vec);


      poly_out_vec = compute_boolean_operation(poly_1,
                                               poly_2,
                                               BooleanOperation::compute_union);
      deallog << "Union: " << poly_out_vec.size() << " polygons" << std::endl;
      write_volumes(poly_out_vec);

      tria_2.clear();
      poly_2.clear();
      deallog << std::endl;
    }
}


int
main()
{
  initlog();
  test(0, 0);
  test(2, 2);
}

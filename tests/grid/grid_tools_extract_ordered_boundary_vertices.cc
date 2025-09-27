// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check extract used_vertices.

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


void
test()
{
  deallog << "dim: " << 2 << ", spacedim: " << 2 << std::endl;

  Triangulation<2, 2> tria_in;
  MappingQ<2, 2>      mapping(1);

  std::vector<std::pair<std::string, std::string>> names_and_args;

  names_and_args = {{"hyper_cube", "0.0 : 1.0 : false"},
                    {"hyper_ball_balanced", "0.,0. : 1. "},
                    {"simplex", "0.0, 0.0 ; 1.0 , 0.0 ; 0.0, 1.0"},
                    {"channel_with_cylinder", "0.03 : 2 : 2.0 : false"},
                    {"cheese", "2 , 2"}};

  for (const auto &[name, args] : names_and_args)
    {
      deallog << "Name: " << name << std::endl;
      GridGenerator::generate_from_name_and_arguments(tria_in, name, args);

      auto boundaries =
        GridTools::extract_ordered_boundary_vertices(tria_in, mapping);

      unsigned int count = 0;
      for (const auto &boundary : boundaries)
        {
          deallog << "Closed boundary number: " << count++ << std::endl;
          for (const auto &v : boundary)
            {
              deallog << "Vertex: " << v.first << ": " << v.second << std::endl;
            }
        }

      tria_in.clear();
      deallog << std::endl;
    }
};


int
main()
{
  initlog();
  test();
  return 0;
}

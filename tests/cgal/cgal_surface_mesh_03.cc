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

// Perform corefinement and boolean operations between surface_meshes.

#include <deal.II/base/config.h>

#include <deal.II/cgal/utilities.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <CGAL/IO/io.h>

#include "../tests.h"


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
using CGALPoint = CGAL::Point_3<K>;
using namespace CGALWrappers;
void
test()
{
  const std::vector<std::string> fnames{SOURCE_DIR "/input_grids/cube.off",
                                        SOURCE_DIR "/input_grids/hedra.off"};
  CGAL::Surface_mesh<CGALPoint>  sm0, sm1, outsm;
  std::ifstream                  input0(fnames[0]);
  std::ifstream                  input1(fnames[1]);
  input0 >> sm0;
  input1 >> sm1;
  CGAL::Polygon_mesh_processing::triangulate_faces(sm0);
  CGAL::Polygon_mesh_processing::triangulate_faces(sm1);
  const std::vector<BooleanOperation> bool_operations{
    BooleanOperation::compute_union,
    BooleanOperation::compute_intersection,
    BooleanOperation::compute_difference,
    BooleanOperation::compute_corefinement};

  for (const auto &bool_op : bool_operations)
    {
      compute_boolean_operation(sm0, sm1, bool_op, outsm);
      Assert((outsm.is_valid() && sm0.is_valid() && sm1.is_valid()),
             ExcMessage("Result is not valid."));
      outsm.clear(); // reset surface
    }
}

int
main()
{
  initlog();
  test();
}

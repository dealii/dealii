// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check that the number of total children of a given level is correct.

#include <deal.II/base/patterns.h>

#include <deal.II/boost_adaptors/bounding_box.h>
#include <deal.II/boost_adaptors/point.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/rtree.h>

#include <algorithm>

#include "../tests.h"

int
main()
{
  initlog(true);

  Triangulation<2> tria;
  MappingQ1<2>     mapping;
  GridGenerator::hyper_ball(tria);
  tria.refine_global(6);

  namespace bgi = boost::geometry::index;
  std::vector<
    std::pair<BoundingBox<2>, typename Triangulation<2>::active_cell_iterator>>
               boxes(tria.n_active_cells());
  unsigned int i = 0;
  for (const auto &cell : tria.active_cell_iterators())
    boxes[i++] = std::make_pair(mapping.get_bounding_box(cell), cell);

  const auto tree = pack_rtree<bgi::linear<16>>(boxes);
  for (const unsigned int level : {0, 1, 2})
    {
      const auto bboxes = extract_rtree_level(tree, level + 1);
      deallog << "LEVEL " + std::to_string(level + 1) + "  N boxes: "
              << bboxes.size() << std::endl;

      const auto   boxes_in_boxes = extract_children_of_level(tree, level);
      unsigned int total_bboxes   = 0;
      for (unsigned int i = 0; i < boxes_in_boxes.size(); ++i)
        total_bboxes += boxes_in_boxes[i].size();

      if (level == 2)
        {
          Assert(boxes_in_boxes.size() == 0,
                 ExcMessage("Leaves have no children."));
        }
      else
        {
          Assert(total_bboxes == bboxes.size(),
                 ExcMessage(
                   "The number of total children of level " +
                   std::to_string(level) +
                   " should be equal to the number of boxes on level " +
                   std::to_string(level + 1)));
        }

      deallog << "OK" << std::endl;
    }

  // for (unsigned int i = 0; i < boxes_in_boxes.size(); ++i)
  //   for (const auto &b : boxes_in_boxes[i])
  //     deallog << "Box: " <<
  //     Patterns::Tools::to_string(b.get_boundary_points())
  //             << std::endl;
}

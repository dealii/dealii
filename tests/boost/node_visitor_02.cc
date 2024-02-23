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

// Use information from NodeVisitor to print show explicitly the boxes
// associated to each parent node on the previous level

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/bounding_box_data_out.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <algorithm>

#include "../tests.h"

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;



template <int dim, int spacedim, unsigned int max_elem_per_node>
void
test(const unsigned int ref = 6, const unsigned int level = 0)
{
  Triangulation<dim, spacedim> tria;
  MappingQ<dim>                mapping(1);
  GridGenerator::hyper_ball(tria);
  tria.refine_global(ref);

  std::vector<
    std::pair<BoundingBox<spacedim>,
              typename Triangulation<dim, spacedim>::active_cell_iterator>>
               boxes(tria.n_active_cells());
  unsigned int i = 0;
  for (const auto &cell : tria.active_cell_iterators())
    boxes[i++] = std::make_pair(mapping.get_bounding_box(cell), cell);

  const auto tree = pack_rtree<bgi::linear<max_elem_per_node>>(boxes);
  Assert(n_levels(tree) == 2,
         ExcMessage("Two levels are needed for this test."));
  const auto boxes_in_level = extract_rtree_level(tree, level);
  deallog << "N Boxes in level " + std::to_string(level) + " "
          << boxes_in_level.size() << std::endl;
  const auto boxes_in_boxes = extract_children_of_level(tree, level);

  for (unsigned int i = 0; i < boxes_in_boxes.size(); ++i)
    deallog << "boxes_in_boxes[" << i
            << "] has size: " << boxes_in_boxes[i].size() << std::endl;

  // Uncomment to display the partition
  // {
  //   i = 0;
  //   std::vector<
  //     std::pair<BoundingBox<spacedim>,
  //               typename Triangulation<dim, spacedim>::active_cell_iterator>>
  //     cells;

  //   for (unsigned int k = 0; k < boxes_in_boxes.size(); ++k)
  //     {
  //       deallog << "boxes_in_boxes[" << k << "]=" << boxes_in_boxes[k].size()
  //               << std::endl;
  //       for (const auto &box : boxes_in_boxes[k])
  //         {
  //           tree.query(bgi::within(box), std::back_inserter(cells));
  //           deallog << "Number of cells inside " << std::to_string(k)
  //                   << "-th bounding box: " << cells.size() << std::endl;
  //           for (const auto &my_pair : cells)
  //             {
  //               deallog << my_pair.second->active_cell_index() << std::endl;
  //               my_pair.second->set_subdomain_id(k);
  //             }

  //           cells.clear();
  //         }
  //     }


  //   for (unsigned int j = 0; j < i; ++j)
  //     deallog << GridTools::count_cells_with_subdomain_association(tria, j)
  //             << " cells have subdomain " + std::to_string(j) << std::endl;
  //   GridOut           grid_out_svg;
  //   GridOutFlags::Svg svg_flags;
  //   svg_flags.label_subdomain_id = true;
  //   svg_flags.coloring           = GridOutFlags::Svg::subdomain_id;
  //   grid_out_svg.set_flags(svg_flags);
  //   std::ofstream out("grid_ball_agglomerates" +
  //                     std::to_string(level) + ".svg");
  //   grid_out_svg.write_svg(tria, out);
  // }

  // Print hierarchy using boost utilities
  // bgi::detail::rtree::utilities::print(std::cout, tree);
}

int
main()
{
  initlog(true);

  static constexpr unsigned int global_refinements = 4; // 6
  static constexpr unsigned int max_per_node       = 16;
  const unsigned int            level              = 0; // 1
  test<2, 2, max_per_node>(global_refinements, level);
}

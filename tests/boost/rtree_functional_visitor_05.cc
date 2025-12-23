/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2025 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

// Visit an rtree of indices, and print out its structure. All lambdas are
// optional. We specify all 4 of them here.

#include <deal.II/base/patterns.h>
#include <deal.II/base/point.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/rtree.h>

#include "../tests.h"

namespace bgi = boost::geometry::index;

int
main()
{
  initlog();
  std::vector<Point<2>> pts(8);
  for (auto &p : pts)
    p = random_point<2>();

  auto tree = pack_rtree_of_indices<bgi::linear<2>>(pts);

  // internal node visitor
  auto internal = [](std::size_t           level,
                     auto                 &parent,
                     const BoundingBox<2> &box,
                     auto                 &child) {
    deallog << "Internal node (level " << level << ") box "
            << Patterns::Tools::to_string(box.get_boundary_points())
            << std::endl;
  };

  // leaf visitor
  auto leaf = [](std::size_t level, auto &, const BoundingBox<2> &box, auto &) {
    deallog << "    Leaf      (level " << level << ") box "
            << Patterns::Tools::to_string(box.get_boundary_points())
            << std::endl;
  };

  // element visitor (receives the stored value)
  auto element = [](auto &, const auto &value) {
    deallog << "        element " << value;
  };

  // indexable visitor (receives the indexable geometry)
  auto indexable = [](auto &, const auto &indexable) {
    deallog << ",  indexable " << indexable << std::endl;
  };

  deallog << "Visiting R-tree with " << pts.size() << " points, and "
          << n_levels(tree) << " levels." << std::endl;

  visit_rtree(tree, internal, leaf, element, indexable);
}

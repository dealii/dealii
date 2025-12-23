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

// Visit an rtree, without any action associated.

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

  auto tree = pack_rtree<bgi::linear<2>>(pts);

  // Do nothing.
  visit_rtree(tree);

  deallog << "OK" << std::endl;
}

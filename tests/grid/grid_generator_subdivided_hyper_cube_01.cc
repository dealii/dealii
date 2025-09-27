// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test for GridGenerator::subdivided_hyper_cube() with colorized boundary
// conditions

#include <deal.II/base/tensor.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
test()
{
  deallog << "dim = " << dim << std::endl;
  Triangulation<dim> tr;
  GridGenerator::subdivided_hyper_cube(tr, 1, -1, 1, true);
  for (const auto &cell : tr.active_cell_iterators())
    {
      deallog << "cell:" << std::endl;

      for (const auto &face : cell->face_iterators())
        {
          if (face->at_boundary())
            deallog << "boundary id = " << face->boundary_id()
                    << " center = " << face->center()
                    << " faceidx = " << face->index() << std::endl;
        }
    }
}

int
main()
{
  initlog();
  test<2>();
  test<3>();
}

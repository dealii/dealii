// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test grid generation for colorized cylinder_shell

#include <deal.II/base/tensor.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



template <int dim>
void
test(std::ostream &out)
{
  deallog << "cylinder_shell colorized with default aspect ratio" << std::endl;
  {
    Triangulation<dim> tr;
    GridGenerator::cylinder_shell(tr, 2., 5., 6., 0, 0, true);

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

  deallog << "cylinder_shell colorized with 3:2 ratio " << std::endl;
  {
    Triangulation<dim> tr;
    GridGenerator::cylinder_shell(tr, 2., 5., 6., 3, 2, true);

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
}


int
main()
{
  initlog();
  test<3>(deallog.get_file_stream());
}

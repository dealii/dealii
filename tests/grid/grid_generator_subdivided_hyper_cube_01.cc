// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

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

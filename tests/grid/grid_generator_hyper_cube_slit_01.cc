// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test output for GridGenerator::hyper_cube_slit()

#include <deal.II/base/tensor.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



template <int dim>
void
test()
{
  deallog << "dim = " << dim << std::endl;
  Triangulation<dim> tr;
  GridGenerator::hyper_cube_slit(tr, -1, 1, true);
  typename Triangulation<dim>::active_cell_iterator cell = tr.begin_active(),
                                                    endc = tr.end();
  for (; cell != endc; ++cell)
    {
      deallog << "cell:" << std::endl;

      for (const unsigned int face : GeometryInfo<dim>::face_indices())
        {
          if (cell->face(face)->at_boundary() &&
              cell->face(face)->boundary_id() != 0)
            deallog << "boundary id = " << (int)cell->face(face)->boundary_id()
                    << " center = " << cell->face(face)->center()
                    << " faceidx = " << face << std::endl;
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

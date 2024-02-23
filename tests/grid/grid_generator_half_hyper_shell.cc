// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test boundary ids for GridGenerator::half_hyper_shell()

#include <deal.II/base/tensor.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
dim_test()
{
  Triangulation<dim> tr;
  double             r_inner = 1.;
  double             r_outer = 2.;

  GridGenerator::half_hyper_shell(tr, Point<dim>(), r_inner, r_outer, 0, true);

  typename Triangulation<dim>::cell_iterator cell = tr.begin();

  for (; cell != tr.end(); ++cell)
    {
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        {
          if (cell->face(f)->at_boundary())
            deallog << "half_hyper_shell<" << dim << ">:: cell "
                    << cell->index() << ", face " << f << ", id "
                    << cell->face(f)->boundary_id() << std::endl;
        }
    }
}


int
main()
{
  initlog();

  dim_test<2>();
  dim_test<3>();
}

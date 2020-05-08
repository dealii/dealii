// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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

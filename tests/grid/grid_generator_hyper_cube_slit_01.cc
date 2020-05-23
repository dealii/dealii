// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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

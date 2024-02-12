// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2019 by the deal.II authors
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
  deallog << "cylinder_shell colorized" << std::endl;
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


int
main()
{
  initlog();
  test<3>(deallog.get_file_stream());
}

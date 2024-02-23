// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// in 1d, we have to read vertex information to set boundary
// indicators
//
// test case by Jan Strebel
//
// this is a variation of the grid_in_msh_02 testcase, but as in Jan's
// original code snippet, it uses GridIn<1,3> (which wasn't
// instantiated at the time)
//
//  This is the same as grid_in_msh_02_13 but uses an input file in GMSH-4
//  format


#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <string>

#include "../tests.h"

void
check_file()
{
  Triangulation<1, 3> tria;
  GridIn<1, 3>        gi;
  gi.attach_triangulation(tria);
  std::ifstream in(SOURCE_DIR "/../grid/grids/grid_in_msh_02_v4.msh");
  gi.read_msh(in);

  for (Triangulation<1, 3>::active_cell_iterator cell = tria.begin_active();
       cell != tria.end();
       ++cell)
    {
      for (unsigned int face = 0; face < 2; ++face)
        {
          if (cell->at_boundary(face))
            deallog << "vertex " << cell->face_index(face)
                    << " has boundary indicator "
                    << (int)cell->face(face)->boundary_id() << std::endl;
        }
    }
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  check_file();
}

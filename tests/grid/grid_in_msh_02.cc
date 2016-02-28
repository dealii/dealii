// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// in 1d, we have to read vertex information to set boundary
// indicators
//
// test case by Jan Strebel


#include "../tests.h"
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/logstream.h>

#include <fstream>
#include <iomanip>
#include <string>

std::ofstream logfile("output");


void check_file ()
{
  Triangulation<1> tria;
  GridIn<1> gi;
  gi.attach_triangulation (tria);
  std::ifstream in (SOURCE_DIR "/../grid/grids/grid_in_msh_02.msh");
  gi.read_msh(in);

  for (Triangulation<1>::active_cell_iterator cell = tria.begin_active(); cell != tria.end(); ++cell)
    {
      for (unsigned int face = 0; face < 2; ++face)
        {
          if (cell->at_boundary(face))
            deallog << "vertex " << cell->face_index(face)
                    << " has boundary indicator " << (int)cell->face(face)->boundary_id()
                    << std::endl;
        }
    }
}


int main ()
{
  deallog << std::setprecision (2);
  logfile << std::setprecision (2);
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  check_file ();
}


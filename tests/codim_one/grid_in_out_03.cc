// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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



// like grid_in_out but test persistent triangulation

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/persistent_tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

#include <string>

std::ofstream logfile("output");

template <int dim, int spacedim>
void test(std::string filename)
{
  logfile << " Tria<" << dim << "," << spacedim << ">: "
          << std::endl;
  Triangulation<dim, spacedim> tria;
  GridIn<dim, spacedim> gi;
  gi.attach_triangulation (tria);
  std::ifstream in (filename.c_str());
  gi.read_ucd (in);

  PersistentTriangulation<dim, spacedim> ptria(tria);
  typename Triangulation<dim,spacedim>::active_cell_iterator cell;

  ptria.restore();
  for (unsigned int i=0; i<2; ++i)
    {
      for (cell=ptria.begin_active(); cell != ptria.end(); ++cell)
        if (cell->center()[0]<.5)
          cell->set_refine_flag();
      ptria.execute_coarsening_and_refinement();
    }

  ptria.write_flags(logfile);
}

int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<2,3>(SOURCE_DIR "/grids/square.inp");
  test<2,3>(SOURCE_DIR "/grids/sphere_1.inp");

  return 0;
}


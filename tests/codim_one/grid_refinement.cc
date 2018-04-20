// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2017 by the deal.II authors
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



// see what happens when creating a surface mesh and then refining it

#include "../tests.h"

// all include files you need here

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <string>

std::ofstream logfile("output");

template <int dim, int spacedim>
void test(std::string filename)
{
  SphericalManifold<dim, spacedim> boundary;
  Triangulation<dim, spacedim> tria;
  GridIn<dim, spacedim> gi;
  gi.attach_triangulation (tria);
  std::ifstream in (filename.c_str());
  gi.read_ucd (in);
  tria.set_all_manifold_ids(1);
  tria.set_manifold(1, boundary);

  GridOut grid_out;
  grid_out.set_flags (GridOutFlags::Ucd(true));
  for (unsigned int cycle=0; cycle<3; ++cycle)
    {
      tria.refine_global(1);
      grid_out.write_msh (tria, logfile);
    }
}

int main ()
{
  deallog.attach(logfile);

  deallog << "Test<1,2>" << std::endl;
  test<1,2>(SOURCE_DIR "/grids/circle_1.inp");

  deallog << std::endl << "Test<1,2>" << std::endl;
  test<2,3>(SOURCE_DIR "/grids/sphere_1.inp");

  return 0;
}


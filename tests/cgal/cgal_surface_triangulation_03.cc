// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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

// Convert a closed CGAL::Surface_mesh or a closed CGAL::Polyhedron_3 to a
// deal.II triangulation. The input mesh is a tripod.

#include <deal.II/base/config.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <CGAL/IO/io.h>
#include <deal.II/cgal/triangulation.h>

#include <fstream>

#include "../tests.h"

using namespace CGALWrappers;
using Kernel     = CGAL::Simple_cartesian<double>;
using CGALPoint  = CGAL::Point_3<Kernel>;
using Mesh       = CGAL::Surface_mesh<CGALPoint>;
using Polyhedron = CGAL::Polyhedron_3<Kernel>;

int
main()
{
  initlog();
  Triangulation<2, 3> tria;
  GridOut             go;
  {
    deallog << "Surface mesh" << std::endl;
    Mesh          sm;
    std::ifstream input_sm(SOURCE_DIR "/input_grids/tripod.off");
    input_sm >> sm;
    cgal_surface_mesh_to_dealii_triangulation(sm, tria);
    {
      std::ofstream out("temp.msh");
      go.write_msh(tria, out);
    }
    remove("temp.msh");
    // If we got here, everything is fine.
    deallog << "OK" << std::endl;
  }
  tria.clear();

  // Now do the same with a Polyhedron
  deallog << "Polyhedron test" << std::endl;
  {
    Polyhedron    P;
    std::ifstream input_p(SOURCE_DIR "/input_grids/tripod.off");
    input_p >> P;
    cgal_surface_mesh_to_dealii_triangulation(P, tria);
    {
      std::ofstream out("temp.msh");
      go.write_msh(tria, out);
    }
    remove("temp.msh");
    // If we got here, everything is fine.
    deallog << "OK" << std::endl;
  }
}

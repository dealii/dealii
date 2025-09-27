// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Convert a closed CGAL::Surface_mesh or a closed CGAL::Polyhedron_3 to a
// deal.II triangulation. The input mesh is a tripod.

#include <deal.II/base/config.h>

#include <deal.II/cgal/triangulation.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <CGAL/IO/io.h>

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

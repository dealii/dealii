// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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



// like mesh_3d_22, but even further reduced: take a few points on a
// line in reference space and see where they are transformed to real
// space using a MappingQ(3) on a cell with a flipped face. the cell
// is a cube, so the points should also be mapped to a straight line,
// but they are not, unfortunately, when this test was written.

#include "../tests.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_q.h>
#include <fstream>
#include <cmath>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <vector>


int main ()
{
  std::ofstream logfile ("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.0e-10);

  Triangulation<3> triangulation;
  GridIn<3> grid_in;
  grid_in.attach_triangulation(triangulation);
  std::ifstream inputStream(SOURCE_DIR "/grids/mesh.msh");
  grid_in.read_msh (inputStream);

  MappingQ<3> mapping(3);

  Triangulation<3>::active_cell_iterator cell;
  cell = triangulation.begin_active();
  ++cell;

  for (unsigned int f=0; f<6; ++f)
    deallog << "face: " << f
            << std::endl
            << "  Face orientation status: "
            << cell->face_orientation (f)
            << std::endl
            << "  On boundary: " << cell->at_boundary(f)
            << std::endl
            << "  Neighbor: "
            << (cell->at_boundary(f) ? triangulation.end() : cell->neighbor(f))
            << std::endl;

  for (double x=0; x<=1; x+=1.0/3.0)
    {
      const Point<3> p(x, 1./3., 0);
      deallog << p
              << " in unit coordinates maps to real coordinates "
              << mapping.transform_unit_to_real_cell(cell,p)
              << std::endl;
    }
}


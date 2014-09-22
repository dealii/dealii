// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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



// transform_real_to_unit_cell works for MappingQ(1) but fails for MappingQ1.
//
// this is a case where a cell is extremely elongated, and just taking the l2
// norm to measure progress in a Newton search is inadequate. we need to take a
// variable norm that keeps the elongation of the cell in mind.

#include "../tests.h"

#include <stdio.h>
#include <cstdlib>

#include <base/quadrature_lib.h>
#include <fe/mapping_q.h>
#include <base/function.h>
#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_tools.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <grid/grid_in.h>


#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <sstream>
#include <time.h>

using namespace dealii;


void create_coarse_grid(Triangulation<2> &coarse_grid)
{

  static const Point<2> vertices_1[]
    = {  Point<2> (9.6982181981258408e-02, 1.1255621492491609e+03),//0
         Point<2> (6.1219285295807092e-02, 1.1256062663438720e+03),//1
         Point<2> (0.00000, 1.1255179557007179e+03),//2
         Point<2> (0.00000, 1.1255994426210491e+03),//3
      };
  const unsigned int
  n_vertices = sizeof(vertices_1) / sizeof(vertices_1[0]);

  const std::vector<Point<2> > vertices (&vertices_1[0],
                                         &vertices_1[n_vertices]);

  static const int cell_vertices[][GeometryInfo<2>::vertices_per_cell]
  = {{0, 1, 2, 3}
  };
  const unsigned int
  n_cells = sizeof(cell_vertices) / sizeof(cell_vertices[0]);

  std::vector<CellData<2> > cells (n_cells, CellData<2>());
  for (unsigned int i=0; i<n_cells; ++i)
    {
      for (unsigned int j=0;
           j<GeometryInfo<2>::vertices_per_cell;
           ++j)
        cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    }

  coarse_grid.create_triangulation (vertices,
                                    cells,
                                    SubCellData());
}

void test()
{
  deallog << std::scientific;
  deallog.precision(16);
  
  Point< 2 > ePos;
  ePos(0) = 0.0653630060373507487669897386695;
  ePos(1) = 1125.59175030825804242340382189;

  MappingQ<2> mapping(1);
  MappingQ1<2> &mapping2 = StaticMappingQ1< 2 >::mapping;

  Triangulation<2>     triangulation;
  create_coarse_grid(triangulation); // first Tria with just one cell

  Triangulation<2>::active_cell_iterator it = triangulation.begin();
				   
  Point<2> p;
  p = mapping.transform_real_to_unit_cell(it, ePos);
  deallog << "A: " << p << std::endl;

				   // throws:
  p = mapping2.transform_real_to_unit_cell(it, ePos);
  deallog << "B: " << p << std::endl;

  deallog << "done" << std::endl;
}

int main (int argc, char **argv)
{
  initlog();

  test();

  return 0;
}

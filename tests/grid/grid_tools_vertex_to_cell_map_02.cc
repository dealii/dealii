
//
// Copyright (C) 2015, 2016 by the deal.II authors
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


// Test vertex_to_cell_map for 3D problem with hanging nodes.

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/cell_id.h>

#include <vector>

std::ofstream logfile("output");

void test()
{
  Triangulation<3> tria;
  Point<3> a(0.,0.,0.);
  Point<3> b(1.,1.,1.);
  std::vector<unsigned int> repetitions(3);
  repetitions[0] = 2;
  repetitions[1] = 2;
  repetitions[2] = 1;
  GridGenerator::subdivided_hyper_rectangle(tria,repetitions,a,b);
  Triangulation<3>::active_cell_iterator cell = tria.begin_active();
  cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  std::vector<std::set<Triangulation<3>::active_cell_iterator> >
  vertex_to_cell = GridTools::vertex_to_cell_map(tria);

  AssertThrow(tria.n_vertices()==vertex_to_cell.size(),
              ExcMessage("Wrong number of vertices"));

  std::vector<unsigned int> n_cells;
  for (unsigned int i=0; i<vertex_to_cell.size(); ++i)
    n_cells.push_back(vertex_to_cell[i].size());

  std::vector<unsigned int> histogram(9,0);
  for (unsigned int i=0; i<n_cells.size(); ++i)
    histogram[n_cells[i]] += 1;

  AssertThrow(histogram[0]==0, ExcMessage("Wrong cell distribution"));
  AssertThrow(histogram[1]==8, ExcMessage("Wrong cell distribution"));
  AssertThrow(histogram[2]==13, ExcMessage("Wrong cell distribution"));
  AssertThrow(histogram[3]==6, ExcMessage("Wrong cell distribution"));
  AssertThrow(histogram[4]==6, ExcMessage("Wrong cell distribution"));
  AssertThrow(histogram[5]==3, ExcMessage("Wrong cell distribution"));
  AssertThrow(histogram[6]==0, ExcMessage("Wrong cell distribution"));
  AssertThrow(histogram[7]==0, ExcMessage("Wrong cell distribution"));
  AssertThrow(histogram[8]==1, ExcMessage("Wrong cell distribution"));
}

int main (int argc, char *argv[])
{
  deallog << std::setprecision(4);
  logfile << std::setprecision(4);
  deallog.attach(logfile);

  test();

  deallog << "OK" << std::endl;

  return 0;
}

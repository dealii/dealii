// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// read and write a 2d file in the VTU format

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <boost/archive/text_oarchive.hpp>

#include <string>

#include "../tests.h"

int
main()
{
  initlog();

  Triangulation<2> tria;
  Triangulation<2> tria2;
  GridGenerator::hyper_ball(tria);

  for (const auto &cell : tria.active_cell_iterators())
    {
      if (cell->center()[0] < 0)
        cell->set_refine_flag();
    }
  tria.execute_coarsening_and_refinement();

  for (const auto &cell : tria.active_cell_iterators())
    {
      if (cell->center()[1] < 0)
        cell->set_refine_flag();
    }
  tria.execute_coarsening_and_refinement();

  {
    std::ofstream     out("grid.vtu");
    GridOut           go;
    GridOutFlags::Vtu vtu_flags(true);
    go.set_flags(vtu_flags);
    go.write_vtu(tria, out);
  }
  {
    std::ifstream in("grid.vtu");
    GridIn<2>     gi;
    gi.attach_triangulation(tria2);
    gi.read_vtu(in);
  }
  std::stringstream stream1;
  std::stringstream stream2;
  {
    boost::archive::text_oarchive oa(stream1);
    tria.save(oa, 0);
  }
  {
    boost::archive::text_oarchive oa(stream2);
    tria2.save(oa, 0);
  }
  if (stream1.str() == stream2.str())
    deallog << "OK" << std::endl;
  else
    deallog << "Not OK" << std::endl;
  std::remove("grid.vtu");
}

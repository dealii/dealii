// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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



// Test refinement of simplex mesh.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

#include "./simplex_grids.h"

void
test(const unsigned int v)
{
  const unsigned int dim = 2;

  Triangulation<dim> tria;

  if (v < 4)
    GridGenerator::subdivided_hyper_cube_with_simplices(tria, 1);
  else
    GridGenerator::subdivided_hyper_cube_with_simplices_mix(tria, 2);


  if (v == 1 || v == 4)
    {
      tria.begin()->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }
  else if (v == 2 || v == 5)
    {
      tria.refine_global(1);
    }
  else if (v == 3)
    {
      tria.begin()->set_refine_flag();
      tria.execute_coarsening_and_refinement();
      tria.begin_active(1)->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  GridOut grid_out;
#if false
  std::ofstream out("mesh." + std::to_string(v) + ".vtk");
  grid_out.write_vtk(tria, out);
#else
  grid_out.write_vtk(tria, deallog.get_file_stream());
#endif
}

int
main()
{
  initlog();
  test(0); // no refinement
  test(1); // refinement of a single cell
  test(2); // global refinement
  test(3); // refine left bottom corner twice
  test(4); // mixed mesh: refinement of a single cell
  test(5); // mixed mesh: global refinement
}

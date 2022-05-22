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



// Test refinement of 3D simplex mesh.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
print(Triangulation<dim> &tria, const std::string &label)
{
  GridOut grid_out;

  unsigned int counter = 0;

  for (const auto &cell : tria.active_cell_iterators())
    cell->set_material_id(counter++);

#if false
  std::ofstream out(label);
  grid_out.write_vtk(tria, out);
#else
  (void)label;
  grid_out.write_vtk(tria, deallog.get_file_stream());
#endif
}

void
test(const unsigned int version)
{
  Triangulation<3> tria;

  if (version == 0 || version == 1)
    GridGenerator::reference_cell(tria, ReferenceCells::Tetrahedron);
  else if (version == 2)
    GridGenerator::subdivided_hyper_cube_with_simplices(tria, 1);

  print(tria, "tria.0.vtk");
  if (version == 0 || version == 2)
    tria.refine_global(1);
  else if (version == 1)
    tria.refine_global(2);

  print(tria, "tria.1.vtk");
}

int
main()
{
  initlog();
  test(0); // coarse grid: 1 tet, n_refinements: 1
  test(1); // coarse grid: 1 tet, n_refinements: 2
  test(2); // coarse grid: 5 tets (cube), n_refinements: 1
}

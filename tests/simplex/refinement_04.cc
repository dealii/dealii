// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



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
  else if (version == 2 || version == 3)
    GridGenerator::subdivided_hyper_cube_with_simplices(tria, 1);

  print(tria, "tria.0.vtk");
  if (version == 0 || version == 2)
    tria.refine_global(1);
  else if (version == 1 || version == 3)
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
  test(3); // coarse grid: 5 tets (cube), n_refinements: 2
}

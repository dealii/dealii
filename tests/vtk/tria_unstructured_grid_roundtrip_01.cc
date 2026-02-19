// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Test deal.II triangulation <-> vtkUnstructuredGrid round-trip conversion.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools_topology.h>
#include <deal.II/grid/tria.h>

#include <deal.II/vtk/utilities.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> tria_in;
  GridGenerator::subdivided_hyper_cube(tria_in, 2);

  const auto grid =
    VTKWrappers::dealii_triangulation_to_unstructured_grid(tria_in);

  Triangulation<dim, spacedim> tria_out;
  VTKWrappers::unstructured_grid_to_dealii_triangulation(*grid, tria_out);

  AssertThrow(GridTools::have_same_coarse_mesh(tria_in, tria_out),
              ExcInternalError());
  AssertDimension(tria_in.n_active_cells(), tria_out.n_active_cells());
  AssertDimension(tria_in.n_vertices(), tria_out.n_vertices());

  deallog << "Round-trip OK for Triangulation<" << dim << ", " << spacedim
          << ">: " << tria_out.n_active_cells() << " cells, "
          << tria_out.n_vertices() << " vertices" << std::endl;
}

int
main()
{
  initlog();

  test<1, 1>();
  test<2, 2>();
  test<3, 3>();

  return 0;
}

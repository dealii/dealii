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

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools_topology.h>
#include <deal.II/grid/tria.h>

#include <deal.II/vtk/utilities.h>

#include "../tests.h"

using namespace dealii;

// Test: round-trip conversion Triangulation -> VTK UnstructuredGrid ->
// Triangulation, ensuring material, boundary and manifold ids are preserved.

template <int dim, int spacedim>
void
test_roundtrip_ids()
{
  Triangulation<dim, spacedim> tria_in;
  GridGenerator::subdivided_hyper_cube(tria_in, 1);

  // assign material ids and set boundary/manifold ids on boundary faces
  unsigned int cell_index = 0;
  for (auto cell = tria_in.begin_active(); cell != tria_in.end();
       ++cell, ++cell_index)
    {
      cell->set_material_id(100 + static_cast<types::material_id>(cell_index));
      for (unsigned int f = 0; f < cell->n_faces(); ++f)
        if (cell->face(f)->at_boundary())
          {
            if constexpr (dim == 2)
              cell->face(f)->set_boundary_id(13);
            else if constexpr (dim == 3)
              cell->face(f)->set_manifold_id(5);
          }
    }

  auto grid = VTKWrappers::dealii_triangulation_to_unstructured_grid(
    tria_in, "material", "boundary", "manifold");

  Triangulation<dim, spacedim> tria_out;
  VTKWrappers::unstructured_grid_to_dealii_triangulation(
    *grid, tria_out, "material", "boundary", "manifold");

  AssertThrow(GridTools::have_same_coarse_mesh(tria_in, tria_out),
              ExcInternalError());
  AssertDimension(tria_in.n_active_cells(), tria_out.n_active_cells());
  AssertDimension(tria_in.n_vertices(), tria_out.n_vertices());

  // verify ids
  auto it_in  = tria_in.begin_active();
  auto it_out = tria_out.begin_active();
  for (; it_in != tria_in.end(); ++it_in, ++it_out)
    {
      deallog << "cell material in/out: "
              << static_cast<int>(it_in->material_id()) << " "
              << static_cast<int>(it_out->material_id()) << std::endl;
      AssertDimension(it_in->material_id(), it_out->material_id());
      AssertDimension(it_in->manifold_id(), it_out->manifold_id());
      for (unsigned int f = 0; f < it_in->n_faces(); ++f)
        {
          deallog << "face " << f << " boundary in/out: "
                  << static_cast<int>(it_in->face(f)->boundary_id()) << " "
                  << static_cast<int>(it_out->face(f)->boundary_id())
                  << std::endl;
          AssertDimension(it_in->face(f)->boundary_id(),
                          it_out->face(f)->boundary_id());
          deallog << "face " << f << " manifold in/out: "
                  << static_cast<int>(it_in->face(f)->manifold_id()) << " "
                  << static_cast<int>(it_out->face(f)->manifold_id())
                  << std::endl;
          AssertDimension(it_in->face(f)->manifold_id(),
                          it_out->face(f)->manifold_id());
        }
    }

  deallog << "Round-trip-ID OK for Triangulation<" << dim << ", " << spacedim
          << ">: " << tria_out.n_active_cells() << " cells, "
          << tria_out.n_vertices() << " vertices" << std::endl;
}

int
main()
{
  initlog();

  test_roundtrip_ids<2, 2>();
  test_roundtrip_ids<3, 3>();

  return 0;
}

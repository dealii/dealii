/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2020 - 2025 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

// Test function GridGenerator::convert_simplex_to_hypercube_mesh() in
// 2D and 3D for a mesh that consists of only a single simplex that is then
// once refined. The function can deal with refined meshes, just not
// adaptively refined ones, which it simply flattens into a single
// coarse mesh before converting into a quad/hex mesh.


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <fstream>
#include <iostream>

#include "../tests.h"



template <int dim, int spacedim>
void
check()
{
  Triangulation<dim, spacedim> in_tria, out_tria;
  GridGenerator::reference_cell(in_tria, ReferenceCells::get_simplex<dim>());
  in_tria.refine_global();

  // make each cell a different material id
  unsigned int m_id = 1;
  for (const auto &cell : in_tria)
    {
      cell.set_material_id(m_id++);
    }

  // set different boundary ids and output
  unsigned int b_id = 42;
  for (const auto &cell : in_tria)
    {
      for (const auto f : cell.face_indices())
        {
          if (cell.face(f)->at_boundary())
            {
              cell.face(f)->set_boundary_id(b_id);
              b_id++;
            }
        }
    }

  // Now convert
  GridGenerator::convert_simplex_to_hypercube_mesh(in_tria, out_tria);

  // And write the mesh:
  GridOut grid_out;
  grid_out.write_vtk(out_tria, deallog.get_file_stream());

  // Then also make sure we get the right material and boundary ids
  for (const auto &cell : out_tria.active_cell_iterators())
    deallog << cell << ", material id=" << cell->material_id() << std::endl;
  for (const auto &cell : out_tria.active_cell_iterators())
    for (const auto &face : cell->face_iterators())
      if (face->at_boundary())
        deallog << face << ", boundary id=" << face->boundary_id() << std::endl;
}


int
main()
{
  initlog();
  check<2, 2>();
  check<3, 3>();
}

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test that cell manifold ids are also set for the
// "apply_all_indicators_to_manifolds" option in GridIn::read_ucd

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <sstream>
#include <string>

#include "../tests.h"

template <int dim, int spacedim>
void
test(const std::string &filename)
{
  std::ifstream                tmp_in(filename);
  Triangulation<dim, spacedim> tria;
  GridIn<dim, spacedim>        gi;
  gi.attach_triangulation(tria);
  gi.read_ucd(tmp_in, true);
  deallog << "Testing dim=" << dim << " spacedim=" << spacedim << std::endl;
  for (const auto &cell : tria.active_cell_iterators())
    {
      deallog << "cell->manifold_id: " << cell->manifold_id() << std::endl;
      if (dim > 1)
        for (const unsigned int face_no : GeometryInfo<dim>::face_indices())
          deallog << "cell->face(" << face_no
                  << ")->manifold_id: " << cell->face(face_no)->manifold_id()
                  << std::endl;
      if (dim > 1)
        for (unsigned int line_no = 0;
             line_no < GeometryInfo<dim>::lines_per_cell;
             ++line_no)
          deallog << "cell->line(" << line_no
                  << ")->manifold_id: " << cell->line(line_no)->manifold_id()
                  << std::endl;
    }
  deallog << std::endl;
}

int
main()
{
  initlog();
  deallog.get_file_stream() << std::setprecision(2);

  test<1, 1>(SOURCE_DIR "/grid_in_ucd_02_grids/grid_1.inp");
  test<1, 2>(SOURCE_DIR "/grid_in_ucd_02_grids/grid_1.inp");
  test<1, 3>(SOURCE_DIR "/grid_in_ucd_02_grids/grid_1.inp");
  test<2, 2>(SOURCE_DIR "/grid_in_ucd_02_grids/grid_2.inp");
  test<2, 3>(SOURCE_DIR "/grid_in_ucd_02_grids/grid_2.inp");
  test<3, 3>(SOURCE_DIR "/grid_in_ucd_02_grids/grid_3.inp");
}

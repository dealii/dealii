// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// verify the distortion in cells of a hyper shell with 6 cells upon
// refinement

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <iostream>

#include "../tests.h"

template <int dim>
void
check(double r1, double r2, unsigned int n)
{
  deallog << "dim=" << dim << std::endl;

  Point<dim>         center;
  Triangulation<dim> tria(Triangulation<dim>::none, true);

  // create a hyper shell. before a bug fix in early 2015, the coloring would
  // only set the boundary indicators of the faces, not of the edges in
  // 3d. the test was written in a way that it tests with this bug present, so
  // undo the big fix here
  GridGenerator::hyper_shell(tria, center, r1, r2, n, true);
  tria.reset_manifold(0);
  tria.set_all_manifold_ids(numbers::flat_manifold_id);
  GridTools::copy_boundary_to_manifold_id(tria);

  for (const auto bid : tria.get_boundary_ids())
    tria.set_manifold(bid, FlatManifold<dim>());

  if (dim == 3)
    for (typename Triangulation<dim>::active_cell_iterator c =
           tria.begin_active();
         c != tria.end();
         ++c)
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        if (c->face(f)->at_boundary())
          for (unsigned int e = 0; e < GeometryInfo<dim>::lines_per_face; ++e)
            c->face(f)->line(e)->set_manifold_id(0);

  static const SphericalManifold<dim> boundary(center);
  tria.set_manifold(0, boundary);

  for (unsigned int i = 0; i < 2; ++i)
    {
      try
        {
          tria.refine_global(1);
        }
      catch (typename Triangulation<dim>::DistortedCellList &dcv)
        {
          deallog << "Found " << dcv.distorted_cells.size()
                  << " distorted cells" << std::endl;

          typename Triangulation<dim>::DistortedCellList subset =
            GridTools::fix_up_distorted_child_cells(dcv, tria);
          deallog << subset.distorted_cells.size()
                  << " distorted cells remaining" << std::endl;
        }
    }

  GridOut          grid_out;
  GridOutFlags::DX flags;
  flags.write_faces = true;
  grid_out.set_flags(flags);
  grid_out.write_dx(tria, deallog.get_file_stream());
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  check<2>(4., 5., 10);
  check<3>(3., 5., 6);
}

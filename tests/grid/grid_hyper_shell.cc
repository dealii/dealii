// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2018 by the deal.II authors
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



// verify the distortion in cells of a hyper shell with 6 cells upon
// refinement

#include "../tests.h"
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <iostream>

std::ofstream logfile("output");


template <int dim>
void check (double r1, double r2, unsigned int n)
{
  deallog << "dim=" << dim << std::endl;

  Point<dim> center;
  Triangulation<dim> tria (Triangulation<dim>::none, true);

  // create a hyper shell. before a bug fix in early 2015, the coloring would
  // only set the boundary indicators of the faces, not of the edges in
  // 3d. the test was written in a way that it tests with this bug present, so
  // undo the big fix here
  GridGenerator::hyper_shell (tria, center, r1, r2, n, true);
  tria.reset_manifold(0);
  tria.set_all_manifold_ids(numbers::flat_manifold_id);
  GridTools::copy_boundary_to_manifold_id(tria);
  if (dim==3)
    for (typename Triangulation<dim>::active_cell_iterator c=tria.begin_active(); c!=tria.end(); ++c)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (c->face(f)->at_boundary())
          for (unsigned int e=0; e<GeometryInfo<dim>::lines_per_face; ++e)
            c->face(f)->line(e)->set_manifold_id(0);

  static const SphericalManifold<dim> boundary(center);
  tria.set_manifold(0, boundary);

  for (unsigned int i=0; i<2; ++i)
    {
      try
        {
          tria.refine_global(1);
        }
      catch (typename Triangulation<dim>::DistortedCellList &dcv)
        {
          deallog << "Found " << dcv.distorted_cells.size()
                  << " distorted cells" << std::endl;

          typename Triangulation<dim>::DistortedCellList
          subset = GridTools::fix_up_distorted_child_cells (dcv,
                                                            tria);
          deallog << subset.distorted_cells.size()
                  << " distorted cells remaining" << std::endl;
        }
    }

  GridOut grid_out;
  GridOutFlags::DX flags;
  flags.write_faces = true;
  grid_out.set_flags(flags);
  grid_out.write_dx (tria, logfile);
}


int main()
{
  deallog << std::setprecision(3);
  deallog.attach(logfile);

  check<2> (4., 5., 10);
  check<3> (3., 5., 6);
}

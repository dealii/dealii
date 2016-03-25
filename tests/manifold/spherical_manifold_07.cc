// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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

// test get_normals_at_vertices for a SphericalManifold.

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/manifold_lib.h>

#include <fstream>
#include <iomanip>


template <int dim, int spacedim>
void test ()
{
  Triangulation<dim, spacedim> triangulation;

  GridGenerator::hyper_ball(triangulation);

  static const SphericalManifold<dim, spacedim> manifold;
  triangulation.set_all_manifold_ids_on_boundary(0);
  triangulation.set_manifold (0, manifold);

  triangulation.refine_global(1);

  unsigned int c = 0;
  deallog.get_file_stream()
      << "set view equal xyz" << std::endl
      << "set size ratio -1" << std::endl
      << "set multiplot" << std::endl
      << (dim == 3 ? "s" : "")
      << "plot '-' with vectors " << std::endl;

  for (typename Triangulation<dim,spacedim>::active_cell_iterator
       cell = triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    {
      typename Manifold<dim,spacedim>::FaceVertexNormals n;

      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->face(f)->at_boundary())
          {
            manifold.get_normals_at_vertices(cell->face(f), n);

            for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_face; ++v)
              {
                Tensor<1,spacedim> dn = n[v]-manifold.normal_vector
                                        (cell->face(f), cell->face(f)->vertex(v));

                if ( dn.norm() > 1e-10)
                  deallog << "Error on vertex " << f << ", face "
                          << cell->face(f) << ": " << dn << std::endl;

                deallog.get_file_stream() << cell->face(f)->vertex(v)
                                          << " " << n[v]*.1 << std::endl;
              }
          }
    }
  deallog.get_file_stream() << "e" << std::endl;
}


int main ()
{
  initlog ();

  test<2,2> ();
  test<3,3> ();

  return 0;
}

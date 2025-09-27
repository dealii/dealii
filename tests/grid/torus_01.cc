// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test GridTools::torus() the manifolds attached to the torus (TorusManifold,
// ToroidalManifold, TransfiniteInterpolationManifold), output visually
// checked

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



template <int dim, int spacedim>
void
test();

template <>
void
test<3, 3>()
{
  const int                    dim      = 3;
  const int                    spacedim = 3;
  Triangulation<dim, spacedim> triangulation;

  GridGenerator::torus(triangulation, 1.0, 0.4);

  triangulation.begin_active()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();

  const TorusManifold<3> desc_torus(1.0, 0.4);
  unsigned int           c = 0;
  for (Triangulation<dim, spacedim>::active_vertex_iterator v =
         triangulation.begin_active_vertex();
       v != triangulation.end_vertex();
       ++v, ++c)
    {
      if (c % 3 != 0)
        continue;
      Point<3> p = v->vertex(0);
      Point<3> x(numbers::PI / 2.5, numbers::PI / 3.5, 1.0);
      x              = desc_torus.push_forward(x);
      Tensor<1, 3> t = desc_torus.get_tangent_vector(p, x);

      deallog.get_file_stream()
        << "set arrow from " << p[0] << ", " << p[1] << ", " << p[2] << " to "
        << t[0] << ", " << t[1] << ", " << t[2] << std::endl;
    }
  deallog.get_file_stream() << "set view equal xyz" << std::endl
                            << "splot '-' w l" << std::endl;
  GridOut().write_gnuplot(triangulation, deallog.get_file_stream());
  deallog.get_file_stream() << 'e' << std::endl << "pause -1" << std::endl;
}


int
main()
{
  initlog();

  test<3, 3>();

  return 0;
}

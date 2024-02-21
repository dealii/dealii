// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test mapping surfaces in higher dimensions. when we use the
// MappingQ(1) class, each 1d cell in 2d space is mapped to a straight
// line and so all cell normals should be parallel. likewise, if the four
// vertices of a 2d cell in 3d space are in a plane, then the cell normal
// vectors at all quadrature points of the same cell should be parallel, even
// though the true surface is curved

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test()
{
  deallog << "Testing hyper_ball in dim: " << dim << "..." << std::endl;

  Triangulation<dim> volume_mesh;
  GridGenerator::hyper_ball(volume_mesh);

  const SphericalManifold<dim - 1, dim> surface_description;
  Triangulation<dim - 1, dim>           boundary_mesh;
  boundary_mesh.set_manifold(0, surface_description);

  GridGenerator::extract_boundary_mesh(volume_mesh, boundary_mesh);

  QGauss<dim - 1>        quadrature(2);
  MappingQ<dim - 1, dim> mapping(1);
  FE_Q<dim - 1, dim>     fe(1);

  FEValues<dim - 1, dim> fe_values(mapping,
                                   fe,
                                   quadrature,
                                   update_normal_vectors);

  for (typename Triangulation<dim - 1, dim>::active_cell_iterator cell =
         boundary_mesh.begin_active();
       cell != boundary_mesh.end();
       ++cell)
    {
      deallog << "Cell = " << cell << ", with center at " << cell->center()
              << std::endl;
      fe_values.reinit(cell);

      for (unsigned int q = 0; q < quadrature.size(); ++q)
        deallog << "  cell_normal[" << q << "] = " << fe_values.normal_vector(q)
                << std::endl;
    }
}



int
main()
{
  initlog();

  test<2>();
  test<3>();

  return 0;
}

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// since early 2009, the FEValues objects try to be more efficient by only
// recomputing things like gradients of shape functions if the cell on which
// we are is not a translation of the previous one. in this series of tests we
// make sure that this actually works the way it's supposed to be. As opposed
// to the tests 01-10, this test uses a more complicated geometry that forces
// a MappingQ to switch between a Q1 while passing through the mesh
//
// this tests outputs the Jacobian of the transformation and the derivatives
// in real space in two locations inside the cell.

// To make sure the cell similarity code is used, only run the program with
// one thread.

#include <deal.II/base/function.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include "../tests.h"



template <int dim>
void
test(const Triangulation<dim> &tr)
{
  FE_Q<dim> fe(2);
  deallog << "FE=" << fe.get_name() << std::endl;

  MappingQ<dim> mapping(2);
  deallog << "Mapping=MappingQ2" << std::endl;


  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  std::vector<Point<dim>> points(2);
  for (unsigned int d = 0; d < dim; ++d)
    points[0][d] = 0.1;
  for (unsigned int d = 0; d < dim; ++d)
    points[1][d] = 0.85;

  const Quadrature<dim> quadrature(points);
  FEValues<dim>         fe_values(mapping,
                          fe,
                          quadrature,
                          update_gradients | update_jacobians);

  for (typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
       cell != dof.end();
       ++cell)
    {
      fe_values.reinit(cell);

      deallog << "Jacobians: ";
      for (const auto q : fe_values.quadrature_point_indices())
        {
          deallog << "[ ";
          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = 0; e < dim; ++e)
              deallog << fe_values.jacobian(q)[d][e] << ' ';
          deallog << " ] ";
        }
      deallog << std::endl;
      deallog << "Derivatives of shape function: ";
      for (const auto q : fe_values.quadrature_point_indices())
        {
          deallog << "[ ";
          for (unsigned int d = 0; d < dim; ++d)
            deallog << fe_values.shape_grad(fe.dofs_per_cell / 2, q)[d] << ' ';
          deallog << " ] ";
        }
      deallog << std::endl;
    }
}



void
test()
{
  // Create a mesh where the cell in the center uses a Q1 mapping but one cell
  // on the boundary uses a Q3 mapping from a hyper ball. The radius is
  // sqrt(1/9+1) which makes sure the mapping connects to the linear boundary
  // for the other cells.
  const int          dim = 2;
  Triangulation<dim> tr;
  GridGenerator::subdivided_hyper_cube(tr, 3, -1, 1);

  for (Triangulation<dim>::cell_iterator cell = tr.begin(); cell != tr.end();
       ++cell)
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      if (cell->face(f)->at_boundary() &&
          std::abs(cell->face(f)->center()[0] + 1.) < 1e-12 &&
          std::abs(cell->face(f)->center()[1]) < 1e-12)
        cell->face(f)->set_manifold_id(1);

  static const SphericalManifold<dim> boundary;
  tr.set_manifold(1, boundary);

  test<dim>(tr);
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(8);
  MultithreadInfo::set_thread_limit(1);

  deallog.attach(logfile);

  test();
}

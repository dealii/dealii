// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2017 by the deal.II authors
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

#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>



template <int dim>
void test (const Triangulation<dim> &tr)
{
  FE_Q<dim> fe(2);
  deallog << "FE=" << fe.get_name() << std::endl;

  MappingQ<dim> mapping(2);
  deallog << "Mapping=MappingQ2" << std::endl;


  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  std::vector<Point<dim> > points(2);
  for (unsigned int d=0; d<dim; ++d)
    points[0][d] = 0.1;
  for (unsigned int d=0; d<dim; ++d)
    points[1][d] = 0.85;

  const Quadrature<dim> quadrature(points);
  FEValues<dim> fe_values (mapping, fe, quadrature,
                           update_gradients | update_jacobians);

  for (typename DoFHandler<dim>::active_cell_iterator
       cell = dof.begin_active();
       cell != dof.end(); ++cell)
    {
      fe_values.reinit (cell);

      deallog << "Jacobians: ";
      for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
        {
          deallog << "[ ";
          for (unsigned int d=0; d<dim; ++d)
            for (unsigned int e=0; e<dim; ++e)
              deallog << fe_values.jacobian(q)[d][e] << " ";
          deallog << " ] ";
        }
      deallog << std::endl;
      deallog << "Derivatives of shape function: ";
      for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
        {
          deallog << "[ ";
          for (unsigned int d=0; d<dim; ++d)
            deallog << fe_values.shape_grad(fe.dofs_per_cell/2,q)[d] << " ";
          deallog << " ] ";
        }
      deallog << std::endl;
    }
}



void test()
{
  // Create a mesh where the cell in the center uses a Q1 mapping but one cell
  // on the boundary uses a Q3 mapping from a hyper ball. The radius is
  // sqrt(1/9+1) which makes sure the mapping connects to the linear boundary
  // for the other cells.
  const int dim = 2;
  Triangulation<dim> tr;
  GridGenerator::subdivided_hyper_cube(tr, 3, -1, 1);

  for (Triangulation<dim>::cell_iterator cell = tr.begin();
       cell != tr.end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if (cell->face(f)->at_boundary() &&
          std::abs(cell->face(f)->center()[0]+1.) < 1e-12 &&
          std::abs(cell->face(f)->center()[1]) < 1e-12)
        cell->face(f)->set_manifold_id(1);

  static const SphericalManifold<dim> boundary;
  tr.set_manifold (1, boundary);

  test<dim>(tr);
}


int main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (8);
  MultithreadInfo::set_thread_limit(1);

  deallog.attach(logfile);

  test();
}

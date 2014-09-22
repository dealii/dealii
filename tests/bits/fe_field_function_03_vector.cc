// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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



// Test the functionality of the laplacian in the FEFieldFunction class.

#include "../tests.h"
#include <fstream>

// all include files you need here
#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/numerics/vector_tools.h>

template <int dim>
void test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(9/dim);

  FE_Q<dim> fe(2);
  DoFHandler<dim> dh(tria);

  dh.distribute_dofs(fe);
  deallog << "Ndofs :" << dh.n_dofs() << std::endl;

  Functions::CosineFunction<dim> ff;

  Vector<double> v1(dh.n_dofs()), v2(dh.n_dofs());

  VectorTools::interpolate(dh, ff, v1);
  deallog << "V norm: " << v1.l2_norm() << std::endl;

  Functions::FEFieldFunction<dim, DoFHandler<dim>, Vector<double> >
  fef(dh, v1);

  //create the origin
  Point<dim> p;
  //compute the error of the laplacian in this point
  deallog << "Value of the laplacian in 0:" << std::endl;
  deallog << "correct value: " << ff.laplacian(p) <<", approximation: "<< fef.laplacian(p) << std::endl;

  //now we want to test the list version
  Point<dim> p1 = Point<dim>::unit_vector(0);
  p1 = p1 * 0.5;
  Point<dim> p2 = p1 * 0.5;
  std::vector<Point<dim> > vec;
  vec.push_back(p1);
  vec.push_back(p2);
  std::vector<double> values_c(2);
  std::vector<double> values_a(2);

  //get the laplacians at these two points
  ff.laplacian_list(vec, values_c);
  fef.laplacian_list(vec, values_a);
  deallog << "Value of the laplacian in 0.5*e1 and 0.25 * e1:" << std::endl;
  deallog << " correct values: " <<values_c[0] << " " << values_c[1]  <<", approximations: "<<values_a[0] << " " << values_a[1] << std::endl;
}

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1>();
  test<2>();
  test<3>();

  return 0;
}


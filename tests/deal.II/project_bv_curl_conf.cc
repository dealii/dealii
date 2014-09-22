// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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

#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/numerics/vector_tools.h>

std::ofstream logfile ("output");

template <int dim>
class BoundaryFunction: public Function<dim>
{
public:
  BoundaryFunction ();
  virtual void vector_value (const Point<dim> &p, Vector<double> &values) const;
};

template <int dim>
BoundaryFunction<dim>::BoundaryFunction (): Function<dim> (dim)
{
}

template <int dim>
void BoundaryFunction<dim>::vector_value (const Point<dim> &,
                                          Vector<double> &values) const
{
  for (unsigned int d = 0; d < dim; ++d)
    values (d) = d + 1.0;
}

template <int dim>
void test_boundary_values (const FiniteElement<dim> &fe)
{
  Triangulation<dim> triangulation;

  GridGenerator::subdivided_hyper_cube (triangulation, 2);

  DoFHandler<dim> dof_handler (triangulation);

  dof_handler.distribute_dofs (fe);

  BoundaryFunction<dim> boundary_function;
  ConstraintMatrix constraints;

  constraints.clear ();
  VectorTools::project_boundary_values_curl_conforming (dof_handler, 0, boundary_function, 0, constraints);
  constraints.close ();
  constraints.print (logfile);
}

int main ()
{
  deallog << std::setprecision (2);
  deallog.attach (logfile);
  deallog.depth_console (0);
  deallog.threshold_double (1e-12);

  FE_Nedelec<2> fe_2 (1);

  deallog << "dim=2:" << std::endl;
  test_boundary_values (fe_2);

  FE_Nedelec<3> fe_3 (1);

  deallog << "dim=3:" << std::endl;
  test_boundary_values (fe_3);
}

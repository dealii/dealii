// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------
//
// Same as project_bv_curl_conf.cc for hp

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_nedelec_sz.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_boundary.h>

#include "../tests.h"

template <int dim>
class BoundaryFunction : public Function<dim>
{
public:
  BoundaryFunction();
  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const;
};

template <int dim>
BoundaryFunction<dim>::BoundaryFunction()
  : Function<dim>(dim)
{}

template <int dim>
void
BoundaryFunction<dim>::vector_value(const Point<dim> &,
                                    Vector<double> &values) const
{
  for (unsigned int d = 0; d < dim; ++d)
    values(d) = d + 1.0;
}

template <int dim>
void
test_boundary_values(const FiniteElement<dim> &fe)
{
  Triangulation<dim> triangulation;

  GridGenerator::subdivided_hyper_cube(triangulation, 2);

  DoFHandler<dim> dof_handler(triangulation);

  dof_handler.distribute_dofs(fe);

  BoundaryFunction<dim>     boundary_function;
  AffineConstraints<double> constraints;

  constraints.clear();
  VectorTools::project_boundary_values_curl_conforming_l2(
    dof_handler,
    0,
    boundary_function,
    0,
    constraints,
    hp::StaticMappingQ1<dim>::mapping_collection);
  constraints.close();
  constraints.print(deallog.get_file_stream());
}

int
main()
{
  initlog();

  {
    FE_Nedelec<2> fe_2(1);

    deallog << "dim=2:" << std::endl;
    test_boundary_values(fe_2);

    FE_Nedelec<3> fe_3(1);

    deallog << "dim=3:" << std::endl;
    test_boundary_values(fe_3);
  }

  {
    FE_NedelecSZ<2> fe_2(1);

    deallog << "dim=2:" << std::endl;
    test_boundary_values(fe_2);

    FE_NedelecSZ<3> fe_3(1);

    deallog << "dim=3:" << std::endl;
    test_boundary_values(fe_3);
  }
}

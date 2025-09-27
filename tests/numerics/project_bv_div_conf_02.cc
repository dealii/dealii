// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test VectorTools::project_boundary_values_div_conforming for the
// case that the DoFHandler contains is more than one FE_RaviartThomas element.

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "../tests.h"


template <int dim>
class BoundaryFunctionDisp : public dealii::Function<dim>
{
public:
  BoundaryFunctionDisp()
    : Function<dim>(dim)
  {}

  virtual double
  value(const Point<dim> & /*point*/,
        const unsigned int /*component*/) const override
  {
    return 1.;
  }
};

template <int dim>
class BoundaryFunctionVelo : public dealii::Function<dim>
{
public:
  BoundaryFunctionVelo()
    : Function<dim>(dim)
  {}

  virtual double
  value(const Point<dim> & /*point*/,
        const unsigned int /* component */) const override
  {
    return -1.;
  }
};

template <int dim>
void
test_boundary_values(const FiniteElement<dim> &fe)
{
  Triangulation<dim> triangulation;

  GridGenerator::hyper_cube(triangulation, -1., 1., true);
  triangulation.refine_global(1);

  DoFHandler<dim> dof_handler(triangulation);

  dof_handler.distribute_dofs(fe);

  BoundaryFunctionDisp<dim> boundary_function_disp;
  BoundaryFunctionVelo<dim> boundary_function_velo;

  AffineConstraints<double> constraints;

  {
    constraints.clear();
    VectorTools::project_boundary_values_div_conforming(
      dof_handler,
      0, /*first_vector_component*/
      boundary_function_disp,
      0, /*bdry_id*/
      constraints,
      StaticMappingQ1<dim>::mapping);
    VectorTools::project_boundary_values_div_conforming(
      dof_handler,
      dim, /*first_vector_component*/
      boundary_function_velo,
      0, /*bdry_id*/
      constraints,
      StaticMappingQ1<dim>::mapping);
    constraints.close();
  }

  constraints.print(deallog.get_file_stream());
}

int
main()
{
  initlog();
  {
    constexpr unsigned int dim = 2;
    FE_RaviartThomas<dim>  u(1);
    FE_DGQ<dim>            p(1);
    FESystem<dim>          fesys(u, 2, p, 1);
    test_boundary_values(fesys);
  }

  {
    constexpr unsigned int dim = 3;
    FE_RaviartThomas<dim>  u(1);
    FE_DGQ<dim>            p(1);
    FESystem<dim>          fesys(u, 2, p, 1);
    test_boundary_values(fesys);
  }
}

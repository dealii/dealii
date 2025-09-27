// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Same as project_bv_curl_conf.cc for std::complex<double>

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_nedelec_sz.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim>
class BoundaryFunction : public Function<dim, std::complex<double>>
{
public:
  BoundaryFunction();
  virtual void
  vector_value(const Point<dim> &p, Vector<std::complex<double>> &values) const;
};

template <int dim>
BoundaryFunction<dim>::BoundaryFunction()
  : Function<dim, std::complex<double>>(dim)
{}

template <int dim>
void
BoundaryFunction<dim>::vector_value(const Point<dim> &,
                                    Vector<std::complex<double>> &values) const
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

  BoundaryFunction<dim>                   boundary_function;
  AffineConstraints<std::complex<double>> constraints;

  constraints.clear();
  VectorTools::project_boundary_values_curl_conforming_l2(
    dof_handler,
    0,
    boundary_function,
    0,
    constraints,
    StaticMappingQ1<dim>::mapping);
  constraints.close();
  constraints.print(deallog.get_file_stream());
}

int
main()
{
  initlog();

  {
    deallog << "Tests with FE_NEdelec:" << std::endl;

    FE_Nedelec<2> fe_2(1);

    deallog << "dim=2:" << std::endl;
    test_boundary_values(fe_2);

    FE_Nedelec<3> fe_3(1);

    deallog << "dim=3:" << std::endl;
    test_boundary_values(fe_3);
  }

  {
    deallog << "Tests with FE_NEdelecSZ:" << std::endl;

    FE_NedelecSZ<2> fe_2(1);

    deallog << "dim=2:" << std::endl;
    test_boundary_values(fe_2);

    FE_NedelecSZ<3> fe_3(1);

    deallog << "dim=3:" << std::endl;
    test_boundary_values(fe_3);
  }
}

// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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

// Test the functionality of the CutOffTensorProductFunction

#include <deal.II/base/function_lib.h>
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <map>

#include "../tests.h"

using namespace Functions;

template <int dim>
struct Domain
{
  Domain()
    : fe(1)
    , dh(tria)
    , quad(5)
  {
    deallog << "Testing dim = " << dim << std::endl;
    GridGenerator::hyper_cube(tria, -2, 2);
    tria.refine_global(5 - dim);
    dh.distribute_dofs(fe);
    zero.reinit(dh.n_dofs());
    cell_integral.reinit(tria.n_active_cells());
  }

  template <class Function>
  void
  integrate(const Function &fun)
  {
    deallog << "Integrating " << Utilities::type_to_string(fun)
            << " -2,2 cube: ";

    VectorTools::integrate_difference(
      dh, zero, fun, cell_integral, quad, VectorTools::L1_norm);
    const double integral =
      VectorTools::compute_global_error(tria,
                                        cell_integral,
                                        VectorTools::L1_norm);

    deallog << integral << std::endl;
  }

  Triangulation<dim> tria;
  FE_Q<dim>          fe;
  DoFHandler<dim>    dh;
  Vector<double>     zero;
  Vector<double>     cell_integral;
  QGauss<dim>        quad;
};

template <int dim, template <int> class TestCutOffFunction>
void
test()
{
  static Domain<dim> domain;

  CutOffFunctionTensorProduct<dim> fun(
    1, Point<dim>(), 1, CutOffFunctionBase<dim>::no_component, true);
  fun.template set_base<TestCutOffFunction>();

  deallog << "Center: " << fun.get_center() << std::endl
          << "Radius: " << fun.get_radius() << std::endl;

  domain.integrate(fun);

  Point<dim> new_center;
  for (unsigned int i = 0; i < dim; ++i)
    new_center[i] = .5;

  fun.set_center(new_center);
  fun.set_radius(.5);

  deallog << "Center: " << fun.get_center() << std::endl
          << "Radius: " << fun.get_radius() << std::endl;

  domain.integrate(fun);
}

int
main()
{
  initlog(1);

  test<1, Functions::CutOffFunctionLinfty>();
  test<2, Functions::CutOffFunctionLinfty>();
  test<3, Functions::CutOffFunctionLinfty>();

  test<1, Functions::CutOffFunctionW1>();
  test<2, Functions::CutOffFunctionW1>();
  test<3, Functions::CutOffFunctionW1>();

  test<1, Functions::CutOffFunctionCinfty>();
  test<2, Functions::CutOffFunctionCinfty>();
  test<3, Functions::CutOffFunctionCinfty>();

  test<1, Functions::CutOffFunctionC1>();
  test<2, Functions::CutOffFunctionC1>();
  test<3, Functions::CutOffFunctionC1>();
}

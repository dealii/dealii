// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test integrate_difference for complex-valued vectors. This test is
// like integrate_difference_01_complex_04, but compares the interpolated
// solution to the exact solution.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/signaling_nan.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



// x+y(+z), x^2+y^2 (, z+xy) times (1+i)/sqrt(2)
template <int dim>
class Ref : public Function<dim, std::complex<double>>
{
public:
  Ref()
    : Function<dim, std::complex<double>>(dim)
  {}

  std::complex<double>
  value(const Point<dim> &p, const unsigned int c) const
  {
    if (c == 0)
      return (p[0] + p[1] + ((dim == 3) ? p[2] : 0.0)) *
             std::complex<double>(1. / std::sqrt(2.), 1. / std::sqrt(2.));
    if (c == 1)
      return (p[0] * p[0] + p[1] * p[1]) *
             std::complex<double>(1. / std::sqrt(2.), 1. / std::sqrt(2.));
    if (c == 2)
      return (p[2] + p[0] * p[1]) *
             std::complex<double>(1. / std::sqrt(2.), 1. / std::sqrt(2.));
    else
      return numbers::signaling_nan<double>();
  }
};



template <int dim>
void
test(VectorTools::NormType norm, double value)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);

  FESystem<dim>   fe(FE_Q<dim>(4), dim);
  DoFHandler<dim> dofh(tria);
  dofh.distribute_dofs(fe);

  Vector<std::complex<double>> solution(dofh.n_dofs());
  VectorTools::interpolate(dofh, Ref<dim>(), solution);

  Vector<double> cellwise_errors(tria.n_active_cells());
  VectorTools::integrate_difference(
    dofh, solution, Ref<dim>(), cellwise_errors, QGauss<dim>(5), norm);

  const double error = cellwise_errors.l2_norm();

  const double difference = std::abs(error - value);
  deallog << "computed: " << error << " expected: " << value
          << " difference: " << difference << std::endl;
  Assert(difference < 1e-10, ExcMessage("Error in integrate_difference"));
}

template <int dim>
void
test()
{
  deallog << "L2_norm:" << std::endl;
  // sqrt(\int_\Omega f^2) = sqrt(\int (x+y)^2+(x^2+y^2)^2)
  test<dim>(VectorTools::L2_norm, 0.);

  deallog << "OK" << std::endl;
}


int
main(int argc, char **argv)
{
  initlog();
  test<2>();
}

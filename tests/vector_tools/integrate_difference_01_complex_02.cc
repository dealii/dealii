// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2018 by the deal.II authors
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


// Test integrate_difference for complex-valued vectors. This test is
// like integrate_difference_01, but multiplies the real-valued
// function by sqrt(-1) and stores everything using complex data types

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



// x+y(+z), x^2+y^2 (, z+xy) times i
// div = (1+2y (+1))i
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
      return std::complex<double>(0, p[0] + p[1] + ((dim == 3) ? p[2] : 0.0));
    if (c == 1)
      return std::complex<double>(0, p[0] * p[0] + p[1] * p[1]);
    if (c == 2)
      return std::complex<double>(0, p[2] + p[0] * p[1]);
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
    dofh,
    solution,
    Functions::ZeroFunction<dim, std::complex<double>>(dim),
    cellwise_errors,
    QGauss<dim>(5),
    norm);

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
  deallog << "Hdiv_seminorm:" << std::endl;
  // sqrt(\int (div f)^2 = sqrt(\int (1+2y)^2)
  test<dim>(VectorTools::Hdiv_seminorm, std::sqrt(13.0 / 3.0));
  deallog << "L2_norm:" << std::endl;
  // sqrt(\int_\Omega f^2) = sqrt(\int (x+y)^2+(x^2+y^2)^2)
  test<dim>(VectorTools::L2_norm, std::sqrt(161.0 / 90.0));
  deallog << "H1_seminorm:" << std::endl;
  // sqrt( sum | d/dxi f |_0^2 ) = sqrt( sum \int   )
  test<dim>(VectorTools::H1_seminorm, std::sqrt(14.0 / 3.0));

  deallog << "OK" << std::endl;
}


int
main(int argc, char **argv)
{
  initlog();
  test<2>();
}

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test integrate_difference with Hdiv_seminorm. Specifically this test
// deals with the case when the ComponentSelectFunction is used and
// the components selected for the Hdiv seminorm computations are not
// the first dim components in the Finite Element.

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



// First dim components:
// f_x = x^2+y(+z), f_y = x^2+y^2, f_z = z+xy
// div(f) = 2x+2y   in 2d
// div(f) = 2x+2y+1 in 3d
// Second dim components:
// g(x,y,z) = 2 * f(x,y,z)
template <int dim>
class Ref : public Function<dim>
{
public:
  Ref()
    : Function<dim>(2 * dim)
  {}

  void
  vector_value(const Point<dim> &p, Vector<double> &values) const
  {
    switch (dim)
      {
        case 2:
          values[0] = p[0] * p[0] + p[1];
          values[1] = p[0] * p[0] + p[1] * p[1];
          values[2] = 2 * (p[0] * p[0] + p[1]);
          values[3] = 2 * (p[0] * p[0] + p[1] * p[1]);
          break;
        case 3:
          values[0] = p[0] * p[0] + p[1] + p[2];
          values[1] = p[0] * p[0] + p[1] * p[1];
          values[2] = p[2] + p[0] * p[1];
          values[3] = 2 * (p[0] * p[0] + p[1] + p[2]);
          values[4] = 2 * (p[0] * p[0] + p[1] * p[1]);
          values[5] = 2 * (p[2] + p[0] * p[1]);
          break;
        default:
          DEAL_II_NOT_IMPLEMENTED();
      }
  }
};


template <int dim>
void
test(VectorTools::NormType norm, double value)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);

  FESystem<dim>   fe(FE_Q<dim>(3), dim, FE_Q<dim>(3), dim);
  DoFHandler<dim> dofh(tria);
  dofh.distribute_dofs(fe);

  Vector<double> solution(dofh.n_dofs());
  VectorTools::interpolate(dofh, Ref<dim>(), solution);

  Vector<double> cellwise_errors(tria.n_active_cells());

  const ComponentSelectFunction<dim> mask_2(std::make_pair(dim, 2 * dim),
                                            2 * dim);
  VectorTools::integrate_difference(dofh,
                                    solution,
                                    Functions::ZeroFunction<dim>(2 * dim),
                                    cellwise_errors,
                                    QGauss<dim>(5),
                                    norm);

  double error = cellwise_errors.l2_norm();

  const double difference_1 = std::abs(error - value);
  deallog << "computed: " << error << " expected: " << value
          << " difference: " << difference_1 << std::endl;
  Assert(difference_1 < 1e-10,
         ExcMessage("Error in integrate_difference, first components"));

  VectorTools::integrate_difference(dofh,
                                    solution,
                                    Functions::ZeroFunction<dim>(2 * dim),
                                    cellwise_errors,
                                    QGauss<dim>(5),
                                    norm,
                                    &mask_2);

  error                     = cellwise_errors.l2_norm();
  const double difference_2 = std::abs(error - 2.0 * value);
  deallog << "computed: " << error << " expected: " << 2.0 * value
          << " difference: " << difference_2 << std::endl;
  Assert(difference_2 < 1e-10,
         ExcMessage("Error in integrate_difference, second components"));
}


template <int dim>
void
test()
{
  deallog << dim << " dimensions, Hdiv_seminorm:" << std::endl;
  double true_value = 0;
  switch (dim)
    {
      case 2:
        true_value = std::sqrt(14.0 / 3.0);
        break;
      case 3:
        true_value = std::sqrt(29.0 / 3.0);
        break;
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }

  test<dim>(VectorTools::Hdiv_seminorm, true_value);
  deallog << "OK" << std::endl;
}


int
main(int argc, char **argv)
{
  initlog();
  test<2>();
  test<3>();
}

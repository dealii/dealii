// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test integrate_difference
//
// The VectorTools::integrate_difference() function allows users
// to provide a weight function that can also serve as a component mask
// to select individual components of the solution vector for error
// computation. For components not selected, such a mask would then
// simply be zero.
//
// In some cases, the solution vector contains NaN numbers, for example
// when one uses the FE_FaceQ element for certain components of the
// solution vector and uses a quadrature formula for error evaluation
// that has quadrature points in the interior of the cell. For any
// "regular" solution component for which the component mask has a zero
// weight, the value of that component will be multiplied by zero and
// consequently does not add anything to the error computation. However,
// if the NaNs of a FE_FaceQ are multiplied with zero weights, the result
// is still a NaN, and adding it to the values times weights of the other
// components results in NaNs -- in effect rendering it impossible to get
// any information out of the VectorTools::integrate_difference()
// function if one of the finite elements involved is FE_FaceQ.
//
// This is now fixed by simply skipping vector components for which the
// weight vector is zero. This has the same result as before for all
// "normal" situations, but also properly skips the NaN case outlined
// above.


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/signaling_nan.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
void
test(VectorTools::NormType norm, double exp = 2.0)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  FESystem<dim>   fe(FE_Q<dim>(1), 1, FE_FaceQ<dim>(1), 1);
  DoFHandler<dim> dofh(tria);
  dofh.distribute_dofs(fe);

  // Create a zero vector. (The values in there are not important,
  // just that the difference between this and later the exact
  // solution happens to have signaling_nans in it. It will have these
  // values because FE_FaceQ refuses to fill FEValues values --
  // because it's a face element, so this refusal actually makes sense
  // -- and thus leaves the signaling NaNs from the original invalid
  // initialization.)
  Vector<double> solution(dofh.n_dofs());

  Vector<double> cellwise_errors(tria.n_active_cells());
  QIterated<dim> quadrature(QTrapezoid<1>(), 2);

  const ComponentSelectFunction<dim> component_mask(1, 2);
  VectorTools::integrate_difference(dofh,
                                    solution,
                                    Functions::ZeroFunction<dim>(2),
                                    cellwise_errors,
                                    quadrature,
                                    norm,
                                    &component_mask,
                                    exp);

  const double error =
    VectorTools::compute_global_error(tria, cellwise_errors, norm, exp);

  deallog << "computed: " << error << std::endl;
}



template <int dim>
void
test()
{
  deallog << "L2_norm:" << std::endl;
  test<dim>(VectorTools::L2_norm);

  deallog << "H1_seminorm:" << std::endl;
  test<dim>(VectorTools::H1_seminorm);

  deallog << "H1_norm:" << std::endl;
  test<dim>(VectorTools::H1_norm);

  deallog << "L1_norm:" << std::endl;
  test<dim>(VectorTools::L1_norm);

  deallog << "Linfty_norm:" << std::endl;
  test<dim>(VectorTools::Linfty_norm);

  deallog << "mean:" << std::endl;
  test<dim>(VectorTools::mean);

  deallog << "Lp_norm:" << std::endl;
  test<dim>(VectorTools::Lp_norm, 3.0);

  deallog << "W1p_seminorm:" << std::endl;
  test<dim>(VectorTools::W1p_seminorm, 3.0);

  deallog << "OK" << std::endl;
}


int
main(int argc, char **argv)
{
  initlog();
  test<3>();
}

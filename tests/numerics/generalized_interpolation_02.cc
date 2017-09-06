// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

// Check that VectorTools::interpolate correctly recovers a
// constant vector field for H1, Hdiv and Hcurl conforming elements.

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/fe_field_function.h>

#include "../tests.h"

template <int dim>
void test(const FiniteElement<dim> &fe,
          const unsigned int n_comp,
          const unsigned int order_mapping,
          bool distort_mesh)
{
  deallog << "dim " << dim << " " << fe.get_name() << std::endl;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -0.3, 0.7);
  triangulation.refine_global(dim == 2 ? 2 : 1);
  if (distort_mesh)
    GridTools::distort_random(0.03, triangulation);

  ConstantFunction<dim> f (1.0, n_comp);

  MappingQ<dim> mapping(order_mapping);

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  Vector<double> interpolant(dof_handler.n_dofs());
  VectorTools::interpolate(mapping, dof_handler, f, interpolant);

  // Check that the interpoland really returns 1.0...

  Functions::FEFieldFunction<dim> f2(dof_handler, interpolant, mapping);
  deallog << "Function value at (0.0,0.0): ";
  for (unsigned int i = 0; i < n_comp; ++i)
    deallog << f2.value(Point<dim>(), i) << " ";
  deallog << std::endl;

  // Check that VectorTools::interpolate is in fact a
  // projection, i.e. applying the interpolation twice results in the same
  // vector:

  Vector<double> interpolant2(dof_handler.n_dofs());
  VectorTools::interpolate(
    mapping, dof_handler, f2, interpolant2);

  interpolant2 -= interpolant;
  deallog << "Check projection property: " << interpolant2.linfty_norm()
          << std::endl;
}


int main ()
{
  deallog.depth_console(3);

  test<2>(FE_Q<2>(1), 1, 1, false);
  test<2>(FE_Q<2>(2), 1, 2, false);
  test<3>(FE_Q<3>(3), 1, 3, false);

  test<2>(FE_Q<2>(1), 1, 1, true);
  test<2>(FE_Q<2>(2), 1, 2, true);
  test<3>(FE_Q<3>(3), 1, 3, true);

  test<2>(FE_RaviartThomas<2>(0), 2, 1, false);
  test<2>(FE_RaviartThomas<2>(1), 2, 2, false);
  test<2>(FE_RaviartThomas<2>(2), 2, 3, false);
  test<3>(FE_RaviartThomas<3>(0), 3, 1, false);
  test<3>(FE_RaviartThomas<3>(1), 3, 2, false);

  test<2>(FE_RaviartThomas<2>(0), 2, 1, true);
  test<2>(FE_RaviartThomas<2>(1), 2, 2, true);
  test<2>(FE_RaviartThomas<2>(2), 2, 3, true);
  // lowest order RT in 3D does not contain constant 1 function on a
  // distorted mesh.
  test<3>(FE_RaviartThomas<3>(1), 3, 2, true);

  // FIXME: Reenable, when FE_Nedelec is fixed :-/
  // test<2>(FE_Nedelec<2>(0), 2, 1, false);
  // test<2>(FE_Nedelec<2>(1), 2, 2, false);
  // test<2>(FE_Nedelec<2>(2), 2, 3, false);
  // test<2>(FE_Nedelec<2>(0), 2, 1, true);
  // test<2>(FE_Nedelec<2>(1), 2, 2, true);
  // test<2>(FE_Nedelec<2>(2), 2, 3, true);
}

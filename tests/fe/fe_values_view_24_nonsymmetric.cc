// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2018 by the deal.II authors
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



// output value of divergence field for Tensor extractors
// for a simple mesh with linear elements.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"



template <int dim>
void
test(const Triangulation<dim> &tr, const FiniteElement<dim> &fe)
{
  deallog << "FE=" << fe.get_name() << std::endl;

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  Vector<double> fe_function(dof.n_dofs());
  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    fe_function(i) = i + 1;

  const QGauss<dim> quadrature(2);
  FEValues<dim>     fe_values(fe, quadrature, update_values | update_gradients);
  fe_values.reinit(dof.begin_active());

  // let the FEValues object compute the
  // divergences at quadrature points
  std::vector<Tensor<1, dim>>   divergences(quadrature.size());
  FEValuesExtractors::Tensor<2> extractor(0);
  fe_values[extractor].get_function_divergences(fe_function, divergences);

  for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
    {
      deallog << "i=" << i << std::endl;

      for (unsigned int q = 0; q < quadrature.size(); ++q)
        deallog << "  q_point=" << q << std::endl
                << "    div= " << fe_values[extractor].divergence(i, q)
                << std::endl;
    }
}



template <int dim>
void
test_hyper_sphere()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_ball(tr);

  static const SphericalManifold<dim> boundary;
  tr.set_manifold(0, boundary);

  FESystem<dim> fe(FE_Q<dim>(1), Tensor<2, dim>::n_independent_components);
  test(tr, fe);
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test_hyper_sphere<2>();
  test_hyper_sphere<3>();
}

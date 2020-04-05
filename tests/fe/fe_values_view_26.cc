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



// a bit like _25, but test for the curl of a function. there was a
// bug in get_function_curls

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



Tensor<1, 1>
curl(const Tensor<2, 2> &grads)
{
  return Point<1>(grads[1][0] - grads[0][1]);
}


Tensor<1, 3>
curl(const Tensor<2, 3> &grads)
{
  return Point<3>(grads[2][1] - grads[1][2],
                  grads[0][2] - grads[2][0],
                  grads[1][0] - grads[0][1]);
}



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
  FEValues<dim>     fe_values(fe,
                          quadrature,
                          update_values | update_gradients |
                            update_quadrature_points);
  fe_values.reinit(dof.begin_active());

  // let the FEValues object compute the
  // divergences at quadrature points
  std::vector<typename dealii::internal::CurlType<dim>::type> curls(
    quadrature.size());
  std::vector<Tensor<2, dim>> grads(quadrature.size());
  FEValuesExtractors::Vector  extractor(0);
  fe_values[extractor].get_function_curls(fe_function, curls);
  fe_values[extractor].get_function_gradients(fe_function, grads);

  // now compare
  for (unsigned int q = 0; q < quadrature.size(); ++q)
    {
      deallog << "  curls[q]= " << curls[q] << std::endl
              << "  grads[q]= " << grads[q] << std::endl;
      Assert((curl(grads[q]) - curls[q]).norm() <= 1e-10, ExcInternalError());
    }
}



template <int dim>
void
test_hyper_cube()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);

  FESystem<dim> fe(FE_Q<dim>(1), dim);
  test(tr, fe);
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test_hyper_cube<2>();
  test_hyper_cube<3>();
}

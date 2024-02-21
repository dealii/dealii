// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

/*
 * Test if FEPointEval::evaluate() works for FESystems with FE_Q and FE_DGQ
 * using the face path.
 *
 * We check FEPointEval with n_components=2 and first_selected_component=0
 * and FEPointEval with n_components=1 and first_selected_component=1.
 */

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/non_matching/mapping_info.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
class AnalyticalFunction : public Function<dim>
{
public:
  AnalyticalFunction(unsigned int n_components)
    : Function<dim>(n_components)
  {}

  double
  value(const Point<dim> &p, const unsigned int component) const final
  {
    return p[component];
  }
};

template <int first_selected_component, int dim, typename Number = double>
void
test(const FiniteElement<dim> &fe)
{
  constexpr unsigned int n_components = dim;

  // setup tria (which is one reference cell)
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);

  MappingQ<dim> mapping(1);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);


  Vector<double> vector(dof_handler.n_dofs());
  VectorTools::interpolate(mapping,
                           dof_handler,
                           AnalyticalFunction<dim>(n_components),
                           vector);

  // setup FEPointEvaluation which uses the face path
  std::vector<std::vector<Quadrature<dim - 1>>> quad_vec;
  quad_vec.emplace_back(
    std::vector<Quadrature<dim - 1>>(4, QGauss<dim - 1>(7)));
  dealii::NonMatching::MappingInfo<dim> mapping_info(mapping, update_values);
  mapping_info.reinit_faces(tria.active_cell_iterators(), quad_vec);
  FEFacePointEvaluation<dim - first_selected_component, dim, dim, Number>
    fe_point_eval(mapping_info, fe, true, first_selected_component);

  std::vector<Number> buffer(fe.dofs_per_cell);

  const auto         cell = *dof_handler.active_cell_iterators().begin();
  const unsigned int face = 3;

  fe_point_eval.reinit(cell->active_cell_index(), face);
  cell->get_dof_values(vector, buffer.begin(), buffer.end());

  fe_point_eval.evaluate(buffer, EvaluationFlags::values);

  for (const unsigned int q : fe_point_eval.quadrature_point_indices())
    {
      deallog << "Value at q " << q << ": " << fe_point_eval.get_value(q)
              << std::endl;
    }
  fe_point_eval.submit_value(fe_point_eval.get_value(0), 0);

  deallog << "Tested DoF values " << std::endl;
  fe_point_eval.test_and_sum(buffer, EvaluationFlags::values);

  for (const auto &b : buffer)
    deallog << "  " << b << std::endl;

  deallog << std::endl;
}


int
main(int argc, char **argv)
{
  initlog();

  const unsigned int degree = 2;
  const unsigned int dim    = 2;

  deallog << "FE_Q" << std::endl;
  FESystem<dim> fe_q(FE_Q<dim>(degree), dim);
  deallog << "All components" << std::endl;
  test<0>(fe_q);
  deallog << "Second component" << std::endl;
  test<1>(fe_q);

  deallog << "FE_DGQ" << std::endl;
  FESystem<dim> fe_dgq(FE_DGQ<dim>(degree), dim);
  deallog << "All components" << std::endl;
  test<0>(fe_dgq);
  deallog << "Second component" << std::endl;
  test<1>(fe_dgq);

  return 0;
}

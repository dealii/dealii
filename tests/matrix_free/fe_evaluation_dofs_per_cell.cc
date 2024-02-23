// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test the basic member variables of FEEvaluation

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"


template <typename FEEval>
void
print_info(const FEEval &eval)
{
  // copy static variables to int to avoid taking references (with possibly
  // undefined references) when inside deallog::operator<<
  unsigned int v = FEEval::dimension;
  deallog << "FEEvaluation::dimension: " << v << std::endl;
  v = FEEval::n_components;
  deallog << "FEEvaluation::n_components: " << v << std::endl;
  v = FEEval::static_n_q_points;
  deallog << "FEEvaluation::static_n_q_points: " << v << std::endl;
  v = FEEval::static_dofs_per_component;
  deallog << "FEEvaluation::static_dofs_per_component: " << v << std::endl;
  v = FEEval::tensor_dofs_per_cell;
  deallog << "FEEvaluation::tensor_dofs_per_cell: " << v << std::endl;
  v = FEEval::static_dofs_per_cell;
  deallog << "FEEvaluation::static_dofs_per_cell: " << v << std::endl;
  deallog << "FEEvaluation::dofs_per_component: " << eval.dofs_per_component
          << std::endl;
  deallog << "FEEvaluation::dofs_per_cell: " << eval.dofs_per_cell << std::endl;
  deallog << "FEEvaluation::n_q_points: " << eval.n_q_points << std::endl;
}



template <int dim>
void
test()
{
  const unsigned int                degree = 1;
  FESystem<dim>                     fe1(FE_Q<dim>(degree + 1), dim);
  FESystem<dim>                     fe2(FE_Q<dim>(degree + 1), 1);
  FESystem<dim>                     fe3(FE_DGP<dim>(degree + 1), 2);
  FESystem<dim>                     fe4(FE_Q<dim>(degree), 1);
  std::vector<FiniteElement<dim> *> fes{&fe1, &fe2, &fe3, &fe4};

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  DoFHandler<dim> dof(tria);
  for (unsigned int i = 0; i < fes.size(); ++i)
    {
      deallog << "Checking " << fes[i]->get_name() << std::endl;
      dof.distribute_dofs(*fes[i]);
      MatrixFree<dim> matrix_free;
      matrix_free.reinit(MappingQ1<dim>{},
                         dof,
                         AffineConstraints<double>(),
                         QGauss<1>(degree + 3),
                         typename MatrixFree<dim>::AdditionalData());
      if (i < 2)
        {
          FEEvaluation<dim, degree + 1, degree + 3, dim> phi(matrix_free);
          print_info(phi);
        }
      if (i != 2)
        {
          FEEvaluation<dim, -1, 0, dim> phi(matrix_free);
          print_info(phi);
        }
      if (i == 2)
        {
          FEEvaluation<dim, degree + 1, degree + 3, 2> phi(matrix_free);
          print_info(phi);
        }
    }
}



int
main()
{
  initlog();
  test<1>();
  test<2>();
  test<3>();
}

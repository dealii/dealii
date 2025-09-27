// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Assemble the sparsity pattern for a Q1 finite element space and an iso Q
// element space with the same number of subdivisions, and check we get the same
// sparsity patterns.

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>

#include "../tests.h"


template <int dim>
void
test(const unsigned int n_subdivisions = 1)
{
  deallog << "dimension: " << dim << ", n_subdivisions = " << n_subdivisions
          << std::endl;
  Triangulation<dim> triangulation;
  FE_Q_iso_Q1<dim>   fe(n_subdivisions);
  DoFHandler<dim>    dof_handler(triangulation);
  SparsityPattern    sparsity_pattern;

  GridGenerator::hyper_cube(triangulation);
  dof_handler.distribute_dofs(fe);

  Triangulation<dim> triangulation_q1;
  DoFHandler<dim>    dof_handler_q1(triangulation_q1);
  SparsityPattern    sparsity_pattern_q1;
  FE_Q<dim>          fe_q1(1);
  GridGenerator::subdivided_hyper_cube(triangulation_q1, n_subdivisions);
  dof_handler_q1.distribute_dofs(fe_q1);

  std::vector<Point<dim>> support_points(dof_handler.n_dofs());
  DoFTools::map_dofs_to_support_points(StaticMappingQ1<dim>::mapping,
                                       dof_handler,
                                       support_points);


  std::vector<Point<dim>> support_points_q1(dof_handler_q1.n_dofs());
  DoFTools::map_dofs_to_support_points(StaticMappingQ1<dim>::mapping,
                                       dof_handler_q1,
                                       support_points_q1);

  AssertDimension(support_points.size(), support_points_q1.size());

  std::vector<types::global_dof_index> new_numbers(dof_handler_q1.n_dofs());

  for (unsigned int i = 0; i < support_points_q1.size(); ++i)
    {
      const Point<dim> &p = support_points_q1[i];
      for (unsigned int j = 0; j < support_points.size(); ++j)
        {
          if (p.distance(support_points[j]) < 1e-10)
            {
              new_numbers[i] = j;
              break;
            }
        }
    }
  dof_handler_q1.renumber_dofs(new_numbers);

  {
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    dsp.print(deallog.get_file_stream());
    sparsity_pattern.copy_from(dsp);
    sparsity_pattern.print(deallog.get_file_stream());
    deallog << std::endl;
  }

  {
    DynamicSparsityPattern dsp(dof_handler_q1.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler_q1, dsp);
    sparsity_pattern_q1.copy_from(dsp);
  }

  // Check that the two sparsity patterns are the same
  auto el = sparsity_pattern.begin();
  for (const auto &el_q1 : sparsity_pattern_q1)
    {
      if ((el->row() != el_q1.row()) || (el->column() != el_q1.column()))
        {
          deallog << "el_q1    : " << el_q1.row() << ", " << el_q1.column()
                  << std::endl;
          deallog << "el_q_iso : " << el->row() << ", " << el->column()
                  << std::endl;
        }
      ++el;
    }
}


int
main()
{
  initlog();

  test<1>(1);
  test<1>(2);
  test<1>(3);

  test<2>(1);
  test<2>(2);
  test<2>(3);

  test<3>(1);
  test<3>(2);
  test<3>(3);
}

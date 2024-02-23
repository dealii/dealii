// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test FETools::get_projection_matrix() for FESystem, FE_RaviartThomas,
// and FE_Nedelec.

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>

#include <iostream>

#include "../tests.h"

template <int dim>
void
test(const FiniteElement<dim> &fe_0,
     const FiniteElement<dim> &fe_1,
     const bool                do_scalar_test)
{
  {
    FullMatrix<double> matrix(fe_1.n_dofs_per_cell(), fe_0.n_dofs_per_cell());
    FETools::get_projection_matrix(fe_0, fe_1, matrix);
    matrix.print_formatted(deallog.get_file_stream(), 8, false, 5, "0");
    deallog << std::endl;

    if (do_scalar_test)
      {
        FullMatrix<double> matrix_1(fe_1.base_element(0).n_dofs_per_cell(),
                                    fe_0.base_element(0).n_dofs_per_cell());
        FETools::get_projection_matrix(fe_0.base_element(0),
                                       fe_1.base_element(0),
                                       matrix_1);

        for (unsigned i = 0; i < matrix.m(); ++i)
          for (unsigned j = 0; j < matrix.n(); ++j)
            {
              const auto component_index_i = fe_1.system_to_component_index(i);
              const auto component_index_j = fe_0.system_to_component_index(j);

              if (component_index_i.first != component_index_j.first)
                {
                  Assert(std::abs(matrix[i][j]) < 1e-8, ExcInternalError());
                }
              else
                {
                  Assert(std::abs(matrix[i][j] -
                                  matrix_1[component_index_i.second]
                                          [component_index_j.second]) < 1e-8,
                         ExcInternalError());
                }
            }
      }
  }
}


template <int dim>
void
test(const unsigned int fe_degree_0, const unsigned int fe_degree_1)
{
  test<dim>(FESystem<dim>(FE_Q<dim>(fe_degree_0), dim),
            FESystem<dim>(FE_Q<dim>(fe_degree_1), dim),
            true);
  test<dim>(FESystem<dim>(FE_Q<dim>(fe_degree_1), dim),
            FESystem<dim>(FE_Q<dim>(fe_degree_0), dim),
            true);

  test<dim>(FE_RaviartThomas<dim>(fe_degree_0),
            FE_RaviartThomas<dim>(fe_degree_1),
            false);
  test<dim>(FE_RaviartThomas<dim>(fe_degree_1),
            FE_RaviartThomas<dim>(fe_degree_0),
            false);

  test<dim>(FE_RaviartThomasNodal<dim>(fe_degree_0),
            FE_RaviartThomas<dim>(fe_degree_1),
            false);
  test<dim>(FE_RaviartThomasNodal<dim>(fe_degree_1),
            FE_RaviartThomas<dim>(fe_degree_0),
            false);

  test<dim>(FE_Nedelec<dim>(fe_degree_0), FE_Nedelec<dim>(fe_degree_1), false);
  test<dim>(FE_Nedelec<dim>(fe_degree_1), FE_Nedelec<dim>(fe_degree_0), false);
}


int
main()
{
  initlog();

  for (unsigned int i = 1; i <= 3; ++i)
    for (unsigned int j = i + i; j <= 3; ++j)
      test<2>(i, j);
}

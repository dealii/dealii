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

// A test that checks make_flux_sparsity_pattern
// with vector-valued hp-objects and coupling masks

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <iostream>
#include <vector>

#include "../tests.h"



template <int dim>
void
check()
{
  Triangulation<dim>        triangulation;
  std::vector<unsigned int> subdivisions(dim, 1U);
  subdivisions[0] = 2;
  Point<dim> p1, p2;
  switch (dim)
    {
      case 2:
        p1[0] = p1[1] = 0.0;
        p2[0]         = 2.0;
        p2[1]         = 1.0;
        break;
      case 3:
        p1[0] = p1[1] = p1[2] = 0.0;
        p2[0]                 = 2.0;
        p2[1] = p2[2] = 1.0;
        break;
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }

  GridGenerator::subdivided_hyper_rectangle(triangulation,
                                            subdivisions,
                                            p1,
                                            p2);

  // Create FE Collection and insert two FE objects
  // (RT0 - DGQ0 - Q1 ^ dim) and (Nothing - Nothing - Q2 ^ dim)
  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FESystem<dim>(
    FE_Nothing<dim>(dim), 1, FE_Nothing<dim>(), 1, FE_Q<dim>(2), dim));
  fe_collection.push_back(FESystem<dim>(FE_RaviartThomasNodal<dim>(0),
                                        1,
                                        FE_DGQ<dim>(0),
                                        1,
                                        FE_Nothing<dim>(),
                                        dim));

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.begin_active()->set_active_fe_index(1);
  dof_handler.distribute_dofs(fe_collection);

  // Create sparsity pattern allowing cell couplings of
  // RT-RT, RT-DGQ, DGQ-RT and Q-Q spaces,
  // and face coupling of
  // RT and Q
  DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());

  Table<2, DoFTools::Coupling> cell_coupling(fe_collection.n_components(),
                                             fe_collection.n_components());
  Table<2, DoFTools::Coupling> face_coupling(fe_collection.n_components(),
                                             fe_collection.n_components());

  for (unsigned int c = 0; c < fe_collection.n_components(); ++c)
    for (unsigned int d = 0; d < fe_collection.n_components(); ++d)
      {
        // RT-RT, RT-DGQ and DGQ-RT
        if ((c < dim + 1 && d < dim + 1) && !((c == dim) && (d == dim)))
          cell_coupling[c][d] = DoFTools::always;

        // Q-Q
        if (c >= dim + 1 && d >= dim + 1)
          cell_coupling[c][d] = DoFTools::always;

        // RT-Q
        if (c < dim && d >= dim + 1)
          face_coupling[c][d] = DoFTools::always;
      }

  DoFTools::make_flux_sparsity_pattern(dof_handler,
                                       dsp,
                                       cell_coupling,
                                       face_coupling);
  dsp.compress();

  // Print sparsity pattern
  dsp.print(deallog.get_file_stream());
}



int
main()
{
  initlog();

  deallog.push("2d");
  check<2>();
  deallog.pop();
  deallog.push("3d");
  check<3>();
  deallog.pop();
}

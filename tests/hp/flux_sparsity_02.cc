// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// A second test that checks make_flux_sparsity_pattern
// with scalar valued hp-objects.

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>

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
  // DGQ0 and DGQ1
  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FE_DGQ<dim>(0));
  fe_collection.push_back(FE_DGQ<dim>(1));

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.begin_active()->set_active_fe_index(1);
  dof_handler.distribute_dofs(fe_collection);

  deallog << "Dimension " << dim << std::endl;

  {
    deallog << "cell and face coupling" << std::endl;
    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());

    Table<2, DoFTools::Coupling> cell_coupling(1, 1);
    Table<2, DoFTools::Coupling> face_coupling(1, 1);
    cell_coupling(0, 0) = DoFTools::always;
    face_coupling(0, 0) = DoFTools::always;

    DoFTools::make_flux_sparsity_pattern(dof_handler,
                                         dsp,
                                         cell_coupling,
                                         face_coupling);
    dsp.compress();

    // Print sparsity pattern, we expect that all dofs couple with each other.
    dsp.print(deallog.get_file_stream());
  }

  {
    deallog << "cell coupling" << std::endl;
    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());

    Table<2, DoFTools::Coupling> cell_coupling(1, 1);
    Table<2, DoFTools::Coupling> face_coupling(1, 1);
    cell_coupling(0, 0) = DoFTools::always;
    face_coupling(0, 0) = DoFTools::none;

    DoFTools::make_flux_sparsity_pattern(dof_handler,
                                         dsp,
                                         cell_coupling,
                                         face_coupling);
    dsp.compress();

    // Print sparsity pattern, we expect no couplings across the face.
    dsp.print(deallog.get_file_stream());
  }

  {
    deallog << "face coupling" << std::endl;
    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());

    Table<2, DoFTools::Coupling> cell_coupling(1, 1);
    Table<2, DoFTools::Coupling> face_coupling(1, 1);
    cell_coupling(0, 0) = DoFTools::none;
    face_coupling(0, 0) = DoFTools::always;

    DoFTools::make_flux_sparsity_pattern(dof_handler,
                                         dsp,
                                         cell_coupling,
                                         face_coupling);
    dsp.compress();

    // Print sparsity pattern, we expect that all dofs couple with each other.
    dsp.print(deallog.get_file_stream());
  }

  {
    deallog << "no coupling" << std::endl;
    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());

    Table<2, DoFTools::Coupling> cell_coupling(1, 1);
    Table<2, DoFTools::Coupling> face_coupling(1, 1);
    cell_coupling(0, 0) = DoFTools::none;
    face_coupling(0, 0) = DoFTools::none;

    DoFTools::make_flux_sparsity_pattern(dof_handler,
                                         dsp,
                                         cell_coupling,
                                         face_coupling);
    dsp.compress();

    // Print sparsity pattern, we expect no couplings.
    dsp.print(deallog.get_file_stream());
  }
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

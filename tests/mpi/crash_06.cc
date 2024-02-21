// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test different elements for hanging_node constraints. This extends
// crash_05.cc to make sure stuff is working.

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_dg0.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <numeric>

#include "../tests.h"


template <int dim>
void
test(FiniteElement<dim> &fe)
{
  deallog << "dim=" << dim << std::endl;

  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_ball(triangulation);
  triangulation.refine_global(1);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      typename Triangulation<dim>::active_cell_iterator cell =
        triangulation.begin_active();
      for (; cell != triangulation.end(); ++cell)
        if (Testing::rand() % 2)
          cell->set_refine_flag();
    }

  triangulation.execute_coarsening_and_refinement();
  deallog << "n_cells: " << triangulation.n_global_active_cells() << std::endl;
  DoFHandler<dim> dof_handler(triangulation);

  dof_handler.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}

template <int dim>
void
testit()
{
  std::vector<std::shared_ptr<FiniteElement<dim>>> fes;
  fes.push_back(
    std::shared_ptr<FiniteElement<dim>>(new FE_RaviartThomas<dim>(0)));
  fes.push_back(
    std::shared_ptr<FiniteElement<dim>>(new FE_RaviartThomas<dim>(1)));
  fes.push_back(std::shared_ptr<FiniteElement<dim>>(new FE_Nedelec<dim>(0)));
  fes.push_back(std::shared_ptr<FiniteElement<dim>>(new FE_Nedelec<dim>(1)));
  fes.push_back(std::shared_ptr<FiniteElement<dim>>(new FE_Q<dim>(3)));
  fes.push_back(std::shared_ptr<FiniteElement<dim>>(new FE_DGQ<dim>(2)));
  fes.push_back(std::shared_ptr<FiniteElement<dim>>(new FE_Q_DG0<dim>(2)));

  for (unsigned int i = 0; i < fes.size(); ++i)
    {
      deallog << fes[i]->get_name() << std::endl;
      test<dim>(*fes[i]);
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  testit<2>();
  testit<3>();
}

// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2022 by the deal.II authors
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


/**
 * Test prolongation and restriction with "non-nested" triangulations where one
 * is obtained from the other by a "random" distortion.
 *
 */

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include "mg_transfer_util.h"

using namespace dealii;

template <int dim, typename Number>
void
do_test(const FiniteElement<dim> &   fe_fine,
        const FiniteElement<dim> &   fe_coarse,
        const Function<dim, Number> &function)
{
  // create coarse grid
  parallel::distributed::Triangulation<dim> tria_coarse(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria_coarse, -1., 1.);
  tria_coarse.refine_global(3);
  const double factor = 0.25;
  GridTools::distort_random(factor,
                            tria_coarse,
                            true,
                            boost::random::mt19937::default_seed);

  // create fine grid
  parallel::distributed::Triangulation<dim> tria_fine(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria_fine, -1., 1.);
  tria_fine.refine_global(4);
  GridTools::distort_random(factor,
                            tria_fine,
                            true,
                            boost::random::mt19937::default_seed);


  {
    // Just print the grids.
    std::ofstream coarse_filename("initial_tria.vtk");
    std::ofstream fine_filename("distorted_tria.vtk");
    GridOut       go;
    go.write_vtk(tria_coarse, coarse_filename);
    go.write_vtk(tria_fine, fine_filename);
  }

  // setup dof-handlers
  DoFHandler<dim> dof_handler_fine(tria_fine);
  dof_handler_fine.distribute_dofs(fe_fine);

  DoFHandler<dim> dof_handler_coarse(tria_coarse);
  dof_handler_coarse.distribute_dofs(fe_coarse);

  // setup constraint matrix
  AffineConstraints<Number> constraint_coarse;
  constraint_coarse.close();

  AffineConstraints<Number> constraint_fine;
  constraint_fine.close();

  // setup transfer operator
  MGTwoLevelTransferNonNested<dim, LinearAlgebra::distributed::Vector<Number>>
                 transfer;
  MappingQ1<dim> mapping_fine, mapping_coarse;
  transfer.reinit(dof_handler_fine,
                  dof_handler_coarse,
                  mapping_fine,
                  mapping_coarse,
                  constraint_fine,
                  constraint_coarse);

  test_non_nested_transfer(transfer,
                           dof_handler_fine,
                           dof_handler_coarse,
                           function);
}

template <int dim, typename Number>
void
test(int fe_degree, const Function<dim, Number> &function)
{
  const auto str_fine   = std::to_string(fe_degree);
  const auto str_coarse = std::to_string(fe_degree);

  if (fe_degree > 0)
    {
      deallog.push("CG<2>(" + str_fine + ")<->CG<2>(" + str_coarse + ")");
      do_test<dim, double>(FE_Q<dim>(fe_degree),
                           FE_Q<dim>(fe_degree),
                           function);
      deallog.pop();
    }


  if (fe_degree > 0)
    {
      deallog.push("DG<2>(" + str_fine + ")<->CG<2>(" + str_coarse + ")");
      do_test<dim, double>(FE_DGQ<dim>(fe_degree),
                           FE_Q<dim>(fe_degree),
                           function);
      deallog.pop();
    }


  {
    deallog.push("DG<2>(" + str_fine + ")<->DG<2>(" + str_coarse + ")");
    do_test<dim, double>(FE_DGQ<dim>(fe_degree),
                         FE_DGQ<dim>(fe_degree),
                         function);
    deallog.pop();
  }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.precision(8);
  Functions::Monomial<2> linear_function(Tensor<1, 2>({1, 0})); // f(x,y)= x
  for (unsigned int i = 0; i < 5; ++i)
    test<2, double>(i, linear_function);
}

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test DataOut::add_mg_data_vector in parallel for a vector valued dof function

#include <deal.II/base/function_lib.h>
#include <deal.II/base/function_parser.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
void
do_test()
{
  parallel::distributed::Triangulation<dim> triangulation(
    MPI_COMM_WORLD,
    dealii::Triangulation<dim>::none,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);
  triangulation.begin_active()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();

  FESystem<dim> fe(FE_DGQ<dim>(1), dim + 1);

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  // Make FE vector
  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  using VectorType = LinearAlgebra::distributed::Vector<double>;
  VectorType global_dof_vector;
  global_dof_vector.reinit(dof_handler.locally_owned_dofs(),
                           locally_relevant_dofs,
                           MPI_COMM_WORLD);

  QGauss<dim> quadrature(3);

  FunctionParser<dim> function(dim == 2 ? "y;-x;x*x+y*y" :
                                          "y;-x;0;x*x+y*y+z*z");

  AffineConstraints<double> constraints(dof_handler.locally_owned_dofs(),
                                        locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();
  VectorTools::project(
    dof_handler, constraints, quadrature, function, global_dof_vector);
  /*  VectorTools::interpolate(dof_handler,
                           Functions::SquareFunction<dim>(),
                           global_dof_vector);
  */

  global_dof_vector.update_ghost_values();

  std::vector<std::string> names(dim, "vec");
  names.emplace_back("my_scalar");
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation(dim,
                   DataComponentInterpretation::component_is_part_of_vector);
  interpretation.emplace_back(DataComponentInterpretation::component_is_scalar);

  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(dof_handler,
                             global_dof_vector,
                             names,
                             interpretation);
    data_out.build_patches(0);
    data_out.write_gnuplot(deallog.get_file_stream());
    data_out.write_vtu_in_parallel("base.vtu", MPI_COMM_WORLD);
  }

  {
    MGLevelObject<VectorType>         dof_vector(0,
                                         triangulation.n_global_levels() - 1);
    MGTransferMatrixFree<dim, double> transfer;

    transfer.build(dof_handler);
    transfer.interpolate_to_mg(dof_handler, dof_vector, global_dof_vector);

    for (unsigned int level = 0; level < triangulation.n_global_levels();
         ++level)
      {
        deallog << "* level " << level << std::endl;
        DataOut<dim> data_out;
        data_out.set_cell_selection(
          [level](const typename Triangulation<dim>::cell_iterator &cell) {
            return (cell->level() == static_cast<int>(level) &&
                    cell->is_locally_owned_on_level());
          });
        data_out.attach_triangulation(triangulation);
        data_out.add_mg_data_vector(dof_handler, dof_vector, "data");
        data_out.add_mg_data_vector(dof_handler,
                                    dof_vector,
                                    names,
                                    interpretation);
        data_out.build_patches(0);
        data_out.write_gnuplot(deallog.get_file_stream());
        data_out.write_vtu_in_parallel(std::string("level") +
                                         Utilities::to_string(level) + ".vtu",
                                       MPI_COMM_WORLD);
      }
  }
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    log;

  do_test<2>();
  //  do_test<3>();
  return 0;
}

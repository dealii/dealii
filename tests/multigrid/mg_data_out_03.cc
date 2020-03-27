// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// Test DataOut::add_mg_data_vector for vector data in serial

#include <deal.II/base/function_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/timer.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

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
  FE_Q<dim>          fe(1);
  Triangulation<dim> triangulation(
    Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);
  triangulation.begin_active()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  Vector<double> global_dof_vector(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler,
                           Functions::SquareFunction<dim>(),
                           global_dof_vector);

  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(dof_handler,
                             global_dof_vector,
                             std::vector<std::string>(1, "data"));
    data_out.build_patches(0);
    data_out.write_gnuplot(deallog.get_file_stream());
    std::ofstream st("base.vtu");
    data_out.write_vtu(st);
  }

  {
    MGLevelObject<Vector<double>> dof_vector(0,
                                             triangulation.n_global_levels() -
                                               1);

    MGTransferPrebuilt<Vector<double>> transfer;

    transfer.build_matrices(dof_handler);

    // This is not quite the correct thing to do, but "interpolate_to_mg"
    // doesn't exist for serial vectors. Here it means that the level 0 is 0
    // instead of having correct values.
    transfer.copy_to_mg(dof_handler, dof_vector, global_dof_vector);

    for (unsigned int level = 0; level < triangulation.n_global_levels();
         ++level)
      {
        deallog << "* level " << level << std::endl;

        DataOut<dim> data_out;
        data_out.set_cell_selection(
          [level](const typename Triangulation<dim>::cell_iterator &cell) {
            return cell->level() == static_cast<int>(level);
          });

        data_out.attach_triangulation(triangulation);
        data_out.add_mg_data_vector(dof_handler, dof_vector, "data");
        data_out.build_patches(0);
        data_out.write_gnuplot(deallog.get_file_stream());
        std::ofstream st(std::string("level") + Utilities::to_string(level) +
                         ".vtu");
        data_out.write_vtu(st);
      }
  }
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  initlog();

  do_test<2>();
  /*
  do_test<3>();*/
  return 0;
}

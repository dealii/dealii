// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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



// check the initialization of MatrixFree with face information on a
// continuous element with hanging nodes


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    log;

  const int dim = 3;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      tria.begin_active()->set_refine_flag();
    }
  tria.execute_coarsening_and_refinement();

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FE_Q<dim>(3));

  MatrixFree<dim>                          matrix_free;
  typename MatrixFree<dim>::AdditionalData data;

  data.mapping_update_flags                = update_values;
  data.mapping_update_flags_inner_faces    = update_values; // problem
  data.mapping_update_flags_boundary_faces = update_values; // problem

  AffineConstraints<double> constraint;
  constraint.clear();
  DoFTools::make_hanging_node_constraints(dof_handler, constraint);
  constraint.close();

  matrix_free.reinit(dof_handler, constraint, QGauss<1>(4), data);

  LinearAlgebra::distributed::Vector<double> src;
  matrix_free.initialize_dof_vector(src);

  deallog << "main partitioner: size="
          << matrix_free.get_dof_info().vector_partitioner->size()
          << " local_size="
          << matrix_free.get_dof_info().vector_partitioner->local_size()
          << " n_ghosts="
          << matrix_free.get_dof_info().vector_partitioner->n_ghost_indices()
          << std::endl;
  for (auto &p : matrix_free.get_dof_info().vector_partitioner_face_variants)
    deallog << "partitioner: size=" << p->size()
            << " local_size=" << p->local_size()
            << " n_ghosts=" << p->n_ghost_indices() << std::endl;
}

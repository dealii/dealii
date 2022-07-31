// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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


// Unit test to check whether VectorTools::project() works for DoFHandler
// objects in hp-mode.

#include <deal.II/base/function_lib.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/vector_tools_project.h>

#include "../tests.h"


template <int dim, int spacedim = dim>
void
test()
{
  MPI_Comm mpi_communicator = MPI_COMM_WORLD;

  parallel::distributed::Triangulation<dim, spacedim> tria(mpi_communicator);
  GridGenerator::hyper_cube(tria);
  tria.refine_global((dim == 2) ? 4 : 2);


  hp::FECollection<dim, spacedim> fes;
  for (unsigned int degree = 1; degree <= 2; ++degree)
    fes.push_back(FE_Q<dim, spacedim>(degree));

  DoFHandler<dim, spacedim> dofh(tria);
  dofh.begin_active()->set_active_fe_index(1);
  dofh.distribute_dofs(fes);


  const IndexSet locally_owned_dofs = dofh.locally_owned_dofs();
  IndexSet       locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dofh, locally_relevant_dofs);


  AffineConstraints<double> constraints;
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dofh, constraints);
  constraints.close();


  LinearAlgebra::distributed::Vector<double> projection(locally_owned_dofs,
                                                        locally_relevant_dofs,
                                                        mpi_communicator);

  projection.zero_out_ghost_values();

  VectorTools::project(dofh,
                       constraints,
                       QGauss<dim>(fes.max_degree() + 2),
                       Functions::SquareFunction<dim>(),
                       projection);

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init_finalize(argc, argv, 1);
  mpi_initlog();
  
  test<2>();
  test<3>();
}

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// tests that the matrix-free implementation works correctly with periodic
// boundary conditions by solving a simple linear system and outputting the
// solution for FE_DGQ(0)

#include <deal.II/base/function.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include "../tests.h"



// We want to use the matrix-vector product provided by this function (which
// also includes a main function)
#include "matrix_vector_faces_common.h"


template <int dim, int fe_degree>
void
test()
{
  if (fe_degree > 1)
    return;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria, 0, 1);

  // set boundary ids on boundaries to the number of the face
  for (unsigned int face = 2; face < GeometryInfo<dim>::faces_per_cell; ++face)
    tria.begin()->face(face)->set_all_boundary_ids(face);

  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodic_faces;
  for (unsigned int d = 1; d < dim; ++d)
    GridTools::collect_periodic_faces(
      tria, 2 * d, 2 * d + 1, d, periodic_faces);
  deallog << "Periodic faces: " << periodic_faces.size() << std::endl;
  tria.add_periodicity(periodic_faces);

  tria.refine_global(8 - 2 * dim);

  FE_DGQ<dim>     fe(0);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  MappingQ<dim> mapping(fe_degree + 1);

  const QGauss<1>                          quad(1);
  typename MatrixFree<dim>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim>::AdditionalData::none;
  data.mapping_update_flags_inner_faces =
    (update_gradients | update_JxW_values);
  data.mapping_update_flags_boundary_faces =
    (update_gradients | update_JxW_values);

  MatrixFree<dim> mf_data;
  mf_data.reinit(MappingQ1<dim>{}, dof, constraints, quad, data);

  LinearAlgebra::distributed::Vector<double> rhs, sol;
  mf_data.initialize_dof_vector(rhs);
  rhs = 1.;
  mf_data.initialize_dof_vector(sol);

  MatrixFreeTest<dim, 0, 1, double, LinearAlgebra::distributed::Vector<double>>
    mf(mf_data);

  SolverControl control(1000, 1e-12 * std::sqrt(rhs.size()));
  SolverCG<LinearAlgebra::distributed::Vector<double>> solver(control);
  solver.solve(mf, sol, rhs, PreconditionIdentity());

  const std::vector<IndexSet> locally_owned_dofs_per_processor =
    Utilities::MPI::all_gather(MPI_COMM_WORLD, dof.locally_owned_dofs());
  // gather all data at root
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      Vector<double> solution_gather0(sol.size());
      double        *sol_gather_ptr = solution_gather0.begin();
      for (unsigned int i = 0; i < sol.locally_owned_size(); ++i)
        *sol_gather_ptr++ = sol.local_element(i);
      for (unsigned int i = 1;
           i < Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
           ++i)
        {
          MPI_Recv(sol_gather_ptr,
                   locally_owned_dofs_per_processor[i].n_elements(),
                   MPI_DOUBLE,
                   i,
                   i,
                   MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
          sol_gather_ptr += locally_owned_dofs_per_processor[i].n_elements();
        }
      solution_gather0.print(deallog.get_file_stream());
    }
  else
    MPI_Send(sol.begin(),
             sol.locally_owned_size(),
             MPI_DOUBLE,
             0,
             Utilities::MPI::this_mpi_process(MPI_COMM_WORLD),
             MPI_COMM_WORLD);

  deallog << "OK" << std::endl;
}

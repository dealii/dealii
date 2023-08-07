// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2022 by the deal.II authors
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


// test MGTransferMatrixFree<dim,Number>::interpolate_to_mg()
// for a scalar field

#include <deal.II/base/function.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/identity_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <vector>

#include "../tests.h"



template <int dim>
class SimpleField : public Function<dim>
{
public:
  SimpleField()
    : Function<dim>(1)
  {}

  double
  value(const Point<dim> &p, const unsigned int component = 0) const
  {
    (void)component;
    return p[0] * 2. + p[1] - 10.;
  }
};


template <int dim,
          int fe_degree            = 2,
          int n_q_points           = fe_degree + 1,
          typename NumberType      = double,
          typename LevelNumberType = NumberType>
void
test(const unsigned int n_glob_ref = 2, const unsigned int n_ref = 0)
{
  SimpleField<dim> function;

  deallog << "dim=" << dim << std::endl;
  MPI_Comm     mpi_communicator(MPI_COMM_WORLD);
  unsigned int myid    = Utilities::MPI::this_mpi_process(mpi_communicator);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(mpi_communicator);

  deallog << "numproc=" << numproc << std::endl;

  parallel::distributed::Triangulation<dim> triangulation(
    mpi_communicator,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::subdivided_hyper_cube(triangulation, 5, -1, 1);
  // now do some refinement
  triangulation.refine_global(n_glob_ref);
  for (unsigned int ref = 0; ref < n_ref; ++ref)
    {
      for (typename Triangulation<dim>::active_cell_iterator cell =
             triangulation.begin_active();
           cell != triangulation.end();
           ++cell)
        if (cell->is_locally_owned() &&
            ((cell->center().norm() < 0.5 &&
              (cell->level() < 5 || cell->center().norm() > 0.45)) ||
             (dim == 2 && cell->center().norm() > 1.2)))
          cell->set_refine_flag();
      triangulation.execute_coarsening_and_refinement();
    }

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof_handler(triangulation);

  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  IndexSet locally_owned_dofs, locally_relevant_dofs;
  locally_owned_dofs = dof_handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  // constraints:
  AffineConstraints<double> constraints;
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();


  // interpolate:
  LinearAlgebra::distributed::Vector<LevelNumberType> fine_projection;
  fine_projection.reinit(locally_owned_dofs,
                         locally_relevant_dofs,
                         mpi_communicator);
  VectorTools::project(dof_handler,
                       constraints,
                       QGauss<dim>(n_q_points + 2),
                       function,
                       fine_projection);
  fine_projection.update_ghost_values();

  // output for debug purposes:
  if (false)
    {
      DataOut<dim> data_out;
      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(fine_projection, "projection");

      Vector<float> subdomain(triangulation.n_active_cells());
      for (unsigned int i = 0; i < subdomain.size(); ++i)
        subdomain(i) = triangulation.locally_owned_subdomain();
      data_out.add_data_vector(subdomain, "subdomain");
      data_out.build_patches();
      const std::string filename =
        "output_" + Utilities::int_to_string(myid) + ".vtu";
      std::ofstream output(filename);
      data_out.write_vtu(output);

      const std::string mg_mesh = "mg_mesh";
      GridOut           grid_out;
      grid_out.write_mesh_per_processor_as_vtu(triangulation,
                                               mg_mesh,
                                               true,
                                               true);
    }

  // MG hanging node constraints
  // we do not add extra zero Dirichlet BC here
  MGConstrainedDoFs mg_constrained_dofs;
  mg_constrained_dofs.initialize(dof_handler);

  // MG transfer:
  MGTransferMF<dim, LevelNumberType> mg_transfer(mg_constrained_dofs);
  mg_transfer.build(dof_handler);

  // now the core of the test:
  const unsigned int max_level =
    dof_handler.get_triangulation().n_global_levels() - 1;
  const unsigned int min_level = 0;
  MGLevelObject<LinearAlgebra::distributed::Vector<LevelNumberType>>
    level_projection(min_level, max_level);
  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      IndexSet set;
      DoFTools::extract_locally_relevant_level_dofs(dof_handler, level, set);
      level_projection[level].reinit(dof_handler.locally_owned_mg_dofs(level),
                                     set,
                                     mpi_communicator);
    }
  mg_transfer.interpolate_to_mg(dof_handler, level_projection, fine_projection);

  // now go through all GMG levels and make sure FE field can represent
  // analytic function exactly:
  QGauss<dim>                  quadrature(n_q_points);
  std::vector<LevelNumberType> q_values(quadrature.size());

  FEValues<dim> fe_values(fe,
                          quadrature,
                          update_values | update_quadrature_points);
  for (unsigned int level = max_level + 1; level != min_level;)
    {
      --level;

      level_projection[level].update_ghost_values();

      std::vector<types::global_dof_index>    dof_indices(fe.dofs_per_cell);
      typename DoFHandler<dim>::cell_iterator cell = dof_handler.begin(level);
      typename DoFHandler<dim>::cell_iterator endc = dof_handler.end(level);
      for (; cell != endc; ++cell)
        if (cell->is_locally_owned_on_level())
          {
            fe_values.reinit(cell);
            cell->get_mg_dof_indices(dof_indices);
            fe_values.get_function_values(level_projection[level],
                                          dof_indices,
                                          q_values);

            const std::vector<Point<dim>> &q_points =
              fe_values.get_quadrature_points();

            for (unsigned int q = 0; q < q_points.size(); ++q)
              {
                const double diff = q_values[q] - function.value(q_points[q]);
                if (std::abs(diff) > 1e-10)
                  {
                    std::cout << "dofs: ";
                    for (const auto i : dof_indices)
                      std::cout << i << ' ';
                    std::cout << std::endl << "values: ";
                    std::vector<LevelNumberType> local_values(
                      dof_indices.size());
                    level_projection[level].extract_subvector_to(dof_indices,
                                                                 local_values);
                    for (const auto v : local_values)
                      std::cout << v << ' ';
                    std::cout << std::endl
                              << "val(q)=" << q_values[q] << std::endl;
                    std::cout << "MGTransfer indices:" << std::endl;
                    mg_transfer.print_indices(std::cout);
                    AssertThrow(false,
                                ExcMessage("Level " + std::to_string(level) +
                                           " Diff " + std::to_string(diff) +
                                           " center_x " +
                                           std::to_string(cell->center()[0]) +
                                           " center_y " +
                                           std::to_string(cell->center()[1])));
                  }
              }
          }
    }

  deallog << "Ok" << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test<2>();
  test<2, 1>(0, 1);
}

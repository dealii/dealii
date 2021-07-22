// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


// Test ConsensusAlgorithms::AnonymousProcess.

#include <deal.II/base/mpi_compute_index_owner_internal.h>
#include <deal.II/base/mpi_consensus_algorithms.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_repartitioner.h>
#include <deal.II/numerics/vector_tools.h>

#include <set>

#include "../tests.h"

using namespace dealii;

template <int dim>
class Solution : public Function<dim>
{
public:
  double
  value(const Point<dim> &p, const unsigned int = 0) const
  {
    return p[0];
  }
};

template <int dim, int spacedim>
void
test(const MPI_Comm &comm)
{
  Triangulation<dim> basetria;
  GridGenerator::subdivided_hyper_cube(basetria, 4);
  basetria.refine_global(2);

  for (unsigned int i = 0; i < 2; ++i)
    {
      for (auto cell : basetria.active_cell_iterators())
        if (cell->center()[0] < 0.5 && cell->center()[1] < 0.5)
          cell->set_refine_flag();
      basetria.execute_coarsening_and_refinement();
    }

  parallel::fullydistributed::Triangulation<dim> tria_1(comm);
  parallel::fullydistributed::Triangulation<dim> tria_2(comm);

  {
    GridTools::partition_triangulation_zorder(
      Utilities::MPI::n_mpi_processes(comm), basetria);

    auto construction_data = TriangulationDescription::Utilities::
      create_description_from_triangulation(basetria, comm);

    tria_1.create_triangulation(construction_data);
  }

  {
    GridTools::partition_triangulation_zorder(
      std::min(Utilities::MPI::n_mpi_processes(comm), 3u), basetria);

    auto construction_data = TriangulationDescription::Utilities::
      create_description_from_triangulation(basetria, comm);

    tria_2.create_triangulation(construction_data);
  }

  FE_Q<dim, spacedim> fe(1);

  DoFHandler<dim, spacedim> dof_handler_1(tria_1);
  dof_handler_1.distribute_dofs(fe);

  DoFHandler<dim, spacedim> dof_handler_2(tria_2);
  dof_handler_2.distribute_dofs(fe);


  // setup first vector
  IndexSet d1 = dof_handler_1.locally_owned_dofs();
  IndexSet dd1;
  DoFTools::extract_locally_relevant_dofs(dof_handler_1, dd1);

  LinearAlgebra::distributed::Vector<double> vec1(d1, dd1, comm);

  VectorTools::interpolate(dof_handler_1, Solution<spacedim>(), vec1);

  // setup second vector
  IndexSet d2 = dof_handler_2.locally_owned_dofs();
  IndexSet dd2;
  DoFTools::extract_locally_relevant_dofs(dof_handler_2, dd2);

  LinearAlgebra::distributed::Vector<double> vec2(d2, dd2, comm);

  VectorRepartitioner<LinearAlgebra::distributed::Vector<double>> vr;

  vr.reinit(dof_handler_2, dof_handler_1);

  vr.update_forwards(vec2, vec1);

  deallog << vec1.l2_norm() << " ";
  deallog << vec2.l2_norm() << " ";

  AssertDimension(vec1.l2_norm(), vec2.l2_norm());

  const auto write_vtu =
    [&comm](const auto &dof_handler_1, const auto &vec1, const auto file_name) {
      return;

      DataOut<dim> data_out;

      data_out.attach_dof_handler(dof_handler_1);

      vec1.update_ghost_values();
      data_out.add_data_vector(vec1, "result");

      Vector<double> ranks(dof_handler_1.get_triangulation().n_active_cells());
      for (const auto &cell : dof_handler_1.active_cell_iterators())
        if (cell->is_locally_owned())
          ranks(cell->active_cell_index()) = cell->subdomain_id();
      data_out.add_data_vector(ranks, "ranks");

      data_out.build_patches(1);
      data_out.write_vtu_in_parallel(file_name, comm);
    };

  write_vtu(dof_handler_1, vec1, "result.0.vtu");
  write_vtu(dof_handler_2, vec2, "result.1.vtu");


  vec1 = 0.0;
  vr.update_backwards(vec1, vec2);
  deallog << vec1.l2_norm() << std::endl;

  AssertDimension(vec1.l2_norm(), vec2.l2_norm());


  write_vtu(dof_handler_1, vec1, "result.2.vtu");
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  const MPI_Comm comm = MPI_COMM_WORLD;

  test<2, 2>(comm);
}

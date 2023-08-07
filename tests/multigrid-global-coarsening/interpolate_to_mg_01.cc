// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2022 by the deal.II authors
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


// test function MGTransferMatrixFree::interpolate_to_mg() for periodic
// boundaries

#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim>
class RightHandSideFunction : public Function<dim>
{
public:
  RightHandSideFunction(const unsigned int n_components)
    : Function<dim>(n_components)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const
  {
    return p[component % dim] * p[component % dim];
  }
};

template <int dim>
void
test()
{
  const unsigned int fe_degree = 1;

  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);

  GridGenerator::hyper_cube(tria, -1., 1., true);

  std::vector<GridTools::PeriodicFacePair<TriaIterator<CellAccessor<dim, dim>>>>
    tria_matched_pairs;
  GridTools::collect_periodic_faces(tria, 0, 1, 0, tria_matched_pairs);
  tria.add_periodicity(tria_matched_pairs);
  tria.refine_global(2);

  const FE_Q<dim> fe(fe_degree);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  AffineConstraints<double> constraints;
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  std::vector<
    GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator>>
    dof_matched_pairs;
  GridTools::collect_periodic_faces(dof_handler, 0, 1, 0, dof_matched_pairs);
  DoFTools::make_periodicity_constraints<dim, dim>(dof_matched_pairs,
                                                   constraints);
  constraints.close();

  LinearAlgebra::distributed::Vector<double> qsol(
    dof_handler.locally_owned_dofs(), locally_relevant_dofs, MPI_COMM_WORLD);

  VectorTools::interpolate(MappingQ1<dim>(),
                           dof_handler,
                           RightHandSideFunction<dim>(1),
                           qsol);

  MGLevelObject<LinearAlgebra::distributed::Vector<double>> mg_qsol;
  MGConstrainedDoFs                                         mg_constrained_dofs;
  MGTransferMF<dim, double>                                 mg_transfer;

  unsigned int n_tria_levels = tria.n_global_levels();
  mg_qsol.resize(0, n_tria_levels - 1);

  mg_constrained_dofs.initialize(dof_handler);
  mg_transfer.initialize_constraints(mg_constrained_dofs);
  mg_transfer.build(dof_handler);

  mg_transfer.interpolate_to_mg(dof_handler, mg_qsol, qsol);

  for (unsigned int i = mg_qsol.min_level(); i <= mg_qsol.max_level(); ++i)
    deallog << mg_qsol[i].l2_norm() << std::endl;

  deallog << "OK" << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test<2>();
}

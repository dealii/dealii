// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Evaluate solution vector along a line.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_evaluate.h>

#include "../tests.h"


template <typename MeshType>
MPI_Comm
get_mpi_comm(const MeshType &mesh)
{
  const auto *tria_parallel = dynamic_cast<
    const parallel::TriangulationBase<MeshType::dimension,
                                      MeshType::space_dimension> *>(
    &(mesh.get_triangulation()));

  return tria_parallel != nullptr ? tria_parallel->get_mpi_communicator() :
                                    MPI_COMM_SELF;
}

template <int dim, int spacedim>
std::shared_ptr<const Utilities::MPI::Partitioner>
create_partitioner(const DoFHandler<dim, spacedim> &dof_handler)
{
  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  return std::make_shared<const Utilities::MPI::Partitioner>(
    dof_handler.locally_owned_dofs(),
    locally_relevant_dofs,
    get_mpi_comm(dof_handler));
}

template <int dim>
void
print(const Mapping<dim>                               &mapping,
      const DoFHandler<dim>                            &dof_handler,
      const LinearAlgebra::distributed::Vector<double> &result,
      const unsigned int                                counter)
{
  result.update_ghost_values();

  DataOutBase::VtkFlags flags;
  flags.write_higher_order_cells = true;

  DataOut<dim> data_out;
  data_out.set_flags(flags);
  data_out.attach_dof_handler(dof_handler);

  const auto &tria = dof_handler.get_triangulation();

  Vector<double> ranks(tria.n_active_cells());
  for (const auto &cell :
       tria.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
    ranks(cell->active_cell_index()) = cell->subdomain_id();
  data_out.add_data_vector(ranks, "rank");
  data_out.add_data_vector(result, "result");

  data_out.build_patches(mapping, dof_handler.get_fe().tensor_degree() + 1);
  data_out.write_vtu_with_pvtu_record(
    "./", "example-7", counter, MPI_COMM_WORLD, 1, 1);

  result.zero_out_ghost_values();
}

template <int dim>
class AnalyticalFunction : public Function<dim>
{
public:
  AnalyticalFunction(const unsigned int n_components)
    : Function<dim>(n_components)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component) const
  {
    return p[component] * p[component];
  }
};

void
test()
{
  const unsigned int dim             = 2;
  const unsigned int degree          = 2;
  const unsigned int n_components    = dim;
  const unsigned int n_refinements_1 = 3;

  parallel::distributed::Triangulation<dim> tria_1(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria_1, -1.0, +1.0);
  tria_1.refine_global(n_refinements_1);
  DoFHandler<dim>       dof_handler_1(tria_1);
  dealii::FESystem<dim> fe_1(dealii::FE_Q<dim>(degree), n_components);
  dof_handler_1.distribute_dofs(fe_1);

  LinearAlgebra::distributed::Vector<double> vector_1(
    create_partitioner(dof_handler_1));

  const MappingQ1<dim> mapping_1;
  VectorTools::interpolate(mapping_1,
                           dof_handler_1,
                           AnalyticalFunction<dim>(n_components),
                           vector_1);

  std::vector<Point<dim>> evaluation_points;

  const unsigned int n_intervals = 100;

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    for (unsigned int i = 0; i <= n_intervals; ++i)
      evaluation_points.emplace_back(-1.0 + 2.0 / n_intervals * i, 0.0);

  vector_1.update_ghost_values();
  Utilities::MPI::RemotePointEvaluation<dim> evaluation_cache;
  const auto                                 evaluation_point_gradient_results =
    VectorTools::point_gradients<n_components>(
      mapping_1, dof_handler_1, vector_1, evaluation_points, evaluation_cache);
  vector_1.zero_out_ghost_values();

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    for (unsigned int i = 0; i <= n_intervals; ++i)
      for (unsigned int c1 = 0; c1 < n_components; ++c1)
        for (unsigned int c2 = 0; c2 < n_components; ++c2)
          deallog << evaluation_point_gradient_results[i][c1][c2] << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  initlog();
  deallog << std::setprecision(8) << std::fixed;

  test();
}

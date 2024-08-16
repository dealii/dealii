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

// Interpolate result between non-matching grids.

#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
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
    "./", "example-3", counter, MPI_COMM_WORLD, 1, 1);

  result.zero_out_ghost_values();
}

void
test()
{
  const unsigned int dim    = 2;
  const unsigned int degree = 2;

  const auto create_triangulation = [](auto &tria) {
    GridGenerator::hyper_ball(tria);
    tria.refine_global(1);
    for (unsigned int i = 0; i < 3; ++i)
      {
        for (const auto &cell : tria.active_cell_iterators())
          if (cell->barycenter()[0] > 0.0)
            cell->set_refine_flag();
        tria.execute_coarsening_and_refinement();
      }
  };

  parallel::shared::Triangulation<dim> tria_1(
    MPI_COMM_WORLD,
    Triangulation<dim>::MeshSmoothing::none,
    false,
    parallel::shared::Triangulation<dim>::Settings::partition_auto);
  create_triangulation(tria_1);

  DoFHandler<dim> dof_handler_1(tria_1);
  FE_Q<dim>       fe_1(degree);
  dof_handler_1.distribute_dofs(fe_1);

  parallel::shared::Triangulation<dim> tria_2(
    MPI_COMM_WORLD,
    Triangulation<dim>::MeshSmoothing::none,
    false,
    parallel::shared::Triangulation<dim>::Settings::partition_custom_signal);
  tria_2.signals.post_refinement.connect([&tria_2]() {
    for (const auto &cell : tria_2.active_cell_iterators())
      {
        const auto point = cell->center();
        if (point[0] < 0.0 && point[1] < 0.0)
          cell->set_subdomain_id(0);
        else if (point[0] >= 0.0 && point[1] < 0.0)
          cell->set_subdomain_id(1);
        else if (point[0] < 0.0 && point[1] >= 0.0)
          cell->set_subdomain_id(2);
        else if (point[0] >= 0.0 && point[1] >= 0.0)
          cell->set_subdomain_id(3);
      }
  });


  create_triangulation(tria_2);

  DoFHandler<dim>           dof_handler_2(tria_2);
  FE_DGQArbitraryNodes<dim> fe_2(QGauss<1>(degree + 1).get_points());
  dof_handler_2.distribute_dofs(fe_2);

  LinearAlgebra::distributed::Vector<double> vector_1(
    create_partitioner(dof_handler_1));

  const MappingQ1<dim> mapping_1;
  VectorTools::interpolate(mapping_1,
                           dof_handler_1,
                           Functions::CosineFunction<2>(),
                           vector_1);

  const MappingQ1<dim> mapping_2;
  const QGauss<dim>    quadrature_2(degree + 1);
  FEValues<dim>        fe_values_2(mapping_2,
                            fe_2,
                            quadrature_2,
                            update_quadrature_points);

  std::vector<Point<dim>> evaluation_points;
  for (const auto &cell : dof_handler_2.active_cell_iterators())
    {
      if (cell->is_locally_owned() == false)
        continue;

      fe_values_2.reinit(cell);
      evaluation_points.insert(evaluation_points.end(),
                               fe_values_2.get_quadrature_points().begin(),
                               fe_values_2.get_quadrature_points().end());
    }

  vector_1.update_ghost_values();
  Utilities::MPI::RemotePointEvaluation<dim> evaluation_cache;
  const auto evaluation_point_results = VectorTools::point_values<1>(
    mapping_1, dof_handler_1, vector_1, evaluation_points, evaluation_cache);
  vector_1.zero_out_ghost_values();

  LinearAlgebra::distributed::Vector<double> vector_2(
    create_partitioner(dof_handler_2));

  auto ptr = evaluation_point_results.begin();
  for (const auto &cell : dof_handler_2.active_cell_iterators())
    {
      if (cell->is_locally_owned() == false)
        continue;

      fe_values_2.reinit(cell);

      Vector<double> temp(ptr, ptr + fe_values_2.n_quadrature_points);

      cell->set_dof_values(temp, vector_2);

      ptr += fe_values_2.n_quadrature_points;
    }

  if (false)
    {
      print(mapping_1, dof_handler_1, vector_1, 0);
      print(mapping_2, dof_handler_2, vector_2, 1);

      vector_1.print(deallog.get_file_stream());
      vector_2.print(deallog.get_file_stream());
    }

  Vector<double> difference(tria_2.n_active_cells());

  VectorTools::integrate_difference(mapping_2,
                                    dof_handler_2,
                                    vector_2,
                                    Functions::CosineFunction<2>(),
                                    difference,
                                    quadrature_2,
                                    VectorTools::L2_norm);

  double error =
    VectorTools::compute_global_error(tria_2, difference, VectorTools::L2_norm);


  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    deallog << error << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  initlog();

  test();
}

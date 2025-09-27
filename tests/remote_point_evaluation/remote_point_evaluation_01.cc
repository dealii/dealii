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

// Integrate surface tension on surface mesh and test the result on background
// mesh.

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_remote_point_evaluation.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/lac/la_parallel_block_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



using VectorType = LinearAlgebra::distributed::Vector<double>;

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

namespace dealii
{
  namespace VectorTools
  {
    template <int dim, int spacedim, typename VectorType>
    void
    get_position_vector(const DoFHandler<dim, spacedim> &dof_handler_dim,
                        VectorType                   &euler_coordinates_vector,
                        const Mapping<dim, spacedim> &mapping)
    {
      FEValues<dim, spacedim> fe_eval(
        mapping,
        dof_handler_dim.get_fe(),
        Quadrature<dim>(dof_handler_dim.get_fe().get_unit_support_points()),
        update_quadrature_points);

      Vector<double> temp;

      for (const auto &cell : dof_handler_dim.active_cell_iterators())
        {
          if (cell->is_locally_owned() == false)
            continue;

          fe_eval.reinit(cell);

          temp.reinit(fe_eval.dofs_per_cell);

          for (const auto q : fe_eval.quadrature_point_indices())
            {
              const auto point = fe_eval.quadrature_point(q);

              const unsigned int comp =
                dof_handler_dim.get_fe().system_to_component_index(q).first;

              temp[q] = point[comp];
            }

          cell->set_dof_values(temp, euler_coordinates_vector);
        }

      euler_coordinates_vector.update_ghost_values();
    }
  } // namespace VectorTools
} // namespace dealii

template <int dim, int spacedim, typename VectorType>
void
compute_normal(const Mapping<dim, spacedim>    &mapping,
               const DoFHandler<dim, spacedim> &dof_handler_dim,
               VectorType                      &normal_vector)
{
  FEValues<dim, spacedim> fe_eval_dim(
    mapping,
    dof_handler_dim.get_fe(),
    dof_handler_dim.get_fe().get_unit_support_points(),
    update_normal_vectors | update_gradients);

  Vector<double> normal_temp;

  for (const auto &cell : dof_handler_dim.active_cell_iterators())
    {
      if (cell->is_locally_owned() == false)
        continue;

      fe_eval_dim.reinit(cell);

      normal_temp.reinit(fe_eval_dim.dofs_per_cell);
      normal_temp = 0.0;

      for (const auto q : fe_eval_dim.quadrature_point_indices())
        {
          const auto normal = fe_eval_dim.normal_vector(q);

          const unsigned int comp =
            dof_handler_dim.get_fe().system_to_component_index(q).first;

          normal_temp[q] = normal[comp];
        }

      cell->set_dof_values(normal_temp, normal_vector);
    }

  normal_vector.update_ghost_values();
}



template <int dim, int spacedim, typename VectorType>
void
compute_curvature(const Mapping<dim, spacedim>    &mapping,
                  const DoFHandler<dim, spacedim> &dof_handler_dim,
                  const DoFHandler<dim, spacedim> &dof_handler,
                  const Quadrature<dim>            quadrature,
                  const VectorType                &normal_vector,
                  VectorType                      &curvature_vector)
{
  FEValues<dim, spacedim> fe_eval(mapping,
                                  dof_handler.get_fe(),
                                  quadrature,
                                  update_gradients);
  FEValues<dim, spacedim> fe_eval_dim(mapping,
                                      dof_handler_dim.get_fe(),
                                      quadrature,
                                      update_gradients);

  Vector<double> curvature_temp;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned() == false)
        continue;

      TriaIterator<DoFCellAccessor<dim, spacedim, false>> dof_cell_dim(
        &dof_handler_dim.get_triangulation(),
        cell->level(),
        cell->index(),
        &dof_handler_dim);

      fe_eval.reinit(cell);
      fe_eval_dim.reinit(dof_cell_dim);

      curvature_temp.reinit(quadrature.size());

      std::vector<std::vector<Tensor<1, spacedim, double>>> normal_gradients(
        quadrature.size(), std::vector<Tensor<1, spacedim, double>>(spacedim));

      fe_eval_dim.get_function_gradients(normal_vector, normal_gradients);

      for (const auto q : fe_eval_dim.quadrature_point_indices())
        {
          double curvature = 0.0;

          for (unsigned c = 0; c < spacedim; ++c)
            curvature += normal_gradients[q][c][c];

          curvature_temp[q] = curvature;
        }

      cell->set_dof_values(curvature_temp, curvature_vector);
    }

  curvature_vector.update_ghost_values();
}



template <int dim, int spacedim>
void
print(std::tuple<
      std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>,
      std::vector<std::vector<Point<dim>>>,
      std::vector<std::vector<unsigned int>>,
      std::vector<std::vector<Point<spacedim>>>,
      std::vector<std::vector<unsigned int>>> result)
{
  for (unsigned int i = 0; i < std::get<0>(result).size(); ++i)
    {
      const unsigned int n_points = std::get<1>(result)[i].size();

      std::cout << std::get<0>(result)[i]->level() << ' '
                << std::get<0>(result)[i]->index() << " x " << n_points
                << std::endl;

      for (unsigned int j = 0; j < n_points; ++j)
        {
          std::cout << std::get<1>(result)[i][j] << " - "
                    << std::get<2>(result)[i][j] << " - "
                    << std::get<3>(result)[i][j] << " - "
                    << std::get<4>(result)[i][j] << std::endl;
        }
    }
}



template <int dim, int spacedim, typename VectorType>
void
compute_force_vector_sharp_interface(
  const Mapping<dim, spacedim>    &surface_mapping,
  const DoFHandler<dim, spacedim> &surface_dofhandler,
  const DoFHandler<dim, spacedim> &surface_dofhandler_dim,
  const Quadrature<dim>           &surface_quadrature,
  const Mapping<spacedim>         &mapping,
  const DoFHandler<spacedim>      &dof_handler,
  const double                     surface_tension,
  const VectorType                &normal_vector,
  const VectorType                &curvature_vector,
  VectorType                      &force_vector)
{
  using T = double;

  const auto integration_points = [&]() {
    std::vector<Point<spacedim>> integration_points;

    FEValues<dim, spacedim> fe_eval(surface_mapping,
                                    surface_dofhandler.get_fe(),
                                    surface_quadrature,
                                    update_values | update_quadrature_points |
                                      update_JxW_values);

    for (const auto &cell :
         surface_dofhandler.get_triangulation().active_cell_iterators())
      {
        if (cell->is_locally_owned() == false)
          continue;

        fe_eval.reinit(cell);

        for (const auto q : fe_eval.quadrature_point_indices())
          integration_points.push_back(fe_eval.quadrature_point(q));
      }

    return integration_points;
  }();

  Utilities::MPI::RemotePointEvaluation<spacedim, spacedim> eval;
  eval.reinit(integration_points, dof_handler.get_triangulation(), mapping);

  const auto integration_values = [&]() {
    std::vector<T> integration_values;

    FEValues<dim, spacedim> fe_eval(surface_mapping,
                                    surface_dofhandler.get_fe(),
                                    surface_quadrature,
                                    update_values | update_quadrature_points |
                                      update_JxW_values);
    FEValues<dim, spacedim> fe_eval_dim(surface_mapping,
                                        surface_dofhandler_dim.get_fe(),
                                        surface_quadrature,
                                        update_values);

    const auto &tria_surface = surface_dofhandler.get_triangulation();

    for (const auto &cell : tria_surface.active_cell_iterators())
      {
        if (cell->is_locally_owned() == false)
          continue;

        TriaIterator<DoFCellAccessor<dim, spacedim, false>> dof_cell(
          &tria_surface, cell->level(), cell->index(), &surface_dofhandler);
        TriaIterator<DoFCellAccessor<dim, spacedim, false>> dof_cell_dim(
          &tria_surface, cell->level(), cell->index(), &surface_dofhandler_dim);

        fe_eval.reinit(dof_cell);
        fe_eval_dim.reinit(dof_cell_dim);

        std::vector<double>         curvature_values(fe_eval.dofs_per_cell);
        std::vector<Vector<double>> normal_values(fe_eval.dofs_per_cell,
                                                  Vector<double>(spacedim));

        fe_eval.get_function_values(curvature_vector, curvature_values);
        fe_eval_dim.get_function_values(normal_vector, normal_values);

        for (const auto q : fe_eval_dim.quadrature_point_indices())
          {
            for (unsigned int i = 0; i < spacedim; ++i)
              integration_values.push_back(-curvature_values[q] *
                                           normal_values[q][i] *
                                           fe_eval.JxW(q) * surface_tension);
          }
      }

    return integration_values;
  }();

  const auto fu = [&](const auto &values, const auto &cell_data) {
    AffineConstraints<double> constraints; // TODO: use the right ones

    FEPointEvaluation<spacedim, spacedim> phi_force(mapping,
                                                    dof_handler.get_fe(),
                                                    update_values);

    std::vector<double>                  buffer;
    std::vector<types::global_dof_index> local_dof_indices;

    for (unsigned int i = 0; i < cell_data.cells.size(); ++i)
      {
        typename DoFHandler<spacedim>::active_cell_iterator cell = {
          &eval.get_triangulation(),
          cell_data.cells[i].first,
          cell_data.cells[i].second,
          &dof_handler};

        local_dof_indices.resize(cell->get_fe().n_dofs_per_cell());
        buffer.resize(cell->get_fe().n_dofs_per_cell());

        cell->get_dof_indices(local_dof_indices);

        const ArrayView<const Point<spacedim>> unit_points(
          cell_data.reference_point_values.data() +
            cell_data.reference_point_ptrs[i],
          cell_data.reference_point_ptrs[i + 1] -
            cell_data.reference_point_ptrs[i]);

        const ArrayView<const Tensor<1, spacedim, T>> force_JxW(
          reinterpret_cast<const Tensor<1, spacedim, T> *>(values.data()) +
            cell_data.reference_point_ptrs[i],
          cell_data.reference_point_ptrs[i + 1] -
            cell_data.reference_point_ptrs[i]);

        phi_force.reinit(cell, unit_points);

        for (unsigned int q = 0; q < unit_points.size(); ++q)
          phi_force.submit_value(force_JxW[q], q);

        phi_force.test_and_sum(buffer, EvaluationFlags::values);

        constraints.distribute_local_to_global(buffer,
                                               local_dof_indices,
                                               force_vector);
      }
  };

  std::vector<T> buffer;

  eval.template process_and_evaluate<T, spacedim>(integration_values,
                                                  buffer,
                                                  fu);
}



template <int dim>
void
test()
{
  const unsigned int spacedim = dim + 1;

  const unsigned int fe_degree      = 3;
  const unsigned int mapping_degree = fe_degree;
  const unsigned int n_refinements  = 5;

  parallel::shared::Triangulation<dim, spacedim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim, spacedim>::none,
    true,
    parallel::shared::Triangulation<dim, spacedim>::Settings::partition_zoltan);
#if false
  GridGenerator::hyper_sphere(tria, Point<spacedim>(), 0.5);
#else
  GridGenerator::hyper_sphere(tria, Point<spacedim>(0.02, 0.03), 0.5);
#endif
  tria.refine_global(n_refinements);

  // quadrature rule and FE for curvature
  FE_Q<dim, spacedim>       fe(fe_degree);
  QGaussLobatto<dim>        quadrature(fe_degree + 1);
  DoFHandler<dim, spacedim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  // FE for normal
  FESystem<dim, spacedim>   fe_dim(fe, spacedim);
  DoFHandler<dim, spacedim> dof_handler_dim(tria);
  dof_handler_dim.distribute_dofs(fe_dim);

  // Set up MappingFEField
  Vector<double> euler_vector(dof_handler_dim.n_dofs());
  VectorTools::get_position_vector(dof_handler_dim,
                                   euler_vector,
                                   MappingQ<dim, spacedim>(mapping_degree));
  MappingFEField<dim, spacedim> mapping(dof_handler_dim, euler_vector);


  // compute normal vector
  VectorType normal_vector(create_partitioner(dof_handler_dim));
  compute_normal(mapping, dof_handler_dim, normal_vector);

  // compute curvature
  VectorType curvature_vector(create_partitioner(dof_handler));
  compute_curvature(mapping,
                    dof_handler_dim,
                    dof_handler,
                    quadrature,
                    normal_vector,
                    curvature_vector);

#if false
  const unsigned int background_n_global_refinements = 6;
#else
  const unsigned int background_n_global_refinements =
    Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1 ? 40 : 80;
#endif
  const unsigned int background_fe_degree = 2;

  parallel::shared::Triangulation<spacedim> background_tria(
    MPI_COMM_WORLD,
    Triangulation<spacedim>::none,
    true,
    parallel::shared::Triangulation<spacedim>::Settings::partition_zoltan);
#if false
  GridGenerator::hyper_cube(background_tria, -1.0, +1.0);
#else
  GridGenerator::subdivided_hyper_cube(background_tria,
                                       background_n_global_refinements,
                                       -2.5,
                                       2.5);
#endif
  if (background_n_global_refinements < 20)
    background_tria.refine_global(background_n_global_refinements);

  FESystem<spacedim>   background_fe(FE_Q<spacedim>{background_fe_degree},
                                   spacedim);
  DoFHandler<spacedim> background_dof_handler(background_tria);
  background_dof_handler.distribute_dofs(background_fe);

  MappingQ1<spacedim> background_mapping;

  VectorType force_vector_sharp_interface(
    create_partitioner(background_dof_handler));

  // write computed vectors to Paraview
  if (false)
    {
      GridOut().write_mesh_per_processor_as_vtu(tria, "grid_surface");
      GridOut().write_mesh_per_processor_as_vtu(background_tria,
                                                "grid_background");
    }

  compute_force_vector_sharp_interface(mapping,
                                       dof_handler,
                                       dof_handler_dim,
                                       QGauss<dim>(fe_degree + 1),
                                       background_mapping,
                                       background_dof_handler,
                                       1.0,
                                       normal_vector,
                                       curvature_vector,
                                       force_vector_sharp_interface);

  force_vector_sharp_interface.update_ghost_values();

  // write computed vectors to Paraview
  if (true)
    {
      DataOutBase::VtkFlags flags;
      // flags.write_higher_order_cells = true;

      DataOut<dim, spacedim> data_out;
      data_out.set_flags(flags);
      data_out.add_data_vector(dof_handler, curvature_vector, "curvature");
      data_out.add_data_vector(dof_handler_dim, normal_vector, "normal");

      data_out.build_patches(
        mapping,
        fe_degree + 1,
        DataOut<dim, spacedim>::CurvedCellRegion::curved_inner_cells);
      data_out.write_vtu_with_pvtu_record("./",
                                          "data_surface",
                                          0,
                                          MPI_COMM_WORLD);
    }

  if (true)
    {
      DataOutBase::VtkFlags flags;
      flags.write_higher_order_cells = true;

      DataOut<spacedim> data_out;
      data_out.set_flags(flags);
      data_out.attach_dof_handler(background_dof_handler);
      data_out.add_data_vector(background_dof_handler,
                               force_vector_sharp_interface,
                               "force");

      data_out.build_patches(background_mapping, background_fe_degree + 1);
      data_out.write_vtu_with_pvtu_record("./",
                                          "data_background",
                                          0,
                                          MPI_COMM_WORLD);
    }

  force_vector_sharp_interface.print(deallog.get_file_stream());
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    all;

  test<1>();
}

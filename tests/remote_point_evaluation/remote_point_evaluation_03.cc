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

// Test Utilities::MPI::RemotePointEvaluation by creating a position vector
// and evaluating the vector at given points and comparing the return value
// with the point positions.

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
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/lac/la_parallel_block_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_evaluate.h>

#include "../tests.h"



using VectorType = LinearAlgebra::distributed::Vector<double>;

template <int dim>
class NormalFunction : public Function<dim>
{
public:
  NormalFunction()
    : Function<dim>(dim)
  {}

  double
  value(const Point<dim> &p, const unsigned int component) const
  {
    return p[component];
  }
};

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
test_1(const Triangulation<dim, spacedim> &surface_mesh,
       const Mapping<dim, spacedim>       &surface_mapping,
       const Quadrature<dim>              &surface_quadrature,
       const Mapping<spacedim>            &mapping,
       const DoFHandler<spacedim>         &dof_handler_dim,
       const VectorType                   &normal_solution)
{
  using T = Point<spacedim>;

  const auto integration_points = [&]() {
    std::vector<Point<spacedim>> integration_points;

    FE_Nothing<dim, spacedim> dummy;

    FEValues<dim, spacedim> fe_eval(surface_mapping,
                                    dummy,
                                    surface_quadrature,
                                    update_quadrature_points);

    for (const auto &cell : surface_mesh.active_cell_iterators())
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
  eval.reinit(integration_points, dof_handler_dim.get_triangulation(), mapping);

  const auto integration_values = [&]() {
    std::vector<T> integration_values;

    FE_Nothing<dim, spacedim> dummy;

    FEValues<dim, spacedim> fe_eval(surface_mapping,
                                    dummy,
                                    surface_quadrature,
                                    update_quadrature_points);

    for (const auto &cell : surface_mesh.active_cell_iterators())
      {
        if (cell->is_locally_owned() == false)
          continue;

        fe_eval.reinit(cell);

        for (const auto q : fe_eval.quadrature_point_indices())
          integration_values.push_back(fe_eval.quadrature_point(q));
      }

    return integration_values;
  }();

  const auto fu = [&](const auto &values, const auto &cell_data) {
    AffineConstraints<double> constraints; // TODO: use the right ones

    FEPointEvaluation<spacedim, spacedim> phi_normal(mapping,
                                                     dof_handler_dim.get_fe(),
                                                     update_values);
    FEPointEvaluation<spacedim, spacedim> phi_force(mapping,
                                                    dof_handler_dim.get_fe(),
                                                    update_values);

    std::vector<double>                  buffer;
    std::vector<double>                  buffer_dim;
    std::vector<types::global_dof_index> local_dof_indices;
    std::vector<types::global_dof_index> local_dof_indices_dim;

    for (unsigned int i = 0; i < cell_data.cells.size(); ++i)
      {
        typename DoFHandler<spacedim>::active_cell_iterator cell_dim = {
          &eval.get_triangulation(),
          cell_data.cells[i].first,
          cell_data.cells[i].second,
          &dof_handler_dim};

        const ArrayView<const Point<spacedim>> unit_points(
          cell_data.reference_point_values.data() +
            cell_data.reference_point_ptrs[i],
          cell_data.reference_point_ptrs[i + 1] -
            cell_data.reference_point_ptrs[i]);

        const ArrayView<const T> JxW(values.data() +
                                       cell_data.reference_point_ptrs[i],
                                     cell_data.reference_point_ptrs[i + 1] -
                                       cell_data.reference_point_ptrs[i]);

        // gather_evaluate normal
        {
          local_dof_indices_dim.resize(cell_dim->get_fe().n_dofs_per_cell());
          buffer_dim.resize(cell_dim->get_fe().n_dofs_per_cell());

          cell_dim->get_dof_indices(local_dof_indices_dim);

          constraints.get_dof_values(normal_solution,
                                     local_dof_indices_dim.begin(),
                                     buffer_dim.begin(),
                                     buffer_dim.end());

          phi_normal.reinit(cell_dim, unit_points);
          phi_normal.evaluate(make_array_view(buffer_dim),
                              EvaluationFlags::values);
        }

        // perform operation on quadrature points
        for (unsigned int q = 0; q < unit_points.size(); ++q)
          {
            if (JxW[q].distance(Point<spacedim>(phi_normal.get_value(q))) >
                1e-7)
              std::cout << JxW[q] << ' ' << phi_normal.get_value(q)
                        << std::endl;
          }
      }
  };

  std::vector<T> buffer;

  eval.template process_and_evaluate<T>(integration_values, buffer, fu);
}



template <int dim, int spacedim, typename VectorType>
void
test_2(const Triangulation<dim, spacedim> &surface_mesh,
       const Mapping<dim, spacedim>       &surface_mapping,
       const Quadrature<dim>              &surface_quadrature,
       const Mapping<spacedim>            &mapping,
       const DoFHandler<spacedim>         &dof_handler_dim,
       const VectorType                   &normal_solution)
{
  using T = Point<spacedim>;

  const auto integration_points = [&]() {
    std::vector<Point<spacedim>> integration_points;

    FE_Nothing<dim, spacedim> dummy;

    FEValues<dim, spacedim> fe_eval(surface_mapping,
                                    dummy,
                                    surface_quadrature,
                                    update_quadrature_points);

    for (const auto &cell : surface_mesh.active_cell_iterators())
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
  eval.reinit(integration_points, dof_handler_dim.get_triangulation(), mapping);

  const auto evaluation_point_results =
    VectorTools::point_values<spacedim>(eval, dof_handler_dim, normal_solution);

  for (unsigned int i = 0; i < evaluation_point_results.size(); ++i)
    {
      if (integration_points[i].distance(
            Point<spacedim>(evaluation_point_results[i])) > 1e-7)
        std::cout << integration_points[i] << ' ' << evaluation_point_results[i]
                  << std::endl;
    }
}



template <int dim>
void
test(const MPI_Comm comm)
{
  const unsigned int spacedim = dim + 1;

  const unsigned int fe_degree      = 3;
  const unsigned int mapping_degree = fe_degree;
  const unsigned int n_refinements  = 5;

  parallel::shared::Triangulation<dim, spacedim> tria(
    comm,
    Triangulation<dim, spacedim>::none,
    true,
    parallel::shared::Triangulation<dim, spacedim>::Settings::partition_zorder);
#if false
  GridGenerator::hyper_sphere(tria, Point<spacedim>(), 0.5);
#else
  GridGenerator::hyper_sphere(tria, Point<spacedim>(0.02, 0.03), 1.0);
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


#if false
  const unsigned int background_n_global_refinements = 6;
#else
  const unsigned int background_n_global_refinements = 40;
#endif
  const unsigned int background_fe_degree = 2;

  parallel::shared::Triangulation<spacedim> background_tria(
    comm,
    Triangulation<spacedim>::none,
    true,
    parallel::shared::Triangulation<spacedim>::Settings::partition_zorder);
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

  FESystem<spacedim>   background_fe_dim(FE_Q<spacedim>{background_fe_degree},
                                       spacedim);
  DoFHandler<spacedim> background_dof_handler_dim(background_tria);
  background_dof_handler_dim.distribute_dofs(background_fe_dim);

  MappingQ1<spacedim> background_mapping;

  VectorType force_vector_sharp_interface(
    create_partitioner(background_dof_handler_dim));
  VectorType normal_vector(create_partitioner(background_dof_handler_dim));

  VectorTools::interpolate(background_mapping,
                           background_dof_handler_dim,
                           NormalFunction<spacedim>(),
                           normal_vector);

  normal_vector.update_ghost_values();

  // write computed vectors to Paraview
  if (true)
    {
      GridOut().write_mesh_per_processor_as_vtu(tria, "grid_surface");
      GridOut().write_mesh_per_processor_as_vtu(background_tria,
                                                "grid_background");
    }

  test_1(tria,
         mapping,
         QGauss<dim>(fe_degree + 1),
         background_mapping,
         background_dof_handler_dim,
         normal_vector);
  deallog << "OK1!" << std::endl;

  test_2(tria,
         mapping,
         QGauss<dim>(fe_degree + 1),
         background_mapping,
         background_dof_handler_dim,
         normal_vector);
  deallog << "OK2!" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    all;

  const int n_ranks = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const int my_rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  for (int i = 1; i <= n_ranks; ++i)
    {
      MPI_Comm comm;
      MPI_Comm_split(MPI_COMM_WORLD, my_rank < i, my_rank, &comm);

      if (my_rank < i)
        {
          deallog.push("n_ranks=" + std::to_string(i));
          test<1>(comm);
          deallog.pop();
        }

      MPI_Comm_free(&comm);
    }
}

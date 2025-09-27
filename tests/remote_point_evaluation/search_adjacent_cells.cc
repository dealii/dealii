// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q_cache.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_evaluate.h>

#include <deal.II/particles/data_out.h>
#include <deal.II/particles/particle_handler.h>

#include "../tests.h"

template <int dim>
double
displacement(const Point<dim>  &point,
             const unsigned int component,
             const double       factor)
{
  if (component == 0)
    return (factor * std::pow(point[1], 2) * std::pow(point[2], 2));
  else
    return 0.0;
}

template <int dim>
void
construct_triangulation(Triangulation<dim> &tria)
{
  const double L_F = 0.7;
  const double B_F = 1.0;
  const double H_F = 0.5;

  const double T_S = 0.05;
  const double B_S = 0.6;
  const double H_S = 0.4;

  const double L_IN = 0.6;

  const unsigned int N_CELLS_X_OUTFLOW = 1;
  const unsigned int N_CELLS_Y_LOWER   = 2;
  const unsigned int N_CELLS_Z_MIDDLE  = 2;

  std::vector<dealii::Triangulation<3>> tria_vec;
  tria_vec.resize(4);

  dealii::GridGenerator::subdivided_hyper_rectangle(
    tria_vec[0],
    std::vector<unsigned int>({N_CELLS_X_OUTFLOW, 1, N_CELLS_Z_MIDDLE}),
    dealii::Point<3>(L_IN + T_S, H_S - H_F / 2.0, -B_S / 2.0),
    dealii::Point<3>(L_F, H_F / 2.0, B_S / 2.0));

  dealii::GridGenerator::subdivided_hyper_rectangle(
    tria_vec[1],
    std::vector<unsigned int>(
      {N_CELLS_X_OUTFLOW, N_CELLS_Y_LOWER, N_CELLS_Z_MIDDLE}),
    dealii::Point<3>(L_IN + T_S, -H_F / 2.0, -B_S / 2.0),
    dealii::Point<3>(L_F, H_S - H_F / 2.0, B_S / 2.0));

  dealii::GridGenerator::subdivided_hyper_rectangle(
    tria_vec[2],
    std::vector<unsigned int>({N_CELLS_X_OUTFLOW, N_CELLS_Y_LOWER, 1}),
    dealii::Point<3>(L_IN + T_S, -H_F / 2.0, -B_F / 2.0),
    dealii::Point<3>(L_F, H_S - H_F / 2.0, -B_S / 2.0));

  dealii::GridGenerator::subdivided_hyper_rectangle(
    tria_vec[3],
    std::vector<unsigned int>({1, N_CELLS_Y_LOWER, 1}),
    dealii::Point<3>(L_IN, -H_F / 2.0, -B_F / 2.0),
    dealii::Point<3>(L_IN + T_S, H_S - H_F / 2.0, -B_S / 2.0));

  std::vector<const dealii::Triangulation<3> *> tria_vec_ptr(tria_vec.size());
  for (unsigned int i = 0; i < tria_vec.size(); ++i)
    tria_vec_ptr[i] = &tria_vec[i];

  dealii::GridGenerator::merge_triangulations(tria_vec_ptr, tria, 1.e-10);
}

template <int dim>
void
test(const unsigned int mapping_degree,
     const double       mapping_factor,
     const double       tolerance,
     const unsigned int n_points)
{
  Triangulation<dim> tria;
  construct_triangulation(tria);

  FESystem<dim>   fe(FE_Q<dim>(mapping_degree), dim);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  Quadrature<dim>    quadrature(fe.base_element(0).get_unit_support_points());
  FEValues<dim>      fe_values(fe,
                          quadrature,
                          update_quadrature_points | update_values);
  MappingQCache<dim> mapping(mapping_degree);
  mapping.initialize(
    tria,
    [&](const typename dealii::Triangulation<dim>::cell_iterator &cell_tria)
      -> std::vector<dealii::Point<dim>> {
      std::vector<dealii::Point<dim>> grid_coordinates(quadrature.size());

      fe_values.reinit(cell_tria);
      // extract displacement and add to original position
      for (unsigned int i = 0; i < grid_coordinates.size(); ++i)
        {
          grid_coordinates[i] = fe_values.quadrature_point(i);
          grid_coordinates[i][0] +=
            displacement(grid_coordinates[i], 0, mapping_factor);
        }
      return grid_coordinates;
    });

  // initialize points with noise
  std::vector<Point<dim>> points(n_points * n_points);
  unsigned int            point_idx = 0;
  for (unsigned int i = 0; i < n_points; ++i)
    {
      for (unsigned int j = 0; j < n_points; ++j)
        {
          points[point_idx][0] = 0.65;
          points[point_idx][1] = -0.25 + 0.4 * (1.0 / (n_points - 1) * j);
          points[point_idx][2] = -0.3 + 0.6 * (1.0 / (n_points - 1) * i);

          points[point_idx][0] +=
            displacement(points[point_idx], 0, mapping_factor) +
            (tolerance * 0.1) * (-0.5 + (rand() / (double)RAND_MAX));

          point_idx += 1;
        }
    }

  // initialize RPE without any marked points
  typename Utilities::MPI::RemotePointEvaluation<dim>::AdditionalData
    additional_data;
  additional_data.tolerance = tolerance;
  dealii::Utilities::MPI::RemotePointEvaluation<dim> rpe(additional_data);
  rpe.reinit(points, tria, mapping);

  unsigned int                    n_points_not_found_rpe = 0;
  std::vector<dealii::Point<dim>> points_not_found;
  if (!rpe.all_points_found())
    {
      // get vector of points not found
      points_not_found.reserve(points.size());

      for (unsigned int i = 0; i < points.size(); ++i)
        {
          if (!rpe.point_found(i))
            {
              n_points_not_found_rpe += 1;
              points_not_found.push_back(points[i]);
            }
        }
    }
  std::cout << "points not found (no marked points)  : "
            << n_points_not_found_rpe << "\n";

  // initialize RPE with all points marked
  std::vector<bool> marked_vertices(tria.n_vertices(), true);

  typename Utilities::MPI::RemotePointEvaluation<dim>::AdditionalData
    additional_data2(tolerance, false, 0, [marked_vertices]() {
      return marked_vertices;
    });

  dealii::Utilities::MPI::RemotePointEvaluation<dim> rpe2(additional_data2);

  rpe2.reinit(points, tria, mapping);

  unsigned int                    n_points_not_found_rpe2 = 0;
  std::vector<dealii::Point<dim>> points_not_found2;
  if (!rpe2.all_points_found())
    {
      // get vector of points not found
      std::vector<dealii::Point<dim>> points_not_found2;
      points_not_found2.reserve(points.size());

      for (unsigned int i = 0; i < points.size(); ++i)
        {
          if (!rpe2.point_found(i))
            {
              n_points_not_found_rpe2 += 1;
              points_not_found2.push_back(points[i]);
            }
        }
    }
  std::cout << "points not found (all points marked) : "
            << n_points_not_found_rpe2 << "\n";

  // output in case of failure
  if (n_points_not_found_rpe > 0 || n_points_not_found_rpe2 > 0)
    {
      // output triangulation
      DataOut<dim>          data_out;
      DataOutBase::VtkFlags flags;
      flags.write_higher_order_cells = true;
      data_out.set_flags(flags);
      data_out.attach_triangulation(tria);
      data_out.build_patches(mapping,
                             mapping_degree,
                             DataOut<dim>::curved_inner_cells);
      std::ofstream stream("grid.vtu");
      data_out.write_vtu(stream);

      // enclosing dummy triangulation for point plots
      dealii::BoundingBox<dim> bounding_box(points);
      const auto               boundary_points =
        bounding_box.create_extended_relative(1e-3).get_boundary_points();

      dealii::Triangulation<dim> particle_dummy_tria;
      dealii::GridGenerator::hyper_rectangle(particle_dummy_tria,
                                             boundary_points.first,
                                             boundary_points.second);

      dealii::MappingQGeneric<dim> particle_dummy_mapping(
        1 /* mapping_degree */);

      // output all points
      {
        dealii::Particles::ParticleHandler<dim, dim> particle_handler(
          particle_dummy_tria, particle_dummy_mapping);
        particle_handler.insert_particles(points);
        dealii::Particles::DataOut<dim, dim> particle_output;
        particle_output.build_patches(particle_handler);
        std::ofstream filestream("all_points.vtu");
        particle_output.write_vtu(filestream);
      }

      // output points not found (no vertices marked)
      {
        dealii::Particles::ParticleHandler<dim, dim> particle_handler(
          particle_dummy_tria, particle_dummy_mapping);
        particle_handler.insert_particles(points_not_found);
        dealii::Particles::DataOut<dim, dim> particle_output;
        particle_output.build_patches(particle_handler);
        std::ofstream filestream("points_not_found.vtu");
        particle_output.write_vtu(filestream);
      }

      // output points not found
      {
        dealii::Particles::ParticleHandler<dim, dim> particle_handler(
          particle_dummy_tria, particle_dummy_mapping);
        particle_handler.insert_particles(points_not_found2);
        dealii::Particles::DataOut<dim, dim> particle_output;
        particle_output.build_patches(particle_handler);
        std::ofstream filestream("points_not_found2.vtu");
        particle_output.write_vtu(filestream);
      }
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  deallog << std::setprecision(8) << std::fixed;

  test<3>(2, 100.0, 1e-5, 100);
  test<3>(3, 100.0, 1e-5, 100);
  test<3>(4, 100.0, 1e-5, 100);
}

// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// Verify that hp::DoFHandler::distribute_dofs runs in linear time for a mesh
// which (prior to commit 23ab0eb0c0) it previously ran in quadratic time. This
// test checks this in three ways:
// 1. The amount of time needed to create the grid and the amount of time
//    needed to distribute dofs should be nearly equal for bilinears. The test
//    below checks that they are within 50% of each-other.
// 2. The linear interpolant of the timing data should be fairly accurate. The test
//    below checks that the least squares error is not too large.
// 3. Finally, this test should time out if the old quadratic time algorithm is
//    used (that numbering algorithm takes roughly 1000 seconds on 2016
//    hardware, but the linear time algorithm takes about 6 seconds)

#include "../tests.h"
#include <deal.II/base/mpi.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/householder.h>

#include <deal.II/base/std_cxx11/array.h>
#include <deal.II/base/std_cxx11/shared_ptr.h>

#include <iostream>

using namespace dealii;

static const types::manifold_id circular_manifold_id = 1;
static const types::manifold_id straight_manifold_id = 3;

/*
 * Declaration has default arguments
 *
 * Slight modification of Konstantin Ladutenko's better conditioned circular
 * grid: see step-6.
 */
template <int dim>
std_cxx11::shared_ptr<dealii::Manifold<dim> >
ladutenko_circle(dealii::Triangulation<dim> &triangulation,
                 const dealii::Point<dim>    center = dealii::Point<dim>(),
                 const double                radius = 1.0);

template <int dim>
std_cxx11::shared_ptr<Manifold<dim > >
ladutenko_circle(Triangulation<dim> &triangulation,
                 const Point<dim>    center,
                 const double        radius)
{
  std_cxx11::shared_ptr<Manifold<dim> > boundary(new SphericalManifold<dim>(center));
  GridGenerator::hyper_ball (triangulation, center, radius);
  triangulation.set_all_manifold_ids(circular_manifold_id);
  triangulation.set_manifold (circular_manifold_id, *boundary);

  const double core_radius  = 1.0/4.8*radius;
  const double inner_radius = 1.0/2.4*radius;

  // Step 1: Shrink the inner cell
  // and
  // Step 2: set the central cell to have a straight manifold
  // and
  // Step 3: Refine all cells except the central one
  typename Triangulation<dim>::active_cell_iterator
  cell = triangulation.begin_active(),
  endc = triangulation.end();
  for (; cell != endc; ++cell)
    {
      if (cell->center().distance(center) < 1e-10)
        {
          for (unsigned int vertex_n = 0;
               vertex_n < GeometryInfo<dim>::vertices_per_cell;
               ++vertex_n)
            {
              cell->vertex(vertex_n) *= core_radius
                                        /center.distance (cell->vertex(vertex_n));
            }
          cell->set_all_manifold_ids(straight_manifold_id);
        }
      else
        {
          cell->set_refine_flag ();
        }
    }
  triangulation.execute_coarsening_and_refinement();

  // Step 4: Resize the inner children of the outer cells
  // and
  // Step 5: Refine the outer loop
  cell = triangulation.begin_active();
  for (; cell != endc; ++cell)
    {
      for (unsigned int v=0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
        {
          const double dist = center.distance (cell->vertex(v));
          if (dist > core_radius*1.0001 && dist < radius - 1.0e-5)
            cell->vertex(v) *= inner_radius/dist;
        }
      if (cell->at_boundary())
        {
          cell->set_refine_flag ();
        }
    }
  triangulation.execute_coarsening_and_refinement();

  return boundary;
}

template <int dim>
class QuadraticTimeCircle
{
public:
  QuadraticTimeCircle(const unsigned int n_global_refines);

  void run();

  Timer create_grid_timer;
  Timer distribute_dofs_timer;

  const unsigned int n_global_refines;

  std_cxx11::shared_ptr<Manifold<dim> > boundary_manifold;
  Triangulation<dim> triangulation;
  hp::FECollection<dim> finite_elements;
  hp::DoFHandler<dim> dof_handler;

  void setup_dofs();
};



template <int dim>
QuadraticTimeCircle<dim>::QuadraticTimeCircle(const unsigned int n_global_refines) :
  n_global_refines (n_global_refines),
  dof_handler(triangulation)
{
  create_grid_timer.start();

  boundary_manifold = ladutenko_circle(triangulation);
  typename Triangulation<dim>::active_cell_iterator
  cell = triangulation.begin_active(),
  endc = triangulation.end();
  for (; cell != endc; ++cell)
    {
      // do not use any curved cells on the interior.
      if (!cell->at_boundary())
        {
          cell->set_all_manifold_ids(straight_manifold_id);
        }
    }
  triangulation.refine_global(n_global_refines);

  finite_elements.push_back(FE_Q<dim>(1));
  finite_elements.push_back(FE_Q<dim>(1));
  create_grid_timer.stop();
}



template <int dim>
void QuadraticTimeCircle<dim>::setup_dofs()
{
  deallog << "Number of cells: " << triangulation.n_active_cells() << std::endl;

  typename hp::DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  {
    cell->set_active_fe_index(0);
  }

  {
    distribute_dofs_timer.start();
    dof_handler.distribute_dofs(finite_elements); // <-- test this!
    distribute_dofs_timer.stop();
  }
  deallog << "Number of DoFs: " << dof_handler.n_dofs() << std::endl;
}



template <int dim>
void QuadraticTimeCircle<dim>::run()
{
  setup_dofs();
}



int main(int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);

  static const int dim = 2;
  static const unsigned int n_cycles = 6;

  unsigned int n_global_refines = 2;

  Vector<double> n_dofs(n_cycles);
  Vector<double> distribute_dofs_run_times(n_cycles);
  Vector<double> create_grid_run_times(n_cycles);
  FullMatrix<double> least_squares_data(n_cycles, 1);
  for (unsigned int i = 0; i < n_cycles; ++i)
    {
      n_global_refines += 1;
      QuadraticTimeCircle<dim> quadratic_time_circle(n_global_refines);
      quadratic_time_circle.run();
      n_dofs[i] = quadratic_time_circle.dof_handler.n_dofs();
      create_grid_run_times[i] = quadratic_time_circle.create_grid_timer.wall_time();
      distribute_dofs_run_times[i] = quadratic_time_circle.distribute_dofs_timer.wall_time();

      // if the times are sufficiently large they should be nearly equal
      if (create_grid_run_times[i] > 0.1)
        {
          AssertThrow(std::abs(create_grid_run_times[i] - distribute_dofs_run_times[i])
                      /create_grid_run_times[i] < 0.5,
                      ExcMessage("The run times for creating the grid and "
                                 "distributing dofs should be within 50% for "
                                 "sufficiently large grids."));
        }
      least_squares_data(i, 0) = n_dofs[i];
    }

  // finally, check that the run time is linear:
  Householder<double> householder(least_squares_data);
  Vector<double> slope(1);
  const double error = householder.least_squares(slope, distribute_dofs_run_times);
  // on my machine the error is about 0.1 so this is a substantial margin for
  // error
  AssertThrow(error < 0.3, ExcMessage("The run time should be close to linear."));

  deallog << "OK" << std::endl;
}

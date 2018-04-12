// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

// Test the distributed version of FEFieldFunction::vector_value_list

#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

using namespace dealii;

template <int dim>
void
test_vector_value_list(unsigned int ref_cube, unsigned int ref_sphere)
{
  MPI_Comm     mpi_communicator = MPI_COMM_WORLD;
  unsigned int n_procs = Utilities::MPI::n_mpi_processes(mpi_communicator);

  deallog << "Testing for dim = " << dim << " on " << n_procs << " processes"
          << std::endl;
  deallog << "Cube refinements: " << ref_cube << std::endl;
  deallog << "Sphere refinements:" << ref_sphere << std::endl;

  // Initializing classes:
  // Serial version of meshes and FeFieldFunction for comparison
  Triangulation<dim> cube_s;
  GridGenerator::hyper_cube(cube_s);
  cube_s.refine_global(ref_cube);
  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler_s(cube_s);
  dof_handler_s.distribute_dofs(fe);
  // Creating a serial data vector
  Functions::CosineFunction<dim> f_cosine;
  Vector<double>                 data_s(dof_handler_s.n_dofs());
  VectorTools::interpolate(dof_handler_s, f_cosine, data_s);

  Functions::FEFieldFunction<dim> fe_function_s(
    dof_handler_s, data_s, StaticMappingQ1<dim, dim>::mapping);

  // We shall use the points on the sphere/circle as shared points among all
  // processes
  Triangulation<dim - 1, dim> sphere;
  Point<dim>                  sphere_center;
  // Defining center and radius
  for (unsigned int i = 0; i < dim; ++i)
    sphere_center[i] = 0.47 - i * 0.05;
  double radius = 0.4 - dim * 0.05;
  GridGenerator::hyper_sphere(sphere, sphere_center, radius);
  static SphericalManifold<dim - 1, dim> surface_description(sphere_center);
  sphere.set_all_manifold_ids(0);
  sphere.set_manifold(0, surface_description);
  sphere.refine_global(ref_sphere);

  deallog << "Sphere center:" << sphere_center << std::endl;
  deallog << "Sphere radius:" << radius << std::endl;
  std::vector<Point<dim>> points;
  for (auto cell : sphere.active_cell_iterators())
    points.emplace_back(cell->center());

  // Calling the serial version of the function
  std::vector<Vector<double>> values_s(points.size());
  fe_function_s.vector_value_list(points, values_s);
  deallog << "Serial part: completed" << std::endl;

  // Distributed cube and FeFieldFunction, on which to test the distributed
  // version of FeFieldFunction

  parallel::distributed::Triangulation<dim> cube_d(mpi_communicator);
  GridGenerator::hyper_cube(cube_d);
  cube_d.refine_global(ref_cube);
  DoFHandler<dim> dof_handler_d(cube_d);
  dof_handler_d.distribute_dofs(fe);


  // Computing bounding boxes describing the locally owned part of the mesh
  IteratorFilters::LocallyOwnedCell locally_owned_cell_predicate;
  std::vector<BoundingBox<dim>>     local_bbox =
    GridTools::compute_mesh_predicate_bounding_box(
      cube_d, locally_owned_cell_predicate);

  // Obtaining the global mesh description through an all to all communication
  std::vector<std::vector<BoundingBox<dim>>> global_bboxes;
  global_bboxes = Utilities::MPI::all_gather(mpi_communicator, local_bbox);

  IndexSet ghost;
  DoFTools::extract_locally_relevant_dofs(dof_handler_d, ghost);

  LinearAlgebra::distributed::Vector<double> data_d_ghosted(
    dof_handler_d.locally_owned_dofs(), ghost, mpi_communicator);
  data_d_ghosted.zero_out_ghosts();

  VectorTools::interpolate(dof_handler_d, f_cosine, data_d_ghosted);
  data_d_ghosted.update_ghost_values();

  Functions::FEFieldFunction<dim,
                             DoFHandler<dim>,
                             LinearAlgebra::distributed::Vector<double>>
    fe_function_d(dof_handler_d,
                  data_d_ghosted,
                  StaticMappingQ1<dim, dim>::mapping,
                  true);


  fe_function_d.set_up_bounding_boxes(global_bboxes);

  // Calling the distributed version of the function
  deallog << "Using the distributed version of the function:" << std::endl;
  std::vector<Vector<double>> values_d(points.size());
  fe_function_d.vector_value_list(points, values_d);

  deallog << "Checking results" << std::endl;
  bool test_passed = true;

  for (unsigned int i = 0; i < points.size(); ++i)
    {
      values_d[i] -= values_s[i];
      if ((values_d[i]).l2_norm() > 1e-12)
        {
          deallog << "ERROR: value " << i << " is wrong!" << std::endl;
          deallog << "L2 norm for the difference is: "
                  << (values_d[i]).l2_norm() << std::endl;
          test_passed = false;
        }
    }

  if (test_passed)
    deallog << "Test passed" << std::endl;
  else
    deallog << "Test FAILED" << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog << "Deal.II FEFieldFunction::vector_value_list" << std::endl;
  deallog << "for the distributed version" << std::endl;
  deallog << "2D tests:" << std::endl;
  test_vector_value_list<2>(3, 3);
  deallog << "3D tests" << std::endl;
  test_vector_value_list<3>(3, 2);
}

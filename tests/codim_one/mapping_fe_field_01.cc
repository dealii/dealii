// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// DataOut::build_patches appeared to have a problem when outputting
// data using MappingFEField in combination with high subdivisions
// and/or high degree. This turned out to be related to multiple
// threads accessing the same global values for local_dof_values and
// local_dof_indices, where they should have been accessing thread
// local data. Running this test on a version of deal.II without the
// fix produces an "exploded" output, where the location of the
// vertices are randomly mixed between cells.

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test(const unsigned int refs,
     const unsigned int degree,
     const unsigned int subdivisions)
{
  const unsigned int id = degree + 10 * refs + 100 * subdivisions;

  Triangulation<dim, spacedim> triangulation;

  FE_Q<dim, spacedim>     fe(degree);
  FESystem<dim, spacedim> fe_euler(FE_Q<dim, spacedim>(degree), spacedim);

  DoFHandler<dim, spacedim> dof_handler(triangulation);
  DoFHandler<dim, spacedim> map_dh(triangulation);

  Vector<double> euler_vec;
  Vector<double> scal_sol;

  GridIn<dim, spacedim> grid_in;
  grid_in.attach_triangulation(triangulation);
  std::ifstream fname(SOURCE_DIR "/grids/sphere_0.inp");
  grid_in.read_ucd(fname);

  SphericalManifold<dim, spacedim> manifold;
  triangulation.set_all_manifold_ids(0);
  triangulation.set_manifold(0, manifold);

  triangulation.refine_global(refs);
  dof_handler.distribute_dofs(fe);
  map_dh.distribute_dofs(fe_euler);

  euler_vec.reinit(map_dh.n_dofs());
  scal_sol.reinit(dof_handler.n_dofs());
  scal_sol = 1;
  VectorTools::get_position_vector(map_dh, euler_vec);

  MappingFEField<dim, spacedim> mapping(map_dh, euler_vec);
  DataOut<dim, spacedim>        data_out_scal;
  data_out_scal.attach_dof_handler(dof_handler);

  data_out_scal.add_data_vector(scal_sol,
                                "scalar_data",
                                DataOut<dim, spacedim>::type_dof_data);

  data_out_scal.build_patches(mapping,
                              subdivisions,
                              DataOut<dim, spacedim>::curved_inner_cells);

  std::string filename_scal =
    ("scal_check_" + Utilities::int_to_string(id) + ".vtu");
  std::ofstream file_scal(filename_scal);
  data_out_scal.write_vtu(file_scal);
  data_out_scal.write_vtk(deallog.get_file_stream());


  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      spacedim, DataComponentInterpretation::component_is_part_of_vector);

  DataOut<dim, spacedim> data_out_euler;
  data_out_euler.attach_dof_handler(map_dh);

  data_out_euler.add_data_vector(euler_vec,
                                 "euler_vec",
                                 DataOut<dim, spacedim>::type_dof_data,
                                 data_component_interpretation);
  data_out_euler.build_patches(mapping,
                               degree,
                               DataOut<dim, spacedim>::curved_inner_cells);

  std::string filename_euler =
    ("euler_check_" + Utilities::int_to_string(id) + ".vtu");
  std::ofstream file_euler(filename_euler);

  data_out_euler.write_vtu(file_euler);
  data_out_euler.write_vtk(deallog.get_file_stream());


  triangulation.reset_manifold(0);
}

int
main()
{
  initlog();

  test<2, 3>(4, 3, 3);
  test<2, 3>(3, 4, 4);
  test<2, 3>(2, 5, 5);

  test<2, 3>(2, 3, 5);
  test<2, 3>(1, 4, 6);
  test<2, 3>(0, 5, 7);

  return 0;
}

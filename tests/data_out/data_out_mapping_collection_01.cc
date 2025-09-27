// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test DataOut::build_patch() that takes MappingCollection.


#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q_cache.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim, int spacedim = dim>
void
test()
{
  // create quarter hyper ball with 3 cells
  Triangulation<dim> triangulation;
  GridGenerator::quarter_hyper_ball(triangulation);

  // create first mapping class, that does NOT preserve position of vertices
  const FE_Q<dim, spacedim>     feq(1);
  const FESystem<dim, spacedim> fesystem(feq, spacedim);
  DoFHandler<dim, spacedim>     dhq(triangulation);
  dhq.distribute_dofs(fesystem);
  const ComponentMask mask(spacedim, true);
  Vector<double>      eulerq(dhq.n_dofs());
  VectorTools::get_position_vector(dhq, eulerq, mask);
  eulerq.add(-0.01);
  MappingFEField<dim, spacedim> mapping_1(dhq, eulerq, mask);

  // create first mapping class, that does preserve position of vertices
  MappingQ<dim, spacedim> mapping_2(1);

  // create mapping collection
  hp::FECollection<dim>      fe_collection(FE_Q<dim>(1), FE_Q<dim>(1));
  hp::MappingCollection<dim> mapping_collection(mapping_1, mapping_2);

  // create dof-handler and assign cells to different fes/manifolds
  DoFHandler<dim> dof_handler(triangulation);

  for (const auto &cell : dof_handler.active_cell_iterators())
    cell->set_active_fe_index(std::min(cell->active_cell_index(), 1u));

  dof_handler.distribute_dofs(fe_collection);

  // output vector
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  Vector<double> solution(dof_handler.n_dofs());

  data_out.add_data_vector(solution, "solution");

  data_out.build_patches(mapping_collection, 3);

#if false
  std::ofstream output("test.vtk");
  data_out.write_vtk(output);
#else
  data_out.write_vtk(deallog.get_file_stream());
#endif
}



int
main()
{
  initlog();
  deallog.get_file_stream().precision(2);

  test<2>();
}

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that DataOut::add_data_vector checks for consistent number of
// components in the names and in the finite element


#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"



void
test()
{
  const int          dim = 2;
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1, 1);
  FESystem<dim>   fe(FE_Q<dim>(1), dim + 1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  Vector<double> vec(dof_handler.n_dofs());
  vec = 1.;

  // "forget" to put a name on the last component -> should trigger assertion
  std::vector<std::string> names(dim, "u");

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    component_interpretation(
      dim, DataComponentInterpretation::component_is_part_of_vector);

  MappingQ<dim> mapping(2);

  // variant 1
  {
    DataOut<dim> data_out;
    try
      {
        data_out.add_data_vector(dof_handler,
                                 vec,
                                 names,
                                 component_interpretation);
        data_out.build_patches(mapping, 2);
      }
    catch (ExceptionBase &exc)
      {
        deallog << exc.get_exc_name() << std::endl;
      }
  }

  // variant 2
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    try
      {
        data_out.add_data_vector(vec,
                                 names,
                                 DataOut<dim>::type_automatic,
                                 component_interpretation);
        data_out.build_patches(mapping, 2);
      }
    catch (ExceptionBase &exc)
      {
        deallog << exc.get_exc_name() << std::endl;
      }
  }

  deallog << "OK" << std::endl;
}


int
main()
{
  deal_II_exceptions::disable_abort_on_exception();
  initlog();
  MultithreadInfo::set_thread_limit(1);
  test();
}

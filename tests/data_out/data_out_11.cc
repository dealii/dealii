// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// Check that DataOut::add_data_vector checks for consistent number of
// components in the names and in the finite element


#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q_generic.h>

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

  MappingQGeneric<dim> mapping(2);

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

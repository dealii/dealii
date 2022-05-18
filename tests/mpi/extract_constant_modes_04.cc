// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2018 by the deal.II authors
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



// test DoFTools::extract_constant_modes for multi-physical problems like in
// step-46

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"



template <int dim>
void
test(const unsigned int fe_degree)
{
  Triangulation<dim> triangulation;

  GridGenerator::hyper_cube(triangulation, 0, 1);
  triangulation.refine_global(1);

  DoFHandler<dim> dof_handler(triangulation);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->center()[0] < 0.5)
        {
          cell->set_active_fe_index(0);
        }
      else if (cell->center()[0] > 0.5)
        {
          cell->set_active_fe_index(1);
        }
      else
        {
          Assert(false, ExcNotImplemented());
        }
    }

  FESystem<dim>         left_fe(FE_Q<dim>(fe_degree), 1, FE_Nothing<dim>(), 1);
  FESystem<dim>         right_fe(FE_Nothing<dim>(), 1, FE_Q<dim>(fe_degree), 1);
  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(left_fe);
  fe_collection.push_back(right_fe);

  dof_handler.distribute_dofs(fe_collection);
  deallog << "FE degree=" << fe_degree << std::endl;
  deallog << "Total dofs=" << dof_handler.n_dofs() << std::endl;



  // extract constant modes and print
  for (unsigned int component = 0; component < fe_collection.n_components();
       ++component)
    {
      const FEValuesExtractors::Scalar component_extractor(component);
      const ComponentMask              component_mask(
        fe_collection.component_mask(component_extractor));

      std::vector<std::vector<bool>> constant_modes(
        fe_collection.n_components(),
        std::vector<bool>(dof_handler.n_locally_owned_dofs(), false));

      DoFTools::extract_constant_modes(dof_handler,
                                       component_mask,
                                       constant_modes);

      for (unsigned int i = 0; i < constant_modes.size(); ++i)
        {
          for (unsigned int j = 0; j < constant_modes[i].size(); ++j)
            {
              deallog << (constant_modes[i][j] ? '1' : '0') << ' ';
            }
          deallog << std::endl;
        }
    }



  // renumber dofs and do the same again
  DoFRenumbering::component_wise(dof_handler);
  for (unsigned int component = 0; component < fe_collection.n_components();
       ++component)
    {
      const FEValuesExtractors::Scalar component_extractor(component);
      const ComponentMask              component_mask(
        fe_collection.component_mask(component_extractor));

      std::vector<std::vector<bool>> constant_modes(
        fe_collection.n_components(),
        std::vector<bool>(dof_handler.n_locally_owned_dofs(), false));


      DoFTools::extract_constant_modes(dof_handler,
                                       component_mask,
                                       constant_modes);

      for (unsigned int i = 0; i < constant_modes.size(); ++i)
        {
          for (unsigned int j = 0; j < constant_modes[i].size(); ++j)
            {
              deallog << (constant_modes[i][j] ? '1' : '0') << ' ';
            }
          deallog << std::endl;
        }
    }
}


int
main(int argc, char *argv[])
{
  initlog();

  test<1>(1);
  test<1>(2);
  test<1>(3);
  test<1>(4);

  test<2>(1);
  test<2>(2);
  test<2>(3);
  test<2>(4);

  test<3>(1);
  test<3>(2);
  test<3>(3);
  test<3>(4);
}

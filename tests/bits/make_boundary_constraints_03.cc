// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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



// bug: https://code.google.com/p/dealii/issues/detail?id=154
// fixed in r31734

#include "../tests.h"
#include <deal.II/base/function_lib.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_values.h>
#include <fstream>
#include <iomanip>


template <int dim>
void test ()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation, -1, 1);

  FESystem<dim> fe(FE_Q<2>(1),1,FE_Q<2>(2),2);
  
  deallog << "Number of cells: "
          << triangulation.n_active_cells() << std::endl;

  DoFHandler<dim> dof_handler (triangulation);
  dof_handler.distribute_dofs (fe);
				   //DoFRenumbering::component_wise(dof_handler);
  deallog << "Number of dofs: "
          << dof_handler.n_dofs() << std::endl;

  ConstraintMatrix constraints;
  FEValuesExtractors::Vector velocities(1);
  ComponentMask mask = fe.component_mask (velocities);

  deallog << "ComponentMask " << mask[0] << mask[1] << mask[2] << std::endl;
  DoFTools::make_zero_boundary_constraints(dof_handler, constraints, mask);
  constraints.close();
  deallog << "Number of Constraints: " << constraints.n_constraints()<<std::endl;

  Quadrature<2> quadrature_formula(fe.get_unit_support_points()); 
  FEValues<2> fe_values (fe, quadrature_formula,
                         update_values | update_gradients | update_JxW_values|update_q_points);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  DoFHandler<2>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);

      std::vector<Point<2> > locations = fe_values.get_quadrature_points();
      cell->get_dof_indices (local_dof_indices);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
        if (constraints.is_constrained(local_dof_indices[i]))
         deallog << "DoF "<< local_dof_indices[i] 
                   << ", copy " << fe.system_to_base_index(i).first.second 
                   << ", base element " << fe.system_to_base_index(i).first.first 
                   << ", index " << fe.system_to_base_index(i).second 
                   << ", position " << locations[i] << std::endl;
      }
    }
}

int main ()
{
  initlog();

  try
    {
      test<2> ();
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}

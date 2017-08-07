// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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



// Check compatibility of Sacado and FEValues and FEFaceValues classes

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/sacado_product_type.h>

#include <fstream>

template<typename NumberType, int dim = 2>
void test ()
{
  FE_Q<dim>   fe (1);
  QGauss<dim> quadrature_formula(2);

  Triangulation<dim> tria;
  DoFHandler<dim>    dof_handler (tria);
  Vector<double>     solution;

  GridGenerator::hyper_cube (tria, -1, 1);
  dof_handler.distribute_dofs (fe);
  solution.reinit (dof_handler.n_dofs());

  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values);
  std::vector<types::global_dof_index> local_dof_indices (fe.dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
  const typename DoFHandler<dim>::active_cell_iterator endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      cell->get_dof_indices (local_dof_indices);

      std::vector<NumberType> local_dof_values(fe.dofs_per_cell);
      std::vector<NumberType> local_values_1(quadrature_formula.size());
      std::vector<double> local_values_2(quadrature_formula.size());

      cell->get_dof_values(solution,local_dof_values.begin(), local_dof_values.end());
      // for(unsigned int i =0; i<fe.dofs_per_cell; ++i) local_values.diff(i,fe.dofs_per_cell);

      fe_values.get_function_values_from_local_vector(local_dof_values, local_values_1);
      fe_values.get_function_values(solution, local_values_2);
      for(unsigned int q=0; q<quadrature_formula.size(); ++q)
          if(NumberType(local_values_2[q]) != local_values_1[q])
              deallog << "NOT OK: " << q << std::endl;
    }


  deallog << "OK" << std::endl;
}

int main (int argc, char **argv)
{
  initlog();
  deallog.threshold_double(1.e-10);

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, testing_max_num_threads());

  deallog.push("Double");
  {
    test<double>();
  }
  deallog.pop();

  deallog.push("Sacado<double>");
  {
    test<Sacado::Fad::DFad<double>>();
  }
  deallog.pop();

  deallog << "OK" << std::endl;
}

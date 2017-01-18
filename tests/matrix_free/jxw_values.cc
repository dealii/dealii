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



// this function tests the correctness of JxW values returned by FEEvaluation
// when compared to FEValues

#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>

#include "create_mesh.h"

std::ofstream logfile("output");


template <int dim>
void test ()
{
  Triangulation<dim> tria;
  create_mesh (tria);
  tria.refine_global(4-dim);

  // refine a few cells
  for (unsigned int i=0; i<10-3*dim; ++i)
    {
      typename Triangulation<dim>::active_cell_iterator
      cell = tria.begin_active (),
      endc = tria.end();
      unsigned int counter = 0;
      for (; cell!=endc; ++cell, ++counter)
        if (counter % (7-i) == 0)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  FE_Q<dim> fe (1);
  DoFHandler<dim> dof (tria);
  dof.distribute_dofs(fe);

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof, constraints);
  constraints.close();

  MatrixFree<dim> mf_data;
  {
    const QGauss<1> quad (2);
    typename MatrixFree<dim>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim>::AdditionalData::none;
    data.mapping_update_flags = update_JxW_values;
    mf_data.reinit (dof, constraints, quad, data);
  }

  double error = 0, error2 = 0, abs = 0;

  QGauss<dim> quad(2);
  FEValues<dim> fe_values(fe, quad, update_JxW_values);
  FEEvaluation<dim,1> fe_eval(mf_data);
  AlignedVector<VectorizedArray<double> > jxw_values_manual(fe_eval.n_q_points);
  for (unsigned int cell=0; cell<mf_data.n_macro_cells(); ++cell)
    {
      fe_eval.reinit(cell);
      fe_eval.fill_JxW_values(jxw_values_manual);
      for (unsigned int v=0; v<mf_data.n_components_filled(cell); ++v)
        {
          fe_values.reinit(mf_data.get_cell_iterator(cell,v));
          for (unsigned int q=0; q<quad.size(); ++q)
            {
              abs += fe_values.JxW(q);
              error += std::abs(fe_values.JxW(q) - fe_eval.JxW(q)[v]);
              error2 += std::abs(fe_values.JxW(q) - jxw_values_manual[q][v]);
            }
        }
    }

  deallog << "Norm of difference: " << error/abs << " " << error2/abs
          << std::endl << std::endl;
}



int main ()
{
  deallog.attach(logfile);
  deallog << std::setprecision (3);

  deallog.threshold_double(1.e-12);
  test<2>();
  test<3>();
}

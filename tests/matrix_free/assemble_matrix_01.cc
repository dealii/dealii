// ---------------------------------------------------------------------
//
// Copyright (C) 2014 by the deal.II authors
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



// test FEEvaluation for assembling the Laplace matrix. It is enough to just
// consider the resulting element matrices element by element

#include "../tests.h"
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>


std::ofstream logfile("output");



template <int dim, int fe_degree>
void do_test (const DoFHandler<dim> &dof)
{

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;


  MappingQ<dim> mapping(fe_degree+1);

  // compute matrix with (\nabla v, \nabla u) + (v, 10 * u)
  {
    QGauss<dim>  quadrature_formula(fe_degree+1);

    FEValues<dim> fe_values (mapping, dof.get_fe(), quadrature_formula,
                             update_values    |  update_gradients |
                             update_JxW_values);

    FEEvaluation<dim,fe_degree,fe_degree+1>
      fe_eval (mapping, dof.get_fe(), QGauss<1>(fe_degree+1),
               update_values | update_gradients | update_JxW_values);

    const unsigned int   dofs_per_cell = dof.get_fe().dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    FullMatrix<double>   test_matrix (dofs_per_cell, dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof.begin_active(),
    endc = dof.end();
    for (; cell!=endc; ++cell)
      {
        cell_matrix = 0;
        test_matrix = 0;
        fe_values.reinit (cell);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                cell_matrix(i,j) += ((fe_values.shape_grad(i,q_point) *
                                      fe_values.shape_grad(j,q_point)
                                      +
                                      10. *
                                      fe_values.shape_value(i,q_point) *
                                      fe_values.shape_value(j,q_point)) *
                                     fe_values.JxW(q_point));
            }

        fe_eval.reinit(cell);
        for (unsigned int i=0; i<dofs_per_cell;
             i += VectorizedArray<double>::n_array_elements)
          {
            const unsigned int n_items =
              i+VectorizedArray<double>::n_array_elements > dofs_per_cell ?
              (dofs_per_cell - i) : VectorizedArray<double>::n_array_elements;
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              fe_eval.begin_dof_values()[j]  = VectorizedArray<double>();
            for (unsigned int v=0; v<n_items; ++v)
              fe_eval.begin_dof_values()[i+v][v] = 1.;

            fe_eval.evaluate(true, true);
            for (unsigned int q=0; q<n_q_points; ++q)
              {
                fe_eval.submit_value(10.*fe_eval.get_value(q), q);
                fe_eval.submit_gradient(fe_eval.get_gradient(q), q);
              }
            fe_eval.integrate(true, true);

            for (unsigned int v=0; v<n_items; ++v)
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                test_matrix(fe_eval.get_internal_dof_numbering()[j],
                            fe_eval.get_internal_dof_numbering()[i+v])
                  = fe_eval.begin_dof_values()[j][v];
          }
        test_matrix.add(-1., cell_matrix);
        deallog << test_matrix.frobenius_norm() << " ";
      }
    deallog << std::endl;
  }
}



template <int dim, int fe_degree>
void test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_ball (tria);
  static const HyperBallBoundary<dim> boundary;
  tria.set_boundary (0, boundary);
  if (dim < 3 || fe_degree < 2)
    tria.refine_global(1);
  tria.begin(tria.n_levels()-1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FE_Q<dim> fe (fe_degree);
  DoFHandler<dim> dof (tria);
  dof.distribute_dofs(fe);
  do_test<dim, fe_degree> (dof);
}



int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  deallog << std::setprecision (3);

  {
    deallog.threshold_double(5.e-11);
    deallog.push("2d");
    test<2,1>();
    test<2,2>();
    test<2,4>();
    deallog.pop();
    deallog.push("3d");
    test<3,1>();
    test<3,2>();
    deallog.pop();
  }
}

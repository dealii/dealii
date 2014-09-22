// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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



// this function tests the correctness of the implementation of matrix free
// operations in integrating functions and gradients on a hypeball mesh with
// adaptive refinement.

#include "../tests.h"

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>

std::ofstream logfile("output");


template <int dim, int fe_degree, typename Number>
class MatrixFreeTest
{
public:
  typedef std::vector<Vector<Number>*> VectorType;

  MatrixFreeTest(const MatrixFree<dim,Number> &data_in):
    data   (data_in),
    fe_val (data.get_dof_handler().get_fe(),
            Quadrature<dim>(data.get_quadrature(0)),
            update_values | update_gradients | update_JxW_values)
  {};

  void operator () (const MatrixFree<dim,Number> &data,
                    VectorType       &dst,
                    const VectorType &src,
                    const std::pair<unsigned int,unsigned int> &cell_range) const;

  void test_functions (Vector<Number> &dst,
                       Vector<Number> &dst_deal) const
  {
    dst = 0;
    dst_deal = 0;
    VectorType dst_data (2);
    dst_data[0] = &dst;
    dst_data[1] = &dst_deal;
    VectorType src_dummy;
    data.cell_loop (&MatrixFreeTest<dim,fe_degree,Number>::operator(), this,
                    dst_data, src_dummy);
  };

private:
  const MatrixFree<dim,Number> &data;
  mutable FEValues<dim> fe_val;
};




template <int dim, int fe_degree, typename Number>
void MatrixFreeTest<dim,fe_degree,Number>::
operator () (const MatrixFree<dim,Number> &data,
             std::vector<Vector<Number>*> &dst,
             const std::vector<Vector<Number>*> &,
             const std::pair<unsigned int,unsigned int> &cell_range) const
{
  FEEvaluation<dim,fe_degree,fe_degree+1,1,Number> fe_eval (data);
  const unsigned int n_q_points = fe_eval.n_q_points;
  const unsigned int dofs_per_cell = fe_eval.dofs_per_cell;
  AlignedVector<VectorizedArray<Number> > values (n_q_points);
  AlignedVector<VectorizedArray<Number> > gradients (dim*n_q_points);
  std::vector<types::global_dof_index> dof_indices (dofs_per_cell);
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
    {
      fe_eval.reinit(cell);
      // compare values with the ones the FEValues
      // gives us. Those are seen as reference
      for (unsigned int j=0; j<data.n_components_filled(cell); ++j)
        {
          // generate random numbers at quadrature
          // points and test them with basis functions
          // and their gradients
          for (unsigned int q=0; q<n_q_points; ++q)
            {
              values[q][j] = Testing::rand()/(double)RAND_MAX;
              for (unsigned int d=0; d<dim; ++d)
                gradients[q*dim+d][j] = -1. + 2. * (Testing::rand()/(double)RAND_MAX);
            }
          fe_val.reinit (data.get_cell_iterator(cell,j));
          data.get_cell_iterator(cell,j)->get_dof_indices(dof_indices);

          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              double sum = 0.;
              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  sum += values[q][j] * fe_val.shape_value(i,q) * fe_val.JxW(q);
                  for (unsigned int d=0; d<dim; ++d)
                    sum += (gradients[q*dim+d][j] * fe_val.shape_grad(i,q)[d] *
                            fe_val.JxW(q));
                }
              (*dst[1])(dof_indices[i]) += sum;
            }
        }
      for (unsigned int q=0; q<n_q_points; ++q)
        {
          fe_eval.submit_value (values[q], q);
          Tensor<1,dim,VectorizedArray<Number> > submit (false);
          for (unsigned int d=0; d<dim; ++d)
            submit[d] = gradients[q*dim+d];
          fe_eval.submit_gradient (submit, q);
        }
      fe_eval.integrate (true,true);
      fe_eval.distribute_local_to_global (*dst[0]);
    }
}



template <int dim, int fe_degree>
void test ()
{
  typedef double number;
  Triangulation<dim> tria;
  GridGenerator::hyper_ball (tria);
  static const HyperBallBoundary<dim> boundary;
  tria.set_boundary (0, boundary);
  typename Triangulation<dim>::active_cell_iterator
  cell = tria.begin_active (),
  endc = tria.end();
  for (; cell!=endc; ++cell)
    if (cell->center().norm()<1e-8)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  cell = tria.begin_active ();
  for (; cell!=endc; ++cell)
    if (cell->center().norm()<0.2)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  if (dim < 3 || fe_degree < 2)
    tria.refine_global(1);
  tria.begin(tria.n_levels()-1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  cell = tria.begin_active ();
  for (unsigned int i=0; i<7-2*dim; ++i)
    {
      cell = tria.begin_active ();
      unsigned int counter = 0;
      for (; cell!=endc; ++cell, ++counter)
        if (counter % (7-i) == 0)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  FE_Q<dim> fe (fe_degree);
  DoFHandler<dim> dof (tria);
  dof.distribute_dofs(fe);
  deallog << "Testing " << fe.get_name() << std::endl;
  //std::cout << "Number of cells: " << tria.n_active_cells() << std::endl;
  //std::cout << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();

  MatrixFree<dim,number> mf_data;
  {
    const QGauss<1> quad (fe_degree+1);
    mf_data.reinit (dof, constraints, quad,
                    typename MatrixFree<dim,number>::AdditionalData(MPI_COMM_SELF,MatrixFree<dim,number>::AdditionalData::none));
  }

  MatrixFreeTest<dim,fe_degree,number> mf (mf_data);
  Vector<number> solution (dof.n_dofs());
  Vector<number> solution_dist (dof.n_dofs());

  mf.test_functions(solution_dist, solution);

  constraints.condense (solution);

  Vector<number> compare (solution_dist);
  compare -= solution;
  const double diff_norm = compare.linfty_norm();

  deallog << "Norm of difference: " << diff_norm << std::endl << std::endl;
}


int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog << std::setprecision (3);

  {
    deallog.threshold_double(1.e-12);
    deallog.push("2d");
    test<2,1>();
    test<2,2>();
    test<2,3>();
    test<2,4>();
    deallog.pop();
    deallog.push("3d");
    test<3,1>();
    test<3,2>();
    test<3,3>();
    deallog.pop();
  }
}


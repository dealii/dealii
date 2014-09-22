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
// adaptive refinement. It uses two different finite elements inside one
// matrix_free structure.

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
  typedef std::vector<Vector<Number> > VectorType;

  MatrixFreeTest(const MatrixFree<dim,Number> &data_in):
    data   (data_in),
    fe_val0 (data.get_dof_handler(0).get_fe(),
             Quadrature<dim>(data.get_quadrature(0)),
             update_values | update_gradients | update_JxW_values),
    fe_val01 (data.get_dof_handler(0).get_fe(),
              Quadrature<dim>(data.get_quadrature(1)),
              update_values | update_gradients | update_JxW_values),
    fe_val1 (data.get_dof_handler(1).get_fe(),
             Quadrature<dim>(data.get_quadrature(1)),
             update_values | update_gradients | update_JxW_values)
  {};

  void operator () (const MatrixFree<dim,Number> &data,
                    VectorType       &dst,
                    const VectorType &src,
                    const std::pair<unsigned int,unsigned int> &cell_range) const;

  void test_functions (VectorType &dst) const
  {
    for (unsigned int comp=0; comp<dst.size(); ++comp)
      dst[comp] = 0;
    VectorType src_dummy;
    data.cell_loop (&MatrixFreeTest<dim,fe_degree,Number>::operator(), this,
                    dst, src_dummy);
  };

private:
  const MatrixFree<dim,Number> &data;
  mutable FEValues<dim> fe_val0;
  mutable FEValues<dim> fe_val01;
  mutable FEValues<dim> fe_val1;
};




template <int dim, int fe_degree, typename Number>
void MatrixFreeTest<dim,fe_degree,Number>::
operator () (const MatrixFree<dim,Number> &data,
             std::vector<Vector<Number> > &dst,
             const std::vector<Vector<Number> > &,
             const std::pair<unsigned int,unsigned int> &cell_range) const
{
  FEEvaluation<dim,fe_degree,fe_degree+1,1,Number>   fe_eval0 (data, 0, 0);
  FEEvaluation<dim,fe_degree+1,fe_degree+2,1,Number> fe_eval1 (data, 1, 1);
  FEEvaluation<dim,fe_degree,fe_degree+2,1,Number>   fe_eval01 (data, 0, 1);
  const unsigned int n_q_points0 = fe_eval0.n_q_points;
  const unsigned int n_q_points1 = fe_eval1.n_q_points;
  const unsigned int dofs_per_cell0 = fe_eval0.dofs_per_cell;
  const unsigned int dofs_per_cell1 = fe_eval1.dofs_per_cell;
  AlignedVector<VectorizedArray<Number> > values0 (n_q_points0);
  AlignedVector<VectorizedArray<Number> > gradients0 (dim*n_q_points0);
  AlignedVector<VectorizedArray<Number> > values1 (n_q_points1);
  AlignedVector<VectorizedArray<Number> > gradients1 (dim*n_q_points1);
  std::vector<types::global_dof_index> dof_indices0 (dofs_per_cell0);
  std::vector<types::global_dof_index> dof_indices1 (dofs_per_cell1);
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
    {
      fe_eval0.reinit(cell);
      fe_eval1.reinit(cell);
      fe_eval01.reinit(cell);

      // compare values with the ones the FEValues
      // gives us. Those are seen as reference
      for (unsigned int j=0; j<data.n_components_filled(cell); ++j)
        {
          // FE 0, Quad 0
          // generate random numbers at quadrature
          // points and test them with basis functions
          // and their gradients
          for (unsigned int q=0; q<n_q_points0; ++q)
            {
              values0[q][j] = Testing::rand()/(double)RAND_MAX;
              for (unsigned int d=0; d<dim; ++d)
                gradients0[q*dim+d][j] = -1. + 2. * (Testing::rand()/(double)RAND_MAX);
            }
          fe_val0.reinit (data.get_cell_iterator(cell,j,0));
          data.get_cell_iterator(cell,j,0)->get_dof_indices(dof_indices0);

          for (unsigned int i=0; i<dofs_per_cell0; ++i)
            {
              double sum = 0.;
              for (unsigned int q=0; q<n_q_points0; ++q)
                {
                  sum += values0[q][j] * fe_val0.shape_value(i,q) * fe_val0.JxW(q);
                  for (unsigned int d=0; d<dim; ++d)
                    sum += (gradients0[q*dim+d][j] * fe_val0.shape_grad(i,q)[d] *
                            fe_val0.JxW(q));
                }
              dst[0+1](dof_indices0[i]) += sum;
            }

          // FE 1, Quad 1
          fe_val1.reinit (data.get_cell_iterator(cell,j,1));
          data.get_cell_iterator(cell,j,1)->get_dof_indices(dof_indices1);

          for (unsigned int q=0; q<n_q_points1; ++q)
            {
              values1[q][j] = Testing::rand()/(double)RAND_MAX;
              for (unsigned int d=0; d<dim; ++d)
                gradients1[q*dim+d][j] = -1. + 2. * (Testing::rand()/(double)RAND_MAX);
            }
          for (unsigned int i=0; i<dofs_per_cell1; ++i)
            {
              double sum = 0.;
              for (unsigned int q=0; q<n_q_points1; ++q)
                {
                  sum += values1[q][j] * fe_val1.shape_value(i,q) * fe_val1.JxW(q);
                  for (unsigned int d=0; d<dim; ++d)
                    sum += (gradients1[q*dim+d][j] * fe_val1.shape_grad(i,q)[d] *
                            fe_val1.JxW(q));
                }
              dst[2+1](dof_indices1[i]) += sum;
            }

          // FE 0, Quad 1
          fe_val01.reinit (data.get_cell_iterator(cell,j,0));
          for (unsigned int i=0; i<dofs_per_cell0; ++i)
            {
              double sum = 0.;
              for (unsigned int q=0; q<n_q_points1; ++q)
                {
                  sum += values1[q][j] * fe_val01.shape_value(i,q) * fe_val01.JxW(q);
                  for (unsigned int d=0; d<dim; ++d)
                    sum += (gradients1[q*dim+d][j] * fe_val01.shape_grad(i,q)[d] *
                            fe_val01.JxW(q));
                }
              dst[4+1](dof_indices0[i]) += sum;
            }
        }

      // FE 0, Quad 0
      for (unsigned int q=0; q<n_q_points0; ++q)
        {
          fe_eval0.submit_value (values0[q], q);
          Tensor<1,dim,VectorizedArray<Number> > submit (false);
          for (unsigned int d=0; d<dim; ++d)
            submit[d] = gradients0[q*dim+d];
          fe_eval0.submit_gradient (submit, q);
        }
      fe_eval0.integrate (true,true);
      fe_eval0.distribute_local_to_global (dst[0]);

      // FE 1, Quad 1
      for (unsigned int q=0; q<n_q_points1; ++q)
        {
          fe_eval1.submit_value (values1[q], q);
          Tensor<1,dim,VectorizedArray<Number> > submit (false);
          for (unsigned int d=0; d<dim; ++d)
            submit[d] = gradients1[q*dim+d];
          fe_eval1.submit_gradient (submit, q);
        }
      fe_eval1.integrate (true,true);
      fe_eval1.distribute_local_to_global (dst[2]);

      // FE 0, Quad 1
      for (unsigned int q=0; q<n_q_points1; ++q)
        {
          fe_eval01.submit_value (values1[q], q);
          Tensor<1,dim,VectorizedArray<Number> > submit (false);
          for (unsigned int d=0; d<dim; ++d)
            submit[d] = gradients1[q*dim+d];
          fe_eval01.submit_gradient (submit, q);
        }
      fe_eval01.integrate (true,true);
      fe_eval01.distribute_local_to_global (dst[4]);
    }
}



template <int dim, int fe_degree, typename number>
void test ()
{
  // create hyper ball geometry and refine some
  // cells
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

  FE_Q<dim> fe0(fe_degree);
  FE_Q<dim> fe1(fe_degree+1);
  DoFHandler<dim> dof0(tria);
  dof0.distribute_dofs(fe0);
  DoFHandler<dim> dof1(tria);
  dof1.distribute_dofs(fe1);

  std::vector<const DoFHandler<dim> *>dof(2);
  dof[0] = &dof0;
  dof[1] = &dof1;

  deallog << "Testing " << fe0.get_name() << " and " << fe1.get_name() << std::endl;
  //std::cout << "Number of cells: " << tria.n_active_cells() << std::endl;

  std::vector<Vector<number> > dst (6);
  dst[0].reinit (dof[0]->n_dofs());
  dst[1].reinit (dst[0]);
  dst[2].reinit (dof[1]->n_dofs());
  dst[3].reinit (dst[2]);
  dst[4].reinit (dst[0]);
  dst[5].reinit (dst[0]);

  std::vector<const ConstraintMatrix *> constraints(2);
  ConstraintMatrix constraint0;
  DoFTools::make_hanging_node_constraints(*dof[0],constraint0);
  constraint0.close();
  constraints[0] = &constraint0;
  ConstraintMatrix constraint1;
  DoFTools::make_hanging_node_constraints(*dof[1],constraint1);
  constraint1.close();
  constraints[1] = &constraint1;

  //std::cout << "Number of degrees of freedom FE 0: " << dof[0]->n_dofs() << std::endl;
  //std::cout << "Number of constraints FE 0: " << constraints[0]->n_constraints() << std::endl;
  //std::cout << "Number of degrees of freedom FE 1: " << dof[1]->n_dofs() << std::endl;
  //std::cout << "Number of constraints FE 1: " << constraints[1]->n_constraints() << std::endl;

  MatrixFree<dim,number> mf_data;
  {
    std::vector<Quadrature<1> > quad;
    for (unsigned int no=0; no<2; ++no)
      quad.push_back(QGauss<1>(fe_degree+1+no));
    mf_data.reinit (dof, constraints, quad,
                    typename MatrixFree<dim,number>::AdditionalData(MPI_COMM_SELF,MatrixFree<dim,number>::AdditionalData::none));
  }

  MatrixFreeTest<dim,fe_degree,number> mf (mf_data);
  mf.test_functions(dst);

  constraints[0]->condense(dst[1]);
  constraints[1]->condense(dst[3]);
  constraints[0]->condense(dst[5]);

  dst[1] -= dst[0];
  double diff_norm = dst[1].linfty_norm();
  deallog << "FE 0, Quad 0; integration difference: " << diff_norm << std::endl;

  dst[3] -= dst[2];
  diff_norm = dst[3].linfty_norm();
  deallog << "FE 1, Quad 1; integration difference: " << diff_norm << std::endl;

  dst[5] -= dst[4];
  diff_norm = dst[5].linfty_norm();
  deallog << "FE 0, Quad 1; integration difference: " << diff_norm << std::endl << std::endl;
}


int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog << std::setprecision (3);

  {
    deallog.threshold_double(1.e-12);
    deallog.push("2d");
    test<2,1,double>();
    test<2,2,double>();
    test<2,3,double>();
    deallog.pop();
    deallog.push("3d");
    test<3,1,double>();
    test<3,2,double>();
    deallog.pop();
  }

  {
    deallog << std::endl << "Test with floats" << std::endl << std::endl;
    deallog.threshold_double(1.e-7);
    deallog.push("2d");
    test<2,1,float>();
    deallog.pop();
    deallog.push("3d");
    test<3,1,float>();
    deallog.pop();
  }
}


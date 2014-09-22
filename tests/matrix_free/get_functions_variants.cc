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
// operations in getting the function values, the function gradients, and the
// function Laplacians on a cartesian mesh (hyper cube) with the different
// variants (functions only, gradients only, some combinations). The reference
// is computed to the function that extracts all data at the same time (the
// correctness of that must be tested elsewhere). No constraints.

#include "../tests.h"

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <iostream>

std::ofstream logfile("output");


template <int dim, int fe_degree, typename Number>
class MatrixFreeTest
{
public:
  typedef Vector<Number> VectorType;

  MatrixFreeTest(const MatrixFree<dim,Number> &data_in):
    data (data_in)
  {};

  void operator () (const MatrixFree<dim,Number> &data,
                    VectorType       &dst,
                    const VectorType &src,
                    const std::pair<unsigned int,unsigned int> &cell_range) const;

  void test_functions (const VectorType &src) const
  {
    for (unsigned int i=0; i<5; ++i)
      errors[i] = 0;
    data.cell_loop (&MatrixFreeTest<dim,fe_degree,Number>::operator(), this,
                    const_cast<VectorType &>(src), src);

    deallog << "Error val, function values alone: "
            << errors[0] << std::endl;
    deallog << "Error grad, function gradients alone: "
            << errors[1] << std::endl;
    deallog << "Error val, function values and gradients alone: "
            << errors[2] << std::endl;
    deallog << "Error grad, function values and gradients alone: "
            << errors[3] << std::endl;
    deallog << "Error Lapl, function Laplacians alone: "
            << errors[4] << std::endl;
  };

private:
  const MatrixFree<dim,Number> &data;
  mutable double errors [5];
};




template <int dim, int fe_degree, typename Number>
void MatrixFreeTest<dim,fe_degree,Number>::
operator () (const MatrixFree<dim,Number> &data,
             VectorType &,
             const VectorType &src,
             const std::pair<unsigned int,unsigned int> &cell_range) const
{
  FEEvaluation<dim,fe_degree,fe_degree+1,1,Number> fe_eval (data);
  FEEvaluation<dim,fe_degree,fe_degree+1,1,Number> fe_eval2 (data);
  FEEvaluation<dim,fe_degree,fe_degree+1,1,Number> fe_eval3 (data);
  FEEvaluation<dim,fe_degree,fe_degree+1,1,Number> fe_eval4 (data);
  FEEvaluation<dim,fe_degree,fe_degree+1,1,Number> fe_eval5 (data);
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
    {
      fe_eval.reinit (cell);
      fe_eval.read_dof_values(src);
      fe_eval.evaluate (true,true,true);

      // only for values (additional test)
      fe_eval2.reinit (cell);
      fe_eval2.read_dof_values(src);
      fe_eval2.evaluate (true,false,false);

      // only gradients
      fe_eval3.reinit (cell);
      fe_eval3.read_dof_values(src);
      fe_eval3.evaluate (false,true,false);

      // only values and gradients
      fe_eval4.reinit (cell);
      fe_eval4.read_dof_values(src);
      fe_eval4.evaluate(true,true,false);

      // only laplacians
      fe_eval5.reinit (cell);
      fe_eval5.read_dof_values(src);
      fe_eval5.evaluate (false,false,true);


      // compare values with the values that we get
      // when expanding the full
      // FEEvaluations. Those are tested in other
      // functions and seen as reference here
      for (unsigned int q=0; q<fe_eval.n_q_points; ++q)
        for (unsigned int j=0; j<VectorizedArray<Number>::n_array_elements; ++j)
          {
            errors[0] += std::fabs(fe_eval.get_value(q)[j]-
                                   fe_eval2.get_value(q)[j]);
            errors[2] += std::fabs(fe_eval.get_value(q)[j]-
                                   fe_eval4.get_value(q)[j]);
            for (unsigned int d=0; d<dim; ++d)
              {
                errors[1] += std::fabs(fe_eval.get_gradient(q)[d][j]-
                                       fe_eval3.get_gradient(q)[d][j]);
                errors[3] += std::fabs(fe_eval.get_gradient(q)[d][j]-
                                       fe_eval4.get_gradient(q)[d][j]);
              }
            errors[4] += std::fabs(fe_eval.get_laplacian(q)[j]-
                                   fe_eval5.get_laplacian(q)[j]);
          }
    }
}


template <int dim, int fe_degree>
void test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria);
  tria.refine_global(1);

  FE_Q<dim> fe (fe_degree);
  DoFHandler<dim> dof (tria);
  dof.distribute_dofs(fe);

  deallog << "Testing " << fe.get_name() << std::endl;
  //std::cout << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;

  Vector<double> solution_dist (dof.n_dofs());

  // create vector with random entries
  for (unsigned int i=0; i<dof.n_dofs(); ++i)
    {
      const double entry = Testing::rand()/(double)RAND_MAX;
      solution_dist(i) = entry;
    }

  ConstraintMatrix constraints;
  MatrixFree<dim> mf_data;
  {
    const QGauss<1> quad (fe_degree+1);
    typename MatrixFree<dim>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim>::AdditionalData::none;
    data.mapping_update_flags = update_gradients | update_second_derivatives;
    mf_data.reinit (dof, constraints, quad, data);
  }

  MatrixFreeTest<dim,fe_degree,double> mf (mf_data);
  mf.test_functions(solution_dist);
  deallog << std::endl;
}


int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  deallog.threshold_double(1.e-14);
  deallog << std::setprecision (3);

  {
    deallog.push("2d");
    test<2,1>();
    test<2,2>();
    test<2,3>();
    test<2,4>();
    test<2,5>();
    deallog.pop();
    deallog.push("3d");
    test<3,1>();
    test<3,2>();
    test<3,3>();
    test<3,4>();
    deallog.pop();
  }
}


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



// check that FEEvaluation can be evaluated without indices initialized (and
// throws an exception when trying to read/write from/to vectors)

#include "../tests.h"


std::ofstream logfile("output");

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

#include <fstream>
#include <iostream>


// forward declare this function
template <int dim, int fe_degree>
void test ();



template <int dim, int fe_degree, int n_q_points_1d=fe_degree+1, typename Number=double>
class MatrixFreeTest
{
public:
  bool read_vector;

  MatrixFreeTest(const MatrixFree<dim,Number> &data_in):
    read_vector(false),
    data       (data_in)
  {};

  void
  operator () (const MatrixFree<dim,Number> &data,
               Vector<Number> &,
               const Vector<Number> &src,
               const std::pair<unsigned int,unsigned int> &cell_range) const
  {
    FEEvaluation<dim,fe_degree,n_q_points_1d,1,Number> fe_eval (data);

    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        fe_eval.reinit (cell);

        // do not read from vector as that is disallowed
        if (read_vector)
          fe_eval.read_dof_values(src);
        else
          for (unsigned int e=0; e<fe_eval.dofs_per_cell; ++e)
            fe_eval.submit_dof_value(make_vectorized_array<Number>(1.), e);
        fe_eval.evaluate (true,true,true);

        // values should evaluate to one, derivatives to zero
        for (unsigned int q=0; q<fe_eval.n_q_points; ++q)
          for (unsigned int d=0; d<VectorizedArray<Number>::n_array_elements; ++d)
            {
              Assert(fe_eval.get_value(q)[d] == 1., ExcInternalError());
              for (unsigned int e=0; e<dim; ++e)
                Assert(fe_eval.get_gradient(q)[e][d] == 0., ExcInternalError());
              Assert(fe_eval.get_laplacian(q)[d] == 0., ExcInternalError());
            }
      }
  }



  void test_functions (const Vector<Number> &src) const
  {
    Vector<Number> dst_dummy;
    data.cell_loop (&MatrixFreeTest::operator(), this, dst_dummy, src);
    deallog << "OK" << std::endl;
  };

protected:
  const MatrixFree<dim,Number> &data;
};




template <int dim, int fe_degree, typename number>
void do_test (const DoFHandler<dim> &dof,
              const ConstraintMatrix &constraints)
{
  // use this for info on problem
  //std::cout << "Number of cells: " << dof.get_tria().n_active_cells()
  //          << std::endl;
  //std::cout << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;
  //std::cout << "Number of constraints: " << constraints.n_constraints() << std::endl;

  Vector<number> solution (dof.n_dofs());

  // create vector with random entries
  for (unsigned int i=0; i<dof.n_dofs(); ++i)
    {
      if (constraints.is_constrained(i))
        continue;
      const double entry = Testing::rand()/(double)RAND_MAX;
      solution(i) = entry;
    }

  constraints.distribute(solution);
  MatrixFree<dim,number> mf_data;
  {
    const QGauss<1> quad (fe_degree+1);
    typename MatrixFree<dim,number>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim,number>::AdditionalData::none;
    data.mapping_update_flags = update_gradients | update_second_derivatives;
    data.initialize_indices = false;
    mf_data.reinit (dof, constraints, quad, data);
  }

  deallog << "Testing " << dof.get_fe().get_name() << " without read" << std::endl;
  MatrixFreeTest<dim,fe_degree,fe_degree+1,number> mf (mf_data);
  mf.test_functions(solution);

  deallog << "Testing " << dof.get_fe().get_name() << " with read" << std::endl;
  try
    {
      mf.read_vector = true;
      mf.test_functions(solution);
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }
}


int main ()
{
  deal_II_exceptions::disable_abort_on_exception();
  deallog.attach(logfile);
  deallog.depth_console(0);

  deallog << std::setprecision (3);
  test<2,1>();
}


template <int dim, int fe_degree>
void test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_ball (tria);
  static const HyperBallBoundary<dim> boundary;
  tria.set_boundary (0, boundary);
  // refine first and last cell
  tria.begin(tria.n_levels()-1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.refine_global (4-dim);

  FE_Q<dim> fe (fe_degree);
  DoFHandler<dim> dof (tria);
  dof.distribute_dofs(fe);

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof, constraints);
  constraints.close();

  do_test <dim, fe_degree, double> (dof, constraints);
}

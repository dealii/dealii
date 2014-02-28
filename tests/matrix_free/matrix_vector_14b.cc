// ---------------------------------------------------------------------
// $Id$
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



// Same as matrix_vector_14, but using more quadrature points than the FE
// degree would suggest

#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/lac/vector.h>


std::ofstream logfile("output");



template <int dim, int fe_degree, typename VECTOR>
void
helmholtz_operator_dgp (const MatrixFree<dim,typename VECTOR::value_type>  &data,
                        VECTOR       &dst,
                        const VECTOR &src,
                        const std::pair<unsigned int,unsigned int> &cell_range)
{
  typedef typename VECTOR::value_type Number;
  FEEvaluationDGP<dim,fe_degree,fe_degree+2,1,Number> fe_eval (data);
  const unsigned int n_q_points = fe_eval.n_q_points;

  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
    {
      fe_eval.reinit (cell);
      fe_eval.read_dof_values (src);
      fe_eval.evaluate (true, true, false);
      for (unsigned int q=0; q<n_q_points; ++q)
        {
          fe_eval.submit_value (Number(10)*fe_eval.get_value(q),q);
          fe_eval.submit_gradient (fe_eval.get_gradient(q),q);
        }
      fe_eval.integrate (true,true);
      fe_eval.distribute_local_to_global (dst);
    }
}



template <int dim, int fe_degree, typename Number, typename VECTOR=Vector<Number> >
class MatrixFreeTest
{
public:
  typedef VectorizedArray<Number> vector_t;

  MatrixFreeTest(const MatrixFree<dim,Number> &data_in):
    data (data_in)
  {};

  void vmult (VECTOR       &dst,
              const VECTOR &src) const
  {
    dst = 0;
    const std_cxx1x::function<void(const MatrixFree<dim,typename VECTOR::value_type> &,
                                   VECTOR &,
                                   const VECTOR &,
                                   const std::pair<unsigned int,unsigned int> &)>
      wrap = helmholtz_operator_dgp<dim,fe_degree,VECTOR>;
    data.cell_loop (wrap, dst, src);
  };

private:
  const MatrixFree<dim,Number> &data;
};



template <int dim, int fe_degree, typename number>
void do_test (const DoFHandler<dim> &dof,
              const ConstraintMatrix &constraints)
{

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;

  MatrixFree<dim,number> mf_data;
  {
    const QGauss<1> quad (fe_degree+2);
    typename MatrixFree<dim,number>::AdditionalData data;
    data.tasks_parallel_scheme =
      MatrixFree<dim,number>::AdditionalData::partition_color;
    data.tasks_block_size = 3;

    mf_data.reinit (dof, constraints, quad, data);
  }

  MatrixFreeTest<dim,fe_degree,number> mf (mf_data);
  Vector<number> in (dof.n_dofs()), out (dof.n_dofs());
  Vector<number> in_dist (dof.n_dofs());
  Vector<number> out_dist (in_dist);

  for (unsigned int i=0; i<dof.n_dofs(); ++i)
    {
      if (constraints.is_constrained(i))
        continue;
      const double entry = Testing::rand()/(double)RAND_MAX;
      in(i) = entry;
      in_dist(i) = entry;
    }

  mf.vmult (out_dist, in_dist);


  // assemble sparse matrix with (\nabla v, \nabla u) + (v, 10 * u)
  SparsityPattern sparsity;
  {
    CompressedSimpleSparsityPattern csp(dof.n_dofs(), dof.n_dofs());
    DoFTools::make_sparsity_pattern (dof, csp, constraints, true);
    sparsity.copy_from(csp);
  }
  SparseMatrix<double> sparse_matrix (sparsity);
  {
    QGauss<dim>  quadrature_formula(fe_degree+2);

    FEValues<dim> fe_values (dof.get_fe(), quadrature_formula,
                             update_values    |  update_gradients |
                             update_JxW_values);

    const unsigned int   dofs_per_cell = dof.get_fe().dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof.begin_active(),
    endc = dof.end();
    for (; cell!=endc; ++cell)
      {
        cell_matrix = 0;
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

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global (cell_matrix,
                                                local_dof_indices,
                                                sparse_matrix);
      }
  }

  sparse_matrix.vmult (out, in);
  out -= out_dist;
  const double diff_norm = out.linfty_norm() / out_dist.linfty_norm();

  deallog << "Norm of difference: " << diff_norm << std::endl << std::endl;
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
  typename Triangulation<dim>::active_cell_iterator
  cell = tria.begin_active (),
  endc = tria.end();
  for (; cell!=endc; ++cell)
    if (cell->center().norm()<1e-8)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FE_DGP<dim> fe (fe_degree);
  DoFHandler<dim> dof (tria);
  dof.distribute_dofs(fe);
  ConstraintMatrix constraints;

  do_test<dim, fe_degree, double> (dof, constraints);
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
    deallog.pop();
    deallog.push("3d");
    test<3,1>();
    test<3,2>();
    deallog.pop();
  }
}



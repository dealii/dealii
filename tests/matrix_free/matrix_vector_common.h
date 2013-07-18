//------------------  matrix_vector_common.h  ------------------------
//    $Id$
//    Version: $Name$
//
//------------------  matrix_vector_common.h  ------------------------


// this is a template for matrix-vector products with the Helmholtz equation
// (zero and first derivatives) on different kinds of meshes (Cartesian,
// general, with and without hanging nodes). It also tests the multithreading
// in case it was enabled

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
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>


// forward declare this function. will be implemented in .cc files
template <int dim, int fe_degree>
void test ();



template <int dim, int fe_degree, typename Number>
void
helmholtz_operator (const MatrixFree<dim,Number>  &data,
                    Vector<Number>       &dst,
                    const Vector<Number> &src,
                    const std::pair<unsigned int,unsigned int> &cell_range)
{
  FEEvaluation<dim,fe_degree,fe_degree+1,1,Number> fe_eval (data);
  const unsigned int n_q_points = fe_eval.n_q_points;

  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
    {
      fe_eval.reinit (cell);
      fe_eval.read_dof_values (src);
      fe_eval.template evaluate (true, true, false);
      for (unsigned int q=0; q<n_q_points; ++q)
        {
          fe_eval.submit_value (Number(10)*fe_eval.get_value(q),q);
          fe_eval.submit_gradient (fe_eval.get_gradient(q),q);
        }
      fe_eval.template integrate (true,true);
      fe_eval.distribute_local_to_global (dst);
    }
}



template <int dim, int fe_degree, typename Number>
class MatrixFreeTest
{
public:
  typedef VectorizedArray<Number> vector_t;

  MatrixFreeTest(const MatrixFree<dim,Number> &data_in):
    data (data_in)
  {};

  void vmult (Vector<Number>       &dst,
              const Vector<Number> &src) const
  {
    dst = 0;
    const std_cxx1x::function<void(const MatrixFree<dim,Number> &,
                                   Vector<Number> &,
                                   const Vector<Number> &,
                                   const std::pair<unsigned int,unsigned int> &)>
    wrap = helmholtz_operator<dim,fe_degree,Number>;
    data.cell_loop (wrap, dst, src);
  };

private:
  const MatrixFree<dim,Number> &data;
};





template <int dim, int fe_degree, typename number>
void do_test (const DoFHandler<dim> &dof,
              const ConstraintMatrix &constraints,
              const unsigned int     parallel_option = 0)
{

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;
  if (parallel_option > 0)
    deallog << "Parallel option: " << parallel_option << std::endl;
  //std::cout << "Number of cells: " << dof.get_tria().n_active_cells() << std::endl;
  //std::cout << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;
  //std::cout << "Number of constraints: " << constraints.n_constraints() << std::endl;

  MatrixFree<dim,number> mf_data;
  {
    const QGauss<1> quad (fe_degree+1);
    typename MatrixFree<dim,number>::AdditionalData data;
    if (parallel_option == 1)
      data.tasks_parallel_scheme =
        MatrixFree<dim,number>::AdditionalData::partition_color;
    else if (parallel_option == 2)
      data.tasks_parallel_scheme =
        MatrixFree<dim,number>::AdditionalData::color;
    else
      {
        Assert (parallel_option == 0, ExcInternalError());
        data.tasks_parallel_scheme =
          MatrixFree<dim,number>::AdditionalData::partition_partition;
      }
    data.tasks_block_size = 7;

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
      const double entry = rand()/(double)RAND_MAX;
      in(i) = entry;
      in_dist(i) = entry;
    }

  mf.vmult (out_dist, in_dist);


  // assemble sparse matrix with (\nabla v,
  // \nabla u) + (v, 10 * u)
  SparsityPattern sparsity;
  {
    CompressedSimpleSparsityPattern csp(dof.n_dofs(), dof.n_dofs());
    DoFTools::make_sparsity_pattern (dof, csp, constraints, true);
    sparsity.copy_from(csp);
  }
  SparseMatrix<double> sparse_matrix (sparsity);
  {
    QGauss<dim>  quadrature_formula(fe_degree+1);

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


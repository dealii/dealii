//------------------  get_functions_common.h  ------------------------
//    Version: $Name$
//
//------------------  get_functions_common.h  ------------------------


// this is a template for getting the function values and comparing them with
// the output of FEValues on different kinds of meshes (Cartesian, general,
// with and without hanging nodes). The tests does not include multithreading
// because FEValues is not thread-safe

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

#include <fstream>
#include <iostream>


// forward declare this function. will be implemented in .cc files
template <int dim, int fe_degree>
void test ();



template <int dim, int fe_degree, int n_q_points_1d=fe_degree+1, typename Number=double>
class MatrixFreeTest
{
public:
  MatrixFreeTest(const MatrixFree<dim,Number> &data_in):
    data   (data_in),
    fe_val (data.get_dof_handler().get_fe(),
            Quadrature<dim>(data.get_quadrature(0)),
            update_values | update_gradients | update_hessians)
  {};

  MatrixFreeTest(const MatrixFree<dim,Number> &data_in,
                 const Mapping<dim>               &mapping):
    data   (data_in),
    fe_val (mapping, data.get_dof_handler().get_fe(),
            Quadrature<dim>(data.get_quadrature(0)),
            update_values | update_gradients | update_hessians)
  {};

  virtual ~MatrixFreeTest ()
  {}

  // make function virtual to allow derived
  // classes to define a different function
  virtual void
  operator () (const MatrixFree<dim,Number> &data,
               Vector<Number> &,
               const Vector<Number> &src,
               const std::pair<unsigned int,unsigned int> &cell_range) const
  {
    FEEvaluation<dim,fe_degree,n_q_points_1d,1,Number> fe_eval (data);

    std::vector<double> reference_values (fe_eval.n_q_points);
    std::vector<Tensor<1,dim> > reference_grads (fe_eval.n_q_points);
    std::vector<Tensor<2,dim> > reference_hess (fe_eval.n_q_points);

    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        fe_eval.reinit (cell);
        fe_eval.read_dof_values(src);
        fe_eval.evaluate (true,true,true);

        // compare values with the ones the FEValues
        // gives us. Those are seen as reference
        for (unsigned int j=0; j<data.n_components_filled(cell); ++j)
          {
            fe_val.reinit (data.get_cell_iterator(cell,j));
            fe_val.get_function_values(src, reference_values);
            fe_val.get_function_gradients(src, reference_grads);
            fe_val.get_function_hessians(src, reference_hess);

            for (int q=0; q<(int)fe_eval.n_q_points; q++)
              {
                errors[0] += std::fabs(fe_eval.get_value(q)[j]-
                                       reference_values[q]);
                for (unsigned int d=0; d<dim; ++d)
                  errors[1] += std::fabs(fe_eval.get_gradient(q)[d][j]-
                                         reference_grads[q][d]);
                errors[2] += std::fabs(fe_eval.get_laplacian(q)[j]-
                                       trace(reference_hess[q]));
                for (unsigned int d=0; d<dim; ++d)
                  {
                    errors[3] += std::fabs(fe_eval.get_hessian_diagonal(q)[d][j]-
                                           reference_hess[q][d][d]);
                    for (unsigned int e=0; e<dim; ++e)
                      errors[4] += std::fabs(fe_eval.get_hessian(q)[d][e][j]-
                                             reference_hess[q][d][e]);
                  }

                total[0] += std::fabs(reference_values[q]);
                for (unsigned int d=0; d<dim; ++d)
                  total[1] += std::fabs(reference_grads[q][d]);

                // reference for second derivatives computed
                // from fe_eval because FEValues is not
                // accurate enough with finite differences
                total[2] += std::fabs(fe_eval.get_laplacian(q)[j]);
                for (unsigned int d=0; d<dim; ++d)
                  {
                    total[3] += std::fabs(fe_eval.get_hessian_diagonal(q)[d][j]);
                    for (unsigned int e=0; e<dim; ++e)
                      total[4] += std::fabs(fe_eval.get_hessian(q)[d][e][j]);
                  }
              }
          }
      }
  }



  void test_functions (const Vector<Number> &src) const
  {
    for (unsigned int i=0; i<5; ++i)
      {
        errors[i] = 0;
        total[i]  = 0;
      }
    Vector<Number> dst_dummy;
    data.cell_loop (&MatrixFreeTest::operator(), this, dst_dummy, src);

    // for doubles, use a stricter condition than
    // for floats for the relative error size
    if (types_are_equal<Number,double>::value == true)
      {
        deallog.threshold_double (5e-14);
        deallog << "Error function values: "
                << errors[0]/total[0] << std::endl;
        deallog << "Error function gradients: "
                << errors[1]/total[1] << std::endl;

        // need to set quite a loose tolerance because
        // FEValues approximates Hessians with finite
        // differences, which are not so
        // accurate. moreover, Hessians are quite
        // large since we chose random numbers. for
        // some elements, it might also be zero
        // (linear elements on quadrilaterals), so
        // need to check for division by 0, too.
        deallog.threshold_double (5e-7);
        const double output2 = total[2] == 0 ? 0. : errors[2] / total[2];
        deallog << "Error function Laplacians: " << output2 << std::endl;
        const double output3 = total[3] == 0 ? 0. : errors[3] / total[3];
        deallog << "Error function diagonal of Hessian: " << output3 << std::endl;
        const double output4 = total[4] == 0 ? 0. : errors[4] / total[4];
        deallog << "Error function Hessians: " << output4 << std::endl;
      }
    else if (types_are_equal<Number,float>::value == true)
      {
        deallog.threshold_double (1e-6);
        deallog << "Error function values: "
                << errors[0]/total[0] << std::endl;
        deallog << "Error function gradients: "
                << errors[1]/total[1] << std::endl;
        const double output2 = total[2] == 0 ? 0. : errors[2] / total[2];
        deallog.threshold_double (1e-5);
        deallog << "Error function Laplacians: " << output2 << std::endl;
        const double output3 = total[3] == 0 ? 0. : errors[3] / total[3];
        deallog << "Error function diagonal of Hessian: " << output3 << std::endl;
        const double output4 = total[4] == 0 ? 0. : errors[4] / total[4];
        deallog << "Error function Hessians: " << output4 << std::endl;
      }
    deallog << std::endl;
  };

protected:
  const MatrixFree<dim,Number> &data;
  mutable FEValues<dim> fe_val;
  mutable double errors[5], total[5];
};



// dummy with empty quadrature formula
template <int dim, int fe_degree,typename Number>
class MatrixFreeTest<dim,fe_degree,0,Number>
{
public:
  MatrixFreeTest(const MatrixFree<dim,Number> &)
  {};

  MatrixFreeTest(const MatrixFree<dim,Number> &,
                 const Mapping<dim> &)
  {};

  void cell_integration (Vector<Number> &,
                         const Vector<Number> &,
                         const std::pair<unsigned int,unsigned int>) const {}

  void test_functions (const Vector<Number> &) const
  {}
};




template <int dim, int fe_degree, typename number>
void do_test (const DoFHandler<dim> &dof,
              const ConstraintMatrix &constraints)
{
  deallog << "Testing " << dof.get_fe().get_name() << std::endl;
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
    mf_data.reinit (dof, constraints, quad, data);
  }

  MatrixFreeTest<dim,fe_degree,fe_degree+1,number> mf (mf_data);
  mf.test_functions(solution);
}


int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  deallog << std::setprecision (3);
  {
    deallog.push("2d");
    test<2,1>();
    test<2,2>();
    test<2,3>();
    test<2,4>();
    deallog.pop();
    deallog.push("3d");
    test<3,1>();
    test<3,2>();
    deallog.pop();
  }
}

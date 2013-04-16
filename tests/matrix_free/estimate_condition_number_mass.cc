//------------------  estimate_condition_number_mass.cc  ---------------------
//    $Id$
//    Version: $Name$
//
//------------------  estimate_condition_number_mass.cc  ---------------------


// this function computes condition number estimates for the mass matrix at
// different polynomial degrees. The mesh uses a hypercube mesh with no
// hanging nodes and no other constraints

#include "../tests.h"

std::ofstream logfile("estimate_condition_number_mass/output");

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/function_lib.h>


template <int dim, int fe_degree, typename Number>
void
mass_operator (const MatrixFree<dim,Number>  &data,
               Vector<Number>       &dst,
               const Vector<Number> &src,
               const std::pair<unsigned int,unsigned int> &cell_range)
{
  FEEvaluationGeneral<dim,fe_degree,fe_degree+1,1,Number> fe_eval (data);
  const unsigned int n_q_points = fe_eval.n_q_points;

  for(unsigned int cell=cell_range.first;cell<cell_range.second;++cell)
    {
      fe_eval.reinit (cell);
      fe_eval.read_dof_values (src);
      fe_eval.template evaluate (true, false, false);
      for (unsigned int q=0; q<n_q_points; ++q)
        {
          fe_eval.submit_value (fe_eval.get_value(q),q);
        }
      fe_eval.template integrate (true,false);
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
    const std_cxx1x::function<void(const MatrixFree<dim,Number>  &,
                                   Vector<Number>       &,
                                   const Vector<Number> &,
                                   const std::pair<unsigned int,unsigned int>&)>
      wrap = mass_operator<dim,fe_degree,Number>;
    data.cell_loop (wrap, dst, src);
  };

private:
  const MatrixFree<dim,Number> &data;
};



template <int dim, int fe_degree>
void test (const FiniteElement<dim> &fe)
{
  typedef double number;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria);
  tria.refine_global(2);

  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  ConstraintMatrix constraints;
  constraints.close();

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;

  MatrixFree<dim,number> mf_data;
  {
    const QGauss<1> quad (fe_degree+1);
    mf_data.reinit (dof, constraints, quad);
  }

  MatrixFreeTest<dim,fe_degree,number> mf (mf_data);
  Vector<number> in (dof.n_dofs()), out (dof.n_dofs());

  VectorTools::create_right_hand_side(dof, QGauss<dim>(fe_degree+1),
                                      Functions::CosineFunction<dim>(), in);

  SolverControl control(10000, 1e-9*in.l2_norm());
  typename SolverCG<>::AdditionalData data;
  data.compute_condition_number = true;
  SolverCG<> solver(control,data);
  solver.solve(mf, out, in, PreconditionIdentity());
  deallog << std::endl;
}


int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  deallog << std::setprecision (2);

  {
    deallog.threshold_double(1.e-9);
    deallog.push("2d");
    test<2,1>(FE_Q<2>(1));
    test<2,1>(FE_DGQ<2>(1));
    test<2,2>(FE_Q<2>(2));
    test<2,2>(FE_DGQ<2>(2));
    test<2,4>(FE_Q<2>(4));
    test<2,4>(FE_Q<2>(QGaussLobatto<1>(5)));
    test<2,4>(FE_DGQ<2>(4));
    test<2,4>(FE_Q_Hierarchical<2>(4));
    test<2,6>(FE_Q<2>(6));
    test<2,6>(FE_Q<2>(QGaussLobatto<1>(7)));
    test<2,6>(FE_DGQ<2>(6));
    test<2,6>(FE_DGQArbitraryNodes<2>(QGaussLobatto<1>(7)));
    test<2,10>(FE_Q<2>(10));
    test<2,10>(FE_Q<2>(QGaussLobatto<1>(11)));
    test<2,16>(FE_Q<2>(QGaussLobatto<1>(17)));
    deallog.pop();
    deallog.push("3d");
    test<3,1>(FE_Q<3>(1));
    test<3,2>(FE_Q<3>(2));
    test<3,5>(FE_Q<3>(5));
    test<3,5>(FE_Q<3>(QGaussLobatto<1>(6)));
    deallog.pop();
  }
}



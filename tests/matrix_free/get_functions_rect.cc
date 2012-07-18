//------------------  get_functions_rect.cc  ------------------------
//    $Id$
//    Version: $Name$
//
//------------------  get_functions_rect.cc  ------------------------


// this function tests the correctness of the implementation of matrix free
// operations in getting the function values, the function gradients, and the
// function Laplacians on a hyperball mesh with different sizes in the number
// of degrees of freedom per cell and quadrature points per cell. The test
// case includes hanging node constraints, but no constraints from boundary
// conditions

#include "../tests.h"

std::ofstream logfile("get_functions_rect/output");

#include "get_functions_common.h"


template <int dim, int fe_degree, int n_q_points_1d, typename number>
void sub_test (const DoFHandler<dim> &dof,
               const ConstraintMatrix &constraints,
               MatrixFree<dim,number> &mf_data,
               Vector<number> &solution)
{
  deallog << "Test with fe_degree " << fe_degree
          << ", n_q_points_1d: " << (n_q_points_1d) << std::endl;
  const QGauss<1> quad (n_q_points_1d);
  MappingQ<dim> mapping (2);
  typename MatrixFree<dim,number>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim,number>::AdditionalData::none;
  data.mapping_update_flags = update_gradients | update_second_derivatives;
  mf_data.reinit (mapping, dof, constraints, quad, data);
  MatrixFreeTest<dim,fe_degree,n_q_points_1d,number> mf (mf_data, mapping);
  mf.test_functions (solution);
}



template <int dim, int fe_degree>
void test ()
{
  typedef double number;
  Triangulation<dim> tria;
  GridGenerator::hyper_ball (tria);
  static const HyperBallBoundary<dim> boundary;
  tria.set_boundary (0, boundary);
                                // refine first and last cell
  tria.begin(tria.n_levels()-1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.refine_global (1);

  FE_Q<dim> fe (fe_degree);
  DoFHandler<dim> dof (tria);
  dof.distribute_dofs(fe);

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof, constraints);
  constraints.close();


                                // in the other functions, use do_test in
                                // get_functions_common, but here we have to
                                // manually choose non-rectangular tests.
  deallog << "Testing " << dof.get_fe().get_name() << std::endl;
  //std::cout << "Number of cells: " << dof.get_tria().n_active_cells()
  //          << std::endl;
  //std::cout << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;
  //std::cout << "Number of constraints: " << constraints.n_constraints() << std::endl;

  Vector<number> solution (dof.n_dofs());

                                // create vector with random entries
  for (unsigned int i=0; i<dof.n_dofs(); ++i)
    {
      if(constraints.is_constrained(i))
        continue;
      const double entry = rand()/(double)RAND_MAX;
      solution(i) = entry;
    }
  constraints.distribute(solution);

  MatrixFree<dim,number> mf_data;
  if (fe_degree > 1)
    sub_test <dim,fe_degree,fe_degree-1,number> (dof, constraints, mf_data,
                                                 solution);
  sub_test <dim,fe_degree,fe_degree,number> (dof, constraints, mf_data,
                                             solution);
  sub_test <dim,fe_degree,fe_degree+2,number> (dof, constraints, mf_data,
                                               solution);
  if (dim == 2)
    sub_test <dim,fe_degree,fe_degree+3,number> (dof, constraints, mf_data,
                                                 solution);
}


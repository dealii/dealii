//------------------  get_functions_gl.cc  ------------------------
//    $Id$
//    Version: $Name$
//
//------------------  get_functions_gl.cc  ------------------------


// this function tests the correctness of the implementation of matrix free
// operations in getting the function values, the function gradients, and the
// function Laplacians on a hyperball mesh for Gauss-Lobatto elements
// (identity values transformation). The test case includes hanging node
// constraints, but no constraints from boundary conditions

#include "../tests.h"


std::ofstream logfile("get_functions_gl/output");

#include "get_functions_common.h"


template <int dim, int n_dofs_1d, typename Number>
class MatrixFreeTestGL : public MatrixFreeTest<dim, n_dofs_1d, n_dofs_1d, Number>
{
 public:
  typedef VectorizedArray<Number> vector_t;
  static const std::size_t n_vectors = VectorizedArray<Number>::n_array_elements;

  MatrixFreeTestGL(const MatrixFree<dim,Number> &data,
	       const Mapping<dim>               &mapping):
    MatrixFreeTest<dim, n_dofs_1d, n_dofs_1d, Number>(data, mapping)
  {};

  void operator() (const MatrixFree<dim,Number> &data,
		   Vector<Number>       &,
		   const Vector<Number> &src,
		   const std::pair<unsigned int,unsigned int> &cell_range) const
  {
    FEEvaluationGL<dim,n_dofs_1d,1,Number> fe_eval (this->data);
    for(unsigned int cell=cell_range.first;cell<cell_range.second;++cell)
      {
	fe_eval.reinit (cell);
	std::vector<double> reference_values (fe_eval.n_q_points);
	std::vector<Tensor<1,dim> > reference_grads (fe_eval.n_q_points);
	std::vector<Tensor<2,dim> > reference_hess (fe_eval.n_q_points);
	fe_eval.read_dof_values(src);
	fe_eval.evaluate (true,true,true);

				// compare values with the ones the FEValues
				// gives us. Those are seen as reference
	for (unsigned int j=0; j<data.n_components_filled(cell); ++j)
	  {
	    this->fe_val.reinit (data.get_cell_iterator(cell,j));
	    this->fe_val.get_function_values(src, reference_values);
	    this->fe_val.get_function_gradients(src, reference_grads);
	    this->fe_val.get_function_hessians(src, reference_hess);

	    for (int q=0; q<(int)fe_eval.n_q_points; q++)
	      {
		this->errors[0] += std::fabs(fe_eval.get_value(q)[j]-
					     reference_values[q]);
		for (unsigned int d=0; d<dim; ++d)
		  this->errors[1] += std::fabs(fe_eval.get_gradient(q)[d][j]-
					       reference_grads[q][d]);
		this->errors[2] += std::fabs(fe_eval.get_laplacian(q)[j]-
					     trace(reference_hess[q]));
		this->total[0] += std::fabs(reference_values[q]);
		for (unsigned int d=0; d<dim; ++d)
		  this->total[1] += std::fabs(reference_grads[q][d]);
		this->total[2] += std::fabs(fe_eval.get_laplacian(q)[j]);
	      }
	  }
      }
  }
};



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
  tria.refine_global (4-dim);

  FE_Q<dim> fe (QGaussLobatto<1>(fe_degree+1));
  DoFHandler<dim> dof (tria);
  dof.distribute_dofs(fe);

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof, constraints);
  constraints.close();

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;
  //std::cout << "Number of cells: " << dof.get_tria().n_active_cells()
  //	  << std::endl;
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
  deallog << "Test with fe_degree " << fe_degree
	  << std::endl;
  const QGaussLobatto<1> quad (fe_degree+1);
  MappingQ<dim> mapping (2);
  typename MatrixFree<dim,number>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim,number>::AdditionalData::none;
  data.mapping_update_flags = update_gradients | update_second_derivatives;
  mf_data.reinit (mapping, dof, constraints, quad, data);
  MatrixFreeTestGL<dim,fe_degree+1,number> mf (mf_data, mapping);
  mf.test_functions (solution);
}


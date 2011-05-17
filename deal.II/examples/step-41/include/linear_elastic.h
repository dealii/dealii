#ifndef LINEARELASTIC
#define LINEARELASTIC

#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/filtered_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/dofs/dof_tools.h>


#include "parsed_symmetric_tensor_function.h"
#include "local_assemble_elastic_matrix.h"
#include "local_assemble_elastic_rhs.h"
#include "vector_space.h"
#include "my_tools.h"

using namespace dealii;
using namespace dealii::Functions;

template <int dim>
class LinearElastic
{
 public:
  
  LinearElastic(); 
  
  ~LinearElastic(); 

  void declare_parameters(ParameterHandler &prm);

  void parse_parameters(ParameterHandler &prm);

  void reinit(VectorSpace<dim> &vspace);

  void build_matrix(ParameterHandler &prm,
		    VectorSpace<dim> &vspace);

  void reinit_step(double &time);

  void build_rhs(VectorSpace<dim> &vspace);

  void solve(VectorSpace<dim> &vspace,
	     double &tolerance);

  SparseMatrix<double> A;

  ParsedSymmetricTensorFunction<4, dim>  C;

  //The solution vectors
  Vector<double> sol_total;
  Vector<double> sol_increment;
			
 private:

  //various functions for the laplace equation
  ParsedFunction<dim> dbc;
  ParsedFunction<dim> nbc;
  ParsedFunction<dim> bf;

  //The system matrix for linear elasticity
  SparsityPattern  sp_A;
  
  //The distributed right hand side      
  Vector<double>               rhs;



};

#endif

#ifndef LINEARELASTIC
#define LINEARELASTIC

#include <base/logstream.h>
#include <base/parameter_handler.h>
#include <base/parsed_function.h>
#include <lac/vector.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>
#include <lac/filtered_matrix.h>
#include <lac/precondition.h>
#include <lac/sparse_ilu.h>
#include <lac/solver_cg.h>
#include <lac/sparse_direct.h>
#include <dofs/dof_tools.h>


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

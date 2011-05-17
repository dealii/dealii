#ifndef CAMCLAY
#define CAMCLAY

//deal.ii packages
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/table.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/matrices.h>

//Trilinos packages
#include <Sacado.hpp>

//LPCM includes
#include "vector_space.h"
#include "local_assemble_plastic_project.h"
#include "local_assemble_scalar_project.h"
#include "my_tools.h"
#include "domain.h"

using namespace dealii;
using namespace dealii::Functions;

template <int dim>
class CamClay
{
 public:
  
  CamClay(); 
 
  ~CamClay(); 
	
  void declare_parameters(ParameterHandler &prm);

  void parse_parameters(ParameterHandler &prm);

  void reinit(VectorSpace<dim> &vspace);

  void reinit_step();

  void project_strain(ParameterHandler &prm,
		      VectorSpace<dim> &vspace,
		      SparseMatrix<double> &A);

  void project_hardening(VectorSpace<dim> &vspace);

  void build_matrix(VectorSpace<dim> &vspace);

  void initial_conditions(VectorSpace<dim> &vspace);

  void update_internal_variables(VectorSpace<dim> &vspace,
				 Vector<double> &elastic_solution,
                                 ParsedSymmetricTensorFunction<4,dim> &C);

  void compute_jacobian(FullMatrix<double> &jac,
                        const Vector<double> &sol,
                        const double p,
                        const double q,
                        const double M,
                        const SymmetricTensor<4,dim> &C_qp,
                        const SymmetricTensor<2,dim> &xi);

  Vector<double> dF_dstress(const SymmetricTensor<2,dim> &stress,
                            const double k);

  //this is the plastic strain projected into the displacement  
  Vector<double> solution;

  //this is the average number of iterations per cell for the current step  
  Vector<double> iterations;
  
  //internal variable vectors
  Table<2, SymmetricTensor<2,dim> > plastic_strain;
  Table<2, double > hardening;
  Table<2, double> iter_table;

  Vector<double> sol_hard_iter;
		
 private:

  //computes the yeild function
  inline double yield_function(const double p, const double q,
			       const double k);

  //does the actual solving for the new internal variables
  double solve(const int &index,
               const unsigned int &qp,
               const SymmetricTensor<2,dim> &trial_strain);
  
  //the potential for the hardening
  double h(const double p,
	   const double k);
  
  //computes the hydrostatic stress invariant
  inline double p(const SymmetricTensor<2,dim> &stress);

  //computes the deviatoric stress invariant
  inline double q(const SymmetricTensor<2,3> &xi);

  //computes the deviatoric stress tensor
  inline SymmetricTensor<2,3> xi(const SymmetricTensor<2,dim> &stress);

  //computer_theresidua and the jacobian
  void compute_res_jac(Vector<double> &res,
		       FullMatrix<double> &jac,
		       const Vector<double> &sol,
		       const Vector<double> &prev_step,
		       const SymmetricTensor<2,dim> &total_strain);
  /*
  //computes the residual
  void compute_residual(Vector<double> &res,
                        const Vector<double> &prev_step,
                        const Vector<double> &sol,
                        const SymmetricTensor<2,dim> &stress);



  //creates the jacobian numerical
  void numerical_jacobian(FullMatrix<double> &jac,
			  const Vector<double> &sol,
			  const Vector<double> &prev_sol,
			  const SymmetricTensor<2,dim> &total_strain);*/

  //translates from sym tensor notation to voigt notation
  inline unsigned int sym2voigt(const unsigned int i,
				const unsigned int j);

  inline unsigned int sac_num(const unsigned int i,
			      const unsigned int j);

  //translates voigt notation to sym tensor indices
  inline std::vector<unsigned int> voigt2sym(const unsigned int i);

  //delta tensor
  inline double delta(const unsigned int i,
		      const unsigned int j);

  //voigt delta tensor
  inline double delta(const unsigned int i);

  //the yeild stress function - not needed for cam clay actually
  ParsedFunction<dim> yield_stress;

  //some constants
  //double shear_mod;
  //double bulk_mod;
  double M;

  SymmetricTensor<4,dim> C_qp;

  SparseMatrix<double> MM;

  SparsityPattern sp_MM;


  
};

#endif

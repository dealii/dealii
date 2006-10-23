/*    $Id: project.cc descends from heat-equation.cc, which       */
/*         descends from step-4.cc (2006/03/01).                  */
/*    Author: Ivan Christov, Texas A&M University, 2006           */

/*    $Id: step-4.cc,v 1.34 2006/02/06 21:33:10 wolf Exp $        */
/*    Version: $Name:  $                                          */
/*                                                                */
/*    Copyright (C) 2006 by the deal.II authors                   */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the text and       */
/*    further information on this license.                        */


 // @sect3{Include files and global variables}

// For an explanation of the include files, the reader should refer to
// the example programs step-1 through step-4. They are in the
// standard order, which is ``base'' -- ``lac'' -- ``grid'' --
// ``dofs'' -- ``fe'' -- ``numerics'' (since each of these categories
// roughly builds upon previous ones), then a few C++ headers for
// file, input/output and string streams.
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_gmres.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_constraints.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <numerics/data_out_stack.h>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace dealii;

// The following global variable is used to determine whether the
// problem being solved is one for which an exact solution is known,
// e.g. we are using the exact solution as the initial condition. It
// is set to zero by default, and modified by ``InitialValues::value''
// (see below). Things such as the computation of the error between
// the numerical and exact solutions depend on the value of this
// variable.
bool exact_solution_known = false;

// @sect3{The ``SineGordonProblem'' class template}

// The entire algorithm for solving the problem is encapsulated in
// this class. Also, note that the class is declared with a template
// parameter, which is the spatial dimension, so that we can solve the
// sine-Gordon equation in one, two or three spatial dimension. For
// more on the dimension-independent class-encapsulation of the
// problem, the reader should consult step-3 and step-4.
template <int dim>
class SineGordonProblem 
{
  public:
    SineGordonProblem ();
    void run ();
    
  private:
    void make_grid_and_dofs ();
    void assemble_system ();
    void compute_nl_term (const Vector<double> &old_data, 
			  const Vector<double> &new_data,
			  Vector<double>       &nl_term) const;
    void compute_nl_matrix (const Vector<double> &old_data, 
			    const Vector<double> &new_data,
			    SparseMatrix<double> &nl_matrix) const;
    void solve ();
    void compute_error  (const unsigned int timestep_number);
    void output_results (const unsigned int timestep_number);

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    SparseMatrix<double> mass_matrix;
    SparseMatrix<double> laplace_matrix;
    
    double time, final_time, time_step;
    double theta;

    Vector<double>       solution, d_solution, old_solution;
    Vector<double>       massmatxvel;
    Vector<double>       system_rhs;
    Vector<double>       fem_errors;

    DataOutStack<dim>    data_out_stack;

    static const unsigned int output_timestep_skip = 1;
    static const int n_global_refinements = 6;
};

// @sect3{Exact solitary wave solutions of the sine-Gordon equation}

// A kink-like solitary wave solution to the (``dim''+1) dimensional
// sine-Gordon equation, which we can test our code against, is given
// by Leibbrandt in \e Phys. \e Rev. \e Lett. \b 41(7), and is
// implemented in the ``ExactSolution'' class.  However, it should be
// noted that a closed-form solution can only be obtained for the
// infinite-line initial-value problem (not the Neumann
// initial-boundary-value problem under consideration here). However,
// given that we impose \e zero Neumann boundary conditions, we expect
// that the solution to our initial-boundary-value problem would be
// close (in fact, equal) to the solution infinite-line initial-value
// problem, if reflections of waves off the boundaries of our domain
// do \e not occur.
// 
// The constants $\vartheta$ (``th'') and $\lambda$ (``ld'') in the 2D
// solution and $\vartheta$ (``th''), $\phi$ (``phi'') and $\tau$
// (``tau'') in the 3D solution are called the B&auml;cklund
// transformation parameters. They control such things as the
// orientation and steepness of the kink. For the purposes of testing
// the code against the exact solution, one should choose the
// parameters so that the kink is aligned with the grid, e.g. $\vartheta
// = \phi = \pi$.
//
// In 1D, more interesting analytical solutions are known. Many of
// them are listed on
// http://mathworld.wolfram.com/Sine-GordonEquation.html . We have
// implemented the one kink, two kink, kink-antikink and stationary
// breather solitary-wave solutions.
template <int dim>
class ExactSolution : public Function<dim>
{
  public:
    ExactSolution (const unsigned int n_components = 1,
		   const double time = 0.) : Function<dim>(n_components, time) {};
    virtual double value (const Point<dim> &p,
			  const unsigned int component = 0) const;
};

template <int dim>
double ExactSolution<dim>::value (const Point<dim> &p,
				  const unsigned int /*component*/) const
{
  double t = this->get_time ();

  switch (dim)
    {
      case 1:
      {
	double m = 0.5;
//      double beta = std::sqrt(m*m-1.)/m;
	double c1 = 0.;
	double c2 = 0.;
//       double s1 = 1.;
//       double s2 = -1.;

					 /* one kink (m>1) */
					 /* return 4.*std::atan(std::exp(s1*(p[0]+s2*beta*t)/std::sqrt(1.-beta*beta))); */

					 /* two kinks (m>1) */
					 /* return 4.*std::atan(beta*std::sinh(beta*m*p[0])/std::cosh(beta*m*t)); */

					 /* kink-antikink (m>1) */
					 /* return -4.*std::atan(m/std::sqrt(m*m-1)*std::sinh(std::sqrt(m*m-1.)*t+c2)/
					    std::cosh(m*p[0]+c1)); */

					 /* stationary breather (m<1), period = 2.*pi*sqrt(1.-m*m) 
					    for m=0.5, -5.4414 <= t <= 2.7207 is a good time interval */
	return -4.*std::atan(m/std::sqrt(1.-m*m)*std::sin(std::sqrt(1.-m*m)*t+c2)
			     /std::cosh(m*p[0]+c1));
      }

      case 2:
      {
	double th  = deal_II_numbers::PI/4.;
	double ld  = 1.;
	double a0  = 1.;
	double s   = 1.;
	double arg = 0.;
	arg = p[0]*std::cos(th) + std::sin(th)*(p[1]*std::cosh(ld)+t*std::sinh(ld));
	return 4.*std::atan(a0*std::exp(s*arg));
      }

      case 3:
      {
	double th  = deal_II_numbers::PI;
	double phi = deal_II_numbers::PI;
	double tau = 1.;
	double c0  = 1.;
	double s   = 1.;
	double arg = 0.;
	arg = (p[0]*std::cos(th) + p[1]*std::sin(th)*std::cos(phi)
	       + std::sin(th)*std::sin(phi)*(p[2]*std::cosh(tau)+t*std::sinh(tau)));
	return 4.*std::atan(c0*std::exp(s*arg));
      }

      default:
	    Assert (false, ExcNotImplemented());
	    return -1e8;
    }
}

// @sect3{Boundary values and initial values}

// For our problem, we do not enforce Dirichlet boundary conditions
// and the Neumann boundary conditions are enforced directly through
// the variational formulation. However, since our problem is time
// dependent, we must specify the value of the independent variable
// $u$ at the initial time $t_0$. We do so via the ``InitialValues''
// class below.
template <int dim>
class InitialValues : public Function<dim>
{
  public:
    InitialValues (const unsigned int n_components = 1, 
		   const double time = 0.) : Function<dim>(n_components, time) {};
  
    virtual double value (const Point<dim> &p,
			  const unsigned int component = 0) const;
};

template <int dim>
double InitialValues<dim>::value (const Point<dim> &p,
				  const unsigned int /*component*/) const 
{   
  // We could also use a localized wave form for our initial
  // condition, and see how it evolves when governed by the
  // sine-Gordon equation. An example of such an initial condition is
  // the following:
  /*
  exact_solution_known = false;
  if ((p[0]>=-M_PI) && (p[0]<=M_PI) && (p[1]>=-M_PI) && (p[1]<=M_PI)) {
    return std::cos(p[0]/2.)*std::cos(p[1]/2.);
  } else {
    return 0.;
  }
  */

  // In 2D, another possibility for a localized-wave initial condition
  // is a separable solution composed of two 1D breathers:
  exact_solution_known = false;
  double m = 0.5;
  double t = this->get_time();
  double argx = m/std::sqrt(1-m*m)*std::sin(std::sqrt(1-m*m)*t)/std::cosh(m*p[0]);
  double argy = m/std::sqrt(1-m*m)*std::sin(std::sqrt(1-m*m)*t)/std::cosh(m*p[1]);
  return 16.*std::atan(argx)*std::atan(argy);

  // For the purposes of validating the program, we can use an exact
  // solution of the sine-Gordon equation, at $t=t_0$, as the initial
  // condition for our problem. Though, perhaps, this is not the most
  // efficient way to implement the exact solution as the initial
  // conditons, it is instuctive.
  /*
  exact_solution_known = false;
  ExactSolution<dim> exact_solution (1, this->get_time());
  return exact_solution.value (p);
  */
}

// @sect3{Implementation of the ``SineGordonProblem'' class}

// \b TO \b DO: present the big picture here?

// @sect4{SineGordonProblem::SineGordonProblem}

// This is the constructor of the ``SineGordonProblem'' class. It
// specifies the desired polynomial degree of the finite elements,
// associates a ``DoFHandler'' to the ``triangulation'' object (just
// as in the example programs step-3 and step-4), initializes the
// current or initial time, the final time, the time step size, and
// the value of $\theta$ for the time stepping scheme.
//
// Note that if we were to chose the explicit Euler time stepping
// scheme ($\theta = 0$), then we must pick a time step $k \le h$,
// otherwise the scheme is not stable and oscillations might arise in
// the solution. The Crank-Nicolson scheme ($\theta = \frac{1}{2}$)
// and the implicit Euler scheme ($\theta=1$) do not suffer from this
// deficiency, since they are unconditionally stable. However, even
// then the time step should be chosen to be on the order of $h$ in
// order to obtain a good solution.
template <int dim>
SineGordonProblem<dim>::SineGordonProblem () :
                fe (1),
		dof_handler (triangulation),
		time (-5.4414/*0.*/),
		final_time (2.7207/*20.*/),
		time_step (10*1./std::pow(2.,n_global_refinements)),
		theta (0.5)
{}

// @sect4{SineGordonProblem::make_grid_and_dofs}

// This function creates a rectangular grid in ``dim'' dimensions and
// refines it several times. Also, all matrix and vector members of
// the ``SineGordonProblem'' class are initialized to their
// approrpiate sizes once the degrees of freedom have been
// assembled. Unlike its analogue in step-3 (and step-4) this function
// uses ``MatrixCreator'' class to generate a mass matrix $M$ and a
// Laplace matrix $A$ and store them in the appropriate variables
// for the remainder of the program's life.
template <int dim>
void SineGordonProblem<dim>::make_grid_and_dofs ()
{
  GridGenerator::hyper_cube (triangulation, -10, 10);
  triangulation.refine_global (n_global_refinements);
  
  std::cout << "   Number of active cells: "
	    << triangulation.n_active_cells()
	    << std::endl
	    << "   Total number of cells: "
	    << triangulation.n_cells()
	    << std::endl;

  dof_handler.distribute_dofs (fe);

  std::cout << "   Number of degrees of freedom: "
	    << dof_handler.n_dofs()
	    << std::endl;

  sparsity_pattern.reinit (dof_handler.n_dofs(), 
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress ();

  system_matrix.reinit  (sparsity_pattern);
  mass_matrix.reinit    (sparsity_pattern);
  laplace_matrix.reinit (sparsity_pattern);

  MatrixCreator::create_mass_matrix (dof_handler, QGauss<dim>(3), 
				     mass_matrix);
  MatrixCreator::create_laplace_matrix (dof_handler, QGauss<dim>(3), 
					laplace_matrix);

  solution.reinit       (dof_handler.n_dofs());
  d_solution.reinit     (dof_handler.n_dofs());
  old_solution.reinit   (dof_handler.n_dofs());
  massmatxvel.reinit    (dof_handler.n_dofs());
  system_rhs.reinit     (dof_handler.n_dofs());

  // We will use the ``fem_errors'' vector, which is of size equal to
  // the number of time steps, to store the errors in the finite
  // element solution after each time step. Note that we must make the
  // first element of the vector equal to zero, since there is no
  // error in the solution after zeroth time step because the solution
  // is just the initial condition.
  const unsigned int n_time_steps
    = static_cast<unsigned int>(std::ceil(std::fabs(final_time-time)/time_step));
  fem_errors.reinit (n_time_steps);
  fem_errors(0) = 0.;
}

// @sect4{SineGordonProblem::assemble_system}

// This functions assembles the system matrix and right-hand side
// vector for each iteration of Newton's method. The reader should
// refer to the last section of the Introduction for the explicit
// formulas for the system matrix and right-hand side.
template <int dim>
void SineGordonProblem<dim>::assemble_system () 
{  
  // First we assemble the Jacobian matrix $F'_h(U^n_l)$, where
  // $U^n_l$ is stored in the vector ``solution'' for convenience.
  system_matrix = 0;
  system_matrix.copy_from (mass_matrix);
  system_matrix.add (std::pow(time_step*theta,2), laplace_matrix);
  SparseMatrix<double> tmp_matrix (sparsity_pattern);
  compute_nl_matrix (old_solution, solution, tmp_matrix);
  system_matrix.add (-std::pow(time_step*theta,2), tmp_matrix);

  // Then, we compute the right-hand side vector $-F_h(U^n_l)$.
  system_rhs = 0;

  tmp_matrix = 0;
  tmp_matrix.copy_from (mass_matrix);
  tmp_matrix.add (std::pow(time_step*theta,2), laplace_matrix);
  Vector<double> tmp_vector (solution.size());
  tmp_matrix.vmult (tmp_vector, solution);
  system_rhs += tmp_vector;

  tmp_matrix = 0;
  tmp_matrix.copy_from (mass_matrix);
  tmp_matrix.add (-std::pow(time_step,2)*theta*(1-theta), laplace_matrix);
  tmp_vector = 0;
  tmp_matrix.vmult (tmp_vector, old_solution);
  system_rhs -= tmp_vector;

  system_rhs.add (-time_step, massmatxvel);

  tmp_vector = 0;
  compute_nl_term (old_solution, solution, tmp_vector);
  system_rhs.add (std::pow(time_step,2)*theta, tmp_vector);

  system_rhs *= -1;
}

// @sect4{SineGordonProblem::compute_nl_term}

// This function computes the vector $S(\cdot,\cdot)$ corresponding to the
// nonlinear term in the auxilliary (second) equation of the split
// formulation. This function not only simplifies the repeated
// computation of this term, but it is also a fundamental part of
// nonlinear iterative solver that we use when the time stepping is
// implicit (i.e. $\theta\ne 0$). Moreover, we must allow the function
// to receive as input an "old" and a "new" solution, which may not be
// the actual solutions of the problem stored in ``old_solution'' and
// ``solution.'' For the purposes of this function, let us call the
// first two arguments $w_{\mathrm{old}}$ and $w_{\mathrm{new}}$,
// respectively.
//
// It is perhaps worth investigating what order quadrature formula is
// best suited for this type of integration, since $\sin(\cdot)$ is an
// oscillatory function.
template <int dim>
void SineGordonProblem<dim>::compute_nl_term (const Vector<double> &old_data,
					      const Vector<double> &new_data,
					      Vector<double>       &nl_term) const
{
  QGauss<dim>   quadrature_formula (3);
  FEValues<dim> fe_values (fe, quadrature_formula, 
			   update_values | update_JxW_values | update_q_points);
  
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.n_quadrature_points;
  
  Vector<double> local_nl_term (dofs_per_cell);      
  std::vector<unsigned int> local_dof_indices (dofs_per_cell); 
  std::vector<double> old_data_values (n_q_points);
  std::vector<double> new_data_values (n_q_points);
  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

  for (; cell!=endc; ++cell)
    { 
      // Once we re-initialize our ``FEValues'' instantiation to the
      // current cell, we make use of the ``get_function_values''
      // routine to get the obtain the values of the "old" data
      // (presumably at $t=t_{n-1}$) and the "new" data (presumably at
      // $t=t_n$) at the nodes of the chosen quadrature formula.
      fe_values.reinit (cell);
      fe_values.get_function_values (old_data, old_data_values);
      fe_values.get_function_values (new_data, new_data_values);
      
      // Now, we can evaluate $\int_K \sin\left[\theta w_{\mathrm{new}} +
      // (1-\theta) w_{\mathrm{old}}\right]\,\varphi_j\,\mathrm{d}x$ using
      // the desired quadrature formula.
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i) 	   
	  local_nl_term(i) += (std::sin(theta*new_data_values.at(q_point) +
					(1-theta)*old_data_values.at(q_point)) *
			       fe_values.shape_value (i, q_point) *
			       fe_values.JxW (q_point));	    
      
      // We conclude by adding up the contributions of the
      // integrals over the cells to the global integral.
      cell->get_dof_indices (local_dof_indices);
      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	nl_term(local_dof_indices[i]) += local_nl_term(i);
      
      local_nl_term = 0;
    }
}

// @sect4{SineGordonProblem::compute_nl_matrix}

// This function computes the matrix $N(\cdot,\cdot)$ corresponding to
// the nonlinear term in the Jacobian of $F(\cdot)$. It is also a
// fundamental part of nonlinear iterative solver. Just as
// ``compute_nl_term'', we must allow this function to receive
// as input an "old" and a "new" solution, which we call the
// $w_{\mathrm{old}}$ and $w_{\mathrm{new}}$, respectively.
template <int dim>
void SineGordonProblem<dim>::compute_nl_matrix (const Vector<double> &old_data, 
						const Vector<double> &new_data,
						SparseMatrix<double> &nl_matrix) const
{
  QGauss<dim>   quadrature_formula (3);
  FEValues<dim> fe_values (fe, quadrature_formula, 
			   update_values | update_JxW_values | update_q_points);
  
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.n_quadrature_points;
  
  FullMatrix<double> local_nl_matrix (dofs_per_cell, dofs_per_cell);
  std::vector<unsigned int> local_dof_indices (dofs_per_cell); 
  std::vector<double> old_data_values (n_q_points);
  std::vector<double> new_data_values (n_q_points);
  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

  for (; cell!=endc; ++cell)
    { 
      // Again, first we re-initialize our ``FEValues'' instantiation
      // to the current cell.
      fe_values.reinit (cell);
      fe_values.get_function_values (old_data, old_data_values);
      fe_values.get_function_values (new_data, new_data_values);
      
      // Then, we evaluate $\int_K \cos\left[\theta w_{\mathrm{new}} +
      // (1-\theta) w_{\mathrm{old}}\right]\, \varphi_i\,
      // \varphi_j\,\mathrm{d}x$ using the desired quadrature formula.	   
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i) 
	  for (unsigned int j=0; j<dofs_per_cell; ++j) 
	    local_nl_matrix(i,j) += (std::cos(theta*new_data_values.at(q_point) +
					      (1-theta)*old_data_values.at(q_point)) *
				     fe_values.shape_value (i, q_point) *
				     fe_values.shape_value (j, q_point) *
				     fe_values.JxW (q_point));
      
      // Finally, we add up the contributions of the integrals over
      // the cells to the global integral.
      cell->get_dof_indices (local_dof_indices);

      for (unsigned int i=0; i<dofs_per_cell; ++i) 
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  nl_matrix.add(local_dof_indices[i], local_dof_indices[j], 
			local_nl_matrix(i,j));

      local_nl_matrix = 0;
    }
}

// @sect4{SineGordonProblem::compute_error}

// This function computes the norm of the difference between the
// computed (i.e., finite element) solution after time step
// ``timestep_number'' and the exact solution to see how well we are
// doing. There are several choices for norms available to us in the
// ``VectorTools'' class. We use the $L^2$ norm because it is a
// natural choice for our problem, since the solutions to the
// sine-Gordon equation have finite energy or, equivalently, are $L^2$
// functions. Given our weak formulation of the sine-Gordon equation,
// we are computing a solution $u\in H^1(\Omega)$, hence we could also
// use the $H^1$ norm to compute the error of the spatial
// discretization. For more information on the details behind this
// computation, the reader should refer to step-7.
template <int dim>
void SineGordonProblem<dim>::compute_error (const unsigned int timestep_number)
{
  ExactSolution<dim> exact_solution (1, time);
  
  Vector<double> difference_per_cell (triangulation.n_active_cells());
  VectorTools::integrate_difference (dof_handler,
				     solution,
				     exact_solution,
				     difference_per_cell,
				     QGauss<dim>(3),
				     VectorTools::L2_norm);
  fem_errors(timestep_number) = difference_per_cell.l2_norm();
  
  std::cout << "   The L^2 error in the solution is " 
	    << fem_errors(timestep_number) << "."
	    << std::endl;
}

// @sect4{SineGordonProblem::solve}

// This function uses the GMRES iterative solver on the linear system
// of equations resulting from the finite element spatial
// discretization of each iteration of Newton's method for the
// (nonlinear) first equation in the split formulation we derived in
// the Introduction. The solution to the system is, in fact, $\delta
// U^n_l$ so it is stored in ``d_solution'' and used to update
// ``solution'' in the ``run'' function. We cannot use the Conjugate
// Gradient solver because the nonlinear term in the Jacobian matrix
// results in a non-positive-definite matrix to invert. Moreover, we
// would like the solver to quit when the \e relative error is
// $10^{-12}$. This function is similar to its analogue in step-3 (and
// step-4); the only difference is the choice of iterative solver and
// the new stopping criterion.
template <int dim>
void SineGordonProblem<dim>::solve () 
{
  SolverControl solver_control (1000, 1e-12*system_rhs.l2_norm());
  SolverGMRES<> gmres (solver_control);
  d_solution = 0;
  gmres.solve (system_matrix, d_solution, system_rhs, PreconditionIdentity());

  std::cout << "   " << solver_control.last_step()
	    << " GMRES iterations needed to obtain convergence."
	    << std::endl;
}

// @sect4{SineGordonProblem::output_results}

// This function outputs the results to a file. It is almost identical
// to its counterpart in step-3 (and step-4). The only new thing is
// that the function now takes a parameter --- the time step number
// --- so that it can append it to the name of the file, which the
// current solution is output to.
template <int dim>
void SineGordonProblem<dim>::output_results (const unsigned int timestep_number)
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "u");
  data_out.build_patches ();

  std::ostringstream filename;
  filename << "solution-" << dim << "d-";

  // Pad the time step number in filename with zeros in the beginning
  // so that the files are ordered correctly in the shell and we can
  // generate a good animation using convert.
  if (timestep_number<10) 
    filename << "0000" << timestep_number;
  else if (timestep_number<100)
    filename << "000" << timestep_number;
  else if (timestep_number<1000)
    filename << "00" << timestep_number;
  else if (timestep_number<10000)
    filename << "0" << timestep_number;
  else 
    filename << timestep_number;
  
  // We output the solution at the desired times in ``vtk'' format, so
  // that we can use VisIt to make plots and/or animations.
  filename << ".vtk";
  std::ofstream output (filename.str().c_str());
  data_out.write_vtk (output);

  // We also store the current solution in our instantiation of a
  // ``DataOutStack'' object, so that we can make a space-time plot of
  // the solution.
  data_out_stack.new_parameter_value (time, time_step*output_timestep_skip);
  data_out_stack.attach_dof_handler (dof_handler);
  data_out_stack.add_data_vector (solution, "solution");
  data_out_stack.build_patches (1);
  data_out_stack.finish_parameter_value ();
}

// @sect4{SineGordonProblem::run}

// This function has the top-level control over everything: it runs
// the (outer) time-stepping loop, the (inner) nonlinear-solver loop,
// outputs the solution after each time step and calls the
// ``compute_error'' routine after each time step if an exact solution
// is known.
template <int dim>
void SineGordonProblem<dim>::run () 
{
  data_out_stack.declare_data_vector ("solution",
				      DataOutStack<dim>::dof_vector);

  std::cout << "Solving problem in " << dim << " space dimensions." 
	    << std::endl;
  
  make_grid_and_dofs ();

  // To aknowledge the initial condition, we must use the function
  // $u_0(x)$ to compute the zeroth time step solution $U^0$. Note
  // that when we create the ``InitialValues'' ``Function'' object, we
  // set its internal time variable to $t_0$, in case our initial
  // condition is a function of space and time evaluated at $t=t_0$.
  InitialValues<dim> initial_condition (1, time);

  // Then, in 2D and 3D, we produce $U^0$ by projecting $u_0(x)$ onto
  // the grid using ``VectorTools::project''. In 1D, however, we
  // obtain the zeroth time step solution by interpolating $u_0(x)$ at
  // the global degrees of freedom using
  // ``VectorTools::interpolate''. We must make an exception for the
  // 1D case because the projection algorithm computes integrals over
  // the boundary of the domain, which do not make sense in 1D, so we
  // cannot use it.
  if (dim == 1) 
    {
      VectorTools::interpolate (dof_handler, initial_condition, solution);
    }
  else 
    {
      ConstraintMatrix constraints;
      constraints.close();
      VectorTools::project (dof_handler, constraints, QGauss<dim>(3),
			    initial_condition, solution);
    }

  // For completeness, we output the zeroth time step to a file just
  // like any other other time step.
  output_results (0);

  // Now we perform the time stepping: at every time step we solve the
  // matrix equation(s) corresponding to the finite element
  // discretization of the problem, and then advance our solution
  // according to the time stepping formulas we discussed in the
  // Introduction.
  unsigned int timestep_number = 1;
  for (time+=time_step; time<=final_time; time+=time_step, ++timestep_number)
    {
      old_solution = solution;

      std::cout << std::endl
		<< " Time step #" << timestep_number << "; "
		<< "advancing to t = " << time << "." 
		<< std::endl;

      // First we must solve the nonlinear equation in the split
      // formulation via Newton's method --- i.e. solve for $\delta
      // U^n_l$ then compute $U^n_{l+1}$ and so on. The stopping
      // criterion is that $\|F_h(U^n_l)\|_2 \le 10^{-6}
      // \|F_h(U^n_0)\|_2$. When the loop below is done, we have (an
      // approximation of) $U^n$.
      double initial_rhs_norm = 0.;
      unsigned int nliter = 1;
      do 
	{
	  assemble_system ();
	  if (nliter == 1) initial_rhs_norm = system_rhs.l2_norm();
	  std::cout << "   [NLITER]"; 
	  solve ();
	  solution += d_solution;
	  d_solution = 0;
	  nliter++;
	} 
      while (system_rhs.l2_norm() > 1e-6 * initial_rhs_norm);
  
      // In the case of the explicit Euler time stepping scheme, we
      // must pick the time step to be quite small in order for the
      // scheme to be stable. Therefore, there are a lot of time steps
      // during which "nothing interesting happens" in the
      // solution. To improve overall efficiency --- in particular,
      // speed up the program and save disk space --- we only output
      // the solution after ``output_timestep_skip'' time steps have
      // been taken.
      if (timestep_number % output_timestep_skip == 0)
	  output_results (timestep_number);
      
      // Upon obtaining the solution to the problem at $t=t_n$, we
      // must update the auxilliary velocity variable $V^n$. However,
      // we do not compute and store $V^n$ since it is not a quantity
      // we use directly in the problem. Hence, for simplicity, we
      // update $MV^n$ directly using the second equation in the last
      // subsection of the Introduction.
      Vector<double> tmp_vector (solution.size());
      laplace_matrix.vmult (tmp_vector, solution);
      massmatxvel.add (-time_step*theta, tmp_vector);

      tmp_vector = 0;
      laplace_matrix.vmult (tmp_vector, old_solution);
      massmatxvel.add (-time_step*(1-theta), tmp_vector);
      
      tmp_vector = 0;
      compute_nl_term (old_solution, solution, tmp_vector);
      massmatxvel.add (-time_step, tmp_vector);

      // Before concluding the $n^{\mathrm{th}}$ time step, we compute
      // the error in the finite element solution at $t=t_n$ if the
      // exact solution to the problem being solved is known.
      if (exact_solution_known)
	compute_error (timestep_number);
    }

  // After the time stepping is complete, we report the maximum (over
  // all time steps) of the errors in the finite element solution if
  // the exact solution of the problem being solved is know.
  if (exact_solution_known)
    std::cout << " The maximum L^2 error in the solution was "
	      << fem_errors.linfty_norm() << "."
	      << std::endl << std::endl;

  // Finally, we output the sequence of solutions stored
  // ``data_out_stack'' to a file of the appropriate format.
  std::ostringstream filename;
  filename << "solution-" << dim << "d-" << "stacked" << ".vtk";
  std::ofstream output (filename.str().c_str());
  data_out_stack.write_vtk (output);
}

// @sect3{The ``main'' function}

// This is the main function of the program. It creates an object of
// top-level class and calls its principal function. Also, we supress
// some of the library output by setting ``deallog.depth_console'' to
// zero. Furthermore, if exceptions are thrown during the execution of
// the run method of the ``SineGordonProblem'' class, we catch and
// report them here. For more information about exceptions the reader
// should consult step-6.
int main () 
{
  try
    {
      deallog.depth_console (0);

      SineGordonProblem<2> sg_problem;
      sg_problem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      
      return 1;
    }
  catch (...) 
    {
       std::cerr << std::endl << std::endl
		 << "----------------------------------------------------"
		 << std::endl;
       std::cerr << "Unknown exception!" << std::endl
		 << "Aborting!" << std::endl
		 << "----------------------------------------------------"
		 << std::endl;
       return 1;
    }
  
  return 0;
}

/* $Id$ */

#include <grid/tria.h>
#include <grid/dof.h>
#include <grid/tria_accessor.h>
#include <grid/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary.h>
#include <grid/dof_constraints.h>
#include <basic/function.h>
#include <basic/data_io.h>
#include <basic/parameter_handler.h>
#include <fe/fe_lib.h>
#include <fe/quadrature_lib.h>
#include <numerics/base.h>
#include <numerics/assembler.h>
#include <numerics/error_estimator.h>

#include <map.h>
#include <fstream.h>
#include <cmath>
#include <string>
#include <cstdlib>




template <int dim>
class PoissonEquation :  public Equation<dim> {
  public:
    PoissonEquation (const Function<dim> &rhs) :
		    Equation<dim>(1),
		    right_hand_side (rhs)  {};

    virtual void assemble (dFMatrix            &cell_matrix,
			   dVector             &rhs,
			   const FEValues<dim> &fe_values,
			   const Triangulation<dim>::cell_iterator &cell) const;
    virtual void assemble (dFMatrix            &cell_matrix,
			   const FEValues<dim> &fe_values,
			   const Triangulation<dim>::cell_iterator &cell) const;
    virtual void assemble (dVector             &rhs,
			   const FEValues<dim> &fe_values,
			   const Triangulation<dim>::cell_iterator &cell) const;
  protected:
    const Function<dim> &right_hand_side;
};






template <int dim>
class PoissonProblem : public ProblemBase<dim>,
		       public MultipleParameterLoop::UserClass {
  public:
    enum RefineMode {
	  global, true_error, error_estimator
    };
    
    PoissonProblem ();

    void clear ();
    void create_new (const unsigned int);
    void declare_parameters (ParameterHandler &prm);
    void run (ParameterHandler &prm);
    void print_history (const ParameterHandler &prm,
			const RefineMode refine_mode) const;
    
  protected:
    Triangulation<dim> *tria;
    DoFHandler<dim>    *dof;
    
    Function<dim>      *rhs;
    Function<dim>      *solution_function;

    Boundary<dim>      *boundary;
    
    vector<double> l2_error, linfty_error;
    vector<double> h1_error, estimated_error;
    vector<int>    n_dofs;

    vector<double> laplacian;
};





template <int dim>
class Solution {
  public:

    class GaussShape : public Function<dim> {
      public:
	virtual double operator () (const Point<dim> &p) const;
	virtual Point<dim> gradient (const Point<dim> &p) const;
    };

    class Singular : public Function<dim> {
      public:
	virtual double operator () (const Point<dim> &p) const;
	virtual Point<dim> gradient (const Point<dim> &p) const;
    };
};




template <int dim>
class RHS {
  public:
    
				     /**
				      *	Right hand side constructed such that
				      * the exact solution is
				      * $x*y*exp(-(x**2+y**2)*40)$.
				      */
    class GaussShape : public Function<dim> {
      public:
	virtual double operator () (const Point<dim> &p) const;
    };

    				     /**
				      *	Right hand side constructed such that
				      * the exact solution is
				      * $r^{2/3} sin(2\phi)$.
				      */
    class Singular : public Function<dim> {
      public:
	virtual double operator () (const Point<dim> &p) const;
    };
};





double Solution<2>::GaussShape::operator () (const Point<2> &p) const {
  return p(0)*p(1)*exp(-40*p.square());
};


Point<2> Solution<2>::GaussShape::gradient (const Point<2> &p) const {
  return Point<2> ((1-80.*p(0)*p(0))*p(1)*exp(-40*p.square()),
		   (1-80.*p(1)*p(1))*p(0)*exp(-40*p.square()));
};



double Solution<2>::Singular::operator () (const Point<2> &p) const {
  return pow(p.square(), 1./3.);
};


Point<2> Solution<2>::Singular::gradient (const Point<2> &p) const {
  return 2./3.*pow(p.square(), -2./3.) * p;
};



double RHS<2>::GaussShape::operator () (const Point<2> &p) const {
  return (480.-6400.*p.square())*p(0)*p(1)*exp(-40.*p.square());
};



double RHS<2>::Singular::operator () (const Point<2> &p) const {
  return -4./9. * pow(p.square(), -2./3.);
};




  




void PoissonEquation<2>::assemble (dFMatrix            &cell_matrix,
				   dVector             &rhs,
				   const FEValues<2>   &fe_values,
				   const Triangulation<2>::cell_iterator &) const {
  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
    for (unsigned int i=0; i<fe_values.total_dofs; ++i) 
      {
	for (unsigned int j=0; j<fe_values.total_dofs; ++j)
	  cell_matrix(i,j) += (fe_values.shape_grad(i,point) *
			       fe_values.shape_grad(j,point)) *
			      fe_values.JxW(point);
	rhs(i) += fe_values.shape_value(i,point) *
		  right_hand_side(fe_values.quadrature_point(point)) *
		  fe_values.JxW(point);
      };
};



template <int dim>
void PoissonEquation<dim>::assemble (dFMatrix            &,
				     const FEValues<dim> &,
				     const Triangulation<dim>::cell_iterator &) const {
  Assert (false, ExcPureVirtualFunctionCalled());
};



template <int dim>
void PoissonEquation<dim>::assemble (dVector             &,
				     const FEValues<dim> &,
				     const Triangulation<dim>::cell_iterator &) const {
  Assert (false, ExcPureVirtualFunctionCalled());
};









template <int dim>
PoissonProblem<dim>::PoissonProblem () :
		tria(0), dof(0), rhs(0),
		solution_function(0), boundary(0) {};




template <int dim>
void PoissonProblem<dim>::clear () {
  if (tria != 0) {
    delete tria;
    tria = 0;
  };
  
  if (dof != 0) {
    delete dof;
    dof = 0;
  };

  if (rhs != 0) 
    {
      delete rhs;
      rhs = 0;
    };

  if (solution_function != 0) 
    {
      delete solution_function;
      solution_function = 0;
    };

  if (boundary != 0)
    {
      delete boundary;
      boundary = 0;
    };
  
  l2_error.clear ();
  linfty_error.clear ();
  h1_error.clear ();
  estimated_error.clear();
  n_dofs.clear ();
  
  ProblemBase<dim>::clear ();
};




template <int dim>
void PoissonProblem<dim>::create_new (const unsigned int) {
  clear ();

  tria = new Triangulation<dim>();
  dof = new DoFHandler<dim> (tria);
  set_tria_and_dof (tria, dof);
  boundary = new HyperBallBoundary<dim> ();
};



template <int dim>
void PoissonProblem<dim>::declare_parameters (ParameterHandler &prm) {
  prm.declare_entry ("Test case", "Gauss shape", "Gauss shape\\|Singular");
  prm.declare_entry ("Initial refinement", "2",
		     ParameterHandler::RegularExpressions::Integer);
  prm.declare_entry ("Refinement criterion", "estimated error",
		     "global\\|true error\\|estimated error");
  prm.declare_entry ("Refinement fraction", "0.3",
		     ParameterHandler::RegularExpressions::Double);
  prm.declare_entry ("Output base filename", "");
  prm.declare_entry ("Output format", "ucd"
		     "ucd\\|gnuplot");
};




template <int dim>
void PoissonProblem<dim>::run (ParameterHandler &prm) {
  cout << "======================================="
       << "=======================================" << endl
       << "===== Test case: " << prm.get ("Test case") << endl
       << "===== Doing computation with refinement criterion: ";
  RefineMode refine_mode;
  if (prm.get("Refinement criterion")=="global")
    refine_mode = global;
  else
    if (prm.get("Refinement criterion")=="true error")
      refine_mode = true_error;
    else
      if (prm.get("Refinement criterion")=="estimated error")
	refine_mode = error_estimator;
      else
	return;

  switch (refine_mode) 
    {
      case global:
	    cout << "global";
	    break;
      case true_error:
	    cout << "true error";
	    break;
      case error_estimator:
	    cout << "error estimator";
	    break;
    };

  const unsigned int start_level = prm.get_integer("Initial refinement");
  cout << endl
       << "======================================="
       << "=======================================" << endl;
  cout << "Making initial grid... " << endl;
  tria->set_boundary (boundary);
  tria->create_hyper_ball ();
  tria->refine_global (start_level);

  if (prm.get("Test case")=="Gauss shape")
    rhs             = new RHS<dim>::GaussShape();
  else
    if (prm.get("Test case")=="Singular")
      rhs             = new RHS<dim>::Singular();
  
  if (prm.get("Test case")=="Gauss shape")
    solution_function = new Solution<dim>::GaussShape ();
  else
    if (prm.get("Test case")=="Singular")
      solution_function = new Solution<dim>::Singular ();
  
  
  FELinear<dim>                   fe;
  PoissonEquation<dim>            equation (*rhs);
  QGauss3<dim>                    quadrature;

  unsigned int refine_step = 0;
  while (tria->n_active_cells() < 2000)
    {
      cout << "Refinement step " << refine_step
	   << ", using " << tria->n_active_cells() << " active cells on "
	   << tria->n_levels() << " levels."
	   << endl;
      cout << "    Distributing dofs... "; 
      dof->distribute_dofs (fe);
      cout << dof->n_dofs() << " degrees of freedom." << endl;
      n_dofs.push_back (dof->n_dofs());

      cout << "    Assembling matrices..." << endl;
      UpdateFlags update_flags = UpdateFlags(update_q_points  | update_gradients |
					     update_jacobians | update_JxW_values);
  
      ProblemBase<dim>::FunctionMap dirichlet_bc;
      dirichlet_bc[0] = solution_function;
      assemble (equation, quadrature, fe, update_flags, dirichlet_bc);

      cout << "    Solving..." << endl;
      solve ();

      dVector       l2_error_per_cell, linfty_error_per_cell, h1_error_per_cell;
      dVector       estimated_error_per_cell;
      dVector       laplacian_per_cell;
      QGauss3<dim>  q;
  
      cout << "    Calculating L2 error... ";
      integrate_difference (*solution_function, l2_error_per_cell, q,
			    fe, L2_norm);
      cout << l2_error_per_cell.l2_norm() << endl;
      l2_error.push_back (l2_error_per_cell.l2_norm());

      cout << "    Calculating L-infinity error... ";
      integrate_difference (*solution_function, linfty_error_per_cell, q,
			    fe, Linfty_norm);
      cout << linfty_error_per_cell.linfty_norm() << endl;
      linfty_error.push_back (linfty_error_per_cell.linfty_norm());

      cout << "    Calculating H1 error... ";
      integrate_difference (*solution_function, h1_error_per_cell, q, fe, H1_norm);
      cout << h1_error_per_cell.l2_norm() << endl;
      h1_error.push_back (h1_error_per_cell.l2_norm());

      cout << "    Estimating H1 error... ";
      KellyErrorEstimator<dim> ee;
      QSimpson<dim-1> eq;
      ee.estimate_error (*dof, eq, fe, *boundary,
			 KellyErrorEstimator<dim>::FunctionMap(),
			 solution,
			 estimated_error_per_cell);
      cout << estimated_error_per_cell.l2_norm() << endl;
      estimated_error.push_back (estimated_error_per_cell.l2_norm());

      laplacian_per_cell = estimated_error_per_cell;
      DoFHandler<dim>::active_cell_iterator cell = dof->begin_active(),
					    endc = dof->end();
      for (unsigned int i=0; cell!=endc; ++cell, ++i)
	laplacian_per_cell(i) /= (cell->diameter() * cell->diameter() / 24);
      cout << "    Computing second derivative maximum... "
	   <<  laplacian_per_cell.linfty_norm()
	   << endl;
      laplacian.push_back (laplacian_per_cell.linfty_norm());
      
      dVector l2_error_per_dof, linfty_error_per_dof;
      dVector h1_error_per_dof, estimated_error_per_dof;
      dVector laplacian_per_dof;
      dof->distribute_cell_to_dof_vector (l2_error_per_cell, l2_error_per_dof);
      dof->distribute_cell_to_dof_vector (linfty_error_per_cell,
					  linfty_error_per_dof);
      dof->distribute_cell_to_dof_vector (h1_error_per_cell, h1_error_per_dof);
      dof->distribute_cell_to_dof_vector (estimated_error_per_cell,
					  estimated_error_per_dof);
      dof->distribute_cell_to_dof_vector (laplacian_per_cell,
					  laplacian_per_dof);
  
  
      DataOut<dim> out;
      fill_data (out);
      out.add_data_vector (l2_error_per_dof, "L2-Error");
      out.add_data_vector (linfty_error_per_dof, "Linfty-Error");
      out.add_data_vector (h1_error_per_dof, "H1-Error");
      out.add_data_vector (estimated_error_per_dof, "Estimated Error");
      out.add_data_vector (laplacian_per_dof, "Second derivative pointwise");
      String filename = prm.get ("Output base filename");
      switch (refine_mode) 
	{
	  case global:
		filename += "global.";
		break;
	  case true_error:
		filename += "true_error.";
		break;
	  case error_estimator:
		filename += "estimated_error.";
		break;
	};
      filename += ('0'+(start_level+refine_step)/10);
      filename += ('0'+(start_level+refine_step)%10);

      if (prm.get("Output format")=="ucd")
	filename += ".inp";
      else
	if (prm.get("Output format")=="gnuplot")
	  filename += ".gnuplot";
      
      cout << "    Writing error plots to <" << filename << ">..." << endl;
      ofstream outfile(filename);
      if (prm.get("Output format")=="ucd")      
	out.write_ucd (outfile);
      else
	if (prm.get("Output format")=="gnuplot")
	  out.write_gnuplot (outfile);
      
      outfile.close();

      cout << "    Refining triangulation...";
      switch (refine_mode) 
	{
	  case global:
		tria->refine_global (1);
		break;
	  case true_error:
		tria->refine_fixed_number (h1_error_per_cell,
					   prm.get_double("Refinement fraction"));
		tria->execute_refinement ();
		break;
	  case error_estimator:
		tria->refine_fixed_number (estimated_error_per_cell,
					   prm.get_double("Refinement fraction"));
		tria->execute_refinement ();
		break;
	};
      cout << endl << endl;
      ++refine_step;
    };
  
  print_history (prm, refine_mode);
  cout << endl << endl << endl;
};


template <int dim>
void PoissonProblem<dim>::print_history (const ParameterHandler &prm,
					 const RefineMode refine_mode) const {
  String filename(prm.get("Output base filename"));
  filename += "history.";
  switch (refine_mode) 
    {
      case global:
	    filename += "global.";
	    break;
      case true_error:
	    filename += "true_error.";
	    break;
      case error_estimator:
	    filename += "estimated_error.";
	    break;
    };
  filename += "gnuplot";
  
  cout << endl << "Printing convergence history to <" << filename << ">..."
       << endl;
  ofstream out(filename);
  out << "# n_dofs    l2_error linfty_error "
      << "h1_error estimated_error "
      << "laplacian"
      << endl;
  for (unsigned int i=0; i<n_dofs.size(); ++i)
    out << n_dofs[i]
	<< "    "
	<< l2_error[i] << "  "
	<< linfty_error[i] << "  "
	<< h1_error[i] << "  "
	<< estimated_error[i] << "  "
	<< laplacian[i]
	<< endl;

  double average_l2=0,
     average_linfty=0,
	 average_h1=0,
	average_est=0;
  
  for (unsigned int i=1; i<n_dofs.size(); ++i) 
    {
      average_l2 += l2_error[i]/l2_error[i-1];
      average_linfty += linfty_error[i]/linfty_error[i-1];
      average_h1 += h1_error[i]/h1_error[i-1];
      average_est += estimated_error[i]/estimated_error[i-1];
    };

  average_l2 /= (l2_error.size()-1);
  average_linfty /= (l2_error.size()-1);
  average_h1 /= (l2_error.size()-1);
  average_est /= (l2_error.size()-1);

  cout << "Average error reduction rates for h->h/2:" << endl;
  cout << "    L2 error         : " << 1./average_l2 << endl
       << "    Linfty error     : " << 1./average_linfty << endl
       << "    H1 error         : " << 1./average_h1 << endl
       << "    Estimated error  : " << 1./average_est << endl;
};




int main (int argc, char **argv) {
  if (argc!=2) 
    {
      cout << "Usage: error-estimation parameterfile" << endl << endl;
      return 1;
    };

  PoissonProblem<2> poisson;
  MultipleParameterLoop input_data;

  poisson.declare_parameters(input_data);
  input_data.read_input (argv[1]);
  input_data.loop (poisson);
  
  return 0;
};




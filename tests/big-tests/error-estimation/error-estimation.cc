/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include <base/function.h>
#include <base/parameter_handler.h>
#include <base/quadrature_lib.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_boundary_lib.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_constraints.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_tools.h>
#include <fe/fe_lib.lagrange.h>
#include <numerics/data_out.h>
#include <numerics/base.h>
#include <numerics/assembler.h>
#include <numerics/vectors.h>
#include <numerics/error_estimator.h>
#include <numerics/solution_transfer.h>
#include <lac/vector.h>

#include <map>
#include <fstream>
#include <cmath>
#include <string>
#include <cstdlib>




template <int dim>
class PoissonEquation :  public Equation<dim> {
  public:
    PoissonEquation (const Function<dim> &rhs) :
		    Equation<dim>(1),
		    use_coefficient(false),
		    right_hand_side (rhs),
		    coefficient (default_coefficient) {};

    PoissonEquation (const Function<dim> &rhs,
		     const Function<dim> &coefficient ) :
		    Equation<dim>(1),
		    use_coefficient(true),
		    right_hand_side (rhs),
		    coefficient (coefficient) {};

    virtual void assemble (FullMatrix<double>  &cell_matrix,
			   Vector<double>      &rhs,
			   const FEValues<dim> &fe_values,
			   const DoFHandler<dim>::cell_iterator &cell) const;
    virtual void assemble (FullMatrix<double>  &cell_matrix,
			   const FEValues<dim> &fe_values,
			   const DoFHandler<dim>::cell_iterator &cell) const;
    virtual void assemble (Vector<double>      &rhs,
			   const FEValues<dim> &fe_values,
			   const DoFHandler<dim>::cell_iterator &cell) const;
  protected:
    const bool           use_coefficient;
    const Function<dim> &right_hand_side;
    const Function<dim> &coefficient;

    static const ConstantFunction<dim> default_coefficient;
};


const ConstantFunction<2> PoissonEquation<2>::default_coefficient(1);





template <int dim>
class PoissonProblem : public ProblemBase<dim>, public MultipleParameterLoop::UserClass {
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
    Function<dim>      *coefficient;
			   
    Boundary<dim>      *boundary;
    
    vector<double> l2_error, linfty_error;
    vector<double> h1_error, estimated_error;
    vector<int>    n_dofs;
};





template <int dim>
class Solution {
  public:

    class GaussShape : public Function<dim> {
      public:
	virtual double value (const Point<dim> &p,
			      const unsigned int component) const;
	virtual Tensor<1,dim> gradient (const Point<dim> &p,
					const unsigned int component) const;
    };

    class Singular : public Function<dim> {
      public:
	virtual double value (const Point<dim> &p,
			      const unsigned int component) const;
	virtual Tensor<1,dim> gradient (const Point<dim> &p,
					const unsigned int component) const;
    };

    class Kink : public Function<dim> {
      public:
	class Coefficient : public Function<dim> {
	  public:
	    virtual double value (const Point<dim> &p,
				  const unsigned int component) const;
	};
	
	virtual double value (const Point<dim> &p,
			      const unsigned int component) const;
	virtual Tensor<1,dim> gradient (const Point<dim> &p,
					const unsigned int component) const;
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
	virtual double value (const Point<dim> &p,
			      const unsigned int component) const;
    };

    				     /**
				      *	Right hand side constructed such that
				      * the exact solution is
				      * $r^{2/3}$.
				      */
    class Singular : public Function<dim> {
      public:
	virtual double value (const Point<dim> &p,
			      const unsigned int component) const;
    };

    				     /**
				      *	Right hand side constructed such that
				      * the exact solution is
				      * $(1+4\theta(f))*f$ with
				      * $f=y-x**2$.
				      */
    class Kink : public Function<dim> {
      public:
	virtual double value (const Point<dim> &p,
			      const unsigned int component) const;
    };
};




template <>
double Solution<2>::GaussShape::value (const Point<2> &p,
				       const unsigned int) const {
  return p(0)*p(1)*exp(-40*p.square());
};


template <>
Tensor<1,2> Solution<2>::GaussShape::gradient (const Point<2> &p,
					       const unsigned int) const {
  return Point<2> ((1-80.*p(0)*p(0))*p(1)*exp(-40*p.square()),
		   (1-80.*p(1)*p(1))*p(0)*exp(-40*p.square()));
};



template <>
double Solution<2>::Singular::value (const Point<2> &p,
				     const unsigned int) const {
  return pow(p.square(), 1./3.);
};


template <>
Tensor<1,2> Solution<2>::Singular::gradient (const Point<2> &p,
					     const unsigned int) const {
  return 2./3.*pow(p.square(), -2./3.) * p;
};




inline double theta(const double x) {
  return (x>0 ? 1 : 0);
};



template <>
double Solution<2>::Kink::value (const Point<2> &p,
				 const unsigned int) const {
  const double s = p(1)-p(0)*p(0);
  return (1+4*theta(s))*s;
};


template <>
Tensor<1,2> Solution<2>::Kink::gradient (const Point<2> &p,
					 const unsigned int) const {
  const double s = p(1)-p(0)*p(0);
  return (1+4*theta(s))*Point<2>(-2*p(0),1);
};


template <>
double Solution<2>::Kink::Coefficient::value (const Point<2> &p,
					      const unsigned int) const {
  const double s = p(1)-p(0)*p(0);
  return 1./(1.+4.*theta(s));
};



template <>
double RHS<2>::GaussShape::value (const Point<2> &p,
				  const unsigned int) const {
  return (480.-6400.*p.square())*p(0)*p(1)*exp(-40.*p.square());
};


template <>
double RHS<2>::Singular::value (const Point<2> &p,
				const unsigned int) const {
  return -4./9. * pow(p.square(), -2./3.);
};


template <>
double RHS<2>::Kink::value (const Point<2> &,
			    const unsigned int) const {
  return 2;
};




  



template <>
void PoissonEquation<2>::assemble (FullMatrix<double>  &cell_matrix,
				   Vector<double>      &rhs,
				   const FEValues<2>   &fe_values,
				   const DoFHandler<2>::cell_iterator &) const {
  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point) 
    {
      const double c = (use_coefficient ?
			coefficient.value(fe_values.quadrature_point(point)) :
			1);
      for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i) 
	{
	  for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
	    cell_matrix(i,j) += (fe_values.shape_grad(i,point) *
				 fe_values.shape_grad(j,point)) *
				fe_values.JxW(point) *
				c;
	  rhs(i) += fe_values.shape_value(i,point) *
		    right_hand_side.value(fe_values.quadrature_point(point)) *
		    fe_values.JxW(point);
	};
    };
};



template <int dim>
void PoissonEquation<dim>::assemble (FullMatrix<double>  &,
				     const FEValues<dim> &,
				     const DoFHandler<dim>::cell_iterator &) const {
  Assert (false, ExcPureVirtualFunctionCalled());
};



template <int dim>
void PoissonEquation<dim>::assemble (Vector<double>      &,
				     const FEValues<dim> &,
				     const DoFHandler<dim>::cell_iterator &) const {
  Assert (false, ExcPureVirtualFunctionCalled());
};









template <int dim>
PoissonProblem<dim>::PoissonProblem () :
		tria(0), dof(0), rhs(0),
		solution_function(0), coefficient(0),
		boundary(0) {};




template <int dim>
void PoissonProblem<dim>::clear () {
  if (dof != 0)               { delete dof;               dof = 0;               };
  if (tria != 0)              { delete tria;              tria = 0;              };
  if (rhs != 0)               { delete rhs;               rhs = 0;               };
  if (solution_function != 0) { delete solution_function; solution_function = 0; };
  if (coefficient != 0)       { delete coefficient;       coefficient = 0;       };
  if (boundary != 0)          { delete boundary;          boundary = 0;          };

  				   // make it known to the underlying
				   // ProblemBase that tria and dof
				   // are already deleted
  set_tria_and_dof (tria, dof);

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
  prm.declare_entry ("Test case", "Gauss shape",
		     Patterns::Selection("Gauss shape|Singular|Kink"));
  prm.declare_entry ("Initial refinement", "2",
		     Patterns::Integer());
  prm.declare_entry ("Refinement criterion", "estimated error",
		     Patterns::Selection("global|true error|estimated error"));
  prm.declare_entry ("Refinement fraction", "0.3",
		     Patterns::Double());
  prm.declare_entry ("Coarsening fraction", "0.1",
		     Patterns::Double());
  prm.declare_entry ("Maximum cells", "3000",
		     Patterns::Integer());
  prm.declare_entry ("Output base filename", "");
  prm.declare_entry ("Output format", "ucd",
		     Patterns::Selection("ucd|gnuplot"));
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

  cout << endl
       << "======================================="
       << "=======================================" << endl;
  cout << "Making initial grid... " << endl;
  const unsigned int start_level(prm.get_integer("Initial refinement"));
  tria->set_boundary (0, *boundary);
  GridGenerator::hyper_ball (*tria);
  tria->refine_global (start_level);

  if (prm.get("Test case")=="Gauss shape")
    rhs             = new RHS<dim>::GaussShape();
  else
    if (prm.get("Test case")=="Singular")
      rhs             = new RHS<dim>::Singular();
    else
      if (prm.get("Test case")=="Kink")
	rhs             = new RHS<dim>::Kink();
  
  if (prm.get("Test case")=="Gauss shape")
    solution_function = new Solution<dim>::GaussShape ();
  else
    if (prm.get("Test case")=="Singular")
      solution_function = new Solution<dim>::Singular ();
    else
      if (prm.get("Test case")=="Kink")
	solution_function = new Solution<dim>::Kink ();
  
  
  FEQ1<dim>         fe;
  QGauss3<dim>          quadrature;
  PoissonEquation<dim> *equation;
  
  static Solution<dim>::Kink::Coefficient kink_coefficient;
  if (prm.get("Test case")=="Kink")
    equation = new PoissonEquation<dim>(*rhs, kink_coefficient);
  else
    equation = new PoissonEquation<dim>(*rhs);

  SolutionTransfer<dim,double> solution_transfer (*dof_handler);

  unsigned int refine_step = 0;
  const unsigned int max_cells = prm.get_integer("Maximum cells");
  while (tria->n_active_cells() < max_cells)
    {
      Vector<double> old_solution = solution;
      cout << "Refinement step " << refine_step
	   << ", using " << tria->n_active_cells() << " active cells on "
	   << tria->n_levels() << " levels."
	   << endl;
      cout << "    Distributing dofs... "; 
      dof->distribute_dofs (fe);
      cout << dof->n_dofs() << " degrees of freedom." << endl;
      n_dofs.push_back (dof->n_dofs());

      cout << "    Assembling matrices..." << endl;
      UpdateFlags update_flags = UpdateFlags(update_values | update_q_points  |
					     update_gradients | update_JxW_values);
  
      ProblemBase<dim>::FunctionMap dirichlet_bc;
      dirichlet_bc[0] = solution_function;
      assemble (*equation, quadrature, update_flags, dirichlet_bc);

				       // if we have an old solution lying
				       // around, use it to preset the solution
				       // vector. this reduced the quired
				       // number of iterations by about
				       // 10 per cent
      if (refine_step != 0)
	{
	  solution.reinit (dof_handler->n_dofs());
	  solution_transfer.interpolate (old_solution, solution);

					   // if you don't want to preset
					   // the solution vector,
					   // uncomment the following
					   // line and comment out the
					   // preceding one
//        solution.reinit (dof_handler->n_dofs());
	  
	  solution_transfer.clear ();
	};

      cout << "    Solving..." << endl;

      solve ();


      Vector<float>       l2_error_per_cell, linfty_error_per_cell, h1_error_per_cell;
      Vector<float>       estimated_error_per_cell;
      QGauss3<dim>  q;
  
      cout << "    Calculating L2 error... ";
      VectorTools::integrate_difference (*dof_handler,
					      solution, *solution_function,
					      l2_error_per_cell, q,
					      L2_norm);
      cout << l2_error_per_cell.l2_norm() << endl;
      l2_error.push_back (l2_error_per_cell.l2_norm());

      cout << "    Calculating L-infinity error... ";
      VectorTools::integrate_difference (*dof_handler,
					      solution, *solution_function,
					      linfty_error_per_cell, q,
					      Linfty_norm);
      cout << linfty_error_per_cell.linfty_norm() << endl;
      linfty_error.push_back (linfty_error_per_cell.linfty_norm());

      cout << "    Calculating H1 error... ";
      VectorTools::integrate_difference (*dof_handler,
					      solution, *solution_function,
					      h1_error_per_cell, q, 
					      H1_norm);
      cout << h1_error_per_cell.l2_norm() << endl;
      h1_error.push_back (h1_error_per_cell.l2_norm());

      cout << "    Estimating H1 error... ";

      QSimpson<dim-1> eq;
      KellyErrorEstimator<dim>::estimate (*dof, eq,
					  KellyErrorEstimator<dim>::FunctionMap(),
					  solution,
					  estimated_error_per_cell,
					  vector<bool>(), // all components
					  ((prm.get("Test case")=="Kink") ?
					   &kink_coefficient : 0 ));
      cout << estimated_error_per_cell.l2_norm() << endl;
      estimated_error.push_back (estimated_error_per_cell.l2_norm());

      Vector<double> l2_error_per_dof(dof->n_dofs()), linfty_error_per_dof(dof->n_dofs());
      Vector<double> h1_error_per_dof(dof->n_dofs()), estimated_error_per_dof(dof->n_dofs());
      Vector<double> error_ratio (dof->n_dofs());
      DoFTools::distribute_cell_to_dof_vector (*dof, l2_error_per_cell, l2_error_per_dof);
      DoFTools::distribute_cell_to_dof_vector (*dof, linfty_error_per_cell,
					  linfty_error_per_dof);
      DoFTools::distribute_cell_to_dof_vector (*dof, h1_error_per_cell, h1_error_per_dof);
      DoFTools::distribute_cell_to_dof_vector (*dof, estimated_error_per_cell,
					  estimated_error_per_dof);
      error_ratio.ratio (h1_error_per_dof, estimated_error_per_dof);
  
      DataOut<dim> out;
      fill_data (out);
      out.add_data_vector (l2_error_per_dof, "L2_Error");
      out.add_data_vector (linfty_error_per_dof, "Linfty_Error");
      out.add_data_vector (h1_error_per_dof, "H1_Error");
      out.add_data_vector (estimated_error_per_dof, "Estimated_Error");
      out.add_data_vector (error_ratio, "Ratio_True_to_Estimated_Error");
      string filename = prm.get ("Output base filename");
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
      ofstream outfile(filename.c_str());
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
		tria->set_all_refine_flags ();
		break;
	  case true_error:
		tria->refine_and_coarsen_fixed_number (h1_error_per_cell,
						       prm.get_double("Refinement fraction"),
						       prm.get_double("Coarsening fraction"));
		break;
	  case error_estimator:
		tria->refine_and_coarsen_fixed_number (estimated_error_per_cell,
						       prm.get_double("Refinement fraction"),
						       prm.get_double("Coarsening fraction"));
		break;
	};

      tria->prepare_coarsening_and_refinement ();
      solution_transfer.prepare_for_coarsening_and_refinement (solution);
      tria->execute_coarsening_and_refinement ();
      
      cout << endl << endl;
      ++refine_step;
    };

  string filename = prm.get ("Output base filename");
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

  cout << endl;
  
  filename += "finest_mesh.gnuplot";
  cout << "    Writing finest grid to <" << filename << ">... " << endl;
  ofstream finest_mesh (filename.c_str());
  GridOut().write_gnuplot (*tria, finest_mesh);
  finest_mesh.close();

  print_history (prm, refine_mode);
  cout << endl << endl << endl;

  dof->clear ();
  delete equation;
};


template <int dim>
void PoissonProblem<dim>::print_history (const ParameterHandler &prm,
					 const RefineMode refine_mode) const {
  string filename(prm.get("Output base filename"));
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
  ofstream out(filename.c_str());
  out << "# n_dofs    l2_error linfty_error "
      << "h1_error estimated_error"
      << endl;
  for (unsigned int i=0; i<n_dofs.size(); ++i)
    out << n_dofs[i]
	<< "    "
	<< l2_error[i] << "  "
	<< linfty_error[i] << "  "
	<< h1_error[i] << "  "
	<< estimated_error[i] << "  "
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




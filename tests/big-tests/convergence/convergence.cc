/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include <grid/tria.h>
#include <grid/dof.h>
#include <grid/tria_accessor.h>
#include <grid/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary.h>
#include <grid/dof_constraints.h>
#include <grid/grid_generator.h>
#include <base/function.h>
#include <basic/data_out.h>
#include <fe/fe_lib.lagrange.h>
#include <fe/fe_lib.criss_cross.h>
#include <base/quadrature_lib.h>
#include <numerics/base.h>
#include <numerics/assembler.h>
#include <numerics/vectors.h>
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
		    right_hand_side (rhs)  {};

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
    const Function<dim> &right_hand_side;
};






template <int dim>
class PoissonProblem : public ProblemBase<dim> {
  public:
    PoissonProblem (unsigned int order);

    void clear ();
    void create_new ();
    int run (unsigned int level);
    void print_history (string filename) const;
    
  protected:
    Triangulation<dim> *tria;
    DoFHandler<dim>    *dof;
    
    Function<dim>      *rhs;
    Function<dim>      *boundary_values;

    vector<double> l1_error, l2_error, linfty_error, h1_seminorm_error, h1_error;
    vector<int>    n_dofs;

    unsigned int        order;
};





/**
  Right hand side constructed such that the exact solution is
  $sin(2 pi x) + sin(2 pi y)$
  */
template <int dim>
class RHSPoly : public Function<dim> {
  public:
    				     /**
				      * Return the value of the function
				      * at the given point.
				      */
    virtual double operator () (const Point<dim> &p) const;
};



template <int dim>
class Solution : public Function<dim> {
  public:
    				     /**
				      * Return the value of the function
				      * at the given point.
				      */
    virtual double operator () (const Point<dim> &p) const;
				     /**
				      * Return the gradient of the function
				      * at the given point.
				      */
    virtual Tensor<1,dim> gradient (const Point<dim> &p) const;
};




template <>
double RHSPoly<2>::operator () (const Point<2> &p) const {
  const double x = p(0),
	       y = p(1);
  const double pi= 3.1415926536;
  return 4*pi*pi*(sin(2*pi*x)+sin(2*pi*y));
};



template <>
double Solution<2>::operator () (const Point<2> &p) const {
  const double x = p(0),
	       y = p(1);
  const double pi= 3.1415926536;
  return sin(2*pi*x)+sin(2*pi*y);
};


template <>
Tensor<1,2> Solution<2>::gradient (const Point<2> &p) const {
  const double x = p(0),
	       y = p(1);
  const double pi= 3.1415926536;
  return Point<2> (2*pi*cos(2*pi*x),
		   2*pi*cos(2*pi*y));
};

  



template <>
void PoissonEquation<2>::assemble (FullMatrix<double>  &cell_matrix,
				   Vector<double>      &rhs,
				   const FEValues<2>   &fe_values,
				   const DoFHandler<2>::cell_iterator &) const {
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
PoissonProblem<dim>::PoissonProblem (unsigned int order) :
		tria(0), dof(0), rhs(0),
		boundary_values(0), order(order) {};




template <int dim>
void PoissonProblem<dim>::clear () {
  if (dof != 0) {
    delete dof;
    dof = 0;
  };

  if (tria != 0) {
    delete tria;
    tria = 0;
  };
  

  				   // make it known to the underlying
				   // ProblemBase that tria and dof
				   // are already deleted
  set_tria_and_dof (tria, dof);

  if (rhs != 0) 
    {
      delete rhs;
      rhs = 0;
    };

  if (boundary_values != 0) 
    {
      delete boundary_values;
      boundary_values = 0;
    };

  ProblemBase<dim>::clear ();
};




template <int dim>
void PoissonProblem<dim>::create_new () {
  clear ();
  
  tria = new Triangulation<dim>();
  dof = new DoFHandler<dim> (tria);
  set_tria_and_dof (tria, dof);
};




template <int dim>
int PoissonProblem<dim>::run (const unsigned int level) {
  create_new ();
  
  cout << "Refinement level = " << level
       << ", using elements of type <";
  switch (order)
    {
      case 0:
	    cout << "criss-cross";
	    break;
      default:
	    cout << "Lagrange-" << order;
	    break;
    };
  cout << ">" << endl;
  
  cout << "    Making grid... ";
  GridGenerator::hyper_ball (*tria);
  HyperBallBoundary<dim> boundary_description;
  tria->set_boundary (0, boundary_description);
  tria->begin_active()->set_refine_flag();
  (++(++(tria->begin_active())))->set_refine_flag();
  tria->execute_coarsening_and_refinement ();
  tria->refine_global (level);
  cout << tria->n_active_cells() << " active cells." << endl;

  rhs             = new RHSPoly<dim>();
  boundary_values = new Solution<dim> ();
  

  FiniteElement<dim>   *fe;
  PoissonEquation<dim>  equation (*rhs);
  Quadrature<dim>      *quadrature;
  Quadrature<dim-1>    *boundary_quadrature;
  switch (order) {
    case 0:
	  fe         = new FECrissCross<dim>();
	  quadrature = new QCrissCross1<dim>();
	  boundary_quadrature = new QGauss2<dim-1>();
	  break;
    case 1:
	  fe         = new FEQ1<dim>();
	  quadrature = new QGauss3<dim>();
	  boundary_quadrature = new QGauss2<dim-1>();
	  break;
    case 2:
	  fe         = new FEQ2<dim>();
	  quadrature = new QGauss4<dim>();
	  boundary_quadrature = new QGauss3<dim-1>();
	  break;
    case 3:
	  fe         = new FEQ3<dim>();
	  quadrature = new QGauss5<dim>();
	  boundary_quadrature = new QGauss4<dim-1>();
	  break;
    case 4:
	  fe         = new FEQ4<dim>();
	  quadrature = new QGauss6<dim>();
	  boundary_quadrature = new QGauss5<dim-1>();
	  break;
    default:
	  return 100000;
  };
  
  cout << "    Distributing dofs... "; 
  dof->distribute_dofs (*fe);
  cout << dof->n_dofs() << " degrees of freedom." << endl;
  n_dofs.push_back (dof->n_dofs());

  cout << "    Assembling matrices..." << endl;
  UpdateFlags update_flags = UpdateFlags(update_q_points  | update_gradients |
					 update_JxW_values);
  
  ProblemBase<dim>::FunctionMap dirichlet_bc;
  dirichlet_bc[0] = boundary_values;
  assemble (equation, *quadrature, update_flags, dirichlet_bc);

  cout << "    Solving..." << endl;
  solve ();

  Solution<dim> sol;
  Vector<float>       l1_error_per_cell, l2_error_per_cell, linfty_error_per_cell;
  Vector<float>       h1_seminorm_error_per_cell, h1_error_per_cell;
  
  cout << "    Calculating L1 error... ";
  VectorTools<dim>::integrate_difference (*dof_handler,
					  solution, sol,
					  l1_error_per_cell,
					  *quadrature, L1_norm);
  cout << l1_error_per_cell.l1_norm() << endl;
  l1_error.push_back (l1_error_per_cell.l1_norm());

  cout << "    Calculating L2 error... ";
  VectorTools<dim>::integrate_difference (*dof_handler,
					  solution, sol,
					  l2_error_per_cell,
					  *quadrature, L2_norm);
  cout << l2_error_per_cell.l2_norm() << endl;
  l2_error.push_back (l2_error_per_cell.l2_norm());

  cout << "    Calculating L-infinity error... ";
  VectorTools<dim>::integrate_difference (*dof_handler,
					  solution, sol,
					  linfty_error_per_cell,
					  *quadrature, Linfty_norm);
  cout << linfty_error_per_cell.linfty_norm() << endl;
  linfty_error.push_back (linfty_error_per_cell.linfty_norm());
  
  cout << "    Calculating H1-seminorm error... ";
  VectorTools<dim>::integrate_difference (*dof_handler,
					  solution, sol,
					  h1_seminorm_error_per_cell,
					  *quadrature, H1_seminorm);
  cout << h1_seminorm_error_per_cell.l2_norm() << endl;
  h1_seminorm_error.push_back (h1_seminorm_error_per_cell.l2_norm());

  cout << "    Calculating H1 error... ";
  VectorTools<dim>::integrate_difference (*dof_handler,
					  solution, sol,
					  h1_error_per_cell,
					  *quadrature, H1_norm);
  cout << h1_error_per_cell.l2_norm() << endl;
  h1_error.push_back (h1_error_per_cell.l2_norm());

  if (dof->n_dofs()<=5000) 
    {
      Vector<double> l1_error_per_dof, l2_error_per_dof, linfty_error_per_dof;
      Vector<double> h1_seminorm_error_per_dof, h1_error_per_dof;
      dof->distribute_cell_to_dof_vector (l1_error_per_cell, l1_error_per_dof);
      dof->distribute_cell_to_dof_vector (l2_error_per_cell, l2_error_per_dof);
      dof->distribute_cell_to_dof_vector (linfty_error_per_cell,
					  linfty_error_per_dof);
      dof->distribute_cell_to_dof_vector (h1_seminorm_error_per_cell,
					  h1_seminorm_error_per_dof);
      dof->distribute_cell_to_dof_vector (h1_error_per_cell, h1_error_per_dof);

//       Vector<double> projected_solution;
//       ConstraintMatrix constraints;
//       constraints.close ();
//       VectorTools<dim>::project (*dof, constraints, *fe,
// 				 StraightBoundary<dim>(), *quadrature, 
// 				 sol, projected_solution, false,
// 				 *boundary_quadrature);
//       cout << "    Calculating L2 error of projected solution... ";
//       VectorTools<dim>::integrate_difference (*dof_handler,
// 					      projected_solution, sol,
// 					      l2_error_per_cell,
// 					      *quadrature, *fe, L2_norm);
//       cout << l2_error_per_cell.l2_norm() << endl;


      string filename;
      filename = ('0'+order);
      filename += ".";
      filename += ('0'+level);
      filename += ".ucd";
      cout << "    Writing error plots to <" << filename << ">..." << endl;
      
      DataOut<dim> out;
      ofstream o(filename.c_str());
      fill_data (out);
      out.add_data_vector (l1_error_per_dof, "L1_Error");
      out.add_data_vector (l2_error_per_dof, "L2_Error");
      out.add_data_vector (linfty_error_per_dof, "Linfty_Error");
      out.add_data_vector (h1_seminorm_error_per_dof, "H1_seminorm_Error");
      out.add_data_vector (h1_error_per_dof, "H1_Error");
      out.write_ucd (o);
      o.close ();
    }
  else
    cout << "    Not writing error as grid." << endl;
  
  cout << endl;

  const unsigned int n_dofs = dof->n_dofs();
				   // release the lock that the dof object
				   // has to the finite element object
  dof->clear ();
  tria->set_boundary (0);
  
  delete fe;
  delete quadrature;
  delete boundary_quadrature;
  
  return n_dofs;
};


template <int dim>
void PoissonProblem<dim>::print_history (string filename) const {
  ofstream out(filename.c_str());
  out << "# n_dofs    l1_error l2_error linfty_error h1_seminorm_error h1_error"
      << endl;
  for (unsigned int i=0; i<n_dofs.size(); ++i)
    out << n_dofs[i]
	<< "    "
	<< l1_error[i] << "  "
	<< l2_error[i] << "  "
	<< linfty_error[i] << "  "
	<< h1_seminorm_error[i] << "  "
	<< h1_error[i] << endl;

  double average_l1=0,
	 average_l2=0,
     average_linfty=0,
    average_h1_semi=0,
	 average_h1=0;
  for (unsigned int i=1; i<n_dofs.size(); ++i) 
    {
      average_l1 += l1_error[i]/l1_error[i-1];
      average_l2 += l2_error[i]/l2_error[i-1];
      average_linfty += linfty_error[i]/linfty_error[i-1];
      average_h1_semi += h1_seminorm_error[i]/h1_seminorm_error[i-1];
      average_h1 += h1_error[i]/h1_error[i-1];
    };

  average_l1 /= (l1_error.size()-1);
  average_l2 /= (l1_error.size()-1);
  average_linfty /= (l1_error.size()-1);
  average_h1_semi /= (l1_error.size()-1);
  average_h1 /= (l1_error.size()-1);

  cout << "==========================================================\n";
  cout << "Average error reduction rates for h->h/2:" << endl;
  cout << "    L1 error         : " << 1./average_l1 << endl
       << "    L2 error         : " << 1./average_l2 << endl
       << "    Linfty error     : " << 1./average_linfty << endl
       << "    H1 seminorm error: " << 1./average_h1_semi << endl
       << "    H1 error         : " << 1./average_h1 << endl;
  cout << "==========================================================\n";
  cout << "==========================================================\n";
};




int main () {
  for (unsigned int order=0; order<5; ++order) 
    {
      PoissonProblem<2> problem (order);
      
      unsigned int level=0;
      unsigned int n_dofs;
      do
	n_dofs = problem.run (level++);
      while (n_dofs<25000);

      string filename;
      switch (order) 
	{
	  case 0:
		filename = "criss_cross";
		break;
	  case 1:
		filename = "linear";
		break;
	  case 2:
		filename = "quadratic";
		break;
	  case 3:
		filename = "cubic";
		break;
	  case 4:
		filename = "quartic";
		break;
	};
      filename += ".history";
      
      cout << endl << "Printing convergence history to <"
	   << filename << ">..." << endl;
      problem.print_history (filename);
      cout << endl << endl << endl;
    };
  
  return 0;
};

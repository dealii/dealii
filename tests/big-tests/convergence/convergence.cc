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
#include <fe/fe_lib.h>
#include <fe/quadrature_lib.h>
#include <numerics/base.h>
#include <numerics/assembler.h>


#include <map.h>
#include <fstream.h>
#include <cmath>
#include <string>
extern "C" {
#  include <stdlib.h>
}




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
class PoissonProblem : public ProblemBase<dim> {
  public:
    PoissonProblem ();

    void clear ();
    void create_new ();
    void run (unsigned int level);
    void print_history () const;
    
  protected:
    Triangulation<dim> *tria;
    DoFHandler<dim>    *dof;
    
    Function<dim>      *rhs;
    Function<dim>      *boundary_values;

    vector<double> l1_error, l2_error, linfty_error, h1_seminorm_error, h1_error;
    vector<int>    n_dofs;
};





/**
  Right hand side constructed such that the exact solution is
  $x(1-x)$ in 1d, $x(1-x)*y(1-y)$ in 2d, etc.
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
    virtual Point<dim> gradient (const Point<dim> &p) const;
};




double RHSPoly<2>::operator () (const Point<2> &p) const {
  return 2*(p(0)*(1-p(0)) + p(1)*(1-p(1)));
};



double Solution<2>::operator () (const Point<2> &p) const {
  return p(0)*(1-p(0))*p(1)*(1-p(1));
};


Point<2> Solution<2>::gradient (const Point<2> &p) const {
  return Point<2> ((1-2*p(0))*p(1)*(1-p(1)), (1-2*p(1))*p(0)*(1-p(0)));
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
		tria(0), dof(0), rhs(0), boundary_values(0) {};




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
void PoissonProblem<dim>::run (const unsigned int level) {
  create_new ();
  
  cout << "Refinement level = " << level
       << endl;
  
  cout << "    Making grid... ";
  tria->create_hypercube ();
  tria->refine_global (level);
  cout << tria->n_active_cells() << " active cells." << endl;

  rhs             = new RHSPoly<dim>();
  boundary_values = new ZeroFunction<dim> ();
  
  
  FELinear<dim>                   fe;
  PoissonEquation<dim>            equation (*rhs);
  QGauss3<dim>                    quadrature;
  
  cout << "    Distributing dofs... "; 
  dof->distribute_dofs (fe);
  cout << dof->n_dofs() << " degrees of freedom." << endl;
  n_dofs.push_back (dof->n_dofs());

  cout << "    Assembling matrices..." << endl;
  UpdateFields update_flags = UpdateFields(update_q_points  | update_gradients |
					   update_jacobians | update_JxW_values);
  
  ProblemBase<dim>::DirichletBC dirichlet_bc;
  dirichlet_bc[0] = boundary_values;
  assemble (equation, quadrature, fe, update_flags, dirichlet_bc);

  cout << "    Solving..." << endl;
  solve ();

  Solution<dim> sol;
  dVector       l1_error_per_cell, l2_error_per_cell, linfty_error_per_cell;
  dVector       h1_seminorm_error_per_cell, h1_error_per_cell;
  QGauss3<dim>  q;
  
  cout << "    Calculating L1 error... ";
  integrate_difference (sol, l1_error_per_cell, q, fe, L1_norm);
  cout << l1_error_per_cell.l1_norm() << endl;
  l1_error.push_back (l1_error_per_cell.l1_norm());

  cout << "    Calculating L2 error... ";
  integrate_difference (sol, l2_error_per_cell, q, fe, L2_norm);
  cout << l2_error_per_cell.l2_norm() << endl;
  l2_error.push_back (l2_error_per_cell.l2_norm());

  cout << "    Calculating L-infinity error... ";
  integrate_difference (sol, linfty_error_per_cell, q, fe, Linfty_norm);
  cout << linfty_error_per_cell.linfty_norm() << endl;
  linfty_error.push_back (linfty_error_per_cell.linfty_norm());
  
  cout << "    Calculating H1-seminorm error... ";
  integrate_difference (sol, h1_seminorm_error_per_cell, q, fe, H1_seminorm);
  cout << h1_seminorm_error_per_cell.l2_norm() << endl;
  h1_seminorm_error.push_back (h1_seminorm_error_per_cell.l2_norm());

  cout << "    Calculating H1 error... ";
  integrate_difference (sol, h1_error_per_cell, q, fe, H1_norm);
  cout << h1_error_per_cell.l2_norm() << endl;
  h1_error.push_back (h1_error_per_cell.l2_norm());

  if (level<=5) 
    {
      dVector l1_error_per_dof, l2_error_per_dof, linfty_error_per_dof;
      dVector h1_seminorm_error_per_dof, h1_error_per_dof;
      dof->distribute_cell_to_dof_vector (l1_error_per_cell, l1_error_per_dof);
      dof->distribute_cell_to_dof_vector (l2_error_per_cell, l2_error_per_dof);
      dof->distribute_cell_to_dof_vector (linfty_error_per_cell, linfty_error_per_dof);
      dof->distribute_cell_to_dof_vector (h1_seminorm_error_per_cell, h1_seminorm_error_per_dof);
      dof->distribute_cell_to_dof_vector (h1_error_per_cell, h1_error_per_dof);

      string filename = "gnuplot.";
      filename += ('0'+level);
      cout << "    Writing error plots to <" << filename << ">..." << endl;
      
      DataOut<dim> out;
      ofstream gnuplot(filename.c_str());
      fill_data (out);
      out.add_data_vector (l1_error_per_dof, "L1-Error");
      out.add_data_vector (l2_error_per_dof, "L2-Error");
      out.add_data_vector (linfty_error_per_dof, "Linfty-Error");
      out.add_data_vector (h1_seminorm_error_per_dof, "H1-seminorm-Error");
      out.add_data_vector (h1_error_per_dof, "H1-Error");
      out.write_gnuplot (gnuplot);
      gnuplot.close ();
    }
  else
    cout << "    Not writing error as grid." << endl;
  
  cout << endl;
};


template <int dim>
void PoissonProblem<dim>::print_history () const {
  ofstream out("gnuplot.history");
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

  cout << "Average error reduction rates for h->h/2:" << endl;
  cout << "    L1 error         : " << 1./average_l1 << endl
       << "    L2 error         : " << 1./average_l2 << endl
       << "    Linfty error     : " << 1./average_linfty << endl
       << "    H1 seminorm error: " << 1./average_h1_semi << endl
       << "    H1 error         : " << 1./average_h1 << endl;
};




int main () {
  PoissonProblem<2> problem;

  for (unsigned int level=1; level<9; ++level)
    problem.run (level);

  cout << endl << "Printing convergence history to <gnuplot.history>..." << endl;
  problem.print_history ();

  return 0;
};

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
    virtual void create_new ();
    virtual void run (unsigned int level);
    
  protected:
    Triangulation<dim> *tria;
    DoFHandler<dim>    *dof;
    
    Function<dim>      *rhs;
    Function<dim>      *boundary_values;
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
};




template <int dim>
double RHSPoly<dim>::operator () (const Point<dim> &p) const {
  double ret_val = 0;
  for (unsigned int i=0; i<dim; ++i)
    ret_val += 2*p(i)*(1.-p(i));
  return ret_val;
};


template <int dim>
double Solution<dim>::operator () (const Point<dim> &p) const {
  double ret_val = 1;
  for (unsigned int i=0; i<dim; ++i)
    ret_val *= p(i)*(1.-p(i));
  return ret_val;
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

  cout << "    Assembling matrices..." << endl;
  FEValues<dim>::UpdateStruct update_flags;
  update_flags.q_points  = update_flags.gradients  = true;
  update_flags.jacobians = update_flags.JxW_values = true;
  
  ProblemBase<dim>::DirichletBC dirichlet_bc;
  dirichlet_bc[0] = boundary_values;
  assemble (equation, quadrature, fe, update_flags, dirichlet_bc);

  cout << "    Solving..." << endl;
  solve ();

  Solution<dim> sol;
  dVector       l1_error_per_cell, l2_error_per_cell, linfty_error_per_cell;
  QGauss4<dim>  q;
  
  cout << "    Calculating L1 error... ";
  integrate_difference (sol, l1_error_per_cell, q, fe, L1_norm);
  cout << l1_error_per_cell.l1_norm() << endl;

  cout << "    Calculating L2 error... ";
  integrate_difference (sol, l2_error_per_cell, q, fe, L2_norm);
  cout << l2_error_per_cell.l2_norm() << endl;

  cout << "    Calculating L-infinity error... ";
  integrate_difference (sol, linfty_error_per_cell, q, fe, Linfty_norm);
  cout << linfty_error_per_cell.linfty_norm() << endl;

  dVector l1_error_per_dof, l2_error_per_dof, linfty_error_per_dof;
  dof->distribute_cell_to_dof_vector (l1_error_per_cell, l1_error_per_dof);
  dof->distribute_cell_to_dof_vector (l2_error_per_cell, l2_error_per_dof);
  dof->distribute_cell_to_dof_vector (linfty_error_per_cell, linfty_error_per_dof);

  string filename = "gnuplot.";
  filename += ('0'+level);
  cout << "    Writing error plots to <" << filename << ">..." << endl;

  DataOut<dim> out;
  ofstream gnuplot(filename.c_str());
  fill_data (out);
  out.add_data_vector (l1_error_per_dof, "L1-Error");
  out.add_data_vector (l2_error_per_dof, "L2-Error");
  out.add_data_vector (linfty_error_per_dof, "L3-Error");
  out.write_gnuplot (gnuplot);
  gnuplot.close ();
  
  cout << endl;
};





int main () {
  PoissonProblem<2> problem;

  for (unsigned int level=1; level<5; ++level)
    problem.run (level);

  return 0;
};

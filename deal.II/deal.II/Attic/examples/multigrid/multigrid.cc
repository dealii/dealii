/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include <grid/tria.h>
#include <grid/mg_dof.h>
#include <grid/tria_accessor.h>
#include <grid/mg_dof_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary.h>
#include <grid/dof_constraints.h>
#include <base/function.h>
#include <basic/data_io.h>
#include <fe/fe_lib.lagrange.h>
#include <fe/fe_lib.criss_cross.h>
#include <fe/fe_update_flags.h>
#include <base/quadrature_lib.h>
#include <numerics/assembler.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>

#include <lac/dvector.h>
#include <lac/dsmatrix.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>

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
class PoissonProblem {
  public:
				     /**
				      *	Declare a data type which denotes a
				      *	mapping between a boundary indicator
				      *	and the function denoting the boundary
				      *	values on this part of the boundary.
				      *	Only one boundary function may be given
				      *	for each boundary indicator, which is
				      *	guaranteed by the #map# data type.
				      *	
				      *	See the general documentation of this
				      *	class for more detail.
				      */
    typedef map<unsigned char,const Function<dim>*> FunctionMap;
				     /**
				      * Typdedef an iterator which assembles
				      * matrices and vectors.
				      */
    typedef TriaActiveIterator<dim, Assembler<dim> > active_assemble_iterator;

    PoissonProblem (unsigned int order);

    void clear ();
    void create_new ();
    void solve ();
    
				     /**
				      * Initiate the process of assemblage of
				      * vectors and system matrix. Use the
				      * given equation object and the given
				      * quadrature formula. Also use the list
				      * of dirichlet boundary value functions
				      * (by default, no dirichlet bc are assumed
				      * which means that all bc are included
				      * into the weak formulation).
				      *
				      * For what exactly happens here, refer to
				      * the general doc of this class.
				      */
    virtual void assemble (const Equation<dim>      &equation,
			   const Quadrature<dim>    &q,
			   const FiniteElement<dim> &fe,
			   const UpdateFlags         update_flags,
			   const FunctionMap        &dirichlet_bc = FunctionMap(),
			   const Boundary<dim>      &boundary = StraightBoundary<dim>());

    int run (unsigned int level);
    void print_history (string filename) const;
    
  protected:
    Triangulation<dim> *tria;
    MGDoFHandler<dim>  *dof;
    

				     /**
				      * Sparsity pattern of the system matrix.
				      */
    dSMatrixStruct      system_sparsity;

				     /**
				      * System matrix.
				      */
    dSMatrix            system_matrix;

				     /**
				      * Vector storing the right hand side.
				      */
    dVector             right_hand_side;

				     /**
				      * Solution vector.
				      */
    dVector             solution;

				     /**
				      * List of constraints introduced by
				      * hanging nodes.
				      */
    ConstraintMatrix    constraints;

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
    virtual Point<dim> gradient (const Point<dim> &p) const;
};




double RHSPoly<2>::operator () (const Point<2> &p) const {
  const double x = p(0),
	       y = p(1);
  const double pi= 3.1415926536;
  return 4*pi*pi*(sin(2*pi*x)+sin(2*pi*y));
};



double Solution<2>::operator () (const Point<2> &p) const {
  const double x = p(0),
	       y = p(1);
  const double pi= 3.1415926536;
  return sin(2*pi*x)+sin(2*pi*y);
};


Point<2> Solution<2>::gradient (const Point<2> &p) const {
  const double x = p(0),
	       y = p(1);
  const double pi= 3.1415926536;
  return Point<2> (2*pi*cos(2*pi*x),
		   2*pi*cos(2*pi*y));
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
PoissonProblem<dim>::PoissonProblem (unsigned int order) :
		tria(0), dof(0), rhs(0),
		boundary_values(0), order(order) {};




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

  system_sparsity.reinit (0,0,1);
  system_matrix.clear ();
  right_hand_side.reinit (1);
  solution.reinit (1);
  constraints.clear ();
};




template <int dim>
void PoissonProblem<dim>::create_new () {
  clear ();
  
  tria = new Triangulation<dim>();
  dof = new MGDoFHandler<dim> (tria);
};



template <int dim>
void PoissonProblem<dim>::assemble (const Equation<dim>      &equation,
				    const Quadrature<dim>    &quadrature,
				    const FiniteElement<dim> &fe,
				    const UpdateFlags         update_flags,
				    const FunctionMap        &dirichlet_bc,
				    const Boundary<dim>      &boundary) {
  Assert ((tria!=0) && (dof!=0), ExcNoTriaSelected());
  
  system_sparsity.reinit (dof->DoFHandler<dim>::n_dofs(),
			  dof->DoFHandler<dim>::n_dofs(),
			  dof->max_couplings_between_dofs());
  right_hand_side.reinit (dof->DoFHandler<dim>::n_dofs());
  
				   // make up sparsity pattern and
				   // compress with constraints
  constraints.clear ();
  dof->DoFHandler<dim>::make_constraint_matrix (constraints);
  dof->DoFHandler<dim>::make_sparsity_pattern (system_sparsity);
  constraints.condense (system_sparsity);

				   // reinite system matrix
  system_matrix.reinit (system_sparsity);
				   // reinit solution vector, preset
				   // with zeroes.
  solution.reinit (dof->DoFHandler<dim>::n_dofs());
  
				   // create assembler
  AssemblerData<dim> data (*dof,
			   true, true, //assemble matrix and rhs
			   system_matrix,
			   right_hand_side,
			   quadrature,
			   fe,
			   update_flags,
			   boundary);
  active_assemble_iterator assembler (tria,
				      tria->begin_active()->level(),
				      tria->begin_active()->index(),
				      &data);
				   // loop over all cells, fill matrix and rhs
  do 
    {
      assembler->assemble (equation);
    }
  while ((++assembler).state() == valid);

				   // condense system matrix in-place
  constraints.condense (system_matrix);

				   // condense right hand side in-place
  constraints.condense (right_hand_side);

				   // apply Dirichlet bc as described
				   // in the docs
  map<int, double> boundary_value_list;
  VectorTools<dim>::interpolate_boundary_values (*dof,
						 dirichlet_bc, fe, boundary,
						 boundary_value_list);
  MatrixTools<dim>::apply_boundary_values (boundary_value_list,
					   system_matrix, solution,
					   right_hand_side);  
};



template <int dim>
void PoissonProblem<dim>::solve () {
  Assert ((tria!=0) && (dof!=0), ExcNoTriaSelected());
  
  SolverControl                    control(4000, 1e-16);
  PrimitiveVectorMemory<dVector>   memory;
  SolverCG<dSMatrix,dVector>       cg(control,memory);

				   // solve
  cg.solve (system_matrix, solution, right_hand_side);
				   // distribute solution
  constraints.distribute (solution);
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
  tria->create_hyper_ball ();
  HyperBallBoundary<dim> boundary_description;
  tria->set_boundary (&boundary_description);
  tria->begin_active()->set_refine_flag();
  (++(++(tria->begin_active())))->set_refine_flag();
  tria->execute_refinement ();
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
	  fe         = new FELinear<dim>();
	  quadrature = new QGauss3<dim>();
	  boundary_quadrature = new QGauss2<dim-1>();
	  break;
    case 2:
	  fe         = new FEQuadraticSub<dim>();
	  quadrature = new QGauss4<dim>();
	  boundary_quadrature = new QGauss3<dim-1>();
	  break;
    case 3:
	  fe         = new FECubicSub<dim>();
	  quadrature = new QGauss5<dim>();
	  boundary_quadrature = new QGauss4<dim-1>();
	  break;
    case 4:
	  fe         = new FEQuarticSub<dim>();
	  quadrature = new QGauss6<dim>();
	  boundary_quadrature = new QGauss5<dim-1>();
	  break;
    default:
	  return 100000;
  };
  
  cout << "    Distributing dofs... "; 
  dof->distribute_dofs (*fe);
  cout << dof->DoFHandler<dim>::n_dofs() << " degrees of freedom." << endl;
  n_dofs.push_back (dof->DoFHandler<dim>::n_dofs());

  cout << "    Assembling matrices..." << endl;
  UpdateFlags update_flags = UpdateFlags(update_q_points  | update_gradients |
					 update_JxW_values);
  
  map<unsigned char,const Function<dim>*> dirichlet_bc;
  dirichlet_bc[0] = boundary_values;
  assemble (equation, *quadrature, *fe, update_flags, dirichlet_bc);

  cout << "    Solving..." << endl;
  solve ();

  Solution<dim> sol;
  dVector       l1_error_per_cell, l2_error_per_cell, linfty_error_per_cell;
  dVector       h1_seminorm_error_per_cell, h1_error_per_cell;
  
  cout << "    Calculating L1 error... ";
  VectorTools<dim>::integrate_difference (*dof,
					  solution, sol,
					  l1_error_per_cell,
					  *quadrature, *fe, L1_norm);
  cout << l1_error_per_cell.l1_norm() << endl;
  l1_error.push_back (l1_error_per_cell.l1_norm());

  cout << "    Calculating L2 error... ";
  VectorTools<dim>::integrate_difference (*dof,
					  solution, sol,
					  l2_error_per_cell,
					  *quadrature, *fe, L2_norm);
  cout << l2_error_per_cell.l2_norm() << endl;
  l2_error.push_back (l2_error_per_cell.l2_norm());

  cout << "    Calculating L-infinity error... ";
  VectorTools<dim>::integrate_difference (*dof,
					  solution, sol,
					  linfty_error_per_cell,
					  *quadrature, *fe, Linfty_norm);
  cout << linfty_error_per_cell.linfty_norm() << endl;
  linfty_error.push_back (linfty_error_per_cell.linfty_norm());
  
  cout << "    Calculating H1-seminorm error... ";
  VectorTools<dim>::integrate_difference (*dof,
					  solution, sol,
					  h1_seminorm_error_per_cell,
					  *quadrature, *fe, H1_seminorm);
  cout << h1_seminorm_error_per_cell.l2_norm() << endl;
  h1_seminorm_error.push_back (h1_seminorm_error_per_cell.l2_norm());

  cout << "    Calculating H1 error... ";
  VectorTools<dim>::integrate_difference (*dof,
					  solution, sol,
					  h1_error_per_cell,
					  *quadrature, *fe, H1_norm);
  cout << h1_error_per_cell.l2_norm() << endl;
  h1_error.push_back (h1_error_per_cell.l2_norm());

  if (dof->DoFHandler<dim>::n_dofs()<=5000) 
    {
      dVector l1_error_per_dof, l2_error_per_dof, linfty_error_per_dof;
      dVector h1_seminorm_error_per_dof, h1_error_per_dof;
      dof->distribute_cell_to_dof_vector (l1_error_per_cell, l1_error_per_dof);
      dof->distribute_cell_to_dof_vector (l2_error_per_cell, l2_error_per_dof);
      dof->distribute_cell_to_dof_vector (linfty_error_per_cell,
					  linfty_error_per_dof);
      dof->distribute_cell_to_dof_vector (h1_seminorm_error_per_cell,
					  h1_seminorm_error_per_dof);
      dof->distribute_cell_to_dof_vector (h1_error_per_cell, h1_error_per_dof);

//       dVector projected_solution;
//       ConstraintMatrix constraints;
//       constraints.close ();
//       VectorTools<dim>::project (*dof, constraints, *fe,
// 				 StraightBoundary<dim>(), *quadrature, 
// 				 sol, projected_solution, false,
// 				 *boundary_quadrature);
//       cout << "    Calculating L2 error of projected solution... ";
//       VectorTools<dim>::integrate_difference (*dof,
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
      out.attach_dof_handler (*dof);

      out.add_data_vector (solution, "u");
      out.add_data_vector (l1_error_per_dof, "L1-Error");
      out.add_data_vector (l2_error_per_dof, "L2-Error");
      out.add_data_vector (linfty_error_per_dof, "Linfty-Error");
      out.add_data_vector (h1_seminorm_error_per_dof, "H1-seminorm-Error");
      out.add_data_vector (h1_error_per_dof, "H1-Error");
      out.write_ucd (o);
      o.close ();
    }
  else
    cout << "    Not writing error as grid." << endl;
  
  cout << endl;

  delete fe;
  delete quadrature;
  delete boundary_quadrature;
  
  return dof->DoFHandler<dim>::n_dofs();
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

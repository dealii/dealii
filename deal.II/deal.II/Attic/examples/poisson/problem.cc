/* $Id$ */


#include "poisson.h"




template <int dim>
class CurvedLine :
  public StraightBoundary<dim> {
  public:
      virtual Point<dim> in_between (const PointArray &neighbors) const {
	Point<dim> middle = StraightBoundary<dim>::in_between(neighbors);
	double x=middle(0),
	       y=middle(1);
	
	if (y<x)
	  if (y<1-x)
	    middle(1) = 0.03*sin(6*3.141592*middle(0));
	  else
	    middle(0) = 1+0.03*sin(6*3.141592*middle(1));
	
	else
	  if (y<1-x)
	    middle(0) = 0.03*sin(6*3.141592*middle(1));
	  else
	    middle(1) = 1+0.03*sin(6*3.141592*middle(0));

	return middle;
      };
};




template <int dim>
PoissonProblem<dim>::PoissonProblem () :
		tria(0), dof(0) {};




template <int dim>
void PoissonProblem<dim>::clear () {
  if (tria != 0) delete tria;
  if (dof  != 0) delete dof;

  ProblemBase<dim>::clear ();
};




template <int dim>
void PoissonProblem<dim>::create_new (const unsigned int) {
  clear ();
  
  tria = new Triangulation<dim>();
  dof = new DoFHandler<dim> (tria);
  set_tria_and_dof (tria, dof);
};




template <int dim>
void PoissonProblem<dim>::declare_parameters (ParameterHandler &prm) {
  if (dim>=2)
    prm.declare_entry ("Test run", "zoom in", "zoom in\\|ball\\|curved line\\|random");
  else
    prm.declare_entry ("Test run", "zoom in", "zoom in\\|random");

  prm.declare_entry ("Global refinement", "0",
		     ParameterHandler::RegularExpressions::Integer);
  prm.declare_entry ("Output file", "gnuplot.1");
};




template <int dim>
bool PoissonProblem<dim>::make_grid (ParameterHandler &prm) {
  String test = prm.get ("Test run");
  unsigned int test_case;
  if (test=="zoom in") test_case = 1;
  else
    if (test=="ball") test_case = 2;
    else
      if (test=="curved line") test_case = 3;
      else
	if (test=="random") test_case = 4;
	else 
	  {
	    cerr << "This test seems not to be implemented!" << endl;
	    return false;
	  };

  switch (test_case) 
    {
      case 1:
	    make_zoom_in_grid ();
	    break;
      case 2:
	    make_ball_grid ();
	    break;
      case 3:
	    make_curved_line_grid ();
	    break;
      case 4:
	    make_random_grid ();
	    break;
      default:
	    return false;
    };

  int refine_global = prm.get_integer ("Global refinement");
  if ((refine_global < 0) || (refine_global>10))
    return false;
  else
    tria->refine_global (refine_global);

  return true;
};

	  


template <int dim>
void PoissonProblem<dim>::make_zoom_in_grid () {
  tria->create_hypercube ();
				   // refine first cell
  tria->begin_active()->set_refine_flag();
  tria->execute_refinement ();
				   // refine first active cell
				   // on coarsest level
  tria->begin_active()->set_refine_flag ();
  tria->execute_refinement ();
  
  Triangulation<dim>::active_cell_iterator cell;
  for (int i=0; i<17; ++i) 
    {
				       // refine the presently
				       // second last cell 17
				       // times
      cell = tria->last_active(tria->n_levels()-1);
      --cell;
      cell->set_refine_flag ();
      tria->execute_refinement ();
    };
};



template <int dim>
void PoissonProblem<dim>::make_ball_grid () {
				   // make ball grid around origin with
				   // unit radius
  const Point<dim> origin;
  static const HyperBallBoundary<dim> boundary(origin, 1.);
  
  tria->create_hyper_ball (origin, 1.);
  tria->set_boundary (&boundary);
};




template <int dim>
void PoissonProblem<dim>::make_curved_line_grid () {
				   // set the boundary function
  static const CurvedLine<dim> boundary;

  tria->create_hypercube ();
  tria->set_boundary (&boundary);
};



template <int dim>
void PoissonProblem<dim>::make_random_grid () {
  tria->create_hypercube ();
  tria->refine_global (1);
	
  Triangulation<dim>::active_cell_iterator cell, endc;
  for (int i=0; i<12; ++i) 
    {
      int n_levels = tria->n_levels();
      cell = tria->begin_active();
      endc = tria->end();
      
      for (; cell!=endc; ++cell) 
	{
	  double r      = rand()*1.0/RAND_MAX,
		 weight = 1.*
			  (cell->level()*cell->level()) /
			  (n_levels*n_levels);
	  
	  if (r <= 0.5*weight)
	    cell->set_refine_flag ();
	};
      
      tria->execute_refinement ();
    };
};





template <int dim>
void PoissonProblem<dim>::run (ParameterHandler &prm) {
  cout << "Test case = " << prm.get ("Test run") << endl
       << endl;
  
  cout << "    Making grid..." << endl;
  if (!make_grid (prm))
    return;
  

  FELinear<dim>                   fe;
  static const RightHandSide<dim> rhs;
  PoissonEquation<dim>            equation (rhs);
  QGauss4<dim>                    quadrature;
  
  cout << "    Distributing dofs... "; 
  dof->distribute_dofs (fe);
  cout << dof->n_dofs() << " degrees of freedom." << endl;

  cout << "    Assembling matrices..." << endl;
  FEValues<dim>::UpdateStruct update_flags;
  update_flags.q_points  = update_flags.gradients  = true;
  update_flags.jacobians = update_flags.JxW_values = true;
  
  ProblemBase<dim>::DirichletBC dirichlet_bc;
  ZeroFunction<dim> zero;
  dirichlet_bc[0] = &zero;
  assemble (equation, quadrature, fe, update_flags, dirichlet_bc);

  cout << "    Solving..." << endl;
  solve ();

  cout << "    Writing to file <" << prm.get("Output file") << ">..." << endl;
  DataOut<dim> out;

  ofstream gnuplot(prm.get("Output file"));
  fill_data (out); 
  out.write_gnuplot (gnuplot);
  gnuplot.close ();
};





template class PoissonProblem<1>;
template class PoissonProblem<2>;

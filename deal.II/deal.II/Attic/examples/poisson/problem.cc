/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */



#include "poisson.h"
#include <lac/vector.h>



template <int dim>
class BoundaryValuesSine : public Function<dim> {
  public:
    				     /**
				      * Return the value of the function
				      * at the given point.
				      */
    virtual double operator () (const Point<dim> &p) const {
      double x = 1;
      
      for (unsigned int i=0; i<dim; ++i)
	x *= cos(2*3.1415926536*p(i));
      return x;
    };
    

				     /**
				      * Set #values# to the point values
				      * of the function at the #points#.
				      * It is assumed that #values# be
				      * empty.
				      */
    virtual void value_list (const vector<Point<dim> > &points,
			     vector<double>            &values) const {
      Assert (values.size() == points.size(),
	      ExcVectorHasWrongSize(values.size(), points.size()));
      for (unsigned int i=0; i<points.size(); ++i) 
	values[i] = BoundaryValuesSine<dim>::operator() (points[i]);
    };
};



template <int dim>
class BoundaryValuesJump : public Function<dim> {
  public:
    				     /**
				      * Return the value of the function
				      * at the given point.
				      */
    virtual double operator () (const Point<dim> &p) const {
      switch (dim) 
	{
	  case 1:
		return 0;
	  default:
		if (p(0) == p(1))
		  return 0.5;
		else
		  return (p(0)>p(1) ? 0. : 1.);
	};
    };
};




template <int dim>
class RHSTrigPoly : public Function<dim> {
  public:
    				     /**
				      * Return the value of the function
				      * at the given point.
				      */
    virtual double operator () (const Point<dim> &p) const;
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
class CurvedLine :
  public StraightBoundary<dim> {
  public:
    virtual Point<dim>
    get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const;

    virtual Point<dim>
    get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const;
};



template <int dim>
Point<dim>
CurvedLine<dim>::get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const
{
  Point<dim> middle = StraightBoundary<dim>::get_new_point_on_line (line);

				   // if the line is at the top of bottom
				   // face: do a special treatment on
				   // this line. Note that if the
				   // z-value of the midpoint is either
				   // 0 or 1, then the z-values of all
				   // vertices of the line is like that
  if (dim>=3)
    if (((middle(2) == 0) || (middle(2) == 1))
				       // find out, if the line is in the
				       // interior of the top or bottom face
				       // of the domain, or at the edge.
				       // lines at the edge need to undergo
				       // the usual treatment, while for
				       // interior lines taking the midpoint
				       // is sufficient
				       //
				       // note: the trick with the boundary
				       // id was invented after the above was
				       // written, so we are not very strict
				       // here with using these flags
	&& (line->boundary_indicator() == 1))
      return middle;


  double x=middle(0),
	 y=middle(1);
  
  if (y<x)
    if (y<1-x)
      middle(1) = 0.04*sin(6*3.141592*middle(0));
    else
      middle(0) = 1+0.04*sin(6*3.141592*middle(1));
  
  else
    if (y<1-x)
      middle(0) = 0.04*sin(6*3.141592*middle(1));
    else
      middle(1) = 1+0.04*sin(6*3.141592*middle(0));
  
  return middle;
};



template <int dim>
Point<dim>
CurvedLine<dim>::get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const
{
  Point<dim> middle = StraightBoundary<dim>::get_new_point_on_quad (quad);

				   // if the face is at the top of bottom
				   // face: do not move the midpoint in
				   // x/y direction. Note that if the
				   // z-value of the midpoint is either
				   // 0 or 1, then the z-values of all
				   // vertices of the quad is like that
  if ((middle(2) == 0) || (middle(2) == 1))
    return middle;
  
  double x=middle(0),
	 y=middle(1);
  
  if (y<x)
    if (y<1-x)
      middle(1) = 0.04*sin(6*3.141592*middle(0));
    else
      middle(0) = 1+0.04*sin(6*3.141592*middle(1));
  
  else
    if (y<1-x)
      middle(0) = 0.04*sin(6*3.141592*middle(1));
    else
      middle(1) = 1+0.04*sin(6*3.141592*middle(0));
  
  return middle;
};




template <int dim>
double RHSTrigPoly<dim>::operator () (const Point<dim> &p) const {
  const double pi = 3.1415926536;
  switch (dim) 
    {
      case 1:
	    return p(0)*p(0)*cos(2*pi*p(0));
      case 2:
	    return (-2.0*cos(pi*p(0)/2)*p(1)*sin(pi*p(1)) +
		    2.0*p(0)*sin(pi*p(0)/2)*pi*p(1)*sin(pi*p(1)) +
		    5.0/4.0*p(0)*p(0)*cos(pi*p(0)/2)*pi*pi*p(1)*sin(pi*p(1)) -
		    2.0*p(0)*p(0)*cos(pi*p(0)/2)*cos(pi*p(1))*pi);
      default:
	    return 0;
    };
};



template <int dim>
double RHSPoly<dim>::operator () (const Point<dim> &p) const {
  double ret_val = 0;
  for (unsigned int i=0; i<dim; ++i)
    ret_val += 2*p(i)*(1.-p(i));
  return ret_val;
};





template <int dim>
PoissonProblem<dim>::PoissonProblem () :
		tria(0), dof(0), rhs(0), boundary_values(0), boundary(0) {};




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

  if (boundary != 0)
    {
      delete boundary;
      boundary = 0;
    };
  
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
    prm.declare_entry ("Test run", "zoom in",
		       Patterns::Sequence("tensor|zoom in|ball|curved line|"
					  "random|jump|L-region"));
  else
    prm.declare_entry ("Test run", "zoom in",
		       Patterns::Sequence("tensor|zoom in|random"));

  prm.declare_entry ("Global refinement", "0",
		     Patterns::Integer());
  prm.declare_entry ("Right hand side", "zero",
		     Patterns::Sequence("zero|constant|trigpoly|poly"));
  prm.declare_entry ("Boundary values", "zero",
		     Patterns::Sequence("zero|sine|jump"));
  prm.declare_entry ("Output file", "gnuplot.1");
};




template <int dim>
bool PoissonProblem<dim>::make_grid (ParameterHandler &prm) {
  string test = prm.get ("Test run");
  unsigned int test_case;
  if (test=="zoom in") test_case = 1;
  else
    if (test=="ball") test_case = 2;
    else
      if (test=="curved line") test_case = 3;
      else
	if (test=="random") test_case = 4;
	else
	  if (test=="tensor") test_case = 5;
	  else
	    if (test=="jump") test_case = 6;
	    else
	      if (test=="L-region") test_case = 7;
	      else
		{
		  cerr << "This test seems not to be implemented!" << endl;
		  return false;
		};

  switch (test_case) 
    {
      case 1:
	    boundary = new StraightBoundary<dim>();
	    tria->set_boundary (boundary);
	    make_zoom_in_grid ();
	    break;
      case 2:
					     // make ball grid around origin with
					     // unit radius
      {
	    static const Point<dim> origin;
	    boundary = new HyperBallBoundary<dim>(origin, 1.);
	    tria->create_hyper_ball (origin, 1.);
	    tria->set_boundary (boundary);
	    break;
      };
      case 3:
					     // set the boundary function
      {
	    boundary = new CurvedLine<dim>();
	    tria->create_hypercube ();
	    tria->set_boundary (boundary);
	    break;
      };
      case 4:
	    boundary = new StraightBoundary<dim>();
	    tria->set_boundary (boundary);
	    make_random_grid ();
	    break;
      case 5:
	    boundary = new StraightBoundary<dim>();
	    tria->set_boundary (boundary);
	    tria->create_hypercube ();
	    break;
      case 6:
	    boundary = new StraightBoundary<dim>();
	    tria->set_boundary (boundary);
	    tria->create_hypercube ();
	    tria->refine_global (1);
	    for (unsigned int i=0; i<5; ++i)
	      {
		tria->begin_active(tria->n_levels()-1)->set_refine_flag();
		(--(tria->last_active()))->set_refine_flag();
		tria->execute_coarsening_and_refinement ();
	      };
	    break;
      case 7:
	    boundary = new StraightBoundary<dim>();
	    tria->set_boundary (boundary);
	    tria->create_hyper_L ();
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
  tria->execute_coarsening_and_refinement ();
				   // refine first active cell
				   // on coarsest level
  tria->begin_active()->set_refine_flag ();
  tria->execute_coarsening_and_refinement ();
  
  Triangulation<dim>::active_cell_iterator cell;
  for (int i=0; i<(dim==3 ? 5 : 17); ++i) 
    {
				       // refine the presently
				       // second last cell several
				       // times
      cell = tria->last_active(tria->n_levels()-1);
      --cell;
      cell->set_refine_flag ();
      tria->execute_coarsening_and_refinement ();
    };
};




template <int dim>
void PoissonProblem<dim>::make_random_grid () {
  tria->create_hypercube ();
  tria->refine_global (1);
	
  Triangulation<dim>::active_cell_iterator cell, endc;
  for (int i=0; i<(dim==3 ? 7 : 12); ++i)
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
      
      tria->execute_coarsening_and_refinement ();
    };
};
  



template <int dim>
bool PoissonProblem<dim>::set_right_hand_side (ParameterHandler &prm) {
  string rhs_name = prm.get ("Right hand side");

  if (rhs_name == "zero")
    rhs = new ZeroFunction<dim>();
  else
    if (rhs_name == "constant")
      rhs = new ConstantFunction<dim>(1.);
    else
      if (rhs_name == "trigpoly")
	rhs = new RHSTrigPoly<dim>();
      else
	if (rhs_name == "poly")
	  rhs = new RHSPoly<dim> ();
	else
	  return false;

  if (rhs != 0)
    return true;
  else
    return false;
};



template <int dim>
bool PoissonProblem<dim>::set_boundary_values (ParameterHandler &prm) {
  string bv_name = prm.get ("Boundary values");
  
  if (bv_name == "zero")
    boundary_values = new ZeroFunction<dim> ();
  else
    if (bv_name == "sine")
      boundary_values = new BoundaryValuesSine<dim> ();
    else
      if (bv_name == "jump")
	boundary_values = new BoundaryValuesJump<dim> ();
      else 
	{
	  cout << "Unknown boundary value function " << bv_name << endl;
	  return false;
	};

  if (boundary_values != 0)
    return true;
  else
    return false;
};




template <int dim>
void PoissonProblem<dim>::run (ParameterHandler &prm) {
  cout << "Test case = " << prm.get ("Test run")
       << endl;
  
  cout << "    Making grid... ";
  if (!make_grid (prm))
    return;
  cout << tria->n_active_cells() << " active cells." << endl;
  
  if (!set_right_hand_side (prm))
    return;

  if (!set_boundary_values (prm))
    return;
  
  FELinear<dim>                   fe;
  PoissonEquation<dim>            equation (*rhs);
  QGauss2<dim>                    quadrature;
  
  cout << "    Distributing dofs... "; 
  dof->distribute_dofs (fe);
  cout << dof->n_dofs() << " degrees of freedom." << endl;

  cout << "    Assembling matrices..." << endl;
  ProblemBase<dim>::FunctionMap dirichlet_bc;
  dirichlet_bc[0] = boundary_values;
  assemble (equation, quadrature,
	    UpdateFlags(update_gradients | update_JxW_values |
			update_q_points),
	    dirichlet_bc);

  cout << "    Solving..." << endl;
  solve ();

  cout << "    Writing to file <" << prm.get("Output file") << ">..."
       << endl;

  DataOut<dim> out;
  string o_filename = prm.get ("Output file");
  ofstream gnuplot(o_filename.c_str());
  fill_data (out);
  out.write_gnuplot (gnuplot);
  gnuplot.close ();

				   // release the lock of the DoF object to
				   // the FE object
  dof->clear ();
  
  cout << endl;
};





template class PoissonProblem<3>;

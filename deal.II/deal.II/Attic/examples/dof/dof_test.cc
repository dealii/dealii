/* $Id$ */


#include <grid/dof.h>
#include <grid/tria.h>
#include <fe/fe_lib.lagrange.h>
#include <grid/tria_boundary.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <lac/dsmatrix.h>
#include <basic/parameter_handler.h>
#include <grid/dof_constraints.h>


#include <fstream>
#include <cmath>
#include <cstdlib>



// 1: continuous refinement of the unit square always in the middle
// 2: refinement of the circle at the boundary
// 2: refinement of a wiggled area at the boundary
// 4: random refinement






template <int dim>
class Ball :
  public StraightBoundary<dim> {
  public:
      virtual Point<dim> in_between (const PointArray &neighbors) const {
	Point<dim> middle = StraightBoundary<dim>::in_between(neighbors);

	for (int i=0; i<dim; ++i)
	  middle(i) -= .5;
	middle /= (sqrt(middle.square())*sqrt(2));
	for (int i=0; i<dim; ++i)
	  middle(i) += .5;
	
	return middle;
      };
};


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
};



template <int dim>
class TestCases : public MultipleParameterLoop::UserClass{
  public:
    TestCases ();
    
    virtual void create_new (const unsigned int run_no);
    virtual void declare_parameters (ParameterHandler &prm);
    virtual void run (ParameterHandler &prm);

  private:
    Triangulation<dim> *tria;
    DoFHandler<dim>    *dof;
};



template <int dim>
TestCases<dim>::TestCases () :
		tria(0), dof(0) {};



template <int dim>
void TestCases<dim>::create_new (const unsigned int) {
  if (tria != 0) delete tria;
  if (dof  != 0) delete dof;

  tria = new Triangulation<dim>();
  tria->create_hypercube();

  dof = new DoFHandler<dim> (tria);
};



template <int dim>
void TestCases<dim>::declare_parameters (ParameterHandler &prm) {
  if (dim>=2)
    prm.declare_entry ("Test run", "zoom in",
		       Patterns::Sequence("zoom in|ball|curved line|random"));
  else
    prm.declare_entry ("Test run", "zoom in",
		       Patterns::Sequence("zoom in|random"));
  prm.declare_entry ("Grid file", "grid.1");
  prm.declare_entry ("Sparsity file", "sparsity.1");
  prm.declare_entry ("Condensed sparsity file", "sparsity.c.1");
};



template <int dim>
void TestCases<dim>::run (ParameterHandler &prm) {
  cout << "Test case = " << prm.get ("Test run") << endl
       << endl;
  
  cout << "    Making grid..." << endl;

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
	  cerr << "This test seems not to be implemented!" << endl;
  
  
  switch (test_case) 
    {
      case 1: 
      {
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

//	tria->refine_global (5);
	
	break;
      }
      
      case 2:
      case 3:
      {
					 // set the boundary function
	Boundary<dim> *boundary = (test_case==2 ?
				   static_cast<Boundary<dim>*>(new Ball<dim>()) :
				   static_cast<Boundary<dim>*>(new CurvedLine<dim>()));
	tria->set_boundary (boundary);
	
					 // refine once
	tria->begin_active()->set_refine_flag();
	tria->execute_refinement ();
	
	Triangulation<dim>::active_cell_iterator cell, endc;
	for (int i=0; i<4; ++i) 
	  {
	    cell = tria->begin_active();
	    endc = tria->end();
	    
					     // refine all
					     // boundary cells
	    for (; cell!=endc; ++cell)
	      if (cell->at_boundary())
		cell->set_refine_flag();
	    
	    tria->execute_refinement();
	  };
	
	break;
      }

      case 4:
      {
					 // refine once
	tria->begin_active()->set_refine_flag();
	tria->execute_refinement ();
	
	Triangulation<dim>::active_cell_iterator cell, endc;
	for (int i=0; i<(dim==2 ? 12 : 20); ++i) 
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
	break;	
      }
    };

  				   // output the grid
  cout << "    Writing grid..." << endl;
  ofstream out(prm.get("Grid file").c_str());
  tria->print_gnuplot (out);




  cout << "    Distributing degrees of freedom..." << endl;
  FELinear<dim> fe;
  dof->distribute_dofs (fe);

  cout << "    Renumbering degrees of freedom..." << endl;
  dof->renumber_dofs (Cuthill_McKee, false);
    
  dSMatrixStruct sparsity (dof->n_dofs(),
			   dof->max_couplings_between_dofs());
  
  
  dof->make_sparsity_pattern (sparsity);
  int unconstrained_bandwidth = sparsity.bandwidth();

  cout << "    Writing sparsity pattern... (This may take a while)" << endl;
  ofstream sparsity_out (prm.get("Sparsity file").c_str());
  sparsity.print_gnuplot (sparsity_out);


  
				   // computing constraints
  cout << "    Computing constraints..." << endl;
  ConstraintMatrix constraints;
  dof->make_constraint_matrix (constraints);
  constraints.condense (sparsity);
  
  cout << "    Writing condensed sparsity pattern... (This may take a while)" << endl;
  ofstream c_sparsity_out (prm.get("Condensed sparsity file").c_str());
  sparsity.print_gnuplot (c_sparsity_out);


//  DataOut<dim> dataout;
//  dataout.attach_dof_handler (*dof);
//  ofstream o("o.inp");
//  dataout.write_ucd (o);
  
  cout << endl
       << "    Total number of cells         = " << tria->n_cells() << endl
       << "    Total number of active cells  = " << tria->n_active_cells() << endl
       << "    Number of DoFs                = " << dof->n_dofs() << endl
       << "    Number of constraints         = " << constraints.n_constraints() << endl
       << "    Unconstrained matrix bandwidth= " << unconstrained_bandwidth << endl
       << "    Constrained matrix bandwidth  = " << sparsity.bandwidth()
       << endl << endl;
};



int main (int argc, char **argv) {
  if (argc!=2) 
    {
      cout << "Usage: grid_test parameterfile" << endl << endl;
      return 1;
    };

  TestCases<2> tests;
  class MultipleParameterLoop input_data;   //gcc2.7 does nonsense, wait for gcc2.8

  tests.declare_parameters(input_data);
  input_data.read_input (argv[1]);
  input_data.loop (tests);
  
  return 0;
};


/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */



#include <grid/grid_out.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_handler.h>
#include <grid/tria.h>
#include <fe/fe_lib.lagrange.h>
#include <grid/tria_boundary.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_generator.h>
#include <lac/sparse_matrix.h>
#include <base/parameter_handler.h>
#include <dofs/dof_constraints.h>
#include <numerics/dof_renumbering.h>

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
    virtual Point<dim>
    get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const {
      Point<dim> middle = StraightBoundary<dim>::get_new_point_on_line(line);
      
      for (int i=0; i<dim; ++i)
	middle(i) -= .5;
      middle *= sqrt(dim) / (sqrt(middle.square())*2);
      for (int i=0; i<dim; ++i)
	middle(i) += .5;
      
      return middle;
    };

    
    virtual Point<dim>
    get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const {
      Point<dim> middle = StraightBoundary<dim>::get_new_point_on_quad(quad);
      
      for (int i=0; i<dim; ++i)
	middle(i) -= .5;
      middle *= sqrt(dim) / (sqrt(middle.square())*2);
      for (int i=0; i<dim; ++i)
	middle(i) += .5;
      
      return middle;
    };
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
class TestCases : public MultipleParameterLoop::UserClass{
  public:
    TestCases ();
    virtual ~TestCases ();
    
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
TestCases<dim>::~TestCases () 
{
  if (dof)  delete dof;
  if (tria) delete tria;
};



template <int dim>
void TestCases<dim>::create_new (const unsigned int) {
  if (dof  != 0) delete dof;
  if (tria != 0) delete tria;

  tria = new Triangulation<dim>();
  GridGenerator::hyper_cube(*tria);

  dof = new DoFHandler<dim> (*tria);
};



template <int dim>
void TestCases<dim>::declare_parameters (ParameterHandler &prm) {
  if (dim>=2)
    prm.declare_entry ("Test run", "zoom in",
		       Patterns::Selection("zoom in|ball|curved line|random"));
  else
    prm.declare_entry ("Test run", "zoom in",
		       Patterns::Selection("zoom in|random"));
  prm.declare_entry ("Grid file", "grid.1");
  prm.declare_entry ("Sparsity file", "sparsity.1");
  prm.declare_entry ("Condensed sparsity file", "sparsity.c.1");
};



template <int dim>
void TestCases<dim>::run (ParameterHandler &prm) {
  cout << "Dimension = " << dim
       << ", Test case = " << prm.get ("Test run") << endl
       << endl;
  
  string test = prm.get ("Test run");
  unsigned int test_case = 1;
  if (test=="zoom in") test_case = 1;
  else
    if (test=="ball") test_case = 2;
    else
      if (test=="curved line") test_case = 3;
      else
	if (test=="random") test_case = 4;
	else
	  cerr << "This test seems not to be implemented!" << endl;

  
  cout << "    Making grid..." << endl;  
  Boundary<dim> *boundary = 0;
  
  switch (test_case) 
    {
      case 1: 
      {
					 // refine first cell
	tria->begin_active()->set_refine_flag();
	tria->execute_coarsening_and_refinement ();
					 // refine first active cell
					 // on coarsest level
	tria->begin_active()->set_refine_flag ();
	tria->execute_coarsening_and_refinement ();

	Triangulation<dim>::active_cell_iterator cell;
	for (int i=0; i<17; ++i) 
	  {
					     // refine the presently
					     // second last cell 17
					     // times
	    cell = tria->last_active(tria->n_levels()-1);
	    --cell;
	    cell->set_refine_flag ();
	    tria->execute_coarsening_and_refinement ();
	  };

	break;
      }
      
      case 2:
      case 3:
      {
	if (dim==3)
	  {
	    tria->begin_active()->face(2)->set_boundary_indicator(1);
	    tria->begin_active()->face(4)->set_boundary_indicator(1);
	  };
	
					 // set the boundary function
	boundary = (test_case==2 ?
		    static_cast<Boundary<dim>*>(new Ball<dim>()) :
		    static_cast<Boundary<dim>*>(new CurvedLine<dim>()));
	tria->set_boundary (0, *boundary);
	tria->set_boundary (1, *boundary);
	
					 // refine once
	tria->begin_active()->set_refine_flag();
	tria->execute_coarsening_and_refinement ();
	
	Triangulation<dim>::active_cell_iterator cell, endc;
	for (int i=0; i<6-dim; ++i) 
	  {
	    cell = tria->begin_active();
	    endc = tria->end();
	    
					     // refine all
					     // boundary cells
	    for (; cell!=endc; ++cell)
	      if (cell->at_boundary())
		cell->set_refine_flag();
	    
	    tria->execute_coarsening_and_refinement();
	  };
	
	break;
      }

      case 4:
      {
					 // refine once
	tria->begin_active()->set_refine_flag();
	tria->execute_coarsening_and_refinement ();
	
	Triangulation<dim>::active_cell_iterator cell, endc;
	for (int i=0; i<(dim==2 ? 12 : (dim==3 ? 7 : 20)); ++i) 
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
	break;	
      }
    };

  				   // output the grid
  string file_prefix ("results/");
  file_prefix += ('0'+dim);
  file_prefix += "d.";
  
  cout << "    Writing grid..." << endl;
  ofstream out((file_prefix + prm.get("Grid file")).c_str());
  GridOut().write_gnuplot (*tria, out);




  cout << "    Distributing degrees of freedom..." << endl;
  FEQ1<dim> fe;
  dof->distribute_dofs (fe);

  cout << "    Renumbering degrees of freedom..." << endl;
  DoFRenumbering::Cuthill_McKee (*dof);
    
  SparseMatrixStruct sparsity (dof->n_dofs(),
			       dof->max_couplings_between_dofs());
  
  
  DoFTools::make_sparsity_pattern (*dof, sparsity);
  int unconstrained_bandwidth = sparsity.bandwidth();

  cout << "    Writing sparsity pattern..." << endl;
  ofstream sparsity_out ((file_prefix + prm.get("Sparsity file")).c_str());
  sparsity.print_gnuplot (sparsity_out);


  
				   // computing constraints
  cout << "    Computing constraints..." << endl;
  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (*dof, constraints);
  constraints.close ();
  constraints.condense (sparsity);
  
  cout << "    Writing condensed sparsity pattern..." << endl;
  ofstream c_sparsity_out ((file_prefix +
			    prm.get("Condensed sparsity file")).c_str());
  sparsity.print_gnuplot (c_sparsity_out);


  cout << endl
       << "    Total number of cells         = " << tria->n_cells() << endl
       << "    Total number of active cells  = " << tria->n_active_cells() << endl
       << "    Number of DoFs                = " << dof->n_dofs() << endl
       << "    Number of constraints         = " << constraints.n_constraints() << endl
       << "    Unconstrained matrix bandwidth= " << unconstrained_bandwidth << endl
       << "    Constrained matrix bandwidth  = " << sparsity.bandwidth()
       << endl << endl;

				   // release the lock that dof has to the
				   // finite element object
  dof->clear ();
  tria->set_boundary (0);
  tria->set_boundary (1);
  if (boundary)
    delete boundary;
};



int main (int argc, char **argv) {
  if (argc!=3) 
    {
      cout << "Usage: grid_test dimension parameterfile" << endl << endl;
      return 1;
    };

  unsigned int dim;
  if (argv[1][0] == '2')
    dim = 2;
  else
    dim = 3;

  switch (dim)
    {
      case 2:
      {
	    TestCases<2> tests;
	    MultipleParameterLoop input_data;

	    tests.declare_parameters(input_data);
	    input_data.read_input (argv[2]);
	    input_data.loop (tests);

	    break;
      };
       
      case 3:
      {
	    TestCases<3> tests;
	    MultipleParameterLoop input_data;

	    tests.declare_parameters(input_data);
	    input_data.read_input (argv[2]);
	    input_data.loop (tests);

	    break;
      };
    };
  
  return 0;
};


/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include <grid/tria_boundary.h>
#include <grid/tria_iterator.h>
#include <grid/tria.h>
#include <fstream>
#include <string>
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



// since gcc can't resolve "test<2>()", do template parameter
// passing this way...
template <int dim>
void test (const int test_case, const Point<dim> &) {
  Triangulation<dim> tria;
  tria.create_hypercube();
  
  if ((dim==1) && ((test_case==2) || (test_case==3)))
    {
      cout << "Impossible for this dimension." << endl;
      return;
    };

  
  switch (test_case) 
    {
      case 1: 
      {
					 // we want to log the
					 // refinement history
//	ofstream history ("mesh.history");
	
					 // refine first cell
	tria.begin_active()->set_refine_flag();
//	tria.save_refine_flags (history);
	tria.execute_refinement ();
	
					 // refine first active cell
					 // on coarsest level
	tria.begin_active()->set_refine_flag ();
//	tria.save_refine_flags (history);
	tria.execute_refinement ();

	Triangulation<dim>::active_cell_iterator cell;
	for (int i=0; i<17; ++i) 
	  {
					     // refine the presently
					     // second last cell 17
					     // times
	    cell = tria.last_active(tria.n_levels()-1);
	    --cell;
	    cell->set_refine_flag ();
//	    tria.save_refine_flags (history);
	    tria.execute_refinement ();
	  };

//	tria.refine_global (5);
	
	break;
      }
      
      case 2:
      case 3:
      {
					 // set the boundary function
	Ball<dim>       ball;
	CurvedLine<dim> curved_line;
	if (test_case==2)
	  tria.set_boundary (&ball);
	else
	  tria.set_boundary (&curved_line);
	
					 // refine once
	tria.begin_active()->set_refine_flag();
	tria.execute_refinement ();
	
	Triangulation<dim>::active_cell_iterator cell, endc;
	for (int i=0; i<7; ++i) 
	  {
	    cell = tria.begin_active();
	    endc = tria.end();
	    
					     // refine all
					     // boundary cells
	    for (; cell!=endc; ++cell)
	      if (cell->at_boundary())
		cell->set_refine_flag();
	    
	    tria.execute_refinement();
	  };
	
	break;
      }

      case 4:
      {
					 // refine once
	tria.begin_active()->set_refine_flag();
	tria.execute_refinement ();
	
	Triangulation<dim>::active_cell_iterator cell, endc;
	for (int i=0; i<(dim==2 ? 13 : 30); ++i) 
	  {
	    int n_levels = tria.n_levels();
	    cell = tria.begin_active();
	    endc = tria.end();

	    for (; cell!=endc; ++cell) 
	      {
		double r      = rand()*1.0/RAND_MAX,
		       weight = 1.*
				(cell->level()*cell->level()) /
				(n_levels*n_levels);
		
		if (r <= 0.5*weight)
		  cell->set_refine_flag ();
	      };
	    
	    tria.execute_refinement ();
	  };
	break;	
      }
    };
  
  
	
				   // output the grid
  string filename("grid.");
  filename += ('0'+test_case);
  
  ofstream out(filename.c_str());
  tria.print_gnuplot (out);
    
  cout << "     Total number of cells        = " << tria.n_cells() << endl
       << "     Total number of active cells = " << tria.n_active_cells() << endl;
};



int main (int argc, char **argv) {
  if (argc!=2) 
    {
      cout << "Usage: grid_test testcase" << endl << endl
	   << "Testcases:" << endl
	   << "  1: continuous refinement of the unit square always in the middle" << endl
	   << "  2: refinement of the circle at the boundary" << endl
	   << "  3: refinement of a wiggled area at the boundary" << endl
	   << "  4: random refinement" << endl << endl;
      return 1;
    };

  test (argv[1][0]-'0', Point<2>());

  return 0;
};

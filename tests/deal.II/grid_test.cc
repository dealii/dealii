/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */



// deal_II_libraries.g=-ldeal_II_2d.g -ldeal_II_3d.g
// deal_II_libraries=-ldeal_II_2d -ldeal_II_3d




#include <grid/tria_boundary.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria.h>
#include <cmath>
#include <cstdlib>



// 1: continuous refinement of the unit square always in the middle
// 2: refinement of the circle at the boundary
// 2: refinement of a wiggled area at the boundary






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
void test (const int test_case) {
  cout << "Running testcase " << test_case
       << " in " << dim << " dimensions." << endl;
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
	
					 // refine first cell
	tria.begin_active()->set_refine_flag();
	tria.execute_coarsening_and_refinement ();
	
					 // refine first active cell
					 // on coarsest level
	tria.begin_active()->set_refine_flag ();
	tria.execute_coarsening_and_refinement ();

	Triangulation<dim>::active_cell_iterator cell;
	for (int i=0; i<(dim==2 ? 13 : 7); ++i) 
	  {
					     // refine the presently
					     // second last cell 17
					     // times
	    cell = tria.last_active(tria.n_levels()-1);
	    --cell;
	    cell->set_refine_flag ();
	    tria.execute_coarsening_and_refinement ();
	  };
	
	break;
      }
      
      case 2:
      case 3:
      {
	if (dim==3)
	  {
	    tria.begin_active()->face(2)->set_boundary_indicator(1);
	    tria.begin_active()->face(4)->set_boundary_indicator(1);
	  };
	
	
					 // set the boundary function
	Ball<dim>       ball;
	CurvedLine<dim> curved_line;
	if (test_case==2)
	  tria.set_boundary (&ball);
	else
	  tria.set_boundary (&curved_line);
	
					 // refine once
	tria.begin_active()->set_refine_flag();
	tria.execute_coarsening_and_refinement ();
	
 	Triangulation<dim>::active_cell_iterator cell, endc;
	const unsigned int steps[4] = { 0, 10, 5, 2 };
 	for (unsigned int i=0; i<steps[dim]; ++i) 
 	  {
 	    cell = tria.begin_active();
 	    endc = tria.end();
	    
 					     // refine all
 					     // boundary cells
 	    for (; cell!=endc; ++cell)
 	      if (cell->at_boundary())
 		cell->set_refine_flag();
	    
 	    tria.execute_coarsening_and_refinement();
 	  };

	tria.set_boundary (0);
	break;
      }
    };
  
  
  tria.print_gnuplot (cout);
    
  cout << "     Total number of cells        = " << tria.n_cells() << endl
       << "     Total number of active cells = " << tria.n_active_cells() << endl;
};



int main () {
				   // limit output a bit
  cout.precision (4);
  
  for (unsigned int i=1; i<=3; ++i)
    test<2> (i);
  for (unsigned int i=1; i<=3; ++i)
    test<3> (i);
  
  return 0;
};

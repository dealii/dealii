// 	$Id$	
// test for correctness of gradients on a given cell

// deal_II_libraries.g=-ldeal_II_2d.g
// deal_II_libraries=-ldeal_II_2d

#include <grid/tria.h>
#include <grid/tria_boundary.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <fe/fe_values.h>
#include <fe/fe_lib.lagrange.h>
#include <base/quadrature_lib.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <lac/vector.h>
#include <base/logstream.h>

#include <fstream>


int main ()
{
  ofstream logfile("second_derivatives.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  Triangulation<2> tria;
  GridGenerator::hyper_cube (tria,0,1);

  FEQ1<2> fe;
  DoFHandler<2> dof(tria);
  dof.distribute_dofs(fe);

  StraightBoundary<2> b;
  QTrapez<2> q;
  FEValues<2> fevalues(fe,q,update_second_derivatives);
  
  
  Vector<double> val(4);

  deallog << "Testing transformation of 2nd derivatives of shape function:" << endl;
  
				   // test for each of the four
				   // shape functions. first loop:
				   // unit cell, second loop:
				   // one vertex moved
  for (unsigned int loop=0; loop<2; ++loop)
    {
      deallog << "Test loop: " << loop << endl;
	  
      				   // move one vertex of the only cell
      if (loop==1)
	tria.begin_active()->vertex(2)(0) = 2;
      fevalues.reinit (dof.begin_active());
      
      for (unsigned int vertex=0; vertex<4; ++vertex)
	{
	  val.clear ();
	  val(vertex) = 1;
	  
	  vector<Tensor<2,2> > derivs(4);
	  fevalues.get_function_2nd_derivatives (val, derivs);
	  
	  deallog << "Vertex " << vertex << ": " << endl;
	  for (unsigned int point=0; point<4; ++point)
	    for (unsigned int component=0; component<2; ++component)
	      deallog << derivs[point][component] << endl;
	  
	  deallog << endl;
	};
    };
};

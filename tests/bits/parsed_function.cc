//----------------------------  function_parser.cc  ---------------------------
//    $Id: function_parser.cc 11749 2005-11-09 19:11:20Z wolf $
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors 
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  function_parser.cc  ---------------------------


// This program tests the functionality of the function parser
// wrapper. 

#include "../tests.h"
#include <fstream>
#include <iostream>
#include <map>
#include <base/logstream.h>
#include <base/point.h>
#include <lac/vector.h>
#include <base/parsed_function.h>
#include <base/parameter_handler.h>
#include <base/utilities.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>

template <int dim> 
void Test() {
  // A parameter handler
  ParameterHandler prm;
  
  // A triangulation
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria,0,1);
  tria.refine_global(3);
  
  // Vertices
  const std::vector<Point<dim> > &vertices = tria.get_vertices();
  
  // Test vector declaration
  for(unsigned int i=0; i<dim;++i) {
    std::string id = 
      "Function " + Utilities::int_to_string(dim) +
      " - " + Utilities::int_to_string(i);
    prm.enter_subsection(id);
    
    Functions::ParsedFunction<dim>::declare_parameters(prm, i+1);
    prm.set("Function constants", "f=" + Utilities::int_to_string(i+1));
    
    // It is cos(pi f x_i) *t^i, numbering from zero
    std::string expr = "cos(pi*f*x)*t";
    if(i>0) expr+=  "; cos(pi*f*y)*t";
    if(i>1) expr+=  "; cos(pi*f*z)*t";
    
    prm.set("Function expression", expr);
    
    Functions::ParsedFunction<dim> function(i+1);
    function.parse_parameters(prm);
    
    prm.leave_subsection();
    
    // Now test the difference from t=0 to t=1
    for(double t=0.; t<1; t+= .1) {
      function.set_time(t);
      Point<dim> p;
      std::vector<Vector<double> > values(vertices.size(),
					  Vector<double>(i+1));
      function.vector_value_list(vertices, values);
      for(unsigned int j=0; j<vertices.size(); ++j) {
	double delta=0.;
	for(unsigned int di=0; di<i; ++di) {
	  delta = values[j](di) - std::cos(M_PI*(i+1)*vertices[j][di])*t;
	  if(std::abs(delta) > 1e-10) 
	    deallog << "p(" << di << "): " << vertices[j] << ", delta: " 
		    << delta << std::endl;
	}
      }
    }
  }
  deallog << "Tested on: " << std::endl;
  prm.log_parameters(deallog);
}

int main () 
{
  std::ofstream logfile("parsed_function/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Test<1>();
  Test<2>();
  Test<3>();
}

      

  

//----------------------------  function_parser.cc  ---------------------------
//    $Id$
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
#include <base/function_parser.h>


int main () 
{
  std::ofstream logfile("function_parser.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  // Define some constants that will be used by the function parser
  std::map<std::string, double> constants;
  constants["pi"] = M_PI;
 
  // Define the variables that will be used inside the expressions
  std::string variables = "x,y";
 
  // Define the expressions of the vector_valued function.
  std::vector<std::string> expressions;
  expressions.push_back("sin(2*pi*x)");
  expressions.push_back("cos(2*pi*y*x)");
  expressions.push_back("if(x<.5,y,exp(x))");
  expressions.push_back("if(x<.5,y,exp(x)*exp(x*y))");
  
  // Concatenate the declared expressions, to test the second way of
  // initializing
  std::string concatenated="cos(pi*y)";
  // Now test each possibility
  for(unsigned int i=1; i<=expressions.size(); ++i) {
    try {
      {
	FunctionParser<2> function(i);
	function.initialize(variables,
			    expressions,
			    constants);
	deallog << "Initialize Succeded with dim = 2, "
		<< i << " components, "
		<< expressions.size() << " expressions, " 
		<< variables << " as variables." << std::endl;
      }
    } catch(...) {
	deallog << "Initialize Failed with dim = 2, "
		<< i << " components, "
		<< expressions.size() << " expressions, " 
		<< variables << " as variables." << std::endl;
    }
    
    try {
      {
	FunctionParser<2> function_bis(i);
	function_bis.initialize(variables,
				concatenated,
				constants);
	deallog << "Initialize Succeded with dim = 2, "
		<< i << " components, "
		<< concatenated << " as function and "
		<< variables << " as variables." << std::endl; 
      }
    } catch(...) {
      deallog << "Initialize Failed with dim = 2, "
	      << i << " components, "
	      << concatenated  << " as function and "
	      << variables << " as variables." << std::endl;
    }
    
    concatenated += "; " + expressions[i-1];
  }
}

      

  

/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */



#include "poisson.h"
#include <lac/vector.h>
#include <base/logstream.h>


int main (int argc, char **argv) {
  if (argc!=2) 
    {
      cout << "Usage: poisson parameterfile" << endl << endl;
      return 1;
    };

				   // no additional output to console
  deallog.depth_console (0);

  PoissonProblem<2> poisson;
  MultipleParameterLoop input_data;

  poisson.declare_parameters(input_data);
  input_data.read_input (argv[1]);
  input_data.loop (poisson);
  
  return 0;
};

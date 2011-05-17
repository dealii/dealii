//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2007, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// Check consistency of function implementations

#include "../tests.h"
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/job_identifier.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/auto_derivative_function.h>
#include <deal.II/base/flow_function.h>
#include <deal.II/lac/vector.h>

#include "functions.h"

#include <vector>
#include <iomanip>
#include <fstream>
#include <string>

#define CHECK(F) { deallog << #F << std::endl;	\
  F f;								\
  check_function_value_consistency(f, 5);				\
  check_function_gradient_consistency(f, 5);			\
  check_gradient(f, 5); \
}  
  

#define CHECKN(F,arg) { deallog << #F << '(' << arg << ')' << std::endl;	\
  F f(arg);								\
  check_function_value_consistency(f, arg+1);				\
  check_function_gradient_consistency(f, arg+1);			\
}  
  

int main()
{
  std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);

  CHECK(Functions::SquareFunction<1>);
  CHECK(Functions::SquareFunction<2>);
  CHECK(Functions::SquareFunction<3>);
  CHECK(Functions::CosineFunction<1>);
  CHECK(Functions::CosineFunction<2>);
  CHECK(Functions::CosineFunction<3>);
  CHECK(Functions::CosineGradFunction<1>);
  CHECK(Functions::CosineGradFunction<2>);
  CHECK(Functions::CosineGradFunction<3>);
  CHECK(Functions::ExpFunction<1>);
  CHECK(Functions::ExpFunction<2>);
  CHECK(Functions::ExpFunction<3>);  
}

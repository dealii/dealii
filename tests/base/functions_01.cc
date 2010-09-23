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

// Plot flow functions in library and check their consistency

#include "../tests.h"
#include <base/data_out_base.h>
#include <base/logstream.h>
#include <base/job_identifier.h>
#include <base/quadrature_lib.h>
#include <base/function_lib.h>
#include <base/auto_derivative_function.h>
#include <base/flow_function.h>
#include <lac/vector.h>

#include "functions.h"

#include <vector>
#include <iomanip>
#include <fstream>
#include <string>

int main()
{
  std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
}

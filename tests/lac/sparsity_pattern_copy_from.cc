//----------------------------  sparsity_pattern_copy_from.cc,v  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparsity_pattern_copy_from.cc,v  ---------------------------


// SparsityPattern::copy_from crashed when the number of rows or columns
// was zero

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include "testmatrix.h"
#include <fstream>
#include <iomanip>
#include <list>
#include <set>
#include <cstdio>


int
main ()
{
  std::ofstream logfile("sparsity_pattern_copy_from/output");
  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  CompressedSparsityPattern csp (10, 0);

  SparsityPattern sp;
  sp.copy_from (csp);

  deallog << "OK" << std::endl;
}




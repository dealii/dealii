//----------------------------  compressed_sparsity_pattern_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparsity_pattern_01.cc  ---------------------------


// check TrilinosWrappers::SparsityPattern::row_length

#include "sparsity_pattern_common.h"
#include <deal.II/base/utilities.h>

int main (int argc, char** argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init (argc, argv);
  std::ofstream logfile("sparsity_pattern_01/output");
  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  row_length ();
}

  
  

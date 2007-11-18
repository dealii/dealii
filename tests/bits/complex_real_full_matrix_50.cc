//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


// check FullMatrix::residual. like the full_matrix_* tests, but use
// complex-valued matrices and vectors, even though we only store real values
// in them


#include "../tests.h"
#include "full_matrix_common.h"


std::string output_file_name = "complex_real_full_matrix_50/output";


template <typename number>
void
check ()
{
  FullMatrix<std::complex<number> > m;
  make_matrix (m);
  Vector<std::complex<number> > v, w, x;
  make_vector (v);
  make_vector (w);
  make_vector (x);

  m.residual (v, w, x);
  print_vector (v);
  print_vector (w);
  print_vector (x);
}


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


// check FullMatrix::equ (6). like the full_matrix_* tests, but use
// complex-valued matrices and vectors; this time we actually store complex values
// in them


#include "../tests.h"
#include "full_matrix_common.h"


std::string output_file_name = "complex_complex_full_matrix_39/output";


template <typename number>
void
check ()
{
  FullMatrix<std::complex<number> > m,n,o,p;
  make_complex_matrix (m);
  make_complex_matrix (n);
  make_complex_matrix (o);
  make_complex_matrix (p);
  m.equ (3.1415, n, 2.718, o, 1.414, p);
  print_matrix (m);
}


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


// check FullMatrix::matrix_norm_square. like the full_matrix_* tests, but use
// complex-valued matrices and vectors; this time we actually store complex values
// in them


#include "../tests.h"
#include "full_matrix_common.h"


std::string output_file_name = "complex_complex_full_matrix_13/output";


template <typename number>
void
check ()
{
  FullMatrix<std::complex<number> > m;
  make_complex_matrix (m);
  Vector<std::complex<number> > v;
  make_complex_vector (v);
  
  deallog << m.matrix_norm_square (v) << std::endl;
}


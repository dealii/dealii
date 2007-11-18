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


// check FullMatrix::mmult


#include "../tests.h"
#include "full_matrix_common.h"


std::string output_file_name = "full_matrix_43/output";


template <typename number>
void
check ()
{
  FullMatrix<number> m, n, o;
  make_matrix (m);
  make_matrix (n);
  make_matrix (o);

  m.mmult (n, o, true);
  print_matrix (m);
  print_matrix (n);
  print_matrix (o);

  m.mmult (n, o, false);
  print_matrix (m);
  print_matrix (n);
  print_matrix (o);
}


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


// check FullMatrix::add (6)


#include "../tests.h"
#include "full_matrix_common.h"


std::string output_file_name = "full_matrix_26/output";


template <typename number>
void
check ()
{
  FullMatrix<number> m,n,o,p;
  make_matrix (m);
  make_matrix (n);
  make_matrix (o);
  make_matrix (p);
  m.add (3.1415, n, 2.718, o, 1.414, p);
  print_matrix (m);
}


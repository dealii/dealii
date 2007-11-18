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


// check FullMatrix::vmult_add


#include "../tests.h"
#include "full_matrix_common.h"


std::string output_file_name = "full_matrix_46/output";


template <typename number>
void
check ()
{
  FullMatrix<number> m;
  make_matrix (m);
  Vector<double> v, w;
  make_vector (v);
  make_vector (w);

  m.vmult_add (v, w);
  print_vector (v);
  print_vector (w);
}


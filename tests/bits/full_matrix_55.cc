//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2007, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


// check FullMatrix::swap_col and FullMatrix::swap_row for nonsymmetric
// matrices; there used to be a bug in the implementation


#include "../tests.h"
#include "full_matrix_common.h"


std::string output_file_name = "full_matrix_55/output";


template <typename number>
void
check ()
{
  FullMatrix<number> m (5, 7);
  fill_matrix (m);
  print_matrix (m);
  
  m.swap_col (2, 4);
  print_matrix (m);

  m.swap_row (2, 4);
  print_matrix (m);
}


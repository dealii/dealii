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


// check creation and output of a matrix


#include "../tests.h"
#include "full_matrix_common.h"


std::string output_file_name = "full_matrix_11/output";


template <typename number>
void
check ()
{
  FullMatrix<number> m, n;
  make_matrix (m);
  make_matrix (n);
  deallog << (m==n) << std::endl;
}


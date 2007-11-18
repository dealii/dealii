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


// check creation and output of a matrix using an array


#include "../tests.h"
#include "full_matrix_common.h"


std::string output_file_name = "full_matrix_04/output";


template <typename number>
void
check ()
{
  const number array[] = { 1, 2, 3, 4, 5,
			   6, 7, 8, 9, 0,
			   1, 2, 3, 4, 5,
			   6, 7, 8, 9, 0,
			   1, 2, 3, 4, 5 };
  
  FullMatrix<number> m (5,5,array);

  print_matrix (m);
}


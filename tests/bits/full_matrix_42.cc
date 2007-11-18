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


// check FullMatrix::invert


#include "../tests.h"
#include "full_matrix_common.h"


std::string output_file_name = "full_matrix_42/output";


template <typename number>
void
check ()
{
  for (unsigned int n=1; n<=5; ++n)
    {
      const number array[] = { 50, 2, 3, 4, 5,
			       6, 50, 8, 9, 0,
			       1, 2, 50, 4, 5,
			       6, 7, 8, 50, 0,
			       1, 2, 3, 4, 50 };
  
      FullMatrix<number> m (n,n,array), p(n,n);
      p.invert (m);
      print_matrix (p);
    }
}


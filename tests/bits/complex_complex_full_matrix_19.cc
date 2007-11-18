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


// check FullMatrix::determinant. like the full_matrix_* tests, but use
// complex-valued matrices and vectors; this time we actually store complex values
// in them


#include "../tests.h"
#include "full_matrix_common.h"


std::string output_file_name = "complex_complex_full_matrix_19/output";


template <typename number>
void
check ()
{
  for (unsigned int n=1; n<=3; ++n)
    {
      const std::complex<number> array[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    
      FullMatrix<std::complex<number> > m(n,n, array);
      print_matrix (m);
      deallog << m.determinant () << std::endl;
    }
}



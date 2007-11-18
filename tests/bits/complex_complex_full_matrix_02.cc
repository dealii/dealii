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


// check running over const iterators. like the full_matrix_* tests, but use
// complex-valued matrices and vectors; this time we actually store complex values
// in them


#include "../tests.h"
#include "full_matrix_common.h"


std::string output_file_name = "complex_complex_full_matrix_02/output";


template <typename number>
void
check ()
{
  FullMatrix<std::complex<number> > m;
  make_complex_matrix (m);

  for (typename FullMatrix<std::complex<number> >::const_iterator
	 p = m.begin(); p!=m.end(); ++p)
    deallog << p->row() << ' ' << p->column() << ' '
	    << p->value()
	    << std::endl;
}


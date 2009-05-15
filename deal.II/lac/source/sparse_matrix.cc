//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/sparse_matrix.templates.h>
#include <lac/block_vector.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace SparseMatrix
  {
				     // set this value to 1/5 of the value of
				     // the minimum grain size of
				     // vectors. this rests on the fact that
				     // we have to do a lot more work per row
				     // of a matrix than per element of a
				     // vector. it could possibly be reduced
				     // even further but that doesn't appear
				     // worth it any more for anything but
				     // very small matrices that we don't care
				     // that much about anyway.
    unsigned int minimum_parallel_grain_size = 200;
  }
}


#include "sparse_matrix.inst"

DEAL_II_NAMESPACE_CLOSE

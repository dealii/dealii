//----------------------------  vectors.all_dimensions.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  vectors.all_dimensions.cc  ---------------------------


#include <lac/vector.h>
#include <numerics/vectors.h>



// allocate storage and initialize static variables
const Quadrature<0> * const VectorTools::invalid_face_quadrature = 0;



void
VectorTools::subtract_mean_value(Vector<double>     &v,
				 const std::vector<bool> &p_select)
{
  unsigned int n = v.size();
  Assert(n == p_select.size(), ExcDimensionMismatch(n, p_select.size()));

  double       s       = 0;
  unsigned int counter = 0;
  
  for (unsigned int i=0; i<n; ++i)
    if (p_select[i])
      {
	s += v(i);
	++counter;
      };

  s /= counter;
  
  for (unsigned int i=0; i<n; ++i)
    if (p_select[i])
      v(i) -= s;  
};



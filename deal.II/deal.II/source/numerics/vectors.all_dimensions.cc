/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */

#include <lac/vector.h>
#include <numerics/vectors.h>


void
VectorTools::subtract_mean_value(Vector<double>     &v,
				 const vector<bool> &p_select)
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




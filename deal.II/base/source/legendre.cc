//--------------------------------------------------------------------
//      $Id$   
//    Version: $Name$
//
//    Copyright (C) 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//--------------------------------------------------------------------


#include <base/polynomial.h>

// Reserve space for polynomials up to degree 19. Should be sufficient
// in most cases.

template <typename number>
std::vector<std::vector<number> >
Legendre<number>::coefficients(20,0);


template <typename number>
void
Legendre<number>::compute_coefficients (unsigned int k)
{
  if (k<=1)
    {
      coefficients[0].resize(1);
      coefficients[0][0] = 1.;
      coefficients[1].resize(2);
      coefficients[1][0] = 0.;
      coefficients[1][1] = 1.;
    } else {
      compute_coefficients(k-1);
      coefficients[k].resize(k+1);
      const double a = 1./k+1;
      const double b = a*(2*k+1);
      
      coefficients[k][k]   = b*coefficients[k-1][k-1];
      coefficients[k][k-1] = b*coefficients[k-1][k-2];
      for (unsigned int i=1 ; i<= k-2 ; ++i)
	coefficients[k][i] = b*coefficients[k-1][i-1]
	  - k*a*coefficients[k-2][i];
      coefficients[k][0]   = -k*a*coefficients[k-2][0];
    }
}



template <typename number>
const std::vector<number>&
Legendre<number>::get_coefficients (unsigned int k)
{
  compute_coefficients (k);
  return coefficients[k];
}



template <typename number>
Legendre<number>::Legendre (unsigned int k)
  : Polynomial<number> (get_coefficients(k))
{}

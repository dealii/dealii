//----------------------------  polynomial.cc  -----------------------
//      $Id$   
//    Version: $Name$
//
//    Copyright (C) 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------   polynomial.cc  ----------------------


#include <base/polynomial.h>


Polynomial::Polynomial (const vector<double> &a):
		coefficients(a)
{}



Polynomial::Polynomial ()
  :
  coefficients(0)
{}



double Polynomial::value (const double x) const
{
  Assert (coefficients.size() > 0, ExcVoidPolynomial());
  const unsigned int m=coefficients.size();

				   // Horner scheme
  double value = coefficients.back();
  for (int k=m-2; k>=0; --k)
    value = value*x + coefficients[k];

  return value;
}



void Polynomial::value (const double    x,
			vector<double> &values) const
{
  Assert (coefficients.size() > 0, ExcVoidPolynomial());
  Assert (values.size() > 0, ExcEmptyArray());
  const unsigned int values_size=values.size();
  
  
				   // if we only need the value, then
				   // call the other function since
				   // that is significantly faster
				   // (there is no need to allocate
				   // and free memory, which is really
				   // expensive compared to all the
				   // other operations!)
  if (values_size == 1)
    {
      values[0] = value(x);
      return;
    };

				   // if there are derivatives needed,
				   // then do it properly by the
				   // full Horner scheme
  const unsigned int m=coefficients.size();
  vector<double> a(coefficients);
  unsigned int j_faculty=1;
  for (unsigned int j=0; j<values_size; ++j)
    {      
      for (int k=m-1; k>=static_cast<int>(j); --k)
	a[k]+=x*a[k+1];
      values[j]=j_faculty*a[j];

      j_faculty*=j+1;
    }
}


// ------------------------------------------------------------ //



LagrangeEquidistant::LagrangeEquidistant (const unsigned int n,
					  const unsigned int support_point):
		Polynomial(compute_coefficients(n,support_point))
{}



LagrangeEquidistant::LagrangeEquidistant ()
{}



vector<double> 
LagrangeEquidistant::compute_coefficients (const unsigned int n,
					   const unsigned int support_point)
{
  vector<double> a (n+1);
  Assert(support_point<n+1, ExcIndexRange(support_point, 0, n+1));

  switch (n)
    {
      case 0:
	    switch (support_point)
	      {
		case 0:
		      a[0]=1.;
		      break;
		default:
		      Assert(false, ExcInvalidSupportPoint(support_point));
	      }
	    break;
      case 1:
	    switch (support_point)
	      {
		case 0:
		      a[0]=1.;
		      a[1]=-1.;
		      break;
		case 1:
		      a[0]=0.;
		      a[1]=1.;
		      break;
		default:
		      Assert(false, ExcInvalidSupportPoint(support_point));
	      }
	    break;
      case 2:
	    switch (support_point)
	      {
		case 0:
		      a[0]=1.;
		      a[1]=-3.;
		      a[2]=2.;
		      break;
		case 1:
		      a[0]=0.;
		      a[1]=4.;		      
		      a[2]=-4.;
		      break;
		case 2:
		      a[0]=0.;
		      a[1]=-1.;		      
		      a[2]=2.;
		      break;
		default:
		      Assert(false, ExcInvalidSupportPoint(support_point));
	      }
	    break;
      case 3:
	    switch (support_point)
	      {
		case 0:
		      a[0]=1.0;
		      a[1]=-11.0/2.0;
		      a[2]=9.0;
		      a[3]=-9.0/2.0;
		      break;
		case 1:
		      a[0]=0.;
		      a[1]=9.0;		      
		      a[2]=-45.0/2.0;
		      a[3]=27.0/2.0;
		      break;
		case 2:
		      a[0]=0.;
		      a[1]=-9.0/2.0;		      
		      a[2]=18.0;
		      a[3]=-27.0/2.0;
		      break;
		case 3:
		      a[0]=0.;
		      a[1]=1.;		      
		      a[2]=-9.0/2.0;
		      a[3]=9.0/2.0;
		      break;
		default:
		      Assert(false, ExcInvalidSupportPoint(support_point));
	      }
	    break;
      case 4:
	    switch (support_point)
	      {
		case 0:
		      a[0]=1.;
		      a[1]=-25.0/3.0;
		      a[2]=70.0/3.0;
		      a[3]=-80.0/3.0;
		      a[4]=32.0/3.0;
		      break;
		case 1:
		      a[0]=0.;
		      a[1]=16.0;
		      a[2]=-208.0/3.0;
		      a[3]=96.0;
		      a[4]=-128.0/3.0;
		      break;
		case 2:
		      a[0]=0.;
		      a[1]=-12.0;
		      a[2]=76.0;
		      a[3]=-128.0;
		      a[4]=64.0;
		      break;
		case 3:
		      a[0]=0.;
		      a[1]=16.0/3.0;
		      a[2]=-112.0/3.0;
		      a[3]=224.0/3.0;
		      a[4]=-128.0/3.0;
		      break;
		case 4:
		      a[0]=0.;
		      a[1]=-1.;
		      a[2]=22.0/3.0;
		      a[3]=-16.0;
		      a[4]=32.0/3.0;
		      break;
		default:
		      Assert(false, ExcInvalidSupportPoint(support_point));
	      }
	    break;
      default:
	    Assert(false, ExcNotImplemented());
    }
  return a;
}

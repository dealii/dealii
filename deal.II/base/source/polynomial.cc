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


Polynomial::Polynomial(const vector<double> &a):
		coefficients(a)
{}


void Polynomial::value(double x, vector<double> &values) const
{
  const unsigned int m=coefficients.size();
  vector<double> a(coefficients);
				   // Horner scheme
  unsigned int j_faculty=1;
  for (unsigned int j=0; j<values.size(); ++j)
    {      
      for (int k=m-1; k>=static_cast<int>(j); --k)
	a[k]+=x*a[k+1];
      values[j]=j_faculty*a[j];

      j_faculty*=j+1;
    }
}


// ------------------------------------------------------------ //



LagrangeEquidistant::LagrangeEquidistant(unsigned int n, unsigned int support_point):
		Polynomial(compute_coefficients(n,support_point))
{}


vector<double> LagrangeEquidistant::compute_coefficients(unsigned int n, unsigned int support_point)
{
  vector<double> a;
  a.resize(n+1);
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
		      Assert(false, ExcInternalError());
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
		      Assert(false, ExcInternalError());
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
		      a[1]=-1.;		      
		      a[2]=2.;
		      break;
		case 2:
		      a[0]=0.;
		      a[1]=4.;		      
		      a[2]=-4.;
		      break;
		default:
		      Assert(false, ExcInternalError());
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
		      a[1]=1.;		      
		      a[2]=-9.0/2.0;
		      a[3]=9.0/2.0;
		      break;
		case 2:
		      a[0]=0.;
		      a[1]=9.0;		      
		      a[2]=-45.0/2.0;
		      a[3]=27.0/2.0;
		      break;
		case 3:
		      a[0]=0.;
		      a[1]=-9.0/2.0;		      
		      a[2]=18.0;
		      a[3]=-27.0/2.0;
		      break;
		default:
		      Assert(false, ExcInternalError());
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
		      a[1]=-1.;
		      a[2]=22.0/3.0;
		      a[3]=-16.0;
		      a[4]=32.0/3.0;
		      break;
		case 2:
		      a[0]=0.;
		      a[1]=16.0;
		      a[2]=-208.0/3.0;
		      a[3]=96.0;
		      a[4]=-128.0/3.0;
		      break;
		case 3:
		      a[0]=0.;
		      a[1]=-12.0;
		      a[2]=76.0;
		      a[3]=-128.0;
		      a[4]=64.0;
		      break;
		case 4:
		      a[0]=0.;
		      a[1]=16.0/3.0;
		      a[2]=-112.0/3.0;
		      a[3]=224.0/3.0;
		      a[4]=-128.0/3.0;
		      break;
		default:
		      Assert(false, ExcInternalError());
	      }
	    break;
      default:
	    Assert(false, ExcNotImplemented());
    }
  return a;
}

//----------------------------  polynomial.cc  -----------------------
//      $Id$   
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------   polynomial.cc  ----------------------


#include <base/polynomial.h>
#include <base/exceptions.h>
#include <base/thread_management.h>

#include <cmath>


// have a lock that guarantees that at most one thread is changing and
// accessing the @p{coefficients} arrays of classes implementing
// polynomials with tables. make this lock local to this file.
//
// having only one lock for all of these classes is probably not going
// to be a problem since we only need it on very rare occasions. if
// someone finds this is a bottleneck, feel free to replace it by a
// more fine-grained solution
namespace 
{
  Threads::ThreadMutex coefficients_lock;
};


// -------------------- class Polynomial ---------------- //


template <typename number>
Polynomial<number>::Polynomial (const std::vector<number> &a):
		coefficients(a)
{}



template <typename number>
number
Polynomial<number>::value (const number x) const
{
  Assert (coefficients.size() > 0, ExcVoidPolynomial());
  const unsigned int m=coefficients.size();

				   // Horner scheme
  number value = coefficients.back();
  for (int k=m-2; k>=0; --k)
    value = value*x + coefficients[k];

  return value;
}



template <typename number>
void
Polynomial<number>::value (const number         x,
			std::vector<number> &values) const
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
  std::vector<number> a(coefficients);
  unsigned int j_faculty=1;

				   // loop over all requested
				   // derivatives. note that
				   // derivatives @p{j>m} are
				   // necessarily zero, as they
				   // differentiate the polynomial
				   // more often than the highest
				   // power is
  const unsigned int min_valuessize_m=std::min(values_size, m);
  for (unsigned int j=0; j<min_valuessize_m; ++j)
    {
      for (int k=m-2; k>=static_cast<int>(j); --k)
	a[k]+=x*a[k+1];
      values[j]=j_faculty*a[j];

      j_faculty*=j+1;
    }

				   // fill higher derivatives by zero
  for (unsigned int j=min_valuessize_m; j<values_size; ++j)
    values[j] = 0;
}


template <typename number>
void
Polynomial<number>::scale(std::vector<number>& coefficients,
			   const number factor)
{
  double f = 1.;
  for (typename std::vector<number>::iterator c = coefficients.begin();
       c != coefficients.end(); ++c)
    {
      *c *= f;
      f *= factor;
    }  
}



template <typename number>
void
Polynomial<number>::scale(const number factor)
{
  scale (coefficients, factor);
}



template <typename number>
void
Polynomial<number>::multiply(std::vector<number>& coefficients,
			     const number factor)
{
  for (typename std::vector<number>::iterator c = coefficients.begin();
       c != coefficients.end(); ++c)
    *c *= factor;
}



template <typename number>
template <typename number2>
void
Polynomial<number>::shift(std::vector<number>& coefficients,
			  const number2 offset)
{  
#ifdef DEAL_II_LONG_DOUBLE_LOOP_BUG
  AssertThrow (false,
	       ExcMessage("Sorry, but the compiler you are using has a bug that disallows "
			  "compilation of this function, so you cannot use it. Read more "
			  "about the bug and when it occurs in the aclocal.m4 file in the "
			  "top level directory (watch for the string "
			  "DEAL_II_LONG_DOUBLE_LOOP_BUG)"));
				   // calm down warning for unused
				   // args. note that this code is
				   // actually unreachable
  coefficients[0] = offset;
#else  
				   // Copy coefficients to a vector of
				   // accuracy given by the argument
  std::vector<number2> new_coefficients(coefficients.begin(),
					coefficients.end());
  
				   // Traverse all coefficients from
				   // c_1. c_0 will be modified by
				   // higher degrees, only.
  for (unsigned int d=1; d<new_coefficients.size(); ++d)
    {
      const unsigned int n = d;
				       // Binomial coefficients are
				       // needed for the
				       // computation. The rightmost
				       // value is unity.
      unsigned int binomial_coefficient = 1;

				       // Powers of the offset will be
				       // needed and computed
				       // successively.
      number2 offset_power = offset;
      
				       // Compute (x+offset)^d
				       // and modify all values c_k
				       // with k<d.
				       // The coefficient in front of
				       // x^d is not modified in this step.
      for (unsigned int k=0;k<d;++k)
	{
					   // Recursion from Bronstein
					   // Make sure no remainders
					   // occur in integer
					   // division.
	  binomial_coefficient = (binomial_coefficient*(n-k))/(k+1);

	  new_coefficients[d-k-1] += new_coefficients[d]
				 * binomial_coefficient
				 * offset_power;
	  offset_power *= offset;
	}
				       // The binomial coefficient
				       // should have gone through a
				       // whole row of Pascal's
				       // triangle.
      Assert (binomial_coefficient == 1, ExcInternalError());
    }

				   // copy new elements to old vector
  coefficients.assign(new_coefficients.begin(), new_coefficients.end());
#endif
}


template <typename number>
template <typename number2>
void
Polynomial<number>::shift(const number2 offset)
{
  shift(coefficients, offset);
}


template <typename number>
void
Polynomial<number>::print(std::ostream& out) const
{
  for (int i=degree();i>=0;--i)
    {
      out << coefficients[i] << " x^" << i << std::endl;
    }
}



// ------------------ class LagrangeEquidistant --------------- //

LagrangeEquidistant::LagrangeEquidistant (const unsigned int n,
					  const unsigned int support_point):
		Polynomial<double>(compute_coefficients(n,support_point))
{}



std::vector<double> 
LagrangeEquidistant::compute_coefficients (const unsigned int n,
					   const unsigned int support_point)
{
  std::vector<double> a (n+1);
  Assert(support_point<n+1, ExcIndexRange(support_point, 0, n+1));

  unsigned int n_functions=n+1;
  Assert(support_point<n_functions,
	 ExcIndexRange(support_point, 0, n_functions));
  double const *x=0;
  
  switch (n)
    {
      case 1:
      {
	static const double x1[4]=
	{
	      1.0, -1.0,
	      0.0, 1.0
	};
	x=&x1[0];
	break;	
      }
      case 2:
      {
	static const double x2[9]=
	{
	      1.0, -3.0, 2.0,
	      0.0, 4.0, -4.0,
	      0.0, -1.0, 2.0
	};
	x=&x2[0];
	break;
      }
      case 3:
      {
	static const double x3[16]=
	{
	      1.0, -11.0/2.0, 9.0, -9.0/2.0,
	      0.0, 9.0, -45.0/2.0, 27.0/2.0,
	      0.0, -9.0/2.0, 18.0, -27.0/2.0,
	      0.0, 1.0, -9.0/2.0, 9.0/2.0
	};
	x=&x3[0];
	break;
      }
      case 4:
      {
	static const double x4[25]=
	{
	      1.0, -25.0/3.0, 70.0/3.0, -80.0/3.0, 32.0/3.0,
	      0.0, 16.0, -208.0/3.0, 96.0, -128.0/3.0,
	      0.0, -12.0, 76.0, -128.0, 64.0,
	      0.0, 16.0/3.0, -112.0/3.0, 224.0/3.0, -128.0/3.0,
	      0.0, -1.0, 22.0/3.0, -16.0, 32.0/3.0
	};
	x=&x4[0];
	break;	
      }
      case 5:
      {
	static const double x5[36]=
	{
	      1.0, -137.0/12.0, 375.0/8.0, -2125.0/24.0, 625.0/8.0, -625.0/24.0,
	      0.0, 25.0, -1925.0/12.0, 8875.0/24.0, -4375.0/12.0, 3125.0/24.0,
	      0.0, -25.0, 2675.0/12.0, -7375.0/12.0, 8125.0/12.0, -3125.0/12.0,
	      0.0, 50.0/3.0, -325.0/2.0, 6125.0/12.0, -625.0, 3125.0/12.0,
	      0.0, -25.0/4.0, 1525.0/24.0, -5125.0/24.0, 6875.0/24.0, -3125.0/24.0,
	      0.0, 1.0, -125.0/12.0, 875.0/24.0, -625.0/12.0, 625.0/24.0
	};
	x=&x5[0];
	break;
      }
      case 6:
      {
	static const double x6[49]=
	{
	      1.0, -147.0/10.0, 406.0/5.0, -441.0/2.0, 315.0, -1134.0/5.0,
	      324.0/5.0, 0.0, 36.0, -1566.0/5.0, 1044.0, -1674.0, 1296.0,
	      -1944.0/5.0, 0.0, -45.0, 1053.0/2.0, -4149.0/2.0, 3699.0, -3078.0,
	      972.0, 0.0, 40.0, -508.0, 2232.0, -4356.0, 3888.0, -1296.0, 0.0,
	      -45.0/2.0, 297.0, -2763.0/2.0, 2889.0, -2754.0, 972.0, 0.0,
	      36.0/5.0, -486.0/5.0, 468.0, -1026.0, 5184.0/5.0, -1944.0/5.0, 0.0,
	      -1.0, 137.0/10.0, -135.0/2.0, 153.0, -162.0, 324.0/5.0
	};
	x=&x6[0];
	break;
      }
      case 7:
      {
	static const double x7[64]=
	{
	      1.0, -363.0/20.0, 22981.0/180.0, -331681.0/720.0, 16807.0/18.0,
	      -386561.0/360.0, 117649.0/180.0, -117649.0/720.0, 0.0, 49.0,
	      -10927.0/20.0, 109417.0/45.0, -88837.0/16.0, 991613.0/144.0,
	      -352947.0/80.0, 823543.0/720.0, 0.0, -147.0/2.0, 43071.0/40.0,
	      -1347647.0/240.0, 170471.0/12.0, -151263.0/8.0, 1529437.0/120.0,
	      -823543.0/240.0, 0.0, 245.0/3.0, -46501.0/36.0, 133427.0/18.0,
	      -2926819.0/144.0, 4151329.0/144.0, -2941225.0/144.0,
	      823543.0/144.0, 0.0, -245.0/4.0, 2009.0/2.0, -872935.0/144.0,
	      52822.0/3.0, -1899191.0/72.0, 117649.0/6.0, -823543.0/144.0, 0.0,
	      147.0/5.0, -9849.0/20.0, 45962.0/15.0, -444185.0/48.0,
	      1159683.0/80.0, -2705927.0/240.0, 823543.0/240.0, 0.0, -49.0/6.0,
	      49931.0/360.0, -634207.0/720.0, 98441.0/36.0, -319333.0/72.0,
	      1294139.0/360.0, -823543.0/720.0, 0.0, 1.0, -343.0/20.0,
	      9947.0/90.0, -16807.0/48.0, 84035.0/144.0, -117649.0/240.0,
	      117649.0/720.0
	};
	x=&x7[0];
	break;
      }
      case 8:
      {
	static const double x8[81]=
	{
	      1.0, -761.0/35.0, 59062.0/315.0, -4272.0/5.0, 34208.0/15.0,
	      -18432.0/5.0, 53248.0/15.0, -65536.0/35.0, 131072.0/315.0, 0.0,
	      64.0, -30784.0/35.0, 44672.0/9.0, -673792.0/45.0, 235520.0/9.0,
	      -1196032.0/45.0, 131072.0/9.0, -1048576.0/315.0, 0.0, -112.0,
	      9936.0/5.0, -587296.0/45.0, 1956992.0/45.0, -733184.0/9.0,
	      3915776.0/45.0, -2228224.0/45.0, 524288.0/45.0, 0.0, 448.0/3.0,
	      -128192.0/45.0, 102016.0/5.0, -1097728.0/15.0, 145408.0,
	      -2441216.0/15.0, 1441792.0/15.0, -1048576.0/45.0, 0.0, -140.0,
	      2764.0, -186496.0/9.0, 703552.0/9.0, -1466368.0/9.0, 1712128.0/9.0,
	      -1048576.0/9.0, 262144.0/9.0, 0.0, 448.0/5.0, -9024.0/5.0,
	      626048.0/45.0, -2443264.0/45.0, 5285888.0/45.0, -6406144.0/45.0,
	      4063232.0/45.0, -1048576.0/45.0, 0.0, -112.0/3.0, 34288.0/45.0,
	      -5984.0, 358784.0/15.0, -53248.0, 999424.0/15.0, -131072.0/3.0,
	      524288.0/45.0, 0.0, 64.0/7.0, -6592.0/35.0, 67456.0/45.0,
	      -274432.0/45.0, 124928.0/9.0, -802816.0/45.0, 3801088.0/315.0,
	      -1048576.0/315.0, 0.0, -1.0, 726.0/35.0, -7504.0/45.0, 30944.0/45.0,
	      -14336.0/9.0, 94208.0/45.0, -65536.0/45.0, 131072.0/315.0
	};
	x=&x8[0];
	break;
      }
      case 9:
      {
	static const double x9[100]=
	{
	      1.0, -7129.0/280.0, 58635.0/224.0, -40707.0/28.0, 623295.0/128.0,
	      -6589431.0/640.0, 885735.0/64.0, -5137263.0/448.0, 4782969.0/896.0,
	      -4782969.0/4480.0, 0.0, 81.0, -373329.0/280.0, 10307331.0/1120.0,
	      -5589243.0/160.0, 51221727.0/640.0, -4546773.0/40.0,
	      31355019.0/320.0, -52612659.0/1120.0, 43046721.0/4480.0, 0.0,
	      -162.0, 475389.0/140.0, -15190173.0/560.0, 18152829.0/160.0,
	      -44529507.0/160.0, 33244587.0/80.0, -3720087.0/10.0,
	      205667667.0/1120.0, -43046721.0/1120.0, 0.0, 252.0, -56601.0/10.0,
	      1959363.0/40.0, -8776431.0/40.0, 91020753.0/160.0,
	      -71035947.0/80.0, 16474671.0/20.0, -33480783.0/80.0,
	      14348907.0/160.0, 0.0, -567.0/2.0, 526419.0/80.0, -4752351.0/80.0,
	      89119521.0/320.0, -241241409.0/320.0, 195629337.0/160.0,
	      -187598673.0/160.0, 196101729.0/320.0, -43046721.0/320.0, 0.0,
	      1134.0/5.0, -21465.0/4.0, 795339.0/16.0, -3844017.0/16.0,
	      215023653.0/320.0, -18009945.0/16.0, 35606547.0/32.0,
	      -4782969.0/8.0, 43046721.0/320.0, 0.0, -126.0, 60381.0/20.0,
	      -2276289.0/80.0, 22480173.0/160.0, -64448703.0/160.0,
	      55447011.0/80.0, -28166373.0/40.0, 62178597.0/160.0,
	      -14348907.0/160.0, 0.0, 324.0/7.0, -78327.0/70.0, 2989629.0/280.0,
	      -2142531.0/40.0, 25043337.0/160.0, -22025277.0/80.0,
	      80247591.0/280.0, -90876411.0/560.0, 43046721.0/1120.0, 0.0,
	      -81.0/8.0, 275967.0/1120.0, -1328967.0/560.0, 7712091.0/640.0,
	      -22878207.0/640.0, 20490003.0/320.0, -21789081.0/320.0,
	      176969853.0/4480.0, -43046721.0/4480.0, 0.0, 1.0, -6849.0/280.0,
	      265779.0/1120.0, -194643.0/160.0, 2337903.0/640.0, -531441.0/80.0,
	      2302911.0/320.0, -4782969.0/1120.0, 4782969.0/4480.0
	};
	x=&x9[0];
	break;
      }
      case 10:
      {
	static const double x10[121]=
	{
	      1.0, -7381.0/252.0, 177133.0/504.0, -10511875.0/4536.0,
	      42711625.0/4536.0, -5369375.0/216.0, 4695625.0/108.0,
	      -9453125.0/189.0, 6875000.0/189.0, -8593750.0/567.0,
	      1562500.0/567.0, 0.0, 100.0, -121525.0/63.0, 1997825.0/126.0,
	      -82992625.0/1134.0, 3775625.0/18.0, -20965625.0/54.0,
	      4187500.0/9.0, -65937500.0/189.0, 3125000.0/21.0,
	      -15625000.0/567.0, 0.0, -225.0, 153025.0/28.0, -2898075.0/56.0,
	      33095875.0/126.0, -57981875.0/72.0, 56396875.0/36.0,
	      -17546875.0/9.0, 94843750.0/63.0, -41406250.0/63.0, 7812500.0/63.0,
	      0.0, 400.0, -654100.0/63.0, 20028950.0/189.0, -108434750.0/189.0,
	      16686250.0/9.0, -33868750.0/9.0, 43625000.0/9.0, -242500000.0/63.0,
	      325000000.0/189.0, -62500000.0/189.0, 0.0, -525.0, 168775.0/12.0,
	      -1792225.0/12.0, 91073375.0/108.0, -102070625.0/36.0,
	      107321875.0/18.0, -71281250.0/9.0, 19375000.0/3.0, -26562500.0/9.0,
	      15625000.0/27.0, 0.0, 504.0, -13754.0, 149625.0, -7818625.0/9.0,
	      27074375.0/9.0, -58608125.0/9.0, 80000000.0/9.0, -66875000.0/9.0,
	      31250000.0/9.0, -6250000.0/9.0, 0.0, -350.0, 174025.0/18.0,
	      -11544725.0/108.0, 34178875.0/54.0, -80666875.0/36.0,
	      89384375.0/18.0, -62468750.0/9.0, 5937500.0, -76562500.0/27.0,
	      15625000.0/27.0, 0.0, 1200.0/7.0, -100300.0/21.0, 1121950.0/21.0,
	      -60659750.0/189.0, 10401250.0/9.0, -7831250.0/3.0,
	      234625000.0/63.0, -205000000.0/63.0, 100000000.0/63.0,
	      -62500000.0/189.0, 0.0, -225.0/4.0, 88325.0/56.0, -996675.0/56.0,
	      54486625.0/504.0, -28405625.0/72.0, 32584375.0/36.0,
	      -11828125.0/9.0, 73750000.0/63.0, -36718750.0/63.0, 7812500.0/63.0,
	      0.0, 100.0/9.0, -6575.0/21.0, 4033825.0/1134.0, -24717625.0/1134.0,
	      4341875.0/54.0, -10090625.0/54.0, 7437500.0/27.0,
	      -47187500.0/189.0, 71875000.0/567.0, -15625000.0/567.0, 0.0, -1.0,
	      7129.0/252.0, -162875.0/504.0, 1130750.0/567.0, -59375.0/8.0,
	      1883125.0/108.0, -78125.0/3.0, 4531250.0/189.0, -781250.0/63.0,
	      1562500.0/567.0
	};
	x=&x10[0];
	break;
      }
      default:
	    Assert(false, ExcNotImplemented());
    }

  Assert(x!=0, ExcInternalError());
  for (unsigned int i=0; i<n_functions; ++i)
    a[i]=x[support_point*n_functions+i];
  
  return a;
}


std::vector<Polynomial<double> >
LagrangeEquidistant::
generate_complete_basis (const unsigned int degree)
{
  if (degree==0)
				     // create constant polynomial
    return std::vector<Polynomial<double> >
      (1, Polynomial<double> (std::vector<double> (1,1.)));
  else
    {
				       // create array of Lagrange
				       // polynomials
      std::vector<Polynomial<double> > v;
      for (unsigned int i=0; i<=degree; ++i)
	v.push_back(LagrangeEquidistant(degree,i));
      return v;
    };
};



// ------------------ class Legendre --------------- //


//TODO:[?] This class leaks memory, but only at the very end of a program.
// Since it expands the Legendre<number>::coefficients array, the elements
// of this static variable are not destroyed at the end of the program
// run. While this is not a problem (since the returned memory could
// not be used anyway then), it is a little confusing when looking at
// a memory checker such as "purify". Maybe, this should be handled somehow
// to avoid this confusion in future.

// Reserve space for polynomials up to degree 19. Should be sufficient
// for the start.
template <typename number>
std::vector<const std::vector<number> *>
Legendre<number>::recursive_coefficients(
  20, static_cast<const std::vector<number>*>(0));
template <typename number>
std::vector<const std::vector<number> *>
Legendre<number>::shifted_coefficients(
  20, static_cast<const std::vector<number>*>(0));



//TODO[WB]: Treat the same way as above
#if ((__GNUC__ == 3) && (__GNUC_MINOR__ == 1))
#define SHIFT_TYPE double
#else
#define SHIFT_TYPE long double
#endif

template <typename number>
void
Legendre<number>::compute_coefficients (const unsigned int k_)
{
  unsigned int k = k_;

				   // first make sure that no other
				   // thread intercepts the operation
				   // of this function
  coefficients_lock.acquire ();

				   // The first 2 coefficients are hard-coded
  if (k==0)
    k=1;
				   // check: does the information
				   // already exist?
  if ((recursive_coefficients.size() < k+1) ||
      ((recursive_coefficients.size() >= k+1) &&
       (recursive_coefficients[k] == 0)))
				     // no, then generate the
				     // respective coefficients
    {
      recursive_coefficients.resize (k+1, 0);
      
      if (k<=1)
	{
					   // create coefficients
					   // vectors for k=0 and k=1
					   //
					   // allocate the respective
					   // amount of memory and
					   // later assign it to the
					   // coefficients array to
					   // make it const
	  std::vector<number> *c0 = new std::vector<number>(1);
	  (*c0)[0] = 1.;

	  std::vector<number> *c1 = new std::vector<number>(2);
	  (*c1)[0] = 0.;
	  (*c1)[1] = 1.;

					   // now make these arrays
					   // const
	  recursive_coefficients[0] = c0;
	  recursive_coefficients[1] = c1;
					   // Compute polynomials
					   // orthogonal on [0,1]
  	  c0 = new std::vector<number>(*c0);
  	  c1 = new std::vector<number>(*c1);
	  
    	  Polynomial<number>::shift(*c0, (SHIFT_TYPE) -1.);
    	  Polynomial<number>::scale(*c0, 2.);
    	  Polynomial<number>::shift(*c1, (SHIFT_TYPE) -1.);
    	  Polynomial<number>::scale(*c1, 2.);
    	  Polynomial<number>::multiply(*c1, std::sqrt(3.));
  	  shifted_coefficients[0]=c0;
  	  shifted_coefficients[1]=c1;
	}
      else
	{
					   // for larger numbers,
					   // compute the coefficients
					   // recursively. to do so,
					   // we have to release the
					   // lock temporarily to
					   // allow the called
					   // function to acquire it
					   // itself
	  coefficients_lock.release ();
	  compute_coefficients(k-1);
	  coefficients_lock.acquire ();

	  std::vector<number> *ck = new std::vector<number>(k+1);
	  
	  const number a = 1./(k);
	  const number b = a*(2*k-1);
	  const number c = a*(k-1);
	  
	  (*ck)[k]   = b*(*recursive_coefficients[k-1])[k-1];
	  (*ck)[k-1] = b*(*recursive_coefficients[k-1])[k-2];
	  for (unsigned int i=1 ; i<= k-2 ; ++i)
	    (*ck)[i] = b*(*recursive_coefficients[k-1])[i-1]
		       -c*(*recursive_coefficients[k-2])[i];

	  (*ck)[0]   = -c*(*recursive_coefficients[k-2])[0];

					   // finally assign the newly
					   // created vector to the
					   // const pointer in the
					   // coefficients array
	  recursive_coefficients[k] = ck;
					   // and compute the
					   // coefficients for [0,1]
  	  ck = new std::vector<number>(*ck);
    	  shift(*ck,(SHIFT_TYPE) -1.);
    	  Polynomial<number>::scale(*ck, 2.);
    	  Polynomial<number>::multiply(*ck, std::sqrt(2.*k+1.));
  	  shifted_coefficients[k] = ck;
	};
    };

				   // now, everything is done, so
				   // release the lock again
  coefficients_lock.release ();
}



template <typename number>
const std::vector<number> &
Legendre<number>::get_coefficients (const unsigned int k)
{
				   // first make sure the coefficients
				   // get computed if so necessary
  compute_coefficients (k);

				   // then get a pointer to the array
				   // of coefficients. do that in a MT
				   // safe way
  coefficients_lock.acquire ();
  const std::vector<number> *p = shifted_coefficients[k];
  coefficients_lock.release ();

				   // return the object pointed
				   // to. since this object does not
				   // change any more once computed,
				   // this is MT safe
  return *p;
}



template <typename number>
Legendre<number>::Legendre (const unsigned int k)
		:
		Polynomial<number> (get_coefficients(k))
{}



template <typename number>
std::vector<Polynomial<number> >
Legendre<number>::generate_complete_basis (const unsigned int degree)
{
  std::vector<Polynomial<double> > v;
  v.reserve(degree+1);
  for (unsigned int i=0; i<=degree; ++i)
    v.push_back (Legendre<double>(i));
  return v;
};


// ------------------ explicit instantiations --------------- //

template class Polynomial<float>;
template class Polynomial<double>;
template class Polynomial<long double>;

template void Polynomial<float>::shift(const float offset);
template void Polynomial<float>::shift(const double offset);
template void Polynomial<double>::shift(const double offset);
template void Polynomial<long double>::shift(const long double offset);
template void Polynomial<float>::shift(const long double offset);
template void Polynomial<double>::shift(const long double offset);

template class Legendre<double>;

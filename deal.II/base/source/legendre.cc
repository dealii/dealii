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
#include <base/thread_management.h>


// Reserve space for polynomials up to degree 19. Should be sufficient
// for the start.
template <typename number>
std::vector<const std::vector<number>*>
Legendre<number>::coefficients(20,
			       static_cast<const std::vector<number>*>(0));


// have a lock that guarantees that at most one thread is changing and
// accessing the @p{coefficients} array. make this lock local to this
// file
namespace 
{
  Threads::ThreadMutex coefficients_lock;
};



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
  if ((coefficients.size() < k+1) ||
      ((coefficients.size() >= k+1) &&
       (coefficients[k] == 0)))
				     // no, then generate the
				     // respective coefficients
    {
      coefficients.resize (k+1, 0);
      
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
	  coefficients[0] = c0;
	  coefficients[1] = c1;
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
	  
	  (*ck)[k]   = b*(*coefficients[k-1])[k-1];
	  (*ck)[k-1] = b*(*coefficients[k-1])[k-2];
	  for (unsigned int i=1 ; i<= k-2 ; ++i)
	    (*ck)[i] = b*(*coefficients[k-1])[i-1]
		       -c*(*coefficients[k-2])[i];

	  (*ck)[0]   = -c*(*coefficients[k-2])[0];

					   // finally assign the newly
					   // created vector to the
					   // const pointer in the
					   // coefficients array
	  coefficients[k] = ck;
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
  const std::vector<number> *p = coefficients[k];
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



// explicit instantiations
template class Legendre<double>;

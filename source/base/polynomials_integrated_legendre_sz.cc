// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include <deal.II/base/polynomials_integrated_legendre_sz.h>

DEAL_II_NAMESPACE_OPEN

// Reserve space for polynomials up to degree 19.
std::vector<std::shared_ptr<const std::vector<double>>> IntegratedLegendreSZ::recursive_coefficients(20);

// Define the static mutex member.
Threads::Mutex IntegratedLegendreSZ::coefficients_lock;


IntegratedLegendreSZ::IntegratedLegendreSZ (const unsigned int k)
  :
  Polynomials::Polynomial<double> (get_coefficients(k))
{}



void IntegratedLegendreSZ::compute_coefficients (const unsigned int k_)
{
  unsigned int k = k_;

  // first make sure that no other thread intercepts the operation of this function;
  // for this, acquire the lock until we quit this function
  Threads::Mutex::ScopedLock lock(coefficients_lock);

  // The first 2 coefficients are hard-coded
  if (k==0)   k=1;


  // check: does the information already exist?
  if ((recursive_coefficients.size() < k+1) ||
      ((recursive_coefficients.size() >= k+1) &&
       (recursive_coefficients[k] == std::shared_ptr<const std::vector<double> >())))
    // no, then generate the respective coefficients
    {
      // make sure that there is enough space in the array for the coefficients,
      // so we have to resize it to size k+1

      // but it's more complicated than that: we call this function recursively, so if we simply
      // resize it to k+1 here, then compute the coefficients for degree k-1 by calling this
      // function recursively, then it will reset the size to k -- not enough for what we want to do below. the
      // solution therefore is to only resize the size if we are going to *increase* it
      if (recursive_coefficients.size() < k+1)
        {
          recursive_coefficients.resize (k+1);
        }
      if (k<=1)
        {
          // create coefficients vectors for k=0 and k=1
          //
          // allocate the respective later assign it to the coefficients array to make it const
          std::vector<double> *c0 = new std::vector<double>(1);
          (*c0)[0] = -1.;

          std::vector<double> *c1 = new std::vector<double>(2);
          (*c1)[0] = 0.;
          (*c1)[1] = 1.;

          // now make these arrays const. use shared_ptr for recursive_coefficients because
          // that avoids a memory leak that would appear if we used plain pointers.
          recursive_coefficients[0] = std::shared_ptr<const std::vector<double> >(c0);
          recursive_coefficients[1] = std::shared_ptr<const std::vector<double> >(c1);

        }
      else
        {
          // for larger numbers, compute the coefficients recursively. to do so, we have to release the
          // lock temporarily to allow the called function to acquire it itself
          coefficients_lock.release ();
          compute_coefficients(k-1);
          coefficients_lock.acquire ();

          std::vector<double> *ck = new std::vector<double>(k+1);

          const double a = 1.0 / k;
          const double b = 2.0*k - 3.0;
          const double c = k - 3.0;

          // To maintain stability, delay the division (multiplication by a) until the end.

          (*ck)[k]   = b*(*recursive_coefficients[k-1])[k-1];
          (*ck)[k-1] = b*(*recursive_coefficients[k-1])[k-2];
          for (unsigned int i=1; i<= k-2 ; ++i)
            {
              (*ck)[i] = b*(*recursive_coefficients[k-1])[i-1] - c*(*recursive_coefficients[k-2])[i];
            }

          (*ck)[0] = -c*(*recursive_coefficients[k-2])[0];

          for (unsigned int i=0; i<ck->size(); i++)
            {
              (*ck)[i] *=a;
            }

          // finally assign the newly created vector to the const pointer in the/ coefficients array
          recursive_coefficients[k] = std::shared_ptr<const std::vector<double> >(ck);
        }
    }
}



const std::vector<double> &IntegratedLegendreSZ::get_coefficients (const unsigned int k)
{
  // first make sure the coefficients get computed if so necessary
  compute_coefficients (k);

  // then get a pointer to the array of coefficients. do that in a MT safe way
  Threads::Mutex::ScopedLock lock (coefficients_lock);
  return *recursive_coefficients[k];
}



std::vector<Polynomials::Polynomial<double> >
IntegratedLegendreSZ::generate_complete_basis (const unsigned int degree)
{
  std::vector<Polynomials::Polynomial<double> > v;
  v.reserve(degree + 1);
  for (unsigned int i=0; i<=degree; ++i)
    {
      v.push_back (IntegratedLegendreSZ(i));
    }
  return v;
}



DEAL_II_NAMESPACE_CLOSE

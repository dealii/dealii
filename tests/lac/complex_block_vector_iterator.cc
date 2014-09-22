// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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



#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/block_vector.h>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <numeric>
#include <utility>
#include <cmath>

template <typename number>
bool operator == (const BlockVector<number> &v1,
                  const BlockVector<number> &v2)
{
  if (v1.size() != v2.size())
    return false;
  for (unsigned int i=0; i<v1.size(); ++i)
    if (v1(i) != v2(i))
      return false;
  return true;
}


void test ()
{
  std::vector<types::global_dof_index> ivector(4);
  ivector[0] = 2;
  ivector[1] = 4;
  ivector[2] = 3;
  ivector[3] = 5;

  // Check 1: initialization via
  // iterators
  if (true)
    {
      BlockVector<std::complex<double> > v1(ivector);
      BlockVector<std::complex<double> > v2(ivector);

      // initialize first vector with
      // simple loop
      for (unsigned int i=0; i<v1.size(); ++i)
        v1(i) = std::complex<double>(i,i+1);
      // initialize other vector
      // through iterators
      BlockVector<std::complex<double> >::iterator p2 = v2.begin();
      for (unsigned int i=0; i<v1.size(); ++i, ++p2)
        *p2 = std::complex<double>(i,i+1);
      Assert (p2==v2.end(), ExcInternalError());

      // check that the two vectors are equal
      deallog << "Check 1: " << (v1==v2 ? "true" : "false") << std::endl;
    };

  // Check 1: initialization via
  // iterators
  if (true)
    {
      BlockVector<std::complex<double> > v1(ivector);

      // initialize first vector with
      // simple loop
      for (unsigned int i=0; i<v1.size(); ++i)
        v1(i) = std::complex<double>(i,i+1);
      // initialize other vector
      // through iterators into first
      // vector
      BlockVector<std::complex<double> > v2(ivector, v1.begin(), v1.end());
      // check that the two vectors are equal
      deallog << "Check 2: " << (v1==v2 ? "true" : "false") << std::endl;
    };

  // Check 3: loop forward and back
  // and check that things are the
  // same
  if (true)
    {
      BlockVector<std::complex<double> > v1(ivector);
      // initialize first vector with
      // simple loop
      for (unsigned int i=0; i<v1.size(); ++i)
        v1(i) = std::complex<double>(i,i+1);

      BlockVector<std::complex<double> >::iterator p1 = v1.begin();
      for (unsigned int i=0; i<v1.size(); ++i, ++p1)
        Assert (*p1 == std::complex<double>(i,i+1), ExcInternalError());

      Assert (p1 == v1.end(), ExcInternalError());

      // move back into allowable
      // region
      --p1;

      // check backwards
      for (unsigned int i=0; i<v1.size(); ++i, --p1)
        Assert (*p1 == std::complex<double>(v1.size()-i-1, v1.size()-i), ExcInternalError());

      // if we came thus far,
      // everything is alright
      deallog << "Check 3: true" << std::endl;
    };


  // Check 4: same, but this time
  // with const iterators
  if (true)
    {
      BlockVector<std::complex<double> > v1(ivector);
      // initialize first vector with
      // simple loop
      for (unsigned int i=0; i<v1.size(); ++i)
        v1(i) = i;

      BlockVector<std::complex<double> >::const_iterator p1 = v1.begin();
      for (unsigned int i=0; i<v1.size(); ++i, ++p1)
        Assert (*p1 == std::complex<double>(i,0), ExcInternalError());

      Assert (p1 == v1.end(), ExcInternalError());

      // move back into allowable
      // region
      --p1;

      // check backwards
      for (unsigned int i=0; i<v1.size(); ++i, --p1)
        {
          const std::complex<double> val = *p1;
          const std::complex<double> ref(v1.size()-i-1);
          Assert (val==ref, ExcInternalError());
        };

      // if we came thus far,
      // everything is alright
      deallog << "Check 4: true" << std::endl;
    };

  // Checks 5-14: use some standard
  // algorithms
  if (true)
    {
      BlockVector<std::complex<double> > v1(ivector);
      // initialize first vector with
      // simple loop
      for (unsigned int i=0; i<v1.size(); ++i)
        v1(i) = i;

      // check std::distance
      // algorithm
      deallog << "Check 5: "
              << (std::distance (v1.begin(), v1.end()) ==
                  static_cast<signed int>(v1.size()) ?
                  "true" : "false")
              << std::endl;

      // check std::copy
      BlockVector<std::complex<double> > v2(ivector);
      std::copy (v1.begin(), v1.end(), v2.begin());
      deallog << "Check 6: " << (v1 == v2 ? "true" : "false") << std::endl;

      // check std::transform
      std::transform (v1.begin(), v1.end(), v2.begin(),
                      std::bind2nd (std::multiplies<std::complex<double> >(),
                                    2.0));
      v2.scale (std::complex<double>(1./2.));
      deallog << "Check 7: " << (v1 == v2 ? "true" : "false") << std::endl;


      // check operators +/-, +=/-=
      deallog << "Check 8: "
              << (std::distance(v1.begin(), v1.begin()+3) == 3 ?
                  "true" : "false")
              << std::endl;
      deallog << "Check 9: "
              << (std::distance(v1.end()-6, v1.end()) == 6 ?
                  "true" : "false")
              << std::endl;
      deallog << "Check 10: "
              << (std::distance(v1.begin(), v1.end()) == (signed)v1.size() ?
                  "true" : "false")
              << std::endl;
      deallog << "Check 11: "
              << (std::distance(v1.begin(), (v1.begin()+=7)) == 7 ?
                  "true" : "false")
              << std::endl;
      deallog << "Check 12: "
              << (std::distance((v1.end()-=4), v1.end()) == 4 ?
                  "true" : "false")
              << std::endl;

      // check advance
      BlockVector<std::complex<double> >::iterator p2 = v1.begin();
      std::advance (p2, v1.size());
      deallog << "Check 13: " << (p2 == v1.end() ? "true" : "false") << std::endl;

      BlockVector<std::complex<double> >::const_iterator p3 = v1.begin();
      std::advance (p3, v1.size());
      deallog << "Check 14: " << (p3 == v1.end() ? "true" : "false") << std::endl;
    };

  // Check 15: initialization through
  // iterators
  if (true)
    {
      BlockVector<std::complex<double> > v1(ivector);
      // initialize first vector with
      // simple loop
      for (unsigned int i=0; i<v1.size(); ++i)
        v1(i) = i;

      // initialize a normal vector
      // from it
      Vector<std::complex<double> > v2(v1.begin(), v1.end());

      // and reverse way
      BlockVector<std::complex<double> > v3(ivector, v2.begin(), v2.end());
      deallog << "Check 15: " << (v1==v3 ? "true" : "false") << std::endl;
    };

  // Check 16: initialization through
  // iterators. variant with one
  // constant object
  if (true)
    {
      BlockVector<std::complex<double> > v1(ivector);
      // initialize first vector with
      // simple loop
      for (unsigned int i=0; i<v1.size(); ++i)
        v1(i) = i;

      // initialize a normal vector
      // from it
      const Vector<std::complex<double> > v2(v1.begin(), v1.end());

      // and reverse way
      BlockVector<std::complex<double> > v3(ivector, v2.begin(), v2.end());
      deallog << "Check 16: " << (v1==v3 ? "true" : "false") << std::endl;
    };

  // Check 17: initialization through
  // iterators. variant with two
  // constant object
  if (true)
    {
      BlockVector<std::complex<double> > v1(ivector);
      // initialize first vector with
      // simple loop
      for (unsigned int i=0; i<v1.size(); ++i)
        v1(i) = i;

      // initialize a normal vector
      // from it
      const Vector<std::complex<double> > v2(v1.begin(), v1.end());

      // and reverse way
      const BlockVector<std::complex<double> > v3(ivector, v2.begin(), v2.end());
      deallog << "Check 17: " << (v1==v3 ? "true" : "false") << std::endl;
    };

  // Check 18: initialization through
  // iterators. variant with three
  // constant object
  if (true)
    {
      BlockVector<std::complex<double> > v0(ivector);
      // initialize first vector with
      // simple loop
      for (unsigned int i=0; i<v0.size(); ++i)
        v0(i) = i;

      const BlockVector<std::complex<double> > v1 = v0;

      // initialize a normal vector
      // from it
      const Vector<std::complex<double> > v2(v1.begin(), v1.end());

      // and reverse way
      const BlockVector<std::complex<double> > v3(ivector, v2.begin(), v2.end());
      deallog << "Check 18: " << (v1==v3 ? "true" : "false") << std::endl;
    };

  // Check 19: operator[]
  if (true)
    {
      BlockVector<std::complex<double> > v1(ivector);
      for (unsigned int i=0; i<v1.size(); ++i)
        v1(i) = i;

      for (unsigned int i=0; i<v1.size(); ++i)
        {
          const BlockVector<std::complex<double> >::iterator p = (v1.begin()+i);
          for (unsigned int j=0; j<v1.size(); ++j)
            Assert (p[(signed)j-(signed)i] == std::complex<double>(j), ExcInternalError());
        };

      // if we came thus far,
      // everything is alright
      deallog << "Check 19: true" << std::endl;
    };
}





int main ()
{
  std::ofstream logfile("output");
  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      test ();
    }
  catch (std::exception &e)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << e.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      // abort
      return 2;
    }
  catch (...)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      // abort
      return 3;
    };


  return 0;
}


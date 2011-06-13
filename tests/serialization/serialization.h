//----------------------------  serialization.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  serialization.h  ---------------------------

// common include file for all serialization tests

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <sstream>
#include <fstream>
#include <iomanip>


// compare objects for equality and pointers for equality of the object
// pointed to
template <typename T>
bool compare (T t1,
	      T t2)
{
  return t1 == t2;
}

template <typename T>
bool compare (T *t1,
	      T *t2)
{
  return *t1 == *t2;
}


template <typename T>
void verify (const T &t1,
	     T       &t2)
{
				   // save data to archive
  std::ostringstream oss;
  {
    boost::archive::text_oarchive oa(oss);
    oa << t1;
				     // archive and stream closed when
				     // destructors are called
  }
  deallog << oss.str() << std::endl;

				   // verify correctness of the
				   // serialization
  {
    std::istringstream  iss(oss.str());
    boost::archive::text_iarchive ia(iss);


    ia >> t2;

    Assert (compare (t1, t2), ExcInternalError());
  }
}



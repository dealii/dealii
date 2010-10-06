//----------------------------  tensor_base_scalar.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tensor_base_scalar.cc  ---------------------------

// check serialization for Tensor<1,dim>

#include "../tests.h"
#include <base/tensor.h>
#include <base/logstream.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <sstream>
#include <fstream>
#include <iomanip>


template <typename T>
void verify (T &t1,
	     T &t2)
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

    Assert (t1 == t2, ExcInternalError());
  }
}


void test ()
{
  const unsigned int dim=3;

  Tensor<0,dim> t1(1.);

  Tensor<0,dim> t2(2.);

  verify (t1, t2);
}


int main ()
{
  std::ofstream logfile("tensor_base_scalar/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}

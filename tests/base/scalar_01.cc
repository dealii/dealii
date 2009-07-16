//----------------------------  scalar_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  scalar_01.cc  ---------------------------

// check Tensor<0,dim>

#include "../tests.h"
#include <base/tensor.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <fstream>
#include <iomanip>

template <typename U, typename V>
void compare (const U &u, const V &v)
{
  Assert (static_cast<double>(u)==static_cast<double>(v), ExcInternalError());
}



int main ()
{
  std::ofstream logfile("scalar_01/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  typedef Tensor<0,1> T;
  T t1(13.), t2(42);

  compare (T(), 0.);
  compare (T(13.), 13.);
  compare (T(t1), 13.);
  compare (static_cast<double>(t1), 13.);
  compare (static_cast<double&>(t1), 13.);
  compare ((T() = t1), 13.);
  compare ((T() = 13.), 13.);
  compare ((t1==t1), true);
  compare ((t1==t2), false);
  compare ((t1!=t2), true);
  compare ((t1!=t1), false);
  compare ((T() += t1), t1);
  compare ((T() -= t1), -t1);
  compare (T(13.)*=3., 39.);
  compare (T(39)/=3., 13.);
  compare ((t1*t2), 13*42);
  compare ((t1+t2), 13+42);
  compare ((t1-t2), 13-42);
  compare (-t1, -13.);
  compare (T(-12).norm(), 12.);
  compare (T(-12).norm_square(), 12*12.);

  t1.clear();
  compare (t1, 0.);
  
  deallog << "OK" << std::endl;
}

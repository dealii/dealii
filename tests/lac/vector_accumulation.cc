//-----------------------  vector_accumulation.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------  vector_accumulation.cc  ---------------------------

// check that the l2 norm is exactly the same for many runs on random vector
// size with random vector entries. this is to ensure that the result is
// reproducible also when parallel evaluation is done

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <cmath>
#include <fstream>
#include <iomanip>




template <typename number>
void check_norms ()
{
  for (unsigned int test=0; test<20; ++test)
    {
      const unsigned int size = rand() % 100000;
      Vector<number> vec (size);
      for (unsigned int i=0; i<size; ++i)
        vec(i) = static_cast<number>(rand())/static_cast<number>(RAND_MAX);
      const typename Vector<number>::real_type norm = vec.l2_norm();
      for (unsigned int i=0; i<30; ++i)
        Assert (vec.l2_norm() == norm, ExcInternalError());

      Vector<number> vec2 (vec);
      for (unsigned int i=0; i<10; ++i)
        Assert (vec2.l2_norm() == norm, ExcInternalError());
    }
}


int main()
{
  std::ofstream logfile("vector_accumulation/output");
  deallog << std::fixed;
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check_norms<float>();
  check_norms<double>();
  check_norms<long double>();
  check_norms<std::complex<double> >();
  deallog << "OK" << std::endl;
}



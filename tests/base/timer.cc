// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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
#include <deal.II/base/timer.h>
#include <deal.II/base/logstream.h>
#include <fstream>
#include <cmath>
#include <iomanip>

// compute the ratio of two measurements and compare to
// the expected value.

void compare (double t1, double t2, double ratio)
{
  double r = t2/t1;
  double d = std::fabs(r-ratio) / ratio;

  // relative error < 25%?
  if (d <= .25)
    {
      deallog << "OK" << std::endl;
    }
  else
    {
      deallog << "Ratio " << r << " should be " << ratio << std::endl;
    }
}

// burn computer time

double s = 0.;
void burn (unsigned int n)
{
  for (unsigned int i=0 ; i<n ; ++i)
    {
      for (unsigned int j=1 ; j<100000 ; ++j)
        {
          s += 1./j * i;
        }
    }
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Timer t1,t2;
  burn (50);
  double s01 = t1.stop();
  double s02 = t2();
  burn (50);
  double s11 = t1.stop();
  double s12 = t2();
  t1.start();
  burn (50);
  double s21 = t1();
  double s22 = t2();
  burn (50);
  double s31 = t1();
  double s32 = t2();

  compare (s01,s02,1.);
  compare (s11,s12,2.);
  compare (s21,s22,3./2.);
  compare (s31,s32,4./3.);
}


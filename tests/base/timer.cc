//----------------------------  timer.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  timer.cc  ---------------------------


#include <base/timer.h>
#include <base/logstream.h>
#include <fstream>
#include <cmath>
#include <iostream>

// compute the ratio of two measurements and compare to
// the expected value.

void compare (double t1, double t2, double ratio)
{
  double r = t2/t1;
  double d = std::fabs(r-ratio) / ratio;

				   // relative error < 10%?
  if (d <= .1)
    {
      deallog << "OK" << std::endl;
    } else {
      deallog << "Ratio " << r << " should be " << ratio << std::endl;
    }
}

// burn computer time

void burn (unsigned int n)
{
  double s = 0.;
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
  std::ofstream logfile("timer.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  Timer t1,t2;
  burn (20);
  double s01 = t1.stop();
  double s02 = t2();
  burn (20);
  double s11 = t1.stop();
  double s12 = t2();
  t1.start();
  burn (20);
  double s21 = t1();
  double s22 = t2();
  burn (20);
  double s31 = t1();
  double s32 = t2();

  compare (s01,s02,1.);
  compare (s11,s12,2.);
  compare (s21,s22,3./2.);
  compare (s31,s32,4./3.);
}


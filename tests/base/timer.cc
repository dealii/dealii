// $Id$
// (c) Guido Kanschat, 2000

#include <base/timer.h>
#include <base/logstream.h>
#include <fstream>

// compute the ratio of two measurements and compare to
// the expected value.

void compare (double t1, double t2, double ratio)
{
  double r = t2/t1;
  double d = fabs(r-ratio) / ratio;

				   // relative error < 10%?
  if (d <= .1)
    {
      deallog << "OK" << endl;
    } else {
      deallog << "Ratio " << r << " should be " << ratio << endl;
    }
}

// burn computer time

void burn (unsigned int n)
{
  double s;
  for (unsigned int i=0 ; i<n ; ++i)
    {
      for (unsigned int j=0 ; j<100000 ; ++j)
	{
	  s += 1./j * i;
	}
    }
}


int main ()
{
  ofstream logfile("timer.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  Timer t1,t2;
  burn (100);
  double s01 = t1.stop();
  double s02 = t2();
  burn (100);
  double s11 = t1.stop();
  double s12 = t2();
  t1.start();
  burn (100);
  double s21 = t1.stop();
  double s22 = t2();

  compare (s01,s02,1.);
  compare (s11,s12,2.);
  compare (s21,s22,3./2.);
}


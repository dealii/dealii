//----------------------------  solver.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  solver.cc  ---------------------------


#include <base/logstream.h>
#include <lac/vector.h>
#include <cmath>
#include <fstream>




const unsigned int N=50;
unsigned int check_point = 0;




template <typename number>
void print (const Vector<number> &v) 
{
  deallog << "Check point " << check_point << endl;
  check_point++;
  
  for (unsigned int i=0; i<v.size(); ++i)
    deallog << v(i) << ' ';
  deallog << endl;
};



template <typename number1, typename number2>
void check_vectors (Vector<number1> &d1, Vector<number2> &d2)
{
  for (unsigned int i=0; i<N; ++i)
    {
      d1(i) = 1. * i / 3;
      d2(i) = d1(i)*d1(i) / 2;
    };

  print (d1);
  print (d2);

  swap (d1, d2);
  print (d1);
  
  d1 = d2;
  print (d1);
  
  d1 = 2.871;
  print (d1);
  
  deallog << d1 * d2 << ' ' << d2.norm_sqr() << endl;
  deallog << d1.mean_value() << ' ' << d2.l1_norm() << endl;
  deallog << d1.l2_norm() << ' ' << d1.linfty_norm() << endl;

  d1 += d2;
  print (d1);
  
  d2 -= d1;
  print (d2);
  
  d1.add (2.54);
  print (d1);
  
  d1.add (6.7, d2);
  print (d1);
  
  d1.add (2.3, d2, 3.4, d2);
  print (d1);
  
  d2.sadd (1.1, d1);
  print (d2);
  
  d2.sadd (1.3, 1.7, d1);
  print (d2);
  
  d1.sadd (12, 17, d2, 14, d2);
  print (d1);
  
  d1.scale (3.14154);
  print (d1);
  
  d2.equ (1.569, d1);
  print (d2);
  
  d2.equ (1.876, d1, 1867, d1);
  print (d2);
  
  d1.ratio (d1, d2);
  print (d1);
};


int main()
{
  ofstream logfile("vector-vector.output");
  logfile.setf(ios::fixed);
  logfile.precision(5);
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  Vector<double>      d1(N), d2(N);
  Vector<float>       f1(N), f2(N);
  Vector<long double> l1(N), l2(N);

				   // cross-tests with double/float
				   // vectors at the same time are not
				   // supported at present,
				   // unfortunately, as many functions
				   // don't accept other data types as
				   // arguments
  check_vectors (d1, d2);
  check_vectors (f1, f2);
  check_vectors (l1, l2);
};

  

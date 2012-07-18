//-----------------------------------------------------------------------------
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
//-----------------------------------------------------------------------------

// test for AlignedVector<AlignedVector<unsigned int> >

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/vectorization.h>


typedef AlignedVector<unsigned int> VEC;
typedef AlignedVector<VEC> VECVEC;
void print_vec (VECVEC &v)
{
  for (unsigned int i=0; i<v.size(); ++i)
    {
      deallog << "[";
      for (unsigned int j=0; j<v[i].size(); ++j)
	deallog << v[i][j] << " ";
      deallog << "]";
    }
  deallog << std::endl;
}

void test ()
{
  typedef AlignedVector<unsigned int> VEC;
  VEC a(4);
  a[0] = 2;
  a[1] = 1;
  a[2] = 42;
  VECVEC v (2);
  deallog << "Constructor: ";
  print_vec(v);

  v[0] = a;
  v[1] = a;

  deallog << "Assignment: ";
  print_vec(v);

  VECVEC w (v);
  deallog << "Assignment vector: ";
  print_vec(w);
  deallog << "Data consistency after assignment: ";
  print_vec(v);

  a[1] = 41;
  a.push_back (100);
  v.push_back(a);
  deallog << "Insertion: ";
  print_vec (v);

  v.resize(1);
  deallog << "Shrinking: ";
  print_vec (v);

  v.reserve(100);
  deallog << "Reserve: ";
  print_vec (v);

  v.resize (10);
  deallog << "Resize: ";
  print_vec (v);

  v.resize(0);
  deallog << "Clear: ";
  print_vec (v);
}




int main()
{
  std::ofstream logfile("vectorization_03/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}

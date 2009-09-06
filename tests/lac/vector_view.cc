//----------------------------  solver.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  solver.cc  ---------------------------


#include "../tests.h"
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/vector_view.h>
#include <cmath>
#include <fstream>
#include <iomanip>

const unsigned int N=10;
unsigned int check_point = 0;

template <typename number>
void print (const Vector<number> &v) 
{
  for (unsigned int i=0; i<v.size(); ++i)
    deallog << v(i) << '\t';
  deallog << std::endl;
}

template<typename T>
void fill( T &a) {
    for(unsigned int i=0; i<a.size(); ++i)
	a(i) = i;
}

int main()
{
  std::ofstream logfile("vector_view/output");
  deallog << std::fixed;
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  Vector<double>  v1(N);
  fill(v1);
  
  deallog << "Vector" << std::endl;
  print(v1);
 
  VectorView<double> v2(N, v1.begin() );
  deallog << "Vector View" << std::endl;
  print(v2);
  
  v2(4) = 0;
  deallog << "Modified element 4" << std::endl;
  deallog << "Vector" << std::endl;
  print(v1);
 
  deallog << "Vector View" << std::endl;
  print(v2);
  
  // Const vector. 
  const Vector<double> v3(v1);
 
  deallog << "const Vector" << std::endl;
  print(v3);
 
  VectorView<double> v4(N, v3.begin());
  deallog << "const Vector View" << std::endl;
  print(v4); 
}

  

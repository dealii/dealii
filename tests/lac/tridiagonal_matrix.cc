//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include "../tests.h"
#include <base/logstream.h>
#include <lac/tridiagonal_matrix.h>
#include <lac/vector.h>

#include <fstream>
#include <iostream>
#include <cmath>

template<typename number>
void
check_vmult()
{
  TridiagonalMatrix<number> M(4);
  Vector<number> u(4);
  Vector<number> v(4);

  for (unsigned int i=0;i<4;++i)
    {
      u(i) = i+1;
      M(i,i) = i+1;
      if (i>0)
	M(i,i-1) = 0.-i;
      if (i<3)
	M(i,i+1) = 4.-i;
    }
  M.vmult(v,u);
  for (unsigned int i=0;i<v.size();++i)
    deallog << ' ' << v(i);
  deallog << std::endl;

  M.vmult_add(v,u);
  for (unsigned int i=0;i<v.size();++i)
    deallog << ' ' << v(i);
  deallog << std::endl;

  M.Tvmult(v,u);
  for (unsigned int i=0;i<v.size();++i)
    deallog << ' ' << v(i);
  deallog << std::endl;

  M.Tvmult_add(v,u);
  for (unsigned int i=0;i<v.size();++i)
    deallog << ' ' << v(i);
  deallog << std::endl;
}


int main()
{
  std::ofstream logfile("tridiagonal_matrix.output");
  logfile.setf(std::ios::fixed);
  logfile.precision(0);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  check_vmult<double>();
}

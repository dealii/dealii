//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

// check IdentityMatrix::vmult and friends


#include "../tests.h"
#include <base/logstream.h>
#include <lac/identity_matrix.h>
#include <lac/vector.h>

#include <fstream>
#include <iostream>
#include <cmath>

template<typename number>
void
check_vmult()
{
  IdentityMatrix M(4);
  Vector<number> u(4);
  Vector<number> v(4);

  for (unsigned int i=0;i<4;++i)
    u(i) = i+1;
  
  M.vmult(v,u);
  Assert (v == u, ExcInternalError());
  for (unsigned int i=0;i<v.size();++i)
    deallog << ' ' << v(i);
  deallog << std::endl;

  M.vmult_add(v,u);
  v /= 2;
  Assert (v == u, ExcInternalError());
  for (unsigned int i=0;i<v.size();++i)
    deallog << ' ' << v(i);
  deallog << std::endl;

  M.Tvmult(v,u);
  Assert (v == u, ExcInternalError());
  for (unsigned int i=0;i<v.size();++i)
    deallog << ' ' << v(i);
  deallog << std::endl;

  M.Tvmult_add(v,u);
  v /= 2;
  Assert (v == u, ExcInternalError());
  for (unsigned int i=0;i<v.size();++i)
    deallog << ' ' << v(i);
  deallog << std::endl;
}


int main()
{
  std::ofstream logfile("identity_matrix_01/output");
  logfile.setf(std::ios::fixed);
  logfile.precision(0);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  check_vmult<double>();
  check_vmult<float>();
}

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

// See documentation of Product for documentation of this example

#include <base/logstream.h>
#include <lac/matrix_lib.h>
#include <lac/full_matrix.h>
#include <lac/vector.h>
#include <lac/vector_memory.h>


double Adata[] =
{
      .5, .1,
      .4, .2
};

double Bdata[] =
{
      .866, .5,
      -.5, .866
};


int main()
{
  FullMatrix<float> A(2,2);
  FullMatrix<double> B(2,2);

  A.fill(Adata);
  B.fill(Bdata);
  
  GrowingVectorMemory<Vector<double> > mem;
  
  ProductMatrix<Vector<double> > AB(A,B,mem);

  Vector<double> u(2);
  Vector<double> v(2);

  u(0) = 1.;
  u(1) = 2.;

  AB.vmult(v,u);

  deallog << v(0) << '\t' << v(1) << std::endl;

  AB.Tvmult(v,u);

  deallog << v(0) << '\t' << v(1) << std::endl;
}

//----------------------------  tensor_complex.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tensor_complex.cc  ---------------------------

// Same as tensor.cc, but uses tensors based on std::complex<double> instead
// of double

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iomanip>
#include <complex>

int main ()
{
  std::ofstream logfile("tensor_complex/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  double a_double[3][3][2] = {{{1,-1}, {2,0}, {3,0}},
			      {{3,0}, {4,0}, {5,0}},
			      {{6,0}, {7,0}, {8,3}}};
  double b_double[3][3][2] = {{{24,-2}, {31,-2}, {37,6}},
			      {{45,-3}, {57,0},  {69,15}},
			      {{75,12}, {96,21}, {108,48}}};

  const unsigned int dim=3;
  std::complex<double> a[dim][dim], b[dim][dim];
  for (unsigned int d=0; d<dim; ++d)
    for (unsigned int e=0; e<dim; ++e)
      {
	a[d][e] = std::complex<double>(a_double[d][e][0],a_double[d][e][1]);
	b[d][e] = std::complex<double>(b_double[d][e][0],b_double[d][e][1]);
      }

  Tensor<2,dim,std::complex<double> > t(a);
  Tensor<2,dim,std::complex<double> > tt;
  Tensor<2,dim,std::complex<double> > result(b);
  Assert (transpose(transpose(t)) == t, ExcInternalError());
  Assert (transpose(transpose(result)) == result, ExcInternalError());

  Vector<std::complex<double> > unrolled(9);

				// cast result to double to profit from zero
				// threshold and so on
  t.unroll(unrolled);
  deallog << "unrolled:";
  for (unsigned i=0;i<9;i++)
    deallog << ' ' << unrolled(i);
  deallog << std::endl;

  deallog << "t=" << std::endl;
  for (unsigned int i=0; i<dim; ++i)
    {
      for (unsigned int j=0; j<dim; ++j)
	deallog << t[i][j] << ' ';
      deallog << std::endl;
    };

  deallog << "norm(t)=" << t.norm() << std::endl;

  contract (tt,t,t);

  deallog << "tt=" << std::endl;
  for (unsigned int i=0; i<dim; ++i)
    {
      for (unsigned int j=0; j<dim; ++j)
	deallog << tt[i][j] << ' ';
      deallog << std::endl;
    };

  if (true)
    {
      deallog.push("Cross product");
      Tensor<1,3,std::complex<double> > e1;
      Tensor<1,3,std::complex<double> > e2;
      Tensor<1,3,std::complex<double> > e3;
      e1[0] = 1.;
      e2[1] = 1.;
      e3[2] = 1.;
      Tensor<1,3,std::complex<double> > result;
      cross_product(result,e1,e2);
      deallog << '\t' << result[0]
	      << '\t' << result[1]
	      << '\t' << result[2] << std::endl;

      cross_product(result,e2,e3);
      deallog << '\t' << result[0]
	      << '\t' << result[1] << '\t'
	      << result[2] << std::endl;

      cross_product(result,e3,e1);
      deallog << '\t' << result[0]
	      << '\t' << result[1]
	      << '\t' << result[2] << std::endl;

      deallog.pop();
    }

  if (tt==result)
    {
      deallog << "Result OK." << std::endl;
    }
  else
    {
      deallog << "Result WRONG!" << std::endl;
    };

  t = 0;
  deallog << t << std::endl;
}

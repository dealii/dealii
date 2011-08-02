//----------------------------  tensor_float.cc  ---------------------------
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
//----------------------------  tensor_float.cc  ---------------------------

// Same as tensor.cc, but uses tensors based on floats instead of doubles

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iomanip>

int main ()
{
  std::ofstream logfile("tensor_float/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  float a[3][3] = {{1, 2, 3}, {3, 4, 5}, {6, 7, 8}};
  float b[3][3] = {{25,31,37}, {45,57,69}, {75,96,117}};

  const unsigned int dim=3;
  Tensor<2,dim,float> t(a);
  Tensor<2,dim,float> tt;
  Tensor<2,dim,float> result(b);
  Assert (transpose(transpose(t)) == t, ExcInternalError());
  Assert (transpose(transpose(result)) == result, ExcInternalError());

  Vector<float> unrolled(9);

				// cast result to double to profit from zero
				// threshold and so on
  t.unroll(unrolled);
  deallog << "unrolled:";
  for (unsigned i=0;i<9;i++)
    deallog << ' ' << static_cast<double>(unrolled(i));
  deallog << std::endl;

  deallog << "t=" << std::endl;
  for (unsigned int i=0; i<dim; ++i)
    {
      for (unsigned int j=0; j<dim; ++j)
	deallog << static_cast<double>(t[i][j]) << ' ';
      deallog << std::endl;
    };

  deallog << "norm(t)=" << t.norm() << std::endl;

  contract (tt,t,t);

  deallog << "tt=" << std::endl;
  for (unsigned int i=0; i<dim; ++i)
    {
      for (unsigned int j=0; j<dim; ++j)
	deallog << static_cast<double>(tt[i][j]) << ' ';
      deallog << std::endl;
    };

  if (true)
    {
      deallog.push("Cross product");
      Tensor<1,3,float> e1;
      Tensor<1,3,float> e2;
      Tensor<1,3,float> e3;
      e1[0] = 1.;
      e2[1] = 1.;
      e3[2] = 1.;
      Tensor<1,3,float> result;
      cross_product(result,e1,e2);
      deallog << '\t' << static_cast<double>(result[0])
	      << '\t' << static_cast<double>(result[1])
	      << '\t' << static_cast<double>(result[2]) << std::endl;

      cross_product(result,e2,e3);
      deallog << '\t' << static_cast<double>(result[0])
	      << '\t' << static_cast<double>(result[1]) << '\t'
	      << static_cast<double>(result[2]) << std::endl;

      cross_product(result,e3,e1);
      deallog << '\t' << static_cast<double>(result[0])
	      << '\t' << static_cast<double>(result[1])
	      << '\t' << static_cast<double>(result[2]) << std::endl;

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

}

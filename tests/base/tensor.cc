//----------------------------  tensor.cc  ---------------------------
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
//----------------------------  tensor.cc  ---------------------------


#include <base/tensor.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <fstream>


int main ()
{
  ofstream logfile("tensor.output");
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  double a[3][3] = {{1, 2, 3}, {3, 4, 5}, {6, 7, 8}};
  double b[3][3] = {{25,31,37}, {45,57,69}, {75,96,117}};
    
  const unsigned int dim=3;
  Tensor<2,dim> t(a);
  Tensor<2,dim> tt;
  Tensor<2,dim> result(b);

  Vector<double> unrolled(9);
  
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
  deallog << std::endl;
  
  contract (tt,t,t);

  deallog << "tt=" << std::endl;
  for (unsigned int i=0; i<dim; ++i)
    {
      for (unsigned int j=0; j<dim; ++j)
	deallog << tt[i][j] << ' ';
      deallog << std::endl;
    };
  deallog << std::endl;

  if (true)
    {
      deallog.push("Cross product");
      Tensor<1,3> e1;
      Tensor<1,3> e2;
      Tensor<1,3> e3;
      e1[0] = 1.;
      e2[1] = 1.;
      e3[2] = 1.;
      Tensor<1,3> result;
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

};

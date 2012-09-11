//----------------------------  tensor_05.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tensor_05.cc  ---------------------------

// check double_contract(Tensor<2,dim>,Tensor<2,dim>)

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iomanip>

template <int dim>
void test_constant()
{
  Tensor<2,dim> t;
  for (unsigned int i=0;i<dim;++i)
    for (unsigned int j=0;j<dim;++j)
      t[i][j] = 2.;
  deallog << "Constant dim " << dim << '\t' << double_contract(t,t)
	  << " compare " << 4*dim*dim << std::endl;
}


template <int dim>
void test_equal()
{
  Tensor<2,dim> t;
  unsigned int sum = 0;
  for (unsigned int i=0;i<dim;++i)
    for (unsigned int j=0;j<dim;++j)
      {
	t[i][j] = i+dim*j;
	sum += (i+dim*j)*(i+dim*j);
      }
  
  deallog << "Equal    dim " << dim << '\t' << double_contract(t,t)
	  << " compare " << sum << std::endl;
}


template <int dim>
void test_unequal()
{
  Tensor<2,dim> s;
  Tensor<2,dim> t;
  unsigned int sum = 0;
  for (unsigned int i=0;i<dim;++i)
    for (unsigned int j=0;j<dim;++j)
      {
	s[i][j] = i+dim*j;
	t[i][j] = dim*i+j;
	sum += (i+dim*j)*(dim*i+j);
      }
  
  deallog << "Unequal  dim " << dim << '\t' << double_contract(s,t)
	  << " compare " << sum << std::endl;
}


int main ()
{
  std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);

  test_constant<2>();
  test_constant<3>();
  test_constant<4>();

  test_equal<2>();
  test_equal<3>();
  test_equal<4>();

  test_unequal<2>();
  test_unequal<3>();
  test_unequal<4>();
}

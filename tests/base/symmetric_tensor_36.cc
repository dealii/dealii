// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// Check the tensor copy constructors for different number types

#include "../tests.h"
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>

// These boost functions are bundled with deal.II
#include <boost/type_traits.hpp>

#include <complex>
#include <fstream>
#include <iomanip>

// Some number types cannot be converted to floats/doubles;
// for example, complex numbers cannot be converted to primitive numbers
// So in these cases we disable the test if NumberType1 cannot be constructed
// from NumberType2
template<int rank, int dim,
         template<int,int,typename> class TensorType,
         typename NumberType1,
         typename NumberType2>
typename boost::disable_if<boost::is_constructible<NumberType1,NumberType2>,void>::type
test_tensor_constructor (const std::string &, const std::string &)
{}

template<int rank, int dim,
         template<int,int,typename> class TensorType,
         typename NumberType1,
         typename NumberType2>
typename boost::enable_if<boost::is_constructible<NumberType1,NumberType2>,void>::type
test_tensor_constructor (const std::string &type1, const std::string &type2)
{
  deallog
      << "Rank " << rank << ", "
      << "Dim " << dim << ":"
      << "  From " << type2
      << " To " << type1
      << " ... "
      << std::flush;
  TensorType<rank,dim,NumberType2> tmp2;
  TensorType<rank,dim,NumberType1> tmp1 (tmp2);
  deallog << "OK" << std::endl;
}

template<int rank, int dim,
         template<int,int,typename> class TensorType,
         typename NumberType1>
void
test_fixed_NT_2(const std::string &type1)
{
  test_tensor_constructor<rank,dim,TensorType,NumberType1,float>(type1, "float");
  test_tensor_constructor<rank,dim,TensorType,NumberType1,double>(type1, "double");
  test_tensor_constructor<rank,dim,TensorType,NumberType1,std::complex<float> >(type1, "std::complex<float>");
  test_tensor_constructor<rank,dim,TensorType,NumberType1,std::complex<double> >(type1, "std::complex<double>");
}

template<int rank, int dim,
         template<int,int,typename> class TensorType>
void
test_fixed_NT_1()
{
  test_fixed_NT_2<rank,dim,TensorType,float>("float");
  test_fixed_NT_2<rank,dim,TensorType,double>("double");
  test_fixed_NT_2<rank,dim,TensorType,std::complex<float> >("std::complex<float>");
  test_fixed_NT_2<rank,dim,TensorType,std::complex<double> >("std::complex<double>");
}

template<int rank, int dim>
void
test_fixed_TensorType()
{
  test_fixed_NT_1<rank,dim,Tensor>();
  test_fixed_NT_1<rank,dim,SymmetricTensor>();
}

template<int dim>
void
test_fixed_rank()
{
  test_fixed_TensorType<2,dim>();
  test_fixed_TensorType<4,dim>();
}

int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(5);
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  test_fixed_rank<2>();
  test_fixed_rank<3>();

  deallog << "All OK" << std::endl;
}

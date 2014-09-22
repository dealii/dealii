// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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


// Test ConstantTensorFunction

#include "../tests.h"
#include <deal.II/base/tensor_function.h>



// Fill a tensor with values:
template<int rank, int dim>
struct FillTensor
{
  static void fill_tensor(Tensor<rank, dim> &tensor, int base)
  {
    for (int i = 0; i < dim; ++i)
      FillTensor<rank-1,dim>::fill_tensor(tensor[i], 10 * base + i);
  }
};

template<int dim>
struct FillTensor<1, dim>
{
  static void fill_tensor(Tensor<1, dim> &tensor, int base)
  {
    for (int i = 0; i < dim; ++i)
      tensor[i] = 10 * base + i;
  }
};



// Print a tensor of arbitrary rank to the console:
template<int rank, int dim>
struct PrintTensor
{
  static void print_tensor(const Tensor<rank, dim> &tensor)
  {
    for (int i = 0; i < dim; ++i)
      {
        PrintTensor<rank-1,dim>::print_tensor(tensor[i]);
        deallog << " ";
      }
  }
};

template<int dim>
struct PrintTensor<1, dim>
{
  static void print_tensor(const Tensor<1, dim> &tensor)
  {
    for (int i = 0; i < dim; ++i)
      deallog << tensor[i] << " ";
  }
};



template <int rank, int dim>
void check ()
{
  deallog << "ConstantTensorFunction<"
          << rank << ", " << dim << ">:" << std::endl;

  Tensor<rank, dim> value;
  FillTensor<rank, dim>::fill_tensor(value, 0);

  ConstantTensorFunction<rank, dim> tensor_function(value);
  TensorFunction<rank, dim> *foo = &tensor_function;


  Point<dim> point;
  for (int i = 0; i < dim; ++i)
    point(i) = i;

  deallog << "->value:" << std::endl;
  PrintTensor<rank, dim>::print_tensor( foo->value(point) );
  deallog << std::endl;

  deallog << "->gradient:" << std::endl;
  PrintTensor<rank+1, dim>::print_tensor( foo->gradient(point) );
  deallog << std::endl;


  std::vector<Point<dim> > points;
  points.push_back(point);

  for (int i = 0; i < dim; ++i)
    point(i) = dim-i;
  points.push_back(point);

  std::vector<Tensor<rank, dim> > tensors;
  std::vector<Tensor<rank + 1, dim> > gradients;
  tensors.resize(2);
  gradients.resize(2);

  deallog << "->value_list:" << std::endl;
  foo->value_list(points, tensors);
  PrintTensor<rank, dim>::print_tensor(tensors[0]);
  deallog << std::endl;
  PrintTensor<rank, dim>::print_tensor(tensors[1]);
  deallog << std::endl;

  deallog << "->gradient_list:" << std::endl;
  foo->gradient_list(points, gradients);
  PrintTensor<rank + 1, dim>::print_tensor(gradients[0]);
  deallog << std::endl;
  PrintTensor<rank + 1, dim>::print_tensor(gradients[1]);
  deallog << std::endl;
}



int main()
{
  std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<1,1>();
  check<2,1>();
  check<3,1>();
  check<4,1>();

  check<1,2>();
  check<2,2>();
  check<3,2>();
  check<4,2>();

  check<1,3>();
  check<2,3>();
  check<3,3>();
  check<4,3>();
}


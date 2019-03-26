// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <iostream>

int main()
{
  cudaDeviceProp device_properties;
  const cudaError_t error = cudaGetDeviceProperties(&device_properties,
                                                    /*device*/0);
  if( error != cudaSuccess)
  {
    std::cout << "CUDA error: " << cudaGetErrorString(error) << '\n';
    return error;
  }
  std::cout << device_properties.major << device_properties.minor;
  return 0;
}

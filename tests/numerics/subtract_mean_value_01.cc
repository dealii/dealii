// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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



// check VectorTools::subtract_mean_value() for deal.II serial vectors

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <vector>

#include "../tests.h"

template <class VectorType>
void
test(VectorType &v)
{
  std::vector<bool> filter(v.size(), false);
  // set some elements of the vector
  for (unsigned int i = 0; i < v.size(); i += 1 + i)
    {
      filter[i] = true;
      v(i)      = i;
    }

  // then check the norm
  VectorTools::subtract_mean_value(v, filter);
  AssertThrow(std::fabs(v.mean_value()) < 1e-10 * v.l2_norm(),
              ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  try
    {
      {
        Vector<double> v(10);
        test(v);
      }

      {
        Vector<float> v(10);
        test(v);
      }

      {
        BlockVector<double> v(std::vector<types::global_dof_index>(1, 10));
        test(v);
      }

      {
        BlockVector<float> v(std::vector<types::global_dof_index>(1, 10));
        test(v);
      }

      {
        LinearAlgebra::Vector<double> v(10);
        test(v);
      }

      {
        LinearAlgebra::Vector<float> v(10);
        test(v);
      }
    }
  catch (std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}

// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2017 by the deal.II authors
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

// check is_serial_vector type trait for internal deal.II Vector types

#include "../tests.h"
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_vector.h>
#include <deal.II/lac/vector.h>

void
test()
{
  // make sure that is_serial_vector< dealii::Vector<Number> > is working
  Assert(is_serial_vector<dealii::Vector<double>>::value == true,
         ExcInternalError());

  deallog << is_serial_vector<dealii::Vector<float>>::value << std::endl;
  deallog << is_serial_vector<dealii::Vector<double>>::value << std::endl;
  deallog << is_serial_vector<dealii::Vector<long double>>::value << std::endl;
  deallog << is_serial_vector<dealii::Vector<int>>::value << std::endl;

  deallog << is_serial_vector<dealii::Vector<std::complex<float>>>::value
          << std::endl;
  deallog << is_serial_vector<dealii::Vector<std::complex<double>>>::value
          << std::endl;

  deallog << "OK" << std::endl << std::endl;

  // make sure that is_serial_vector< dealii::BlockVector<Number> > is working
  Assert(is_serial_vector<dealii::BlockVector<double>>::value == true,
         ExcInternalError());

  deallog << is_serial_vector<dealii::BlockVector<float>>::value << std::endl;
  deallog << is_serial_vector<dealii::BlockVector<double>>::value << std::endl;

  deallog << is_serial_vector<dealii::BlockVector<std::complex<float>>>::value
          << std::endl;
  deallog << is_serial_vector<dealii::BlockVector<std::complex<double>>>::value
          << std::endl;

  deallog << "OK" << std::endl << std::endl;

  // make sure that is_serial_vector< dealii::LinearAlgebra::Vector<Number> > is working
  Assert(is_serial_vector<dealii::LinearAlgebra::Vector<double>>::value == true,
         ExcInternalError());

  deallog << is_serial_vector<dealii::LinearAlgebra::Vector<float>>::value
          << std::endl;
  deallog << is_serial_vector<dealii::LinearAlgebra::Vector<double>>::value
          << std::endl;

  deallog << is_serial_vector<
               dealii::LinearAlgebra::Vector<std::complex<float>>>::value
          << std::endl;
  deallog << is_serial_vector<
               dealii::LinearAlgebra::Vector<std::complex<double>>>::value
          << std::endl;

  deallog << "OK" << std::endl << std::endl;

  // make sure that dealii::LinearAlgebra::distributed::Vector<Number> > is working
  Assert(
    is_serial_vector<dealii::LinearAlgebra::distributed::Vector<double>>::value
      == false,
    ExcInternalError());

  deallog << is_serial_vector<
               dealii::LinearAlgebra::distributed::Vector<float>>::value
          << std::endl;
  deallog << is_serial_vector<
               dealii::LinearAlgebra::distributed::Vector<double>>::value
          << std::endl;

  deallog
    << is_serial_vector<
         dealii::LinearAlgebra::distributed::Vector<std::complex<float>>>::value
    << std::endl;
  deallog << is_serial_vector<dealii::LinearAlgebra::distributed::Vector<
               std::complex<double>>>::value
          << std::endl;

  deallog << "OK" << std::endl << std::endl;

  // make sure that dealii::LinearAlgebra::distributed::BlockVector<Number> > is working
  Assert(is_serial_vector<
           dealii::LinearAlgebra::distributed::BlockVector<double>>::value
           == false,
         ExcInternalError());

  deallog << is_serial_vector<
               dealii::LinearAlgebra::distributed::BlockVector<float>>::value
          << std::endl;
  deallog << is_serial_vector<
               dealii::LinearAlgebra::distributed::BlockVector<double>>::value
          << std::endl;

  deallog << is_serial_vector<dealii::LinearAlgebra::distributed::BlockVector<
               std::complex<float>>>::value
          << std::endl;
  deallog << is_serial_vector<dealii::LinearAlgebra::distributed::BlockVector<
               std::complex<double>>>::value
          << std::endl;

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();

  try
    {
      test();
    }
  catch(std::exception& exc)
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
  catch(...)
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

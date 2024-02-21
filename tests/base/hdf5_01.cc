// ------------------------------------------------------------------------
//
// Copyright (C) 2019 - 2023 by the deal.II Authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Attributes tests with serial HDF5

#include <deal.II/base/hdf5.h>

#include "../tests.h"

// all include files you need here

void
test()
{
  std::string filename = "test.h5";

  // Write data
  {
    // Create file
    HDF5::File data_file(filename, HDF5::File::FileAccessMode::create);

    // Create attributes attached to the root
    const float        root_float        = 2.45681934e5;
    const double       root_double       = 6.234542e3;
    const int          root_int          = -56;
    const unsigned int root_unsigned_int = 22;
    data_file.set_attribute("root_float", root_float);
    data_file.set_attribute("root_double", root_double);
    data_file.set_attribute("root_int", root_int);
    data_file.set_attribute("root_unsigned_int", root_unsigned_int);
    data_file.set_attribute("root_total",
                            (root_float + root_double) * root_int *
                              root_unsigned_int);

    // Create attributes attached to a group
    auto test_group = data_file.create_group("test_group");

#ifdef DEAL_II_WITH_COMPLEX_VALUES
    const std::complex<float>  group_complex_float  = {2.45681934e5, 45e2};
    const std::complex<double> group_complex_double = {6.234542e3, 2};
    test_group.set_attribute("group_complex_float", group_complex_float);
    test_group.set_attribute("group_complex_double", group_complex_double);
    test_group.set_attribute("group_complex_total",
                             group_complex_float * group_complex_double);
#endif

    test_group.set_attribute("group_string",
                             std::string("test_string_attribute"));

    // Create attributes attached to a dataset
    const std::vector<hsize_t> dimensions = {50, 30};

    auto test_dataset =
      test_group.create_dataset<double>("test_dataset", dimensions);

    test_dataset.set_attribute("dataset_double", 20.2);
    test_dataset.set_attribute("dataset_string",
                               std::string("test_dataset_attribute"));
  }

  // Read data
  {
    // Read attributes attached to the root
    HDF5::File data_file(filename, HDF5::File::FileAccessMode::open);
    auto       root_float  = data_file.get_attribute<float>("root_float");
    auto       root_double = data_file.get_attribute<double>("root_double");
    auto       root_int    = data_file.get_attribute<int>("root_int");
    auto       root_unsigned_int =
      data_file.get_attribute<unsigned int>("root_unsigned_int");
    // calculated and read should be the same
    deallog << "root_total calculated:"
            << (root_float + root_double) * root_int * root_unsigned_int
            << std::endl;
    deallog << "root_total read:"
            << data_file.get_attribute<double>("root_total") << std::endl;

    // Read attributes attached to a group
    auto test_group = data_file.open_group("test_group");
    deallog << "group_name: " << test_group.get_name() << std::endl;

#ifdef DEAL_II_WITH_COMPLEX_VALUES
    auto group_complex_float =
      test_group.get_attribute<std::complex<float>>("group_complex_float");
    auto group_complex_double =
      test_group.get_attribute<std::complex<double>>("group_complex_double");
    deallog << "group_complex_total calculated:"
            << (group_complex_float * group_complex_double) << std::endl;
    deallog << "group_complex_total read:"
            << test_group.get_attribute<std::complex<double>>(
                 "group_complex_total")
            << std::endl;
#endif

    deallog << "group_string read:"
            << test_group.get_attribute<std::string>("group_string")
            << std::endl;

    // Read attributes attached to a dataset
    auto test_dataset = test_group.open_dataset("test_dataset");
    deallog << "dataset_double read:"
            << test_dataset.get_attribute<double>("dataset_double")
            << std::endl;
    deallog << "dataset_string read:"
            << test_dataset.get_attribute<std::string>("dataset_string")
            << std::endl;
  }
}



int
main(int argc, char **argv)
{
  initlog();

  try
    {
      test();
    }
  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}

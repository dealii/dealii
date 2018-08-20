// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II Authors
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

// Test the attributes of serial HDF5

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
    HDF5::File data_file(filename, HDF5::File::Mode::create);

    // Create attributes attached to the root
    const float        root_float        = 2.45681934e5;
    const double       root_double       = 6.234542e3;
    const int          root_int          = -56;
    const unsigned int root_unsigned_int = 22;
    data_file.write_attr("root_float", root_float);
    data_file.write_attr("root_double", root_double);
    data_file.write_attr("root_int", root_int);
    data_file.write_attr("root_unsigned_int", root_unsigned_int);
    data_file.write_attr("root_total",
                         (root_float + root_double) * root_int *
                           root_unsigned_int);

    // Create attributes attached to a group
    auto                      test_group = data_file.create_group("test_group");
    const std::complex<float> group_complex_float   = {2.45681934e5, 45e2};
    const std::complex<double> group_complex_double = {6.234542e3, 2};
    test_group.write_attr("group_complex_float", group_complex_float);
    test_group.write_attr("group_complex_double", group_complex_double);
    test_group.write_attr("group_complex_total",
                          group_complex_float * group_complex_double);
    test_group.write_attr("group_string", std::string("test_string_attribute"));

    // Create attributes attached to a dataset
    const std::vector<hsize_t> dimensions = {50, 30};

    auto test_dataset =
      test_group.create_dataset<double>("test_dataset", dimensions);

    test_dataset.write_attr("dataset_double", 20.2);
    test_dataset.write_attr("dataset_string",
                            std::string("test_dataset_attribute"));
  }

  // Read data
  {
    // Read attributes attached to the root
    HDF5::File data_file(filename, HDF5::File::Mode::open);
    auto       root_float  = data_file.attr<float>("root_float");
    auto       root_double = data_file.attr<double>("root_double");
    auto       root_int    = data_file.attr<int>("root_int");
    auto root_unsigned_int = data_file.attr<unsigned int>("root_unsigned_int");
    // calculated and read should be the same
    deallog << "root_total calculated:"
            << (root_float + root_double) * root_int * root_unsigned_int
            << std::endl;
    deallog << "root_total read:" << data_file.attr<double>("root_total")
            << std::endl;

    // Read attributes attached to a group
    auto test_group = data_file.group("test_group");
    auto group_complex_float =
      test_group.attr<std::complex<float>>("group_complex_float");
    auto group_complex_double =
      test_group.attr<std::complex<double>>("group_complex_double");
    deallog << "group_complex_total calculated:"
            << (group_complex_float * group_complex_double) << std::endl;
    deallog << "group_complex_total read:"
            << test_group.attr<std::complex<double>>("group_complex_total")
            << std::endl;
    deallog << "group_string read:"
            << test_group.attr<std::string>("group_string") << std::endl;

    // Read attributes attached to a dataset
    auto test_dataset = test_group.dataset("test_dataset");
    deallog << "dataset_double read:"
            << test_dataset.attr<double>("dataset_double") << std::endl;
    deallog << "dataset_string read:"
            << test_dataset.attr<std::string>("dataset_string") << std::endl;
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
  catch (std::exception &exc)
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

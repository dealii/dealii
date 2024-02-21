// ------------------------------------------------------------------------
//
// Copyright (C) 2021 - 2023 by the deal.II Authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test bool with parallel HDF5
// This test is based on tests/base/hdf5_02.cc

#include <deal.II/base/hdf5.h>

#include "../tests.h"

// all include files you need here

void
test()
{
  std::string filename = "test.h5";

  MPI_Comm mpi_communicator = MPI_COMM_WORLD;

  // Write data
  {
    // Create file
    HDF5::File data_file(filename,
                         HDF5::File::FileAccessMode::create,
                         mpi_communicator);

    // Create attributes attached to the root
    data_file.set_attribute("root_true", true);
    data_file.set_attribute("root_false", false);

    // Create attributes attached to a group
    auto group = data_file.create_group("test_group");

    group.set_attribute("group_true", true);
    group.set_attribute("group_false", false);



    std::string dataset_name("test_dataset");


    const std::vector<hsize_t> dataset_dimensions = {2, 5, 3};
    auto                       dataset =
      group.template create_dataset<double>(dataset_name, dataset_dimensions);


    dataset.set_attribute("dataset_true", true);
    dataset.set_attribute("dataset_false", false);
  }

  // Read data
  {
    // Read attributes attached to the root
    HDF5::File data_file(filename,
                         HDF5::File::FileAccessMode::open,
                         mpi_communicator);
    auto       root_true  = data_file.get_attribute<bool>("root_true");
    auto       root_false = data_file.get_attribute<bool>("root_false");
    deallog << "root_true:" << root_true << std::endl;
    deallog << "root_false:" << root_false << std::endl;

    // Read attributes attached to a group
    auto test_group = data_file.open_group("test_group");
    deallog << "group_name: " << test_group.get_name() << std::endl;
    auto group_true  = test_group.get_attribute<bool>("group_true");
    auto group_false = test_group.get_attribute<bool>("group_false");
    deallog << "group_true:" << group_true << std::endl;
    deallog << "group_false:" << group_false << std::endl;

    // Read attributes attached to a dataset
    auto test_dataset = test_group.open_dataset("test_dataset");
    deallog << "dataset_name: " << test_dataset.get_name() << std::endl;
    auto dataset_true  = test_dataset.get_attribute<bool>("dataset_true");
    auto dataset_false = test_dataset.get_attribute<bool>("dataset_false");
    deallog << "dataset_true:" << dataset_true << std::endl;
    deallog << "dataset_false:" << dataset_false << std::endl;
  }
}



int
main(int argc, char **argv)
{
  initlog();

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

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

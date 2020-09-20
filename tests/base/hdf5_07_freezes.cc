// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2019 by the deal.II Authors
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

// Test AssertThrow with hdf5

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/hdf5.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <boost/type_traits.hpp>

#include <numeric>

#include "../tests.h"


// This function tests parallel write and gets the group by reference
template <template <class...> class Container, typename Number>
void
write_test(HDF5::Group &              root_group,
           const std::vector<hsize_t> dataset_dimensions,
           MPI_Comm                   mpi_communicator,
           ConditionalOStream         pcout)
{
  AssertDimension(dataset_dimensions.size(), 2);

  std::string container_name;
  std::string type_name;

  if (std::is_same<Container<Number>, std::vector<Number>>::value)
    {
      container_name = std::string("std::vector");
    }
  else if (std::is_same<Container<Number>, FullMatrix<Number>>::value)
    {
      container_name = std::string("FullMatrix");
    }

  if (std::is_same<Number, double>::value)
    {
      type_name = std::string("double");
    }

  deallog << "Write tests for " << container_name << "<" << type_name << ">"
          << " datasets" << std::endl;

  auto group = root_group.create_group(container_name + "<" + type_name + ">");

  {
    std::string dataset_name("dataset_2");
    auto        dataset =
      group.create_dataset<Number>(dataset_name, dataset_dimensions);
    dataset.set_query_io_mode(true);
    deallog << "Dimensions " + dataset_name << " " << container_name << "<"
            << type_name << ">"
            << " (Write): " << dataset.get_dimensions() << std::endl;
    deallog << "Size " + dataset_name << " " << container_name << "<"
            << type_name << ">"
            << " (Write): " << dataset.get_size() << std::endl;
    deallog << "Rank " + dataset_name << " " << container_name << "<"
            << type_name << ">"
            << " (Write): " << dataset.get_rank() << std::endl;

    std::vector<hsize_t> coordinates_a = {0,
                                          0, // first point
                                          0,
                                          2, // second point
                                          3,
                                          4, // third point
                                          25,
                                          12}; // fourth point
    std::vector<Number>  data_a        = {2, 3, 5, 6};

    std::vector<hsize_t> coordinates_b = {5,
                                          0, // first point
                                          0,
                                          4, // second point
                                          5,
                                          4, // third point
                                          26,
                                          12}; // fourth point
    std::vector<Number>  data_b        = {9, 4, 7, 6};

    std::vector<hsize_t> coordinates_c = {5,
                                          7, // first point
                                          1,
                                          8, // second point
                                          5,
                                          5, // third point
                                          27,
                                          12, // fourth point
                                          28,
                                          29}; // fifth point
    std::vector<Number>  data_c        = {2, 3, 5, 6, 3};

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        dataset.write_selection(data_a, coordinates_a);
      }
    else if (Utilities::MPI::this_mpi_process(mpi_communicator) == 1)
      {
        dataset.write_selection(data_b, coordinates_b);
      }
    else if (Utilities::MPI::this_mpi_process(mpi_communicator) == 2)
      {
        dataset.write_selection(data_c, coordinates_c);
      }
    else
      {
        dataset.template write_none<Number>();
      }
    deallog << "Sum " + dataset_name << " " << container_name << "<"
            << type_name << ">"
            << " (Write): "
            << std::accumulate(data_a.begin(), data_a.end(), 0) +
                 std::accumulate(data_a.begin(), data_a.end(), 0) +
                 std::accumulate(data_a.begin(), data_a.end(), 0)
            << std::endl;
  }
}



int
main(int argc, char **argv)
{
  initlog();

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      MPI_Comm mpi_communicator = MPI_COMM_WORLD;

      ConditionalOStream pcout(
        std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0);

      std::string                filename           = "test.h5";
      const std::vector<hsize_t> dataset_dimensions = {50, 30};

      HDF5::File data_file(filename,
                           HDF5::File::FileAccessMode::create,
                           mpi_communicator);

      write_test<std::vector, double>(data_file,
                                      dataset_dimensions,
                                      mpi_communicator,
                                      pcout);

      AssertThrow(Utilities::MPI::this_mpi_process(mpi_communicator) == 0,
                  ExcInternalError());

      write_test<FullMatrix, double>(data_file,
                                     dataset_dimensions,
                                     mpi_communicator,
                                     pcout);
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

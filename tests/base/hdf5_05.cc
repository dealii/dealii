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

// Read/write tests with parallel HDF5. This tests different hyperslab shapes
// with different ranks.

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/hdf5.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <boost/lexical_cast.hpp>
#include <boost/type_traits.hpp>

#include <numeric>

#include "../tests.h"

// This function initializes a container of Number type
template <template <class...> class Container, typename Number>
typename std::enable_if<
  std::is_same<Container<Number>, std::vector<Number>>::value,
  Container<Number>>::type
initialize_container(std::vector<hsize_t> dimensions)
{
  return Container<Number>(std::accumulate(
    dimensions.begin(), dimensions.end(), 1, std::multiplies<int>()));
}



template <template <class...> class Container, typename Number>
typename std::enable_if<std::is_same<Container<Number>, Vector<Number>>::value,
                        Container<Number>>::type
initialize_container(std::vector<hsize_t> dimensions)
{
  return Container<Number>(std::accumulate(
    dimensions.begin(), dimensions.end(), 1, std::multiplies<int>()));
}



template <template <class...> class Container, typename Number>
typename std::enable_if<
  std::is_same<Container<Number>, FullMatrix<Number>>::value,
  Container<Number>>::type
initialize_container(std::vector<hsize_t> dimensions)
{
  return FullMatrix<Number>(dimensions[0], dimensions[1]);
}



// This function assignes data to the elements of the container
template <template <class...> class Container, typename Number>
typename std::enable_if<
  std::is_same<Container<Number>, std::vector<Number>>::value,
  void>::type
assign_data(Container<Number> &data)
{
  for (unsigned int idx = 0; idx < data.size(); ++idx)
    {
      data[idx] = idx;
    }
}



template <template <class...> class Container, typename Number>
typename std::enable_if<std::is_same<Container<Number>, Vector<Number>>::value,
                        void>::type
assign_data(Container<Number> &data)
{
  for (unsigned int idx = 0; idx < data.size(); ++idx)
    {
      data[idx] = idx;
    }
}



template <template <class...> class Container, typename Number>
typename std::enable_if<
  std::is_same<Container<Number>, FullMatrix<Number>>::value,
  void>::type
assign_data(Container<Number> &data)
{
  for (unsigned int row_idx = 0; row_idx < data.m(); ++row_idx)
    {
      for (unsigned int col_idx = 0; col_idx < data.n(); ++col_idx)
        {
          data[row_idx][col_idx] = (row_idx * data.n() + col_idx);
        }
    }
}


// This function converts the data to a string
template <typename Number>
std::string
container_to_string(std::vector<Number> &data)
{
  std::string data_string;

  auto append_to_data_string = [&data_string](const std::string p) {
    if (data_string.length() > 0)
      data_string += ", ";
    data_string += p;
  };

  for (unsigned int idx = 0; idx < data.size(); ++idx)
    {
      append_to_data_string(boost::lexical_cast<std::string>(data[idx]));
    }

  return "[" + data_string + "]";
}

template <typename Number>
std::string
container_to_string(FullMatrix<Number> &data)
{
  std::string data_string;

  for (unsigned int row_idx = 0; row_idx < data.m(); ++row_idx)
    {
      if (row_idx > 0)
        data_string += ", ";
      data_string += "[";
      for (unsigned int col_idx = 0; col_idx < data.n(); ++col_idx)
        {
          if (col_idx > 0)
            data_string += ", ";
          data_string +=
            boost::lexical_cast<std::string>(data[row_idx][col_idx]);
        }
      data_string += "]";
    }

  return "[" + data_string + "]";
}

// This constexpr function converts a type to a std::string
template <typename Number>
std::string
type_to_string()
{
  std::string type_name;

  if (std::is_same<Number, float>::value)
    {
      type_name = std::string("float");
    }
  else if (std::is_same<Number, double>::value)
    {
      type_name = std::string("double");
    }
  else if (std::is_same<Number, std::complex<float>>::value)
    {
      type_name = std::string("std::complex<float>");
    }
  else if (std::is_same<Number, std::complex<double>>::value)
    {
      type_name = std::string("std::complex<double>");
    }
  else if (std::is_same<Number, int>::value)
    {
      type_name = std::string("int");
    }
  else if (std::is_same<Number, unsigned int>::value)
    {
      type_name = std::string("unsigned int");
    }

  return type_name;
}



// This function tests parallel write and gets the group by reference
template <typename Number>
void
write_test(HDF5::Group &      root_group,
           MPI_Comm           mpi_communicator,
           ConditionalOStream pcout)
{
  const std::string type_name = type_to_string<Number>();

  deallog << "Write tests for " << type_name << " datasets" << std::endl;

  auto group = root_group.create_group(type_name);


  {
    std::string dataset_name("dataset_1");

    // Write matrix of 5x3
    // Output generated by h5dump
    //  (0,0): 0, 1, 2,
    //  (1,0): 3, 4, 5,
    //  (2,0): 6, 7, 8,
    //  (3,0): 9, 10, 11,
    //  (4,0): 12, 13, 14

    const std::vector<hsize_t> dataset_dimensions = {5, 3};
    auto                       dataset =
      group.create_dataset<Number>(dataset_name, dataset_dimensions);
    dataset.set_query_io_mode(true);
    deallog << "Dimensions " << dataset_name << "<" << type_name << ">"
            << " (Write): " << dataset.get_dimensions() << std::endl;
    deallog << "Size " << dataset_name << "<" << type_name << ">"
            << " (Write): " << dataset.get_size() << std::endl;
    deallog << "Rank " << dataset_name << "<" << type_name << ">"
            << " (Write): " << dataset.get_rank() << std::endl;
    auto data = initialize_container<std::vector, Number>(dataset_dimensions);
    assign_data(data);
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        dataset.write(data);
      }
    else
      {
        dataset.template write_none<Number>();
      }
    deallog << "Data " + dataset_name << "<" << type_name << ">"
            << " (Write): " << container_to_string(data) << std::endl;
    pcout << "IO mode " + dataset_name << "<" << type_name << ">"
          << " (Write): " << dataset.get_io_mode() << std::endl;
    pcout << "Local no collective cause " + dataset_name << "<" << type_name
          << ">"
          << " (Write): " << dataset.get_local_no_collective_cause()
          << std::endl;
    pcout << "Global no collective cause " + dataset_name << "<" << type_name
          << ">"
          << " (Write): " << dataset.get_global_no_collective_cause()
          << std::endl;
  }
}


// This function tests parallel read. Unlike write_test, this functions gets a
// copy of the group
template <typename Number>
void
read_test(HDF5::Group        root_group,
          MPI_Comm           mpi_communicator,
          ConditionalOStream pcout)
{
  const std::string type_name = type_to_string<Number>();

  deallog << "Read tests for " << type_name << " datasets" << std::endl;

  auto group = root_group.group(type_name);

  {
    std::string dataset_name("dataset_1");
    auto        dataset = group.dataset(dataset_name);
    dataset.set_query_io_mode(true);
    deallog << "Dimensions " + dataset_name << "<" << type_name << ">"
            << " (Read): " << dataset.get_dimensions() << std::endl;
    deallog << "Size " + dataset_name << "<" << type_name << ">"
            << " (Read): " << dataset.get_size() << std::endl;
    deallog << "Rank " + dataset_name << "<" << type_name << ">"
            << " (Read): " << dataset.get_rank() << std::endl;
    {
      auto data = dataset.read<FullMatrix<Number>>();

      pcout << "IO mode " + dataset_name << "<" << type_name << ">"
            << " (Read): " << dataset.get_io_mode() << std::endl;
      pcout << "Local no collective cause " << dataset_name << "<" << type_name
            << ">"
            << " (Read): " << dataset.get_local_no_collective_cause()
            << std::endl;
      pcout << "Global no collective cause " << dataset_name << "<" << type_name
            << ">"
            << " (Read): " << dataset.get_global_no_collective_cause()
            << std::endl;

      deallog << "Data " + dataset_name << "<" << type_name << ">"
              << " (Read): " << container_to_string(data) << std::endl;
    }

    {
      // The output should be the column vector [7, 10, 13]
      const std::vector<hsize_t> hyperslab_offset = {2, 1};
      const std::vector<hsize_t> hyperslab_count  = {3, 1};
      auto data = dataset.read_hyperslab<std::vector<Number>>(hyperslab_offset,
                                                              hyperslab_count);

      pcout << "Column vector" << std::endl;
      pcout << "IO mode " + dataset_name << "<" << type_name << ">"
            << " (Read): " << dataset.get_io_mode() << std::endl;
      pcout << "Local no collective cause " << dataset_name << "<" << type_name
            << ">"
            << " (Read): " << dataset.get_local_no_collective_cause()
            << std::endl;
      pcout << "Global no collective cause " << dataset_name << "<" << type_name
            << ">"
            << " (Read): " << dataset.get_global_no_collective_cause()
            << std::endl;

      deallog << "Column vector " + dataset_name << "<" << type_name << ">"
              << " (Read): " << container_to_string(data) << std::endl;
    }

    {
      // The output should be the row vector [9, 10, 11]
      const std::vector<hsize_t> hyperslab_offset = {3, 0};
      const std::vector<hsize_t> hyperslab_count  = {1, 3};
      auto data = dataset.read_hyperslab<std::vector<Number>>(hyperslab_offset,
                                                              hyperslab_count);

      pcout << "Row vector" << std::endl;
      pcout << "IO mode " + dataset_name << "<" << type_name << ">"
            << " (Read): " << dataset.get_io_mode() << std::endl;
      pcout << "Local no collective cause " << dataset_name << "<" << type_name
            << ">"
            << " (Read): " << dataset.get_local_no_collective_cause()
            << std::endl;
      pcout << "Global no collective cause " << dataset_name << "<" << type_name
            << ">"
            << " (Read): " << dataset.get_global_no_collective_cause()
            << std::endl;

      deallog << "Row vector " + dataset_name << "<" << type_name << ">"
              << " (Read): " << container_to_string(data) << std::endl;
    }

    {
      // The output should be the matrix [[3, 4, 5], [6, 7, 8]]
      const std::vector<hsize_t> hyperslab_offset = {1, 0};
      const std::vector<hsize_t> hyperslab_count  = {2, 3};
      auto data = dataset.read_hyperslab<FullMatrix<Number>>(hyperslab_offset,
                                                             hyperslab_count);

      pcout << "Sub-matrix" << std::endl;
      pcout << "IO mode " + dataset_name << "<" << type_name << ">"
            << " (Read): " << dataset.get_io_mode() << std::endl;
      pcout << "Local no collective cause " << dataset_name << "<" << type_name
            << ">"
            << " (Read): " << dataset.get_local_no_collective_cause()
            << std::endl;
      pcout << "Global no collective cause " << dataset_name << "<" << type_name
            << ">"
            << " (Read): " << dataset.get_global_no_collective_cause()
            << std::endl;

      deallog << "Sub-matrix " + dataset_name << "<" << type_name << ">"
              << " (Read): " << container_to_string(data) << std::endl;
    }
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

      std::string filename = "test.h5";

      {
        HDF5::File data_file(filename,
                             mpi_communicator,
                             HDF5::File::Mode::create);

        write_test<double>(data_file, mpi_communicator, pcout);
        write_test<std::complex<double>>(data_file, mpi_communicator, pcout);
      }

      {
        HDF5::File data_file(filename,
                             mpi_communicator,
                             HDF5::File::Mode::open);
        read_test<double>(data_file, mpi_communicator, pcout);
        read_test<std::complex<double>>(data_file, mpi_communicator, pcout);
      }
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

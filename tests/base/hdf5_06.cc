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

// Read/write tests with serial HDF5. This tests different hyperslab shapes with
// different ranks.

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
std::enable_if_t<std::is_same_v<Container<Number>, std::vector<Number>>,
                 Container<Number>>
initialize_container(std::vector<hsize_t> dimensions)
{
  return Container<Number>(std::accumulate(
    dimensions.begin(), dimensions.end(), 1, std::multiplies<int>()));
}



template <template <class...> class Container, typename Number>
std::enable_if_t<std::is_same_v<Container<Number>, Vector<Number>>,
                 Container<Number>>
initialize_container(std::vector<hsize_t> dimensions)
{
  return Container<Number>(std::accumulate(
    dimensions.begin(), dimensions.end(), 1, std::multiplies<int>()));
}



template <template <class...> class Container, typename Number>
std::enable_if_t<std::is_same_v<Container<Number>, FullMatrix<Number>>,
                 Container<Number>>
initialize_container(std::vector<hsize_t> dimensions)
{
  return FullMatrix<Number>(dimensions[0], dimensions[1]);
}



// This function assigns data to the elements of the container
template <template <class...> class Container, typename Number>
std::enable_if_t<std::is_same_v<Container<Number>, std::vector<Number>>, void>
assign_data(Container<Number> &data)
{
  for (unsigned int idx = 0; idx < data.size(); ++idx)
    {
      data[idx] = idx;
    }
}



template <template <class...> class Container, typename Number>
std::enable_if_t<std::is_same_v<Container<Number>, Vector<Number>>, void>
assign_data(Container<Number> &data)
{
  for (unsigned int idx = 0; idx < data.size(); ++idx)
    {
      data[idx] = idx;
    }
}



template <template <class...> class Container, typename Number>
std::enable_if_t<std::is_same_v<Container<Number>, FullMatrix<Number>>, void>
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
container_to_string(Vector<Number> &data)
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

  if (std::is_same_v<Number, float>)
    {
      type_name = std::string("float");
    }
  else if (std::is_same_v<Number, double>)
    {
      type_name = std::string("double");
    }
  else if (std::is_same_v<Number, std::complex<float>>)
    {
      type_name = std::string("std::complex<float>");
    }
  else if (std::is_same_v<Number, std::complex<double>>)
    {
      type_name = std::string("std::complex<double>");
    }
  else if (std::is_same_v<Number, int>)
    {
      type_name = std::string("int");
    }
  else if (std::is_same_v<Number, unsigned int>)
    {
      type_name = std::string("unsigned int");
    }

  return type_name;
}



// This function tests parallel write and gets the group by reference
template <typename Number>
void
write_test(HDF5::Group &root_group)
{
  const std::string type_name = type_to_string<Number>();

  deallog << "Write tests for " << type_name << " datasets" << std::endl;

  auto group = root_group.create_group(type_name);


  {
    std::string dataset_name("dataset_1");

    // This dataset is a 5x3 matrix
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
    deallog << "Dimensions " << dataset_name << '<' << type_name << '>'
            << " (Write): " << dataset.get_dimensions() << std::endl;
    deallog << "Size " << dataset_name << '<' << type_name << '>'
            << " (Write): " << dataset.get_size() << std::endl;
    deallog << "Rank " << dataset_name << '<' << type_name << '>'
            << " (Write): " << dataset.get_rank() << std::endl;
    auto data = initialize_container<std::vector, Number>(dataset_dimensions);
    assign_data(data);

    dataset.write(data);

    deallog << "Data " + dataset_name << '<' << type_name << '>'
            << " (Write): " << container_to_string(data) << std::endl;
    std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
              << " (Write): " << dataset.get_io_mode() << std::endl;
    std::cout << "Local no collective cause " + dataset_name << '<' << type_name
              << '>' << " (Write): " << dataset.get_local_no_collective_cause()
              << std::endl;
    std::cout << "Global no collective cause " + dataset_name << '<'
              << type_name << '>'
              << " (Write): " << dataset.get_global_no_collective_cause()
              << std::endl;
  }

  {
    std::string dataset_name("dataset_2");

    // This dataset is a 2x5x3 tensor
    // Output generated by h5dump
    //  (0,0,0): 0, 1, 2,         (1,0,0): 15, 16, 17,
    //  (0,1,0): 3, 4, 5,         (1,1,0): 18, 19, 20,
    //  (0,2,0): 6, 7, 8,         (1,2,0): 21, 22, 23,
    //  (0,3,0): 9, 10, 11,       (1,3,0): 24, 25, 26,
    //  (0,4,0): 12, 13, 14,      (1,4,0): 27, 28, 29

    const std::vector<hsize_t> dataset_dimensions = {2, 5, 3};
    auto                       dataset =
      group.create_dataset<Number>(dataset_name, dataset_dimensions);
    dataset.set_query_io_mode(true);
    deallog << "Dimensions " << dataset_name << '<' << type_name << '>'
            << " (Write): " << dataset.get_dimensions() << std::endl;
    deallog << "Size " << dataset_name << '<' << type_name << '>'
            << " (Write): " << dataset.get_size() << std::endl;
    deallog << "Rank " << dataset_name << '<' << type_name << '>'
            << " (Write): " << dataset.get_rank() << std::endl;
    auto data = initialize_container<std::vector, Number>(dataset_dimensions);
    assign_data(data);

    dataset.write(data);

    deallog << "Data " + dataset_name << '<' << type_name << '>'
            << " (Write): " << container_to_string(data) << std::endl;
    std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
              << " (Write): " << dataset.get_io_mode() << std::endl;
    std::cout << "Local no collective cause " + dataset_name << '<' << type_name
              << '>' << " (Write): " << dataset.get_local_no_collective_cause()
              << std::endl;
    std::cout << "Global no collective cause " + dataset_name << '<'
              << type_name << '>'
              << " (Write): " << dataset.get_global_no_collective_cause()
              << std::endl;
  }

  {
    std::string dataset_name("dataset_3");

    // This dataset is a 2x5x3 tensor
    // Output generated by h5dump
    //  (0,0,0): 0, 1, 2,         (1,0,0): 0, 2, 4,
    //  (0,1,0): 3, 4, 5,         (1,1,0): 6, 8, 10,
    //  (0,2,0): 6, 7, 8,         (1,2,0): 12, 14, 16,
    //  (0,3,0): 9, 10, 11,       (1,3,0): 18, 20, 22,
    //  (0,4,0): 12, 13, 14,      (1,4,0): 24, 26, 28

    const std::vector<hsize_t> dataset_dimensions = {2, 5, 3};
    auto                       dataset =
      group.create_dataset<Number>(dataset_name, dataset_dimensions);
    dataset.set_query_io_mode(true);
    deallog << "Dimensions " << dataset_name << '<' << type_name << '>'
            << " (Write): " << dataset.get_dimensions() << std::endl;
    deallog << "Size " << dataset_name << '<' << type_name << '>'
            << " (Write): " << dataset.get_size() << std::endl;
    deallog << "Rank " << dataset_name << '<' << type_name << '>'
            << " (Write): " << dataset.get_rank() << std::endl;

    {
      const std::vector<hsize_t> data_dimensions = {5, 3};
      auto data = initialize_container<std::vector, Number>(data_dimensions);
      assign_data(data);
      const std::vector<hsize_t> hyperslab_offset = {0, 0, 0};
      const std::vector<hsize_t> hyperslab_count  = {1, 5, 3};
      dataset.write_hyperslab(data, hyperslab_offset, hyperslab_count);

      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Write): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_global_no_collective_cause()
                << std::endl;
    }

    {
      const std::vector<hsize_t> data_dimensions = {5, 3};
      auto data = initialize_container<FullMatrix, Number>(data_dimensions);
      assign_data(data);
      // Modify the data in order to have different data than in mpi_comm = 0
      data *= 2;
      const std::vector<hsize_t> hyperslab_offset = {1, 0, 0};
      const std::vector<hsize_t> hyperslab_count  = {1, 5, 3};
      dataset.write_hyperslab(data, hyperslab_offset, hyperslab_count);

      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Write): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_global_no_collective_cause()
                << std::endl;
    }
  }

  {
    std::string dataset_name("dataset_4");

    // This dataset is a 3x2x4 tensor
    // Output generated by h5dump
    //  (0,0,0): 0, 1, 2, 3,
    //  (0,1,0): 0, 2, 4, 6,
    //  (1,0,0): 4, 5, 6, 7,
    //  (1,1,0): 8, 10, 12, 14,
    //  (2,0,0): 8, 9, 10, 11,
    //  (2,1,0): 16, 18, 20, 22

    const std::vector<hsize_t> dataset_dimensions = {3, 2, 4};
    auto                       dataset =
      group.create_dataset<Number>(dataset_name, dataset_dimensions);
    dataset.set_query_io_mode(true);
    deallog << "Dimensions " << dataset_name << '<' << type_name << '>'
            << " (Write): " << dataset.get_dimensions() << std::endl;
    deallog << "Size " << dataset_name << '<' << type_name << '>'
            << " (Write): " << dataset.get_size() << std::endl;
    deallog << "Rank " << dataset_name << '<' << type_name << '>'
            << " (Write): " << dataset.get_rank() << std::endl;

    {
      const std::vector<hsize_t> data_dimensions = {3, 4};
      auto data = initialize_container<std::vector, Number>(data_dimensions);
      assign_data(data);
      const std::vector<hsize_t> hyperslab_offset = {0, 0, 0};
      const std::vector<hsize_t> hyperslab_count  = {3, 1, 4};
      dataset.write_hyperslab(data, hyperslab_offset, hyperslab_count);

      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Write): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_global_no_collective_cause()
                << std::endl;
    }

    {
      const std::vector<hsize_t> data_dimensions = {3, 4};
      auto data = initialize_container<FullMatrix, Number>(data_dimensions);
      assign_data(data);
      // Modify the data in order to have different data than in mpi_comm = 0
      data *= 2;
      const std::vector<hsize_t> hyperslab_offset = {0, 1, 0};
      const std::vector<hsize_t> hyperslab_count  = {3, 1, 4};
      dataset.write_hyperslab(data, hyperslab_offset, hyperslab_count);

      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Write): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_global_no_collective_cause()
                << std::endl;
    }
  }

  {
    std::string dataset_name("dataset_5");

    // This dataset is a 3x2x4 tensor
    // Output generated by h5dump
    //  (0,0,0): 2, 0, 0, 38,
    //  (0,1,0): 6, 5, 0, 3,
    //  (1,0,0): 12, 0, 0, 0,
    //  (1,1,0): 16, 15, 0, 13,
    //  (2,0,0): 32, 0, 0, 0,
    //  (2,1,0): 36, 35, 0, 33


    const std::vector<hsize_t> dataset_dimensions = {3, 2, 4};
    auto                       dataset =
      group.create_dataset<Number>(dataset_name, dataset_dimensions);
    dataset.set_query_io_mode(true);
    deallog << "Dimensions " << dataset_name << '<' << type_name << '>'
            << " (Write): " << dataset.get_dimensions() << std::endl;
    deallog << "Size " << dataset_name << '<' << type_name << '>'
            << " (Write): " << dataset.get_size() << std::endl;
    deallog << "Rank " << dataset_name << '<' << type_name << '>'
            << " (Write): " << dataset.get_rank() << std::endl;

    {
      std::vector<hsize_t> coordinates = {0,
                                          0,
                                          0, // first point
                                          0,
                                          1,
                                          3, // second point
                                          0,
                                          1,
                                          1, // third point
                                          0,
                                          1,
                                          0}; // fourth point
      std::vector<Number>  data        = {2, 3, 5, 6};

      dataset.write_selection(data, coordinates);

      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Write): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_global_no_collective_cause()
                << std::endl;
    }

    {
      std::vector<hsize_t> coordinates = {1,
                                          0,
                                          0, // first point
                                          1,
                                          1,
                                          3, // second point
                                          1,
                                          1,
                                          1, // third point
                                          1,
                                          1,
                                          0}; // fourth point
      std::vector<Number>  data        = {12, 13, 15, 16};

      dataset.write_selection(data, coordinates);

      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Write): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_global_no_collective_cause()
                << std::endl;
    }

    {
      std::vector<hsize_t> coordinates = {2,
                                          0,
                                          0, // first point
                                          2,
                                          1,
                                          3, // second point
                                          2,
                                          1,
                                          1, // third point
                                          2,
                                          1,
                                          0, // fourth point
                                          0,
                                          0,
                                          3}; // fifth point
      std::vector<Number>  data        = {32, 33, 35, 36, 38};

      dataset.write_selection(data, coordinates);

      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Write): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_global_no_collective_cause()
                << std::endl;
    }
  }

  {
    std::string dataset_name("dataset_6");

    // This dataset is a 7x12 matrix
    // Output generated by h5dump
    //  (0,0): -0, 0, 1, -1, 2, 3, -2, 4, 5, -3, 6, 7,
    //  (1,0): -4, 8, 9, -5, 10, 11, -6, 12, 13, -7, 14, 15,
    //  (2,0): -8, 16, 17, -9, 18, 19, -10, 20, 21, -11, 22, 23,
    //  (3,0): -12, 0, 0, -13, 0, 0, -14, 0, 0, -15, 0, 0,
    //  (4,0): -16, 24, 25, -17, 26, 27, -18, 28, 29, -19, 30, 31,
    //  (5,0): -20, 32, 33, -21, 34, 35, -22, 36, 37, -23, 38, 39,
    //  (6,0): -24, 40, 41, -25, 42, 43, -26, 44, 45, -27, 46, 47

    const std::vector<hsize_t> dataset_dimensions = {7, 12};
    auto                       dataset =
      group.create_dataset<Number>(dataset_name, dataset_dimensions);
    dataset.set_query_io_mode(true);
    deallog << "Dimensions " << dataset_name << '<' << type_name << '>'
            << " (Write): " << dataset.get_dimensions() << std::endl;
    deallog << "Size " << dataset_name << '<' << type_name << '>'
            << " (Write): " << dataset.get_size() << std::endl;
    deallog << "Rank " << dataset_name << '<' << type_name << '>'
            << " (Write): " << dataset.get_rank() << std::endl;

    {
      const std::vector<hsize_t> data_dimensions = {6, 8};
      auto data = initialize_container<FullMatrix, Number>(data_dimensions);
      assign_data(data);
      const std::vector<hsize_t> offset = {0, 1};
      const std::vector<hsize_t> stride = {4, 3};
      const std::vector<hsize_t> count  = {2, 4};
      const std::vector<hsize_t> block  = {3, 2};

      dataset.write_hyperslab(
        data, data_dimensions, offset, stride, count, block);

      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Write): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_global_no_collective_cause()
                << std::endl;
    }

    {
      const std::vector<hsize_t> data_dimensions = {7, 4};
      auto data = initialize_container<FullMatrix, Number>(data_dimensions);
      assign_data(data);
      data *= -1;
      const std::vector<hsize_t> offset = {0, 0};
      const std::vector<hsize_t> stride = {7, 3};
      const std::vector<hsize_t> count  = {1, 4};
      const std::vector<hsize_t> block  = {7, 1};

      dataset.write_hyperslab(
        data, data_dimensions, offset, stride, count, block);

      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Write): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_global_no_collective_cause()
                << std::endl;
    }
  }

  {
    std::string dataset_name("dataset_7");

    // This dataset is a 3x2x4 tensor
    // Output generated by h5dump
    //  (0,0,0): 2, 0, 0, 38,
    //  (0,1,0): 6, 5, 0, 3,
    //  (1,0,0): 12, 0, 0, 0,
    //  (1,1,0): 16, 15, 0, 13,
    //  (2,0,0): 32, 0, 0, 0,
    //  (2,1,0): 36, 35, 0, 33


    const std::vector<hsize_t> dataset_dimensions = {3, 2, 4};
    auto                       dataset =
      group.create_dataset<Number>(dataset_name, dataset_dimensions);
    dataset.set_query_io_mode(true);
    deallog << "Dimensions " << dataset_name << '<' << type_name << '>'
            << " (Write): " << dataset.get_dimensions() << std::endl;
    deallog << "Size " << dataset_name << '<' << type_name << '>'
            << " (Write): " << dataset.get_size() << std::endl;
    deallog << "Rank " << dataset_name << '<' << type_name << '>'
            << " (Write): " << dataset.get_rank() << std::endl;

    {
      std::vector<hsize_t> coordinates = {0,
                                          0,
                                          0, // first point
                                          0,
                                          1,
                                          3, // second point
                                          0,
                                          1,
                                          1, // third point
                                          0,
                                          1,
                                          0}; // fourth point
      Vector<Number>       data(4);
      data[0] = 2;
      data[1] = 3;
      data[2] = 5;
      data[3] = 6;

      dataset.write_selection(data, coordinates);

      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Write): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_global_no_collective_cause()
                << std::endl;
    }

    {
      std::vector<hsize_t> coordinates = {1,
                                          0,
                                          0, // first point
                                          1,
                                          1,
                                          3, // second point
                                          1,
                                          1,
                                          1, // third point
                                          1,
                                          1,
                                          0}; // fourth point
      Vector<Number>       data(4);
      data[0] = 12;
      data[1] = 13;
      data[2] = 15;
      data[3] = 16;

      dataset.write_selection(data, coordinates);

      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Write): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_global_no_collective_cause()
                << std::endl;
    }

    {
      std::vector<hsize_t> coordinates = {2,
                                          0,
                                          0, // first point
                                          2,
                                          1,
                                          3, // second point
                                          2,
                                          1,
                                          1, // third point
                                          2,
                                          1,
                                          0, // fourth point
                                          0,
                                          0,
                                          3}; // fifth point
      Vector<Number>       data(5);
      data[0] = 32;
      data[1] = 33;
      data[2] = 35;
      data[3] = 36;
      data[4] = 38;

      dataset.write_selection(data, coordinates);

      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Write): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " + dataset_name << '<'
                << type_name << '>'
                << " (Write): " << dataset.get_global_no_collective_cause()
                << std::endl;
    }
  }
}


// This function tests parallel read. Unlike write_test, this functions gets a
// copy of the group
template <typename Number>
void
read_test(HDF5::Group root_group)
{
  const std::string type_name = type_to_string<Number>();

  deallog << "Read tests for " << type_name << " datasets" << std::endl;

  auto group = root_group.open_group(type_name);

  {
    std::string dataset_name("dataset_1");

    // This dataset is a 5x3 matrix
    // Output generated by h5dump
    //  (0,0): 0, 1, 2,
    //  (1,0): 3, 4, 5,
    //  (2,0): 6, 7, 8,
    //  (3,0): 9, 10, 11,
    //  (4,0): 12, 13, 14

    auto dataset = group.open_dataset(dataset_name);
    dataset.set_query_io_mode(true);
    deallog << "Dimensions " + dataset_name << '<' << type_name << '>'
            << " (Read): " << dataset.get_dimensions() << std::endl;
    deallog << "Size " + dataset_name << '<' << type_name << '>'
            << " (Read): " << dataset.get_size() << std::endl;
    deallog << "Rank " + dataset_name << '<' << type_name << '>'
            << " (Read): " << dataset.get_rank() << std::endl;
    {
      auto data = dataset.read<FullMatrix<Number>>();

      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Read): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_global_no_collective_cause()
                << std::endl;

      deallog << "Data " + dataset_name << '<' << type_name << '>'
              << " (Read): " << container_to_string(data) << std::endl;
    }

    {
      // The output should be the column vector [7, 10, 13]
      const std::vector<hsize_t> hyperslab_offset = {2, 1};
      const std::vector<hsize_t> hyperslab_count  = {3, 1};
      auto data = dataset.read_hyperslab<std::vector<Number>>(hyperslab_offset,
                                                              hyperslab_count);

      std::cout << "Column vector" << std::endl;
      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Read): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_global_no_collective_cause()
                << std::endl;

      deallog << "Column vector " + dataset_name << '<' << type_name << '>'
              << " (Read): " << container_to_string(data) << std::endl;
    }

    {
      // The output should be the row vector [9, 10, 11]
      const std::vector<hsize_t> hyperslab_offset = {3, 0};
      const std::vector<hsize_t> hyperslab_count  = {1, 3};
      auto data = dataset.read_hyperslab<std::vector<Number>>(hyperslab_offset,
                                                              hyperslab_count);

      std::cout << "Row vector" << std::endl;
      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Read): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_global_no_collective_cause()
                << std::endl;

      deallog << "Row vector " + dataset_name << '<' << type_name << '>'
              << " (Read): " << container_to_string(data) << std::endl;
    }

    {
      // The output should be the matrix [[3, 4, 5], [6, 7, 8]]
      const std::vector<hsize_t> hyperslab_offset = {1, 0};
      const std::vector<hsize_t> hyperslab_count  = {2, 3};
      auto data = dataset.read_hyperslab<FullMatrix<Number>>(hyperslab_offset,
                                                             hyperslab_count);

      std::cout << "Sub-matrix" << std::endl;
      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Read): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_global_no_collective_cause()
                << std::endl;

      deallog << "Sub-matrix " + dataset_name << '<' << type_name << '>'
              << " (Read): " << container_to_string(data) << std::endl;
    }
  }

  {
    std::string dataset_name("dataset_2");

    // This dataset is a 2x5x3 tensor
    // Output generated by h5dump
    //  (0,0,0): 0, 1, 2,         (1,0,0): 15, 16, 17,
    //  (0,1,0): 3, 4, 5,         (1,1,0): 18, 19, 20,
    //  (0,2,0): 6, 7, 8,         (1,2,0): 21, 22, 23,
    //  (0,3,0): 9, 10, 11,       (1,3,0): 24, 25, 26,
    //  (0,4,0): 12, 13, 14,      (1,4,0): 27, 28, 29

    auto dataset = group.open_dataset(dataset_name);
    dataset.set_query_io_mode(true);
    deallog << "Dimensions " + dataset_name << '<' << type_name << '>'
            << " (Read): " << dataset.get_dimensions() << std::endl;
    deallog << "Size " + dataset_name << '<' << type_name << '>'
            << " (Read): " << dataset.get_size() << std::endl;
    deallog << "Rank " + dataset_name << '<' << type_name << '>'
            << " (Read): " << dataset.get_rank() << std::endl;
    {
      auto data = dataset.read<std::vector<Number>>();

      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Read): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_global_no_collective_cause()
                << std::endl;

      deallog << "Data " + dataset_name << '<' << type_name << '>'
              << " (Read): " << container_to_string(data) << std::endl;
    }

    {
      // The output should be the column vector [16, 19, 22, 25]
      const std::vector<hsize_t> hyperslab_offset = {1, 0, 1};
      const std::vector<hsize_t> hyperslab_count  = {1, 4, 1};
      auto data = dataset.read_hyperslab<std::vector<Number>>(hyperslab_offset,
                                                              hyperslab_count);

      std::cout << "Column vector" << std::endl;
      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Read): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_global_no_collective_cause()
                << std::endl;

      deallog << "Column vector " + dataset_name << '<' << type_name << '>'
              << " (Read): " << container_to_string(data) << std::endl;
    }

    {
      // The output should be the row vector [1, 16]
      const std::vector<hsize_t> hyperslab_offset = {0, 0, 1};
      const std::vector<hsize_t> hyperslab_count  = {2, 1, 1};
      auto data = dataset.read_hyperslab<std::vector<Number>>(hyperslab_offset,
                                                              hyperslab_count);

      std::cout << "Row vector" << std::endl;
      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Read): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_global_no_collective_cause()
                << std::endl;

      deallog << "Row vector " + dataset_name << '<' << type_name << '>'
              << " (Read): " << container_to_string(data) << std::endl;
    }

    {
      // The output should be the matrix [[18, 19], [21, 22], [24, 25]]
      const std::vector<hsize_t> hyperslab_offset = {1, 1, 0};
      const std::vector<hsize_t> hyperslab_count  = {1, 3, 2};
      auto data = dataset.read_hyperslab<FullMatrix<Number>>(hyperslab_offset,
                                                             hyperslab_count);

      std::cout << "Sub-matrix" << std::endl;
      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Read): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_global_no_collective_cause()
                << std::endl;

      deallog << "Sub-matrix " + dataset_name << '<' << type_name << '>'
              << " (Read): " << container_to_string(data) << std::endl;
    }

    {
      // The output should be the matrix [[3, 4, 5], [18, 19, 20]]
      const std::vector<hsize_t> hyperslab_offset = {0, 1, 0};
      const std::vector<hsize_t> hyperslab_count  = {2, 1, 3};
      auto data = dataset.read_hyperslab<FullMatrix<Number>>(hyperslab_offset,
                                                             hyperslab_count);

      std::cout << "Sub-matrix" << std::endl;
      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Read): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_global_no_collective_cause()
                << std::endl;

      deallog << "Sub-matrix " + dataset_name << '<' << type_name << '>'
              << " (Read): " << container_to_string(data) << std::endl;
    }

    {
      // The output should be the matrix [[0, 2], [6, 8], [12, 14]]

      const std::vector<hsize_t> data_dimensions = {1, 3, 2};
      const std::vector<hsize_t> offset          = {0, 0, 0};
      const std::vector<hsize_t> stride          = {1, 2, 2};
      const std::vector<hsize_t> count           = {1, 3, 2};
      const std::vector<hsize_t> block           = {1, 1, 1};

      auto data = dataset.read_hyperslab<FullMatrix<Number>>(
        data_dimensions, offset, stride, count, block);

      std::cout << "Sub-matrix" << std::endl;
      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Read): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_global_no_collective_cause()
                << std::endl;

      deallog << "Sub-matrix " + dataset_name << '<' << type_name << '>'
              << " (Read): " << container_to_string(data) << std::endl;
    }

    {
      // The output should be the matrix [0, 20, 22, 13, 11]
      std::vector<hsize_t> coordinates = {0,
                                          0,
                                          0, // first point
                                          1,
                                          1,
                                          2, // second point
                                          1,
                                          2,
                                          1, // third point
                                          0,
                                          4,
                                          1, // fourth point
                                          0,
                                          3,
                                          2}; // fifth point

      auto data = dataset.read_selection<std::vector<Number>>(coordinates);

      std::cout << "Sub-matrix" << std::endl;
      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Read): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_global_no_collective_cause()
                << std::endl;

      deallog << "Selection std::vector " + dataset_name << '<' << type_name
              << '>' << " (Read): " << container_to_string(data) << std::endl;
    }

    {
      // The output should be the matrix [0, 20, 22, 13, 11]
      std::vector<hsize_t> coordinates = {0,
                                          0,
                                          0, // first point
                                          1,
                                          1,
                                          2, // second point
                                          1,
                                          2,
                                          1, // third point
                                          0,
                                          4,
                                          1, // fourth point
                                          0,
                                          3,
                                          2}; // fifth point

      auto data = dataset.read_selection<Vector<Number>>(coordinates);

      std::cout << "Sub-matrix" << std::endl;
      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Read): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_global_no_collective_cause()
                << std::endl;

      deallog << "Selection Vector " + dataset_name << '<' << type_name << '>'
              << " (Read): " << container_to_string(data) << std::endl;
    }
  }

  {
    std::string dataset_name("dataset_3");

    // This dataset is a 2x5x3 tensor
    // Output generated by h5dump
    //  (0,0,0): 0, 1, 2,         (1,0,0): 0, 2, 4,
    //  (0,1,0): 3, 4, 5,         (1,1,0): 6, 8, 10,
    //  (0,2,0): 6, 7, 8,         (1,2,0): 12, 14, 16,
    //  (0,3,0): 9, 10, 11,       (1,3,0): 18, 20, 22,
    //  (0,4,0): 12, 13, 14,      (1,4,0): 24, 26, 28

    auto dataset = group.open_dataset(dataset_name);
    dataset.set_query_io_mode(true);
    deallog << "Dimensions " + dataset_name << '<' << type_name << '>'
            << " (Read): " << dataset.get_dimensions() << std::endl;
    deallog << "Size " + dataset_name << '<' << type_name << '>'
            << " (Read): " << dataset.get_size() << std::endl;
    deallog << "Rank " + dataset_name << '<' << type_name << '>'
            << " (Read): " << dataset.get_rank() << std::endl;
    {
      auto data = dataset.read<std::vector<Number>>();

      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Read): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_global_no_collective_cause()
                << std::endl;

      deallog << "Data " + dataset_name << '<' << type_name << '>'
              << " (Read): " << container_to_string(data) << std::endl;
    }
  }

  {
    std::string dataset_name("dataset_4");

    // This dataset is a 3x2x4 tensor
    // Output generated by h5dump
    //  (0,0,0): 0, 1, 2, 3,
    //  (0,1,0): 0, 2, 4, 6,
    //  (1,0,0): 4, 5, 6, 7,
    //  (1,1,0): 8, 10, 12, 14,
    //  (2,0,0): 8, 9, 10, 11,
    //  (2,1,0): 16, 18, 20, 22

    auto dataset = group.open_dataset(dataset_name);
    dataset.set_query_io_mode(true);
    deallog << "Dimensions " + dataset_name << '<' << type_name << '>'
            << " (Read): " << dataset.get_dimensions() << std::endl;
    deallog << "Size " + dataset_name << '<' << type_name << '>'
            << " (Read): " << dataset.get_size() << std::endl;
    deallog << "Rank " + dataset_name << '<' << type_name << '>'
            << " (Read): " << dataset.get_rank() << std::endl;
    {
      auto data = dataset.read<std::vector<Number>>();

      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Read): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_global_no_collective_cause()
                << std::endl;

      deallog << "Data " + dataset_name << '<' << type_name << '>'
              << " (Read): " << container_to_string(data) << std::endl;
    }
  }

  {
    std::string dataset_name("dataset_5");

    // This dataset is a 3x2x4 tensor
    // Output generated by h5dump
    //  (0,0,0): 2, 0, 0, 38,
    //  (0,1,0): 6, 5, 0, 3,
    //  (1,0,0): 12, 0, 0, 0,
    //  (1,1,0): 16, 15, 0, 13,
    //  (2,0,0): 32, 0, 0, 0,
    //  (2,1,0): 36, 35, 0, 33


    auto dataset = group.open_dataset(dataset_name);
    dataset.set_query_io_mode(true);
    deallog << "Dimensions " + dataset_name << '<' << type_name << '>'
            << " (Read): " << dataset.get_dimensions() << std::endl;
    deallog << "Size " + dataset_name << '<' << type_name << '>'
            << " (Read): " << dataset.get_size() << std::endl;
    deallog << "Rank " + dataset_name << '<' << type_name << '>'
            << " (Read): " << dataset.get_rank() << std::endl;
    {
      auto data = dataset.read<std::vector<Number>>();

      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Read): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_global_no_collective_cause()
                << std::endl;

      deallog << "Data " + dataset_name << '<' << type_name << '>'
              << " (Read): " << container_to_string(data) << std::endl;
    }
  }

  {
    std::string dataset_name("dataset_6");

    // This dataset is a 7x12 matrix
    // Output generated by h5dump
    //  (0,0): -0, 0, 1, -1, 2, 3, -2, 4, 5, -3, 6, 7,
    //  (1,0): -4, 8, 9, -5, 10, 11, -6, 12, 13, -7, 14, 15,
    //  (2,0): -8, 16, 17, -9, 18, 19, -10, 20, 21, -11, 22, 23,
    //  (3,0): -12, 0, 0, -13, 0, 0, -14, 0, 0, -15, 0, 0,
    //  (4,0): -16, 24, 25, -17, 26, 27, -18, 28, 29, -19, 30, 31,
    //  (5,0): -20, 32, 33, -21, 34, 35, -22, 36, 37, -23, 38, 39,
    //  (6,0): -24, 40, 41, -25, 42, 43, -26, 44, 45, -27, 46, 47


    auto dataset = group.open_dataset(dataset_name);
    dataset.set_query_io_mode(true);
    deallog << "Dimensions " + dataset_name << '<' << type_name << '>'
            << " (Read): " << dataset.get_dimensions() << std::endl;
    deallog << "Size " + dataset_name << '<' << type_name << '>'
            << " (Read): " << dataset.get_size() << std::endl;
    deallog << "Rank " + dataset_name << '<' << type_name << '>'
            << " (Read): " << dataset.get_rank() << std::endl;
    {
      auto data = dataset.read<FullMatrix<Number>>();

      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Read): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_global_no_collective_cause()
                << std::endl;

      deallog << "Data " + dataset_name << '<' << type_name << '>'
              << " (Read): " << container_to_string(data) << std::endl;
    }
  }

  {
    std::string dataset_name("dataset_7");

    // This dataset is a 3x2x4 tensor
    // Output generated by h5dump
    //  (0,0,0): 2, 0, 0, 38,
    //  (0,1,0): 6, 5, 0, 3,
    //  (1,0,0): 12, 0, 0, 0,
    //  (1,1,0): 16, 15, 0, 13,
    //  (2,0,0): 32, 0, 0, 0,
    //  (2,1,0): 36, 35, 0, 33


    auto dataset = group.open_dataset(dataset_name);
    dataset.set_query_io_mode(true);
    deallog << "Dimensions " + dataset_name << '<' << type_name << '>'
            << " (Read): " << dataset.get_dimensions() << std::endl;
    deallog << "Size " + dataset_name << '<' << type_name << '>'
            << " (Read): " << dataset.get_size() << std::endl;
    deallog << "Rank " + dataset_name << '<' << type_name << '>'
            << " (Read): " << dataset.get_rank() << std::endl;
    {
      auto data = dataset.read<std::vector<Number>>();

      std::cout << "IO mode " + dataset_name << '<' << type_name << '>'
                << " (Read): " << dataset.get_io_mode() << std::endl;
      std::cout << "Local no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_local_no_collective_cause()
                << std::endl;
      std::cout << "Global no collective cause " << dataset_name << '<'
                << type_name << '>'
                << " (Read): " << dataset.get_global_no_collective_cause()
                << std::endl;

      deallog << "Data " + dataset_name << '<' << type_name << '>'
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
      std::string filename = "test.h5";

      {
        HDF5::File data_file(filename, HDF5::File::FileAccessMode::create);

        write_test<float>(data_file);
        write_test<double>(data_file);

#ifdef DEAL_II_WITH_COMPLEX_VALUES
        write_test<std::complex<float>>(data_file);
        write_test<std::complex<double>>(data_file);
#endif
      }

      {
        HDF5::File data_file(filename, HDF5::File::FileAccessMode::open);
        read_test<float>(data_file);
        read_test<double>(data_file);

#ifdef DEAL_II_WITH_COMPLEX_VALUES
        read_test<std::complex<float>>(data_file);
        read_test<std::complex<double>>(data_file);
#endif
      }
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

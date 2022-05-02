// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
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


// Test ConstantVectorFunction class

#include <deal.II/base/function.h>

#include <deal.II/lac/vector.h>

#include <array>
#include <vector>

#include "../tests.h"

#define MAX_DIM 3
#define MAX_N_COMPONENT 7
#define NUMBER double

#define TESTEE Functions::ConstantFunction

// Test a given TESTEE object f on n_points points
template <int dim, typename Number>
void
test_one_object(const TESTEE<dim> &f,
                const unsigned int n_points,
                const unsigned int n_component)
{
  // Create a vector of random points
  std::vector<Point<dim>> points(n_points);
  for (unsigned int i_point = 0; i_point < n_points; ++i_point)
    {
      points[i_point] = Point<dim>();
      for (unsigned int id = 0; id < dim; ++id)
        {
          points[i_point][id] = (Number)Testing::rand();
        }
    }

  // Do test
  deallog.flags(std::ios::scientific);

  // value
  {
    deallog << f.value(points[2], 0) << std::endl;
    deallog << f.value(points[2], n_component - 1) << std::endl;
  }

  // vector_value
  {
    deallog << std::endl;
    Vector<Number> retune_value(n_component);

    f.vector_value(points[0], retune_value);
    retune_value.print(deallog.get_file_stream(), /*precision =*/6);

    f.vector_value(points[4], retune_value);
    retune_value.print(deallog.get_file_stream(), /*precision =*/6);
  }

  // value_list
  {
    deallog << std::endl;
    std::vector<Number> return_values(n_points);
    f.value_list(points, return_values /*, component = 0*/);
    for (unsigned int i = 0; i < n_points; ++i)
      {
        deallog << return_values[i] << ' ';
      }
    deallog << std::endl;

    f.value_list(points, return_values, std::min(n_component - 1, 3u));
    for (unsigned int i = 0; i < n_points; ++i)
      {
        deallog << return_values[i] << ' ';
      }
    deallog << std::endl;
  }

  // vector_value_list
  {
    deallog << std::endl;
    std::vector<Vector<Number>> return_values(n_points,
                                              Vector<Number>(n_component));
    f.vector_value_list(points, return_values);
    for (unsigned int p = 0; p < n_points; ++p)
      {
        for (unsigned int c = 0; c < n_component; ++c)
          deallog << return_values[p][c] << ' ';
        deallog << std::endl;
      }
  }

  points.clear();

  return;
}

// Construct TESTEE with different constructor and test them
template <int dim, typename Number>
void
test_constructor(const unsigned int         n_component,
                 const std::vector<NUMBER> &component_data)
{
  const unsigned int n_points(5);

  // Construct with value and n_component
  {
    deallog.flags(std::ios::fixed);
    deallog << "\nIn " << dim << " dimension, "
            << " Function with " << n_component << " components, "
            << " Constructed with single value and n_component" << std::endl;
    // Close file to let Test use it.

    TESTEE<dim> constant_vector_function(
      component_data[2 /*of MAX_N_COMPONENT*/], n_component);
    test_one_object<dim, Number>(constant_vector_function,
                                 n_points,
                                 n_component);
  }

  // Construct with std::vector<>
  {
    std::vector<Number> v(n_component);
    for (unsigned int i = 0; i < n_component; ++i)
      v[i] = component_data[i];

    deallog.flags(std::ios::fixed);
    deallog << "\nIn " << dim << " dimension, "
            << " Function with " << n_component << " components, "
            << " Constructed with std::vector:" << std::endl;
    // Close file to let Test use it.

    TESTEE<dim> constant_vector_function(v);
    test_one_object<dim, Number>(constant_vector_function,
                                 n_points,
                                 n_component);
  }

  // Construct with Vector
  {
    Vector<Number> v(n_component);
    for (unsigned int i = 0; i < n_component; ++i)
      v[i] = component_data[i];

    deallog.flags(std::ios::fixed);
    deallog << "\nIn " << dim << " dimension, "
            << " Function with " << n_component << " components, "
            << " Constructed with Vector:" << std::endl;
    // Close file to let Test use it.

    TESTEE<dim> constant_vector_function(v);
    test_one_object<dim, Number>(constant_vector_function,
                                 n_points,
                                 n_component);
  }

  // Construct with pointer and offset
  {
    deallog.flags(std::ios::fixed);
    deallog << "\nIn " << dim << " dimension, "
            << " Function with " << n_component << " components, "
            << " Constructed with pointer and offset:" << std::endl;
    // Close file to let Test use it.

    TESTEE<dim> constant_vector_function(&component_data[0], n_component);
    test_one_object<dim, Number>(constant_vector_function,
                                 n_points,
                                 n_component);
  }

  return;
}

// test different number of vector components
template <int dim, typename Number>
void
test_n_components(const std::vector<NUMBER> &component_data)
{
#define N_COMPONENT_CASE 3
  const unsigned int components[N_COMPONENT_CASE] = {1, 4, 7};

  for (unsigned int c = 0; c < N_COMPONENT_CASE; ++c)
    {
      test_constructor<dim, Number>(components[c], component_data);
    }

  return;
}

int
main()
{
  initlog();
  // Create data for each component

  std::vector<NUMBER> component_data(MAX_N_COMPONENT);
  {
    NUMBER data[MAX_N_COMPONENT] = {1, 255, 32768, 98.2, 75, 7.9e+8, 6.7e-4};
    component_data.assign(data, data + MAX_N_COMPONENT);
  }

  test_n_components<1, NUMBER>(component_data);
  test_n_components<2, NUMBER>(component_data);
  test_n_components<3, NUMBER>(component_data);

  return (0);
}

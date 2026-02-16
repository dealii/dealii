// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Check DataOutPoints class

#include "deal.II/numerics/data_component_interpretation.h"
#include <deal.II/numerics/data_out_points.h>

#include "../tests.h"



int
main()
{
  initlog();

  {
    deallog << "* dim <2> no data:" << std::endl;
    DataOutPoints<2>      data_out;
    std::vector<Point<2>> points;
    points.emplace_back(0., 0.);
    points.emplace_back(1., 2.);
    data_out.build_patches(points, 1000);
    data_out.write_gnuplot(deallog.get_file_stream());
    deallog << "OK" << std::endl;
  }

  {
    deallog << "* dim <2,3> scalar and vector data:" << std::endl;
    DataOutPoints<2, 3>   data_out;
    std::vector<Point<3>> points;
    points.emplace_back(0., 0., 0.);
    points.emplace_back(1., 2., 3.);
    std::vector<std::vector<double>> data;
    data.emplace_back(std::vector<double>{1., 2., 3., 4.});
    data.emplace_back(std::vector<double>{4., 5., 6., 7.});
    std::vector<std::string> data_component_names;
    data_component_names.emplace_back("c1");
    data_component_names.emplace_back("tensor");
    data_component_names.emplace_back("tensor");
    data_component_names.emplace_back("tensor");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretations;
    data_component_interpretations.emplace_back(
      DataComponentInterpretation::component_is_scalar);
    data_component_interpretations.emplace_back(
      DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretations.emplace_back(
      DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretations.emplace_back(
      DataComponentInterpretation::component_is_part_of_vector);
    data_out.build_patches(
      points, 0, data, data_component_names, data_component_interpretations);
    data_out.write_gnuplot(deallog.get_file_stream());
    deallog << "OK" << std::endl;
  }
}

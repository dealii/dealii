// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test reading of data contained in a vtk file into a serial deal.II
// Vector<double> using the VTKWrappers::data_to_dealii_vector function.

#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>

#include <deal.II/vtk/utilities.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test(const std::string &filename)
{
  deallog << "Reading Triangulation<" << dim << ", " << spacedim
          << "> from file: " << filename << std::endl;

  Triangulation<dim, spacedim> triangulation;
  DoFHandler<dim, spacedim>    dh(triangulation);
  Vector<double>               vtk_data;
  Vector<double>               dealii_data;

  VTKWrappers::read_tria(std::string(SOURCE_DIR) + "/" + filename,
                         triangulation);

  VTKWrappers::read_all_data(std::string(SOURCE_DIR) + "/" + filename,
                             vtk_data);

  const auto [fe, data_names] =
    VTKWrappers::vtk_to_finite_element<dim, spacedim>(std::string(SOURCE_DIR) +
                                                      "/" + filename);

  dh.distribute_dofs(*fe);
  dealii_data.reinit(dh.n_dofs());

  VTKWrappers::data_to_dealii_vector(triangulation, vtk_data, dh, dealii_data);

  deallog << "All data size: " << vtk_data.size() << std::endl;
  deallog << "Dealii vector size: " << dealii_data.size() << std::endl;

  Assert(vtk_data.size() == dealii_data.size(), ExcInternalError());

  // Output using data out
  if (std::getenv("OUTPUT_VTK"))
    {
      DataOut<dim, spacedim> data_out;
      data_out.attach_dof_handler(dh);
      data_out.add_data_vector(dealii_data, "data");
      data_out.build_patches();

      std::ofstream out("output_" + std::to_string(dim) + "_" +
                        std::to_string(spacedim) + ".vtk");
      data_out.write_vtk(out);
    }
}

int
main()
{
  initlog();
  test<1, 1>("data/data_1_1.vtk");
  test<1, 2>("data/data_1_1.vtk");
  test<1, 3>("data/data_1_1.vtk");

  test<2, 2>("data/data_2_2_simplex.vtk");
  test<2, 3>("data/data_2_2_simplex.vtk");

  test<2, 2>("data/data_2_2_quad.vtk");
  test<2, 3>("data/data_2_2_quad.vtk");

  test<3, 3>("data/data_3_3_tets.vtk");
  test<3, 3>("data/data_3_3_hexes.vtk");
  return 0;
}

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check that gmsh api correctly reads and writes a mesh with manifold
// information, in all coordinate dimensions

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <boost/algorithm/string.hpp>

#include <fstream>

#include "../tests.h"

/*
 * Test that we can use the gmsh-API version of GridIn::read_msh, when the
 * description of the "Physical name" does not contain MaterialID, BoundaryID,
 * or ManifoldID.
 *
 * Create a hypercube triangulation, write it to file so that we get a valid
 * msh-file. Read the filecontent back in and replace MaterialId with a string
 * that will not be recognized. Write it back to the file and try to read the
 * triangulation from the msh-file.
 */
template <int dim>
void
test()
{
  // Create a mesh, write it to file.
  Triangulation<dim> original_mesh;
  const double       left     = 0.;
  const double       right    = 1.;
  const bool         colorize = true;
  GridGenerator::hyper_cube(original_mesh, left, right, colorize);
  // Set materialID to something nonzero so that we
  // get a physical name "MaterialID: 1" in the msh-file.
  auto cell_on_original = original_mesh.begin_active();
  cell_on_original->set_material_id(17);

  const std::string mesh_filename = "output.msh";
  GridOut           grid_out;
  grid_out.write_msh(original_mesh, mesh_filename);

  // Read the file back and change the physical group description to be
  // unrecognizable.
  std::ifstream     input(mesh_filename, std::ifstream::in);
  std::stringstream buffer;
  buffer << input.rdbuf();
  std::string filecontent = buffer.str();
  input.close();

  boost::replace_all(filecontent, "MaterialID", "Unrecognizable");
  boost::replace_all(filecontent, "BoundaryID", "Unrecognizable");
  boost::replace_all(filecontent, "ManifoldID", "Unrecognizable");

  // Write the mesh to file and read the triangulation back in.
  std::ofstream output(mesh_filename);
  output << filecontent;
  output.close();

  Triangulation<dim> triangulation_from_file;
  GridIn<dim>        grid_in(triangulation_from_file);
  grid_in.read_msh(mesh_filename);

  // The function should set material_id and boundary_id based on the physical
  // tag. Write these to log to check them.
  const auto cell = triangulation_from_file.begin_active();
  deallog << "MaterialID = " << cell->material_id() << std::endl;
  for (unsigned int f : cell->face_indices())
    deallog << "BoundaryID = " << cell->face(f)->boundary_id() << std::endl;
}

int
main(int argc, char **argv)
{
  // gmsh might be build with mpi support enabled.
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  initlog();

  test<2>();
}

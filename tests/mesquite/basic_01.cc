/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2018 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 */

/**
 * This demonstrates the application of a high-level API, namely the Laplace
 * smoother, to a distorted grid.
 *
 * The tested mesh is courtesy of the Mesquite project.
 * This simple program is part of a tutorial found in the user guide.
 *
 * Knupp, P.; Freitag-Diachin, L. & Tidwell, B.
 * Mesquite Mesh Quality Improvement Toolkit User's Guide
 * Sandia National Laboratories, Sandia National Laboratories, 2013
 *
 * See sec 3.1.3 Improving the Mesh with a Wrapper Class
 */


#include <Mesquite_all_headers.hpp>

#include <fstream>

#include "../tests.h"

#include "mesquite_utilities.h"


int
main()
{
  initlog();

  // Error tracker
  Mesquite2::MsqError err;

  // Read in a mesh file
  const std::string   filename_input = SOURCE_DIR "/mesh/2d/tangled_quad.vtk";
  Mesquite2::MeshImpl mesh;
  mesh.read_vtk(filename_input.c_str(), err);
  Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

  // Write it out again for reference...
  const std::string filename_output_1 = "tangled_quad.01.vtk";
  mesh.write_vtk(filename_output_1.c_str(), err);
  Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

  // Make a plane on which the mesh lies
  const Mesquite2::Vector3D normal(0, 0, 1);
  const Mesquite2::Vector3D point(0, 0, 5);
  Mesquite2::PlanarDomain   mesh_plane(normal, point);

  // Build a view of the domain
  Mesquite2::MeshDomainAssoc mesh_and_domain(&mesh, &mesh_plane);

  // Perform Laplace smoothing
  Mesquite2::LaplaceWrapper mesh_quality_algorithm;
  mesh_quality_algorithm.run_instructions(&mesh_and_domain, err);
  Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

  // Write out the smoothed mesh
  const std::string filename_output_2 = "tangled_quad.02.vtk";
  mesh.write_vtk(filename_output_2.c_str(), err);
  Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

  deallog << "OK" << std::endl;
}

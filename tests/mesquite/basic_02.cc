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
 * This demonstrates the application of a low-level API to a distorted grid.
 *
 * The tested mesh is courtesy of the Mesquite project.
 * This simple program is part of a tutorial found in the user guide.
 *
 * Knupp, P.; Freitag-Diachin, L. & Tidwell, B.
 * Mesquite Mesh Quality Improvement Toolkit User's Guide
 * Sandia National Laboratories, Sandia National Laboratories, 2013
 *
 * See sec 3.1.4 Improving the Mesh with the Low Level API
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

  // ---------------------

  // Creates a mean ratio quality metric
  Mesquite2::IdealWeightInverseMeanRatio inverse_mean_ratio(err);
  Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

  // Prescribe the objective function template
  Mesquite2::LPtoPTemplate obj_func(&inverse_mean_ratio, 2, err);
  Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

  // Creates the optimisation procedures
  Mesquite2::TrustRegion t_region(&obj_func);

  // Performs global optimisation
  t_region.use_global_patch();

  // Create a termination criterion and add it to
  // the optimisation procedure
  // outer loop: default behavior: 1 iteration
  // inner loop: stop if gradient norm < eps
  Mesquite2::TerminationCriterion tc_inner;
  tc_inner.add_absolute_gradient_L2_norm(1.0e-4);
  t_region.set_inner_termination_criterion(&tc_inner);

  // Create a quality assessor
  Mesquite2::QualityAssessor m_ratio_qa(&inverse_mean_ratio);

  // Creates an instruction queue
  Mesquite2::InstructionQueue queue;
  queue.add_quality_assessor(&m_ratio_qa, err);
  queue.set_master_quality_improver(&t_region, err);
  queue.add_quality_assessor(&m_ratio_qa, err);

  // ---------------------

  // Make a plane on which the mesh lies
  const Mesquite2::Vector3D normal(0, 0, -1);
  const Mesquite2::Vector3D point(0, 0, -5);
  Mesquite2::PlanarDomain   mesh_plane(normal, point);

  // Read in a mesh file
  const std::string   filename_input = SOURCE_DIR "/mesh/2d/hole_in_square.vtk";
  Mesquite2::MeshImpl mesh;
  mesh.read_vtk(filename_input.c_str(), err);
  Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

  // Write it out again for reference...
  const std::string filename_output_1 = "hole_in_square.01.vtk";
  mesh.write_vtk(filename_output_1.c_str(), err);
  Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

  // Build a view of the domain
  Mesquite2::MeshDomainAssoc mesh_and_domain(&mesh, &mesh_plane);

  // Perform optimisation of the mesh
  queue.run_instructions(&mesh_and_domain, err);
  Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

  // Write out the smoothed mesh
  const std::string filename_output_2 = "hole_in_square.02.vtk";
  mesh.write_vtk(filename_output_2.c_str(), err);
  Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

  deallog << "OK" << std::endl;
}

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
 * A 2d example of smoothing of a distorted grid using the Mesquite interface.
 * The grid is in its coarsest state, with no hanging nodes.
 *
 * - The vertices of the triangulation are moved to a distorted position
 * - Mesh adaption is performed
 *     - A customised smoother is used
 *     - The boundary vertices are considered fixed
 * - The vertices of the triangulation are updated to their optimial locations
 */


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_mesquite.h>
#include <deal.II/grid/tria.h>

#include <fstream>

#include "../tests.h"
#include "mesquite_utilities.h"


int
main()
{
  initlog();

  const unsigned int dim                 = 2;
  const unsigned int n_elements_per_side = 4;

  // Make a distored grid
  Triangulation<dim>        tria;
  std::vector<unsigned int> repetitions(dim, n_elements_per_side);
  GridGenerator::subdivided_hyper_rectangle(tria,
                                            repetitions,
                                            Point<2>(0.0, 0.0),
                                            Point<2>(1.0, 1.0));
  GridTools::distort_random(0.3, tria, true);

  // Output the distorted grid
  std::cout << "Distorted mesh" << std::endl;
  output_mesh(tria, "grid.01");

  // Initialize the mesh smoothing interface
  GridTools::MesquiteMeshInterface<dim> mesquite_interface(
    tria, true /*fix_all_boundary_vertices*/);

  // Error tracker
  Mesquite2::MsqError err;

  // Create a mesh smoothing algorithm:
  // - Create a mean ratio quality metric
  Mesquite2::IdealWeightInverseMeanRatio inverse_mean_ratio(err);
  Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
  // - Set the objective function template
  Mesquite2::LPtoPTemplate obj_func(&inverse_mean_ratio, 2, err);
  Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
  // - Create the optimization procedures
  Mesquite2::TrustRegion t_region(&obj_func);
  // - Perform optimization globally
  t_region.use_global_patch();

  // - Create a termination criterion and add it to the optimization procedure
  //    - Outer loop: default behavior: 1 iteration
  //    - Inner loop: stop if gradient norm < eps
  Mesquite2::TerminationCriterion tc_inner;
  tc_inner.add_absolute_gradient_L2_norm(1.0e-4);
  t_region.set_inner_termination_criterion(&tc_inner);

  // - Create a quality assessor
  Mesquite2::QualityAssessor m_ratio_qa(&inverse_mean_ratio);

  // - Create an instruction queue
  Mesquite2::InstructionQueue queue;
  queue.add_quality_assessor(&m_ratio_qa, err);
  Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
  queue.set_master_quality_improver(&t_region, err);
  Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
  queue.add_quality_assessor(&m_ratio_qa, err);
  Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

  // Perform the mesh optimization
  mesquite_interface.execute(queue);

  // Move the triangulation's vertices
  mesquite_interface.move_triangulation_vertices(tria);

  // Output the smoothed grid
  std::cout << "Smoothed mesh" << std::endl;
  output_mesh(tria, "grid.02", true /*write_to_deallog*/);

  std::cout << "OK" << std::endl;
}

/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2019 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

// Compute integral of known function using mesh_loop and ScratchData

#include <deal.II/base/function_parser.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/meshworker/copy_data.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/meshworker/scratch_data.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <unordered_map>

#include "../tests.h"

using namespace MeshWorker;

template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> tria;
  FE_Q<dim, spacedim>          fe(1);
  DoFHandler<dim, spacedim>    dh(tria);

  FunctionParser<spacedim> integral_function("x");
  GridGenerator::hyper_cube(tria);
  tria.refine_global(5);
  tria.execute_coarsening_and_refinement();
  dh.distribute_dofs(fe);

  Vector<double> solution(dh.n_dofs());

  VectorTools::interpolate(dh, integral_function, solution);

  QGauss<dim> quad(3);

  UpdateFlags cell_flags = update_values | update_gradients |
                           update_quadrature_points | update_JxW_values;

  using ScratchData = ScratchData<dim, spacedim>;
  struct CopyData
  {};


  double                           H1_norm = 0;
  const FEValuesExtractors::Scalar scalar(0);

  ScratchData scratch(fe, quad, cell_flags);
  CopyData    copy;

  auto cell = dh.begin_active();
  auto endc = dh.end();

  using Iterator = decltype(cell);


  auto cell_integrator = [&H1_norm, &solution, &scalar](const Iterator &cell,
                                                        ScratchData    &s,
                                                        CopyData       &c) {
    const auto &fev = s.reinit(cell);
    const auto &JxW = s.get_JxW_values();

    s.extract_local_dof_values("solution", solution);

    const auto &u      = s.get_values("solution", scalar);
    const auto &grad_u = s.get_gradients("solution", scalar);

    for (unsigned int q = 0; q < u.size(); ++q)
      H1_norm += (u[q] * u[q] + grad_u[q] * grad_u[q]) * JxW[q];
  };

  const auto empty_copyer = [](const CopyData &) {};

  mesh_loop(cell,
            endc,
            cell_integrator,
            empty_copyer,
            scratch,
            copy,
            assemble_own_cells);

  deallog << "H1 norm: " << std::sqrt(H1_norm) << std::endl;
}


int
main()
{
  initlog(true);
  MultithreadInfo::set_thread_limit(1); // to make output deterministic

  test<2, 2>();
}

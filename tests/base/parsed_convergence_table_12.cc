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

// Check errors and differences on simplex meshes with explicit and default
// mappings. On the unit cube, x_0 has L1 norm 1/2, L2 norm sqrt(1/3), and
// H1 norm sqrt(4/3), independently of the dimension.

#include <deal.II/base/function_lib.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_convergence_table.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools_interpolate.h>

#include "../tests.h"


template <int dim>
void
check()
{
  deallog << "dim = " << dim << std::endl;

  Triangulation<dim> cube, triangulation;
  GridGenerator::hyper_cube(cube);
  GridGenerator::convert_hypercube_to_simplex_mesh(cube, triangulation);

  FE_SimplexP<dim> fe(1);
  MappingFE<dim>   mapping(fe);
  DoFHandler<dim>  dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  Tensor<1, dim> exponents;
  exponents[0] = 1;
  Functions::Monomial<dim> exact(exponents);
  Vector<double>           zero(dof_handler.n_dofs());
  Vector<double>           solution(dof_handler.n_dofs());
  VectorTools::interpolate(mapping, dof_handler, exact, solution);

  ParsedConvergenceTable table(
    {"u"},
    {{VectorTools::L1_norm, VectorTools::L2_norm, VectorTools::H1_norm}});
  ParameterHandler prm;
  table.add_parameters(prm);
  prm.set("Extra columns", "");
  prm.set("Rate key", "");
  prm.set("Rate mode", "none");

  table.error_from_exact(mapping, dof_handler, zero, exact);
  table.error_from_exact(dof_handler, zero, exact);
  table.difference(mapping, dof_handler, solution, zero);
  table.difference(dof_handler, solution, zero);
  table.output_table(deallog.get_file_stream());
}



int
main()
{
  initlog();
  check<2>();
  check<3>();
}

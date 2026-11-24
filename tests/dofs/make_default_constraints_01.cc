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

// Test DoFTools::make_default_constraints() for hanging nodes and different
// scalar types.

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/vector_tools.h>

#include <sstream>

#include "../tests.h"


template <int dim, typename number>
void
check_equal(const AffineConstraints<number> &a,
            const AffineConstraints<number> &b)
{
  std::ostringstream sa, sb;
  a.print(sa);
  b.print(sb);
  AssertThrow(sa.str() == sb.str(),
              ExcMessage("Constraints differ:\n" + sa.str() + "\nvs\n" +
                         sb.str()));
}


template <int dim>
void
test()
{
  deallog << dim << "D" << std::endl;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FE_Q<dim>       fe(2);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  // Default scalar type (double) with additional Dirichlet boundary values
  {
    auto constraints = DoFTools::make_default_constraints(dof_handler);
    AssertThrow(!constraints.is_closed(), ExcInternalError());
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             constraints);
    constraints.close();

    AffineConstraints<double> ref(dof_handler.locally_owned_dofs(),
                                  DoFTools::extract_locally_relevant_dofs(
                                    dof_handler));
    DoFTools::make_hanging_node_constraints(dof_handler, ref);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             ref);
    ref.close();
    check_equal<dim>(constraints, ref);
  }

  // Non-default scalar type (float)
  {
    auto constraints_float =
      DoFTools::make_default_constraints<float>(dof_handler);
    AssertThrow(!constraints_float.is_closed(), ExcInternalError());
    constraints_float.close();

    AffineConstraints<float> ref_float(dof_handler.locally_owned_dofs(),
                                       DoFTools::extract_locally_relevant_dofs(
                                         dof_handler));
    DoFTools::make_hanging_node_constraints(dof_handler, ref_float);
    ref_float.close();
    check_equal<dim>(constraints_float, ref_float);

    deallog << "  hanging nodes OK (" << constraints_float.n_constraints()
            << " constraints)" << std::endl;
  }
}


int
main()
{
  initlog();

  test<2>();
  test<3>();
}

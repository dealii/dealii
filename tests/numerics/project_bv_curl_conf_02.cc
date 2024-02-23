// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// The test checks that project_boundary_values_curl_conforming
// works correctly for high-order FE_Nedelec elements used via
// FESystem. This requires the produced constraints to be the same
// for FE_Nedelec and FESystem(FE_Nedelec, 1).

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim>
class BoundaryFunction : public Function<dim>
{
public:
  BoundaryFunction();
  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const;
};

template <int dim>
BoundaryFunction<dim>::BoundaryFunction()
  : Function<dim>(dim)
{}

template <int dim>
void
BoundaryFunction<dim>::vector_value(const Point<dim> &,
                                    Vector<double> &values) const
{
  for (unsigned int d = 0; d < dim; ++d)
    values(d) = d + 1.0;
}

template <int dim>
void
test_boundary_values(const FiniteElement<dim>  &fe,
                     AffineConstraints<double> &constraints)
{
  Triangulation<dim> triangulation;
  GridGenerator::subdivided_hyper_cube(triangulation, 2);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  BoundaryFunction<dim> boundary_function;
  constraints.clear();
  VectorTools::project_boundary_values_curl_conforming_l2(
    dof_handler,
    0,
    boundary_function,
    0,
    constraints,
    StaticMappingQ1<dim>::mapping);
  constraints.close();
}

template <int dim>
void
test(unsigned order)
{
  deallog << "dim:" << dim << " order:" << order << "\t";

  AffineConstraints<double> constraints_fe, constraints_fes;

  {
    FE_Nedelec<3> fe(order);
    test_boundary_values(fe, constraints_fe);
  }

  {
    FESystem<3> fe(FE_Nedelec<3>(order), 1);
    test_boundary_values(fe, constraints_fes);
  }

  if (constraints_fes.n_constraints() == constraints_fe.n_constraints())
    {
      const IndexSet &lines = constraints_fes.get_local_lines();

      for (unsigned i = 0; i < lines.n_elements(); ++i)
        {
          if (!constraints_fe.is_constrained(lines.nth_index_in_set(i)))
            {
              deallog << "Failed" << std::endl;
              return;
            }

          const std::vector<std::pair<types::global_dof_index, double>> &c1 =
            *constraints_fes.get_constraint_entries(lines.nth_index_in_set(i));
          const std::vector<std::pair<types::global_dof_index, double>> &c2 =
            *constraints_fe.get_constraint_entries(lines.nth_index_in_set(i));

          for (std::size_t j = 0; j < c1.size(); ++j)
            if ((c1[j].first != c2[j].first) ||
                (fabs(c1[j].second - c2[j].second) > 1e-14))
              {
                deallog << "Failed" << std::endl;
                return;
              }
        }
    }
  else
    {
      deallog << "Failed" << std::endl;
      return;
    }

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();
  deallog << std::setprecision(2);

  test<2>(0);
  test<2>(1);
  test<2>(2);

  test<3>(0);
  test<3>(1);
  // test<3>(2);
}

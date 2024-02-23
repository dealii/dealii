// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test the different FEValuesBase::get_function_values

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"


// Call this function with a system consisting of several copies of
// the SAME element
template <int dim>
void
vector_values(const FiniteElement<dim> &fe)
{
  Assert(fe.n_base_elements() == 1, ExcNotImplemented());
  deallog.push(fe.get_name());

  QTrapezoid<dim>                        quadrature;
  std::vector<unsigned int>              renumbering(fe.dofs_per_cell);
  std::vector<std::vector<unsigned int>> component_start;
  FETools::compute_component_wise(fe, renumbering, component_start);

  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);
  DoFRenumbering::component_wise(dof);

  Vector<float> v(dof.n_dofs());
  for (unsigned int i = 0; i < v.size(); ++i)
    v(i) = i;

  FEValues<dim>              feval(fe, quadrature, update_values);
  std::vector<Vector<float>> local(quadrature.size(),
                                   Vector<float>(fe.n_components()));

  typename DoFHandler<dim>::active_cell_iterator cell      = dof.begin_active();
  const typename DoFHandler<dim>::active_cell_iterator end = dof.end();

  unsigned int cell_no = 0;
  while (cell != end)
    {
      deallog << "Cell " << cell_no++ << std::endl;
      feval.reinit(cell);
      feval.get_function_values(v, local);
      for (unsigned int c = 0; c < fe.n_components(); ++c)
        {
          deallog << "Component " << c;
          for (unsigned int k = 0; k < quadrature.size(); ++k)
            deallog << '\t' << (int)local[k](c);
          deallog << std::endl;
        }
      ++cell;
    }
  deallog.pop();
}

template <int dim>
void
test_vectors()
{
  FE_Q<dim> q1(1);

  FESystem<dim> q1_3(q1, 3);
  vector_values(q1_3);
}


int
main()
{
  initlog();

  test_vectors<2>();

  return 0;
}

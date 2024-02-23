// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <string>
#include <vector>

#include "../tests.h"


template <int dim>
inline void
plot_derivatives(Mapping<dim>       &mapping,
                 FiniteElement<dim> &finel,
                 const char         *name)
{
  deallog.push(name);

  Triangulation<dim> tr;
  DoFHandler<dim>    dof(tr);
  GridGenerator::hyper_cube(tr);
  typename DoFHandler<dim>::cell_iterator c = dof.begin();
  dof.distribute_dofs(finel);

  QGauss<dim - 1> q(1);
  //  QIterated<dim> q(q_trapez, div);
  FEFaceValues<dim> fe(mapping,
                       finel,
                       q,
                       UpdateFlags(update_gradients | update_hessians));
  for (const unsigned int face : GeometryInfo<dim>::face_indices())
    {
      fe.reinit(c, face);

      for (unsigned int k = 0; k < q.size(); ++k)
        {
          deallog << "Face " << face << " Point " << q.point(k) << std::endl;
          for (unsigned int i = 0; i < finel.dofs_per_cell; ++i)
            {
              deallog << "\tGrad " << fe.shape_grad(i, k);
              deallog << "\t2nd " << fe.shape_hessian(i, k);
              deallog << std::endl;
            }
        }
    }
  deallog.pop();
}



template <int dim>
void
plot_FE_Q_shape_functions()
{
  MappingQ<dim> m(1);
  //  FE_Q<dim> q1(1);
  //  plot_derivatives(m, q1, "Q1");
  //  plot_face_shape_functions(m, q1, "Q1");
  FE_Q<dim> q2(2);
  plot_derivatives(m, q2, "Q2");
  FE_Q<dim> q3(3);
  plot_derivatives(m, q3, "Q3");
  FE_Q<dim> q4(4);
  plot_derivatives(m, q4, "Q4");
}


template <int dim>
void
plot_FE_DGQ_shape_functions()
{
  MappingQ<dim> m(1);
  FE_DGQ<dim>   q1(1);
  plot_derivatives(m, q1, "DGQ1");
  FE_DGQ<dim> q2(2);
  plot_derivatives(m, q2, "DGQ2");
  FE_DGQ<dim> q3(3);
  plot_derivatives(m, q3, "DGQ3");
  FE_DGQ<dim> q4(4);
  plot_derivatives(m, q4, "DGQ4");
}


int
main()
{
  initlog();
  deallog << std::setprecision(8) << std::fixed;

  deallog.push("2d");
  plot_FE_Q_shape_functions<2>();
  //  plot_FE_DGQ_shape_functions<2>();
  deallog.pop();

  deallog.push("3d");
  //  plot_FE_Q_shape_functions<3>();
  deallog.pop();
  return 0;
}

// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_bernstein.h>
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
//#include "../../include/fe_bernstein.h"

template <int dim>
inline void
plot_derivatives(Mapping<dim> &      mapping,
                 FiniteElement<dim> &finel,
                 const char *        name)
{
  deallog.push(name);

  Triangulation<dim> tr;
  DoFHandler<dim>    dof(tr);
  GridGenerator::hyper_cube(tr, 2., 5.);
  typename DoFHandler<dim>::cell_iterator c = dof.begin();
  dof.distribute_dofs(finel);

  const unsigned int div = 1;

  QTrapez<dim> q;
  //  QIterated<dim> q(q_trapez, div);
  FEValues<dim> fe(mapping,
                   finel,
                   q,
                   UpdateFlags(update_gradients | update_hessians));
  fe.reinit(c);

  unsigned int k = 0;
  for (unsigned int mz = 0; mz <= ((dim > 2) ? div : 0); ++mz)
    {
      for (unsigned int my = 0; my <= ((dim > 1) ? div : 0); ++my)
        {
          for (unsigned int mx = 0; mx <= div; ++mx)
            {
              deallog << q.point(k) << std::endl;

              for (unsigned int i = 0; i < finel.dofs_per_cell; ++i)
                {
                  deallog << "\tGrad " << fe.shape_grad(i, k);
                  deallog << "\t2nd " << fe.shape_hessian(i, k);
                  deallog << std::endl;
                }
              k++;
            }
        }
    }
  deallog.pop();
}



template <int dim>
void
plot_FE_Bernstein_shape_functions()
{
  MappingQGeneric<dim> m(1);
  FE_Bernstein<dim>    b1(1);
  plot_derivatives(m, b1, "B1");

  FE_Bernstein<dim> b2(2);
  plot_derivatives(m, b2, "B2");

  FE_Bernstein<dim> b3(3);
  plot_derivatives(m, b3, "B3");

  FE_Bernstein<dim> b4(4);
  plot_derivatives(m, b4, "B4");
}



int
main()
{
  initlog();
  deallog << std::setprecision(8) << std::fixed;

  deallog.push("1d");
  plot_FE_Bernstein_shape_functions<1>();
  deallog.pop();

  deallog.push("2d");
  plot_FE_Bernstein_shape_functions<2>();

  deallog.pop();

  deallog.push("3d");
  //  plot_FE_Bernstein_shape_functions<3>();
  deallog.pop();


  // FESystem test.
  MappingQGeneric<2> m(1);
  FESystem<2>        q2_q3(FE_Bernstein<2>(2), 1, FE_Bernstein<2>(3), 1);
  //  plot_derivatives(m, q2_q3, "B2_Q3");

  return 0;
}

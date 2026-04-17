// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2006 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Like rt_13, but use FESubfaceValues
//
// the test used to fail because of the issue with computing the
// normals using FEFaceValue, where FEFaceValue by accident uses the
// *face* mapping, not the *cell* mapping to compute the Piola
// transform (leading to a missing power of h in the determinant)

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include <string>
#include <vector>

#include "../tests.h"

#define PRECISION 8


template <int dim>
void
test(const unsigned int degree)
{
  FE_RaviartThomas<dim> fe_rt(degree);

  deallog << "Degree=" << degree << std::endl;

  for (double h = 1; h > 1. / 128; h /= 2)
    {
      deallog << "  h=" << h << std::endl;

      Triangulation<dim> tr;

      // Construct refinement situation where the first active cell is
      // adjacent to a refined cell with face 0
      Point<dim> l, r;
      l[0] = -h;
      for (unsigned int d = 0; d < dim; ++d)
        r[d] = h;
      std::vector<unsigned int> refine(dim, 1);
      refine[0] = 2;
      GridGenerator::subdivided_hyper_rectangle(tr, refine, l, r);

      // Create refinement setting by refining second to last element one more
      // time
      tr.begin_active()->set_refine_flag();
      tr.execute_coarsening_and_refinement();

      DoFHandler<dim> dof(tr);
      dof.distribute_dofs(fe_rt);

      QTrapezoid<dim - 1> quadrature;

      FESubfaceValues<dim> fe_values(fe_rt, quadrature, update_gradients);
      fe_values.reinit(tr.begin_active(), 0, 0);
      for (unsigned int q = 0; q < quadrature.size(); ++q)
        {
          deallog << "    Quadrature point " << q << ": ";
          for (unsigned int i = 0; i < fe_rt.dofs_per_cell; ++i)
            {
              deallog << "[ ";
              for (unsigned int c = 0; c < fe_rt.n_components(); ++c)
                deallog << fe_values.shape_grad_component(i, q, c) << ' ';
              deallog << ']';
            }
          deallog << std::endl;
        }
    }
}


int
main()
{
  initlog();
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;

  for (unsigned int i = 0; i < 4; ++i)
    test<2>(i);

  return 0;
}

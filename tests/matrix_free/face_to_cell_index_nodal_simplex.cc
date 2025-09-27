// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test the correctness of face_to_cell_index_nodal
// in internal::MatrixFreeFunctions::ShapeInfo
// for simplex elements

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/shape_info.h>

#include <iostream>

#include "../tests.h"



template <int dim>
void
test(const FiniteElement<dim> &fe, const Quadrature<dim> &quad)
{
  for (unsigned int i = 0; i < fe.n_base_elements(); ++i)
    {
      internal::MatrixFreeFunctions::ShapeInfo<double> shape_info;
      shape_info.reinit(quad, fe, i);
      deallog << "Testing: " << fe.get_name() << std::endl;

      Triangulation<dim> tria;
      GridGenerator::reference_cell(tria, fe.reference_cell());

      FE_SimplexP<dim> fe_continuous(fe.degree);

      DoFHandler<dim> dof_handler(tria);
      dof_handler.distribute_dofs(fe_continuous);

      std::vector<types::global_dof_index> dof_indices(
        fe_continuous.n_dofs_per_face());
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          unsigned int face_counter = 0;
          for (const auto &face : cell->face_iterators())
            {
              face->get_dof_indices(dof_indices);

              unsigned int dof_counter = 0;
              for (const auto i : dof_indices)
                {
                  if (i == shape_info.face_to_cell_index_nodal[face_counter]
                                                              [dof_counter])
                    deallog << "Ok ";
                  else
                    deallog << "Not ok";
                  ++dof_counter;
                }
              deallog << std::endl;
              ++face_counter;
            }
        }
    }
}


int
main()
{
  initlog();
  for (unsigned int degree = 1; degree < 4; ++degree)
    {
      FE_SimplexP<2>   fe(degree);
      FE_SimplexDGP<2> fe_dg(degree);
      QGaussSimplex<2> quad(1);
      test(fe, quad);
      test(fe_dg, quad);


      FE_SimplexP<3>   fe_3(degree);
      FE_SimplexDGP<3> fe_dg_3(degree);
      QGaussSimplex<3> quad_3(1);
      test(fe_3, quad_3);
      test(fe_dg_3, quad_3);
    }
}

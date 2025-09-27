// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test FESubfaceValues for different subface quadrature rules on triangles.


#include <deal.II/base/function_lib.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/reference_cell.h>

#include <deal.II/hp/q_collection.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

#include "./simplex_grids.h"

template <int dim>
void
test()
{}

template <>
void
test<2>()
{
  const unsigned int dim = 2;

  // test FESubFaceValues for FE_SimplexP
  {
    deallog.push("FESubFaceValues");
    const hp::QCollection<dim - 1> quad_ref(QGaussSimplex<dim - 1>(1),
                                            QGaussSimplex<dim - 1>(2),
                                            QGaussSimplex<dim - 1>(3),
                                            QGaussSimplex<dim - 1>(4));

    MappingFE<dim>   mapping(FE_SimplexP<dim>(1));
    FE_SimplexP<dim> fe(3);

    const UpdateFlags flags =
      update_values | update_quadrature_points | update_JxW_values;

    Triangulation<dim> tria;
    GridGenerator::subdivided_hyper_cube_with_simplices(tria, 1);
    tria.begin_active()->set_refine_flag(
      RefinementCase<dim>::isotropic_refinement);
    tria.execute_coarsening_and_refinement();

    DoFHandler<dim> dof_handler(tria);
    dof_handler.distribute_dofs(fe);

    Vector<double> vector_0(dof_handler.n_dofs());
    Vector<double> vector_1(dof_handler.n_dofs());

    VectorTools::interpolate(mapping,
                             dof_handler,
                             Functions::Monomial<dim, double>(
                               Tensor<1, dim, double>({1., 0.})),
                             vector_0);
    VectorTools::interpolate(mapping,
                             dof_handler,
                             Functions::Monomial<dim, double>(
                               Tensor<1, dim, double>({0., 1.})),
                             vector_1);

    std::vector<double> values_0, values_1;

    for (unsigned int current_quad = 0; current_quad < quad_ref.size();
         ++current_quad)
      {
        FESubfaceValues<dim> fe_face_values(mapping,
                                            fe,
                                            quad_ref[current_quad],
                                            flags);

        for (const auto &cell : dof_handler.active_cell_iterators())
          for (const auto face_no : cell->face_indices())
            for (unsigned int subface_no = 0; subface_no < 2; ++subface_no)
              {
                fe_face_values.reinit(cell, face_no, subface_no);

                values_0.resize(fe_face_values.n_quadrature_points);
                values_1.resize(fe_face_values.n_quadrature_points);

                fe_face_values.get_function_values(vector_0, values_0);
                fe_face_values.get_function_values(vector_1, values_1);

                for (unsigned int q = 0; q < values_0.size(); ++q)
                  {
                    if (std::abs(values_0[q] -
                                 fe_face_values.quadrature_point(q)[0]) < 1e-12)
                      deallog << "Ok"
                              << " ";
                    else
                      deallog << "False"
                              << " ";

                    if (std::abs(values_1[q] -
                                 fe_face_values.quadrature_point(q)[1]) < 1e-12)
                      deallog << "Ok"
                              << " ";
                    else
                      deallog << "False"
                              << " ";
                    deallog << std::endl;
                  }
              }
      }
    deallog << std::endl;
  }
}


int
main()
{
  initlog();

  deallog.push("dim=2");
  test<2>();
  deallog.pop();
}

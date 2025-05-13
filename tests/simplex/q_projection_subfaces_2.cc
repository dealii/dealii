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


// Test FESubfaceValues for different subface quadrature rules on a mixed mesh.


#include <deal.II/base/function_lib.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/reference_cell.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

#include "./simplex_grids.h"

template <int dim>
void
test(const unsigned int n_q)
{}

template <>
void
test<2>(const unsigned int n_q)
{
  const unsigned int dim = 2;

  // test FESubFaceValues for mixed meshes
  {
    deallog.push("FESubFaceValues");
    QGaussSimplex<dim - 1>         quad_simplex(n_q);
    QGauss<dim - 1>                quad_q(n_q);
    const hp::QCollection<dim - 1> quad_ref(quad_simplex, quad_q);

    hp::MappingCollection<dim> mapping(MappingFE<dim>(FE_SimplexP<dim>(1)),
                                       MappingQ1<dim>());
    hp::FECollection<dim>      fe(FE_SimplexP<dim>(3), FE_Q<dim>(3));

    const UpdateFlags flags =
      update_values | update_quadrature_points | update_JxW_values;

    Triangulation<dim> tria;
    {
      std::vector<Point<dim>>    vertices;
      std::vector<CellData<dim>> cells;
      // determine cell sizes
      const Point<dim> dx(0.25, 0.25);

      // create vertices
      for (unsigned int j = 0; j <= 4; ++j)
        for (unsigned int i = 0; i <= 4; ++i)
          vertices.push_back(Point<dim>(dx[0] * i, dx[1] * j));

      // create cells
      for (unsigned int j = 0; j < 4; ++j)
        for (unsigned int i = 0; i < 4; ++i)
          {
            // create reference QUAD cell
            std::array<unsigned int, 4> quad{{
              (j + 0) * (4 + 1) + i + 0, //
              (j + 0) * (4 + 1) + i + 1, //
              (j + 1) * (4 + 1) + i + 0, //
              (j + 1) * (4 + 1) + i + 1  //
            }};                          //

            if (j < 4 / 2 && i < 4 / 2)
              {
                CellData<dim> quad_;
                quad_.vertices = {quad[0], quad[3], quad[1], quad[2]};
                cells.push_back(quad_);

                continue;
              }

            // TRI cell 0
            {
              CellData<dim> tri;
              tri.vertices = {quad[0], quad[1], quad[2]};
              cells.push_back(tri);
            }

            // TRI cell 1
            {
              CellData<dim> tri;
              tri.vertices = {quad[3], quad[2], quad[1]};
              cells.push_back(tri);
            }
          }
      tria.create_triangulation(vertices, cells, SubCellData());
    }

    DoFHandler<dim> dof_handler(tria);
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->reference_cell() == ReferenceCells::Triangle)
          cell->set_active_fe_index(0);
        else
          cell->set_active_fe_index(1);
      }
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

    std::vector<double>      values_0, values_1;
    hp::FESubfaceValues<dim> hp_fe_subface_values(mapping, fe, quad_ref, flags);

    for (const auto &cell : dof_handler.active_cell_iterators())
      for (const auto face_no : cell->face_indices())
        for (unsigned int subface_no = 0; subface_no < 2; ++subface_no)
          {
            hp_fe_subface_values.reinit(cell, face_no, subface_no);
            const FESubfaceValues<dim> &fe_face_values =
              hp_fe_subface_values.get_present_fe_values();

            values_0.resize(
              hp_fe_subface_values.get_quadrature_collection()[0].size());
            values_1.resize(
              hp_fe_subface_values.get_quadrature_collection()[1].size());

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
    deallog << std::endl;
    deallog.pop();
  }
}


int
main()
{
  initlog();

  deallog.push("dim=2");
  for (unsigned int q = 1; q < 5; ++q)
    test<2>(q);
  deallog.pop();
}

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


// Test FEFaceValues for different face quadrature rules.


#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/reference_cell.h>

#include <deal.II/hp/q_collection.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
class Fu : public Function<dim>
{
public:
  Fu(const unsigned int component = numbers::invalid_unsigned_int)
    : Function<dim>(component == numbers::invalid_unsigned_int ? dim : 1)
    , component(component)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int i = 0) const
  {
    if (component == numbers::invalid_unsigned_int)
      return p[i];
    else
      return p[component];
  }

  const unsigned int component;
};

template <int dim>
void
test()
{}

template <>
void
test<2>()
{
  const unsigned int dim = 2;

  // test: QProjector::project_to_all_faces
  {
    const hp::QCollection<dim - 1> quad_ref(QGauss<dim - 1>(1),
                                            QGauss<dim - 1>(2),
                                            QGauss<dim - 1>(3),
                                            QGauss<dim - 1>(4));

    const auto quad =
      QProjector<dim>::project_to_all_faces(ReferenceCells::Quadrilateral,
                                            quad_ref);

    const auto print = [&](const unsigned int face_no) {
      deallog << "face_no=" << face_no << ':' << std::endl;
      for (unsigned int q = 0,
                        i = QProjector<dim>::DataSetDescriptor::face(
                          ReferenceCells::Quadrilateral,
                          face_no,
                          numbers::default_geometric_orientation,
                          quad_ref);
           q < quad_ref[face_no].size();
           ++q, ++i)
        {
          deallog << quad.point(i) << ' ';
          deallog << quad.weight(i) << ' ';
          deallog << std::endl;
        }
      deallog << std::endl;
    };

    for (unsigned int i = 0; i < 4 /*TODO*/; ++i)
      print(i);

    deallog << std::endl;
  }

  // test: Mapping
  {
    const hp::QCollection<dim - 1> quad_ref(QGauss<dim - 1>(1),
                                            QGauss<dim - 1>(2),
                                            QGauss<dim - 1>(3),
                                            QGauss<dim - 1>(4));

    MappingFE<dim> mapping(FE_Q<dim>(1));
    FE_Q<dim>      fe(3);

    const UpdateFlags flags = mapping.requires_update_flags(
      update_values | update_quadrature_points | update_JxW_values);

    auto data_ref = mapping.get_face_data(flags, quad_ref);

    internal::FEValuesImplementation::MappingRelatedData<dim> data;
    data.initialize(quad_ref.max_n_quadrature_points(), flags);

    Triangulation<dim> tria;
    GridGenerator::hyper_cube(tria);

    for (const auto &cell : tria.active_cell_iterators())
      for (const auto face_no : cell->face_indices())
        {
          mapping.fill_fe_face_values(cell, face_no, quad_ref, *data_ref, data);

          deallog << "face_no=" << face_no << ':' << std::endl;
          for (unsigned int q = 0; q < quad_ref[face_no].size(); ++q)
            {
              deallog << data.quadrature_points[q] << ' ';
              deallog << data.JxW_values[q] << ' ';
              deallog << std::endl;
            }
          deallog << std::endl;
        }

#if false
  // not possible to test since functions are protected
  internal::FEValuesImplementation::FiniteElementRelatedData<dim> data_fe;
  data_fe.initialize(quad_ref.max_n_quadrature_points(), fe, flags);

  auto data_fe_ref =fe.get_face_data(flags, mapping, quad_ref, data_fe);

  for (const auto &cell : tria.active_cell_iterators())
    for (const auto face_no : cell->face_indices())
      {
        fe.fill_fe_face_values(cell, face_no, quad_ref, mapping, *data_ref, data, *data_fe_ref, data_fe);
      }
#endif
  }

  // test FEFaceValues for FE_Q
  {
    const hp::QCollection<dim - 1> quad_ref(QGauss<dim - 1>(1),
                                            QGauss<dim - 1>(2),
                                            QGauss<dim - 1>(3),
                                            QGauss<dim - 1>(4));

    MappingFE<dim> mapping(FE_Q<dim>(1));
    FE_Q<dim>      fe(3);

    const UpdateFlags flags = mapping.requires_update_flags(
      update_values | update_quadrature_points | update_JxW_values);


    FEFaceValues<dim> fe_face_values(mapping, fe, quad_ref, flags);

    Triangulation<dim> tria;
    GridGenerator::hyper_cube(tria);

    DoFHandler<dim> dof_handler(tria);
    dof_handler.distribute_dofs(fe);

    Vector<double> vector_0(dof_handler.n_dofs());
    Vector<double> vector_1(dof_handler.n_dofs());

    VectorTools::interpolate(mapping, dof_handler, Fu<dim>(0), vector_0);
    VectorTools::interpolate(mapping, dof_handler, Fu<dim>(1), vector_1);

    std::vector<double> values_0, values_1;

    values_0.reserve(fe_face_values.n_quadrature_points);
    values_1.reserve(fe_face_values.n_quadrature_points);

    for (const auto &cell : dof_handler.active_cell_iterators())
      for (const auto face_no : cell->face_indices())
        {
          fe_face_values.reinit(cell, face_no);

          values_0.resize(fe_face_values.n_quadrature_points);
          values_1.resize(fe_face_values.n_quadrature_points);

          fe_face_values.get_function_values(vector_0, values_0);
          fe_face_values.get_function_values(vector_1, values_1);

          for (unsigned int q = 0; q < values_0.size(); ++q)
            deallog << values_0[q] << ' ' << values_1[q] << ' ' << std::endl;

          deallog << std::endl;
        }
  }

  // test FEFaceValues for FE_System(FE_Q)
  {
    const hp::QCollection<dim - 1> quad_ref(QGauss<dim - 1>(1),
                                            QGauss<dim - 1>(2),
                                            QGauss<dim - 1>(3),
                                            QGauss<dim - 1>(4));

    MappingFE<dim> mapping(FE_Q<dim>(1));
    FESystem<dim>  fe(FE_Q<dim>{3}, dim);

    const UpdateFlags flags = mapping.requires_update_flags(
      update_values | update_quadrature_points | update_JxW_values);


    FEFaceValues<dim> fe_face_values(mapping, fe, quad_ref, flags);

    Triangulation<dim> tria;
    GridGenerator::hyper_cube(tria);

    DoFHandler<dim> dof_handler(tria);
    dof_handler.distribute_dofs(fe);

    Vector<double> vector_0(dof_handler.n_dofs());

    VectorTools::interpolate(mapping, dof_handler, Fu<dim>(), vector_0);

    for (const auto &cell : dof_handler.active_cell_iterators())
      for (const auto face_no : cell->face_indices())
        {
          fe_face_values.reinit(cell, face_no);

          std::vector<Vector<double>> values_0(quad_ref[face_no].size(),
                                               Vector<double>(dim));
          fe_face_values.get_function_values(vector_0, values_0);

          deallog << "face_no=" << face_no << ':' << std::endl;

          for (unsigned int q = 0; q < values_0.size(); ++q)
            deallog << values_0[q][0] << ' ' << values_0[q][1] << ' '
                    << std::endl;

          deallog << std::endl;
        }
  }
}

template <>
void
test<3>()
{
  const unsigned int dim = 3;

  // test FEFaceValues for FE_System(FE_Q)
  {
    const hp::QCollection<dim - 1> quad_ref(QGauss<dim - 1>(1),
                                            QGauss<dim - 1>(2),
                                            QGauss<dim - 1>(3),
                                            QGauss<dim - 1>(4),
                                            QGauss<dim - 1>(1),
                                            QGauss<dim - 1>(2));

    MappingFE<dim> mapping(FE_Q<dim>(1));
    FESystem<dim>  fe(FE_Q<dim>{3}, dim);

    const UpdateFlags flags =
      update_values | update_quadrature_points | update_JxW_values;


    FEFaceValues<dim> fe_face_values(mapping, fe, quad_ref, flags);

    Triangulation<dim> tria;
    GridGenerator::hyper_cube(tria);

    DoFHandler<dim> dof_handler(tria);
    dof_handler.distribute_dofs(fe);

    Vector<double> vector_0(dof_handler.n_dofs());

    VectorTools::interpolate(mapping, dof_handler, Fu<dim>(), vector_0);

    for (const auto &cell : dof_handler.active_cell_iterators())
      for (const auto face_no : cell->face_indices())
        {
          fe_face_values.reinit(cell, face_no);

          std::vector<Vector<double>> values_0(quad_ref[face_no].size(),
                                               Vector<double>(dim));
          fe_face_values.get_function_values(vector_0, values_0);

          deallog << "face_no=" << face_no << ':' << std::endl;

          for (unsigned int q = 0; q < values_0.size(); ++q)
            deallog << values_0[q][0] << ' ' << values_0[q][1] << ' '
                    << values_0[q][2] << ' ' << std::endl;

          deallog << std::endl;
        }
  }
}

int
main()
{
  initlog();

  deallog.push("dim=2");
  test<2>();
  deallog.pop();

  deallog.push("dim=3");
  test<3>();
  deallog.pop();
}

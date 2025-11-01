// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_tests_fe_conformity_test_h
#define dealii_tests_fe_conformity_test_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_data.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>

// STL
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

// My test headers
#include "fe_conformity_fill_vector_random.h"

namespace FEConforimityTest
{
  using namespace dealii;

  template <int dim>
  class FEConformityTest
  {
  public:
    FEConformityTest(const FiniteElement<dim> &fe,
                     const unsigned int        config_switch);

    void
    run();

  private:
    void
    make_grid();

    void
    make_dofs();

    void
    run_test();

    void
    get_function_jump(const FEInterfaceValues<dim> &fe_interface_values,
                      const Vector<double>         &dof_vector,
                      std::vector<double>          &jumps);

    void
    get_normal_jump(const FEInterfaceValues<dim> &fe_interface_values,
                    const Vector<double>         &dof_vector,
                    std::vector<double>          &jumps);

    void
    get_tangential_jump(const FEInterfaceValues<dim> &fe_interface_values,
                        const Vector<double>         &dof_vector,
                        std::vector<double>          &jumps);

    ObserverPointer<const FiniteElement<dim>> fe_ptr;
    Triangulation<dim>                        triangulation;
    DoFHandler<dim>                           dof_handler;
    AffineConstraints<double>                 constraints;

    Vector<double> random_fe_function;

    const unsigned int config_switch;

    const bool plot_mesh_properties = false;
  };

  template <int dim>
  FEConformityTest<dim>::FEConformityTest(const FiniteElement<dim> &_fe,
                                          const unsigned int config_switch)
    : fe_ptr(&_fe)
    , dof_handler(triangulation)
    , config_switch(config_switch)
  {
    deallog
      << "Element Info:  " << std::endl
      << "   name                 : " << fe_ptr->get_name() << std::endl
      << "   is_primitive         : " << fe_ptr->is_primitive() << std::endl
      << "   n_dofs_per_cell      : " << fe_ptr->dofs_per_cell << std::endl
      << "   n_dofs_per_face      : " << fe_ptr->dofs_per_face << std::endl
      << "   n_dofs_per_quad      : " << fe_ptr->dofs_per_quad << std::endl
      << "   n_dofs_per_line      : " << fe_ptr->dofs_per_line << std::endl
      << "   n_dofs_per_vertex    : " << fe_ptr->dofs_per_vertex << std::endl
      << std::endl
      << "   first_line_index     : " << fe_ptr->first_line_index << std::endl
      << "   first_quad_index     : " << fe_ptr->first_quad_index << std::endl
      << "   first_face_line_index: " << fe_ptr->first_face_line_index
      << std::endl
      << "   first_face_quad_index: " << fe_ptr->first_face_quad_index
      << std::endl
      << std::endl
      << "   n_components         : " << fe_ptr->n_components() << std::endl
      << "   n_blocks             : " << fe_ptr->n_blocks() << std::endl
      << "   n_base_elements      : " << fe_ptr->n_base_elements() << std::endl
      << std::endl;

    if (dim == 2)
      AssertThrow(config_switch < 4,
                  ExcMessage("If dim=2 the config witch must be < 4."));

    if (dim == 3)
      AssertThrow(config_switch < 8,
                  ExcMessage("If dim=3 the config witch must be < 8."));
  }



  template <>
  void
  FEConformityTest<3>::make_grid()
  {
    triangulation.clear();

    bool face_orientation = (((config_switch / 4) % 2) == 1);
    bool face_flip        = (((config_switch / 2) % 2) == 1);
    bool face_rotation    = ((config_switch % 2) == 1);

    bool manipulate_first_cube = true;

    GridGenerator::non_standard_orientation_mesh(triangulation,
                                                 face_orientation,
                                                 face_flip,
                                                 face_rotation,
                                                 manipulate_first_cube);

    //    GridTools::distort_random(/* factor */ 0.15,
    //                              triangulation_coarse,
    //                              /* keep_boundary */ false);

    triangulation.refine_global(0);

    if (plot_mesh_properties)
      {
        deallog << std::endl
                << "3D Mesh properties: " << std::endl
                << std::endl;
        for (const auto &cell : triangulation.active_cell_iterators())
          {
            CellId current_cell_id(cell->id());

            deallog
              << "CellId = " << current_cell_id << std::endl
              << "   {face_index -> face_orientation | face_flip | face_rotation}: "
              << std::endl;
            for (unsigned int face_index = 0;
                 face_index < GeometryInfo<3>::faces_per_cell;
                 ++face_index)
              {
                deallog << "      {" << face_index << " -> "
                        << cell->face_orientation(face_index) << " | "
                        << cell->face_flip(face_index) << " | "
                        << cell->face_rotation(face_index) << "}" << std::endl;
              } // face_index

            deallog << "   {line_index -> line_orientation}: " << std::endl;
            for (unsigned int line_index = 0;
                 line_index < GeometryInfo<3>::lines_per_cell;
                 ++line_index)
              {
                deallog << "      {" << line_index << " -> "
                        << (cell->line_orientation(line_index) ==
                            numbers::default_geometric_orientation)
                        << "}" << std::endl;
              } // line_index
          }     // cell
      }         // plot_mesh_properties
  }             // make_grid<3>()


  template <>
  void
  FEConformityTest<2>::make_grid()
  {
    triangulation.clear();

    // alias for better readability
    const unsigned int n_rotate_central_square = config_switch;

    GridGenerator::non_standard_orientation_mesh(triangulation,
                                                 n_rotate_central_square);

    //    GridTools::distort_random(/* factor */ 0.15,
    //                              triangulation_coarse,
    //                              /* keep_boundary */ false);

    triangulation.refine_global(0);

    if (plot_mesh_properties)
      {
        deallog << std::endl
                << "2D Mesh properties: " << std::endl
                << std::endl;
        for (const auto &cell : triangulation.active_cell_iterators())
          {
            CellId current_cell_id(cell->id());

            deallog << "CellId = " << current_cell_id << std::endl
                    << "   {face_index -> face_orientation}: " << std::endl;
            for (unsigned int face_index = 0;
                 face_index < GeometryInfo<2>::faces_per_cell;
                 ++face_index)
              {
                deallog << "      {" << face_index << " -> "
                        << cell->face_orientation(face_index) << "}"
                        << std::endl;
              } // face_index

            deallog << "   {line_index -> line_orientation}: " << std::endl;
            for (unsigned int line_index = 0;
                 line_index < GeometryInfo<2>::lines_per_cell;
                 ++line_index)
              {
                deallog << "      {" << line_index << " -> "
                        << (cell->line_orientation(line_index) ==
                            numbers::default_geometric_orientation)
                        << "}" << std::endl;
              } // line_index
          }     // cell
      }         // plot_mesh_properties
  }             // make_grid<2>()


  template <int dim>
  void
  FEConformityTest<dim>::make_dofs()
  {
    dof_handler.clear();
    constraints.clear();

    dof_handler.distribute_dofs(*fe_ptr);

    { // Constraints
      constraints.clear();
      DoFTools::make_hanging_node_constraints(dof_handler, constraints);
      constraints.close();
    }

    random_fe_function.reinit(dof_handler.n_dofs());

    // Fill vector with pseudo-random values
    fill_vector_randomly(random_fe_function, /* min */ -5.0, /* max */ 5.0);
  } // make_dofs()



  template <int dim>
  void
  FEConformityTest<dim>::get_function_jump(
    const FEInterfaceValues<dim> &fe_interface_values,
    const Vector<double>         &dof_vector,
    std::vector<double>          &jumps)
  {
    const unsigned n_q = fe_interface_values.n_quadrature_points;

    Assert(jumps.size() == n_q,
           ExcMessage(
             "Length of \"jumps\" and number of quad points must coincide"));

    std::array<std::vector<double>, 2> face_values;
    for (unsigned i = 0; i < 2; ++i)
      {
        face_values[i].resize(n_q);
        fe_interface_values.get_fe_face_values(i).get_function_values(
          dof_vector, face_values[i]);
      }

    for (unsigned int q = 0; q < n_q; ++q)
      jumps[q] = face_values[0][q] - face_values[1][q];
  } // get_function_jump()



  template <int dim>
  void
  FEConformityTest<dim>::get_normal_jump(
    const FEInterfaceValues<dim> &fe_interface_values,
    const Vector<double>         &dof_vector,
    std::vector<double>          &jumps)
  {
    const unsigned n_q = fe_interface_values.n_quadrature_points;

    Assert(jumps.size() == n_q,
           ExcMessage(
             "Length of \"jumps\" and number of quad points must coincide"));

    const FEValuesExtractors::Vector normal_flux(0);

    std::array<std::vector<Tensor<1, dim>>, 2> face_values;
    for (unsigned i = 0; i < 2; ++i)
      {
        face_values[i].resize(n_q);

        const auto &fe_face_values = fe_interface_values.get_fe_face_values(i);

        fe_face_values[normal_flux].get_function_values(dof_vector,
                                                        face_values[i]);
      }

    const std::vector<Tensor<1, dim>> interface_normals =
      fe_interface_values.get_normal_vectors();

    for (unsigned int q = 0; q < n_q; ++q)
      {
        jumps[q] =
          (face_values[0][q] - face_values[1][q]) * interface_normals[q];
      }
  } // get_normal_jump()



  template <>
  void
  FEConformityTest<2>::get_tangential_jump(
    const FEInterfaceValues<2> &fe_interface_values,
    const Vector<double>       &dof_vector,
    std::vector<double>        &jumps)
  {
    const unsigned n_q = fe_interface_values.n_quadrature_points;

    Assert(jumps.size() == n_q,
           ExcMessage(
             "Length of \"jumps\" and number of quad points must coincide"));

    const FEValuesExtractors::Vector tangent_flux(0);

    std::array<std::vector<Tensor<1, 2>>, 2> face_values;
    for (unsigned i = 0; i < 2; ++i)
      {
        face_values[i].resize(n_q);

        const auto &fe_face_values = fe_interface_values.get_fe_face_values(i);

        fe_face_values[tangent_flux].get_function_values(dof_vector,
                                                         face_values[i]);
      }

    const std::vector<Tensor<1, 2>> interface_normals =
      fe_interface_values.get_normal_vectors();

    for (unsigned int q = 0; q < n_q; ++q)
      jumps[q] = cross_product_2d(face_values[0][q] - face_values[1][q]) *
                 interface_normals[q];
  } // get_tangential_jump()



  template <>
  void
  FEConformityTest<3>::get_tangential_jump(
    const FEInterfaceValues<3> &fe_interface_values,
    const Vector<double>       &dof_vector,
    std::vector<double>        &jumps)
  {
    const unsigned n_q = fe_interface_values.n_quadrature_points;

    Assert(jumps.size() == n_q,
           ExcMessage(
             "Length of \"jumps\" and number of quad points must coincide"));

    const FEValuesExtractors::Vector tangent_flux(0);

    std::array<std::vector<Tensor<1, 3>>, 2> face_values;
    for (unsigned i = 0; i < 2; ++i)
      {
        face_values[i].resize(n_q);

        const auto &fe_face_values = fe_interface_values.get_fe_face_values(i);

        fe_face_values[tangent_flux].get_function_values(dof_vector,
                                                         face_values[i]);
      }

    const std::vector<Tensor<1, 3>> interface_normals =
      fe_interface_values.get_normal_vectors();

    for (unsigned int q = 0; q < n_q; ++q)
      jumps[q] = cross_product_3d(face_values[0][q] - face_values[1][q],
                                  interface_normals[q])
                   .norm();
  } // get_tangential_jump()



  template <int dim>
  void
  FEConformityTest<dim>::run_test()
  {
    QGauss<dim - 1>    interface_quad_rule(fe_ptr->degree + 2);
    const unsigned int n_interface_q_points = interface_quad_rule.size();

    const UpdateFlags interface_update_flags =
      (update_values | update_quadrature_points | update_JxW_values |
       update_normal_vectors);

    FEInterfaceValues<dim> fe_interface_values(*fe_ptr,
                                               interface_quad_rule,
                                               interface_update_flags);

    std::vector<double> interface_jumps(n_interface_q_points);

    // Loop over all (two) cells and faces. If at inner face check for
    // appropriate jumps.
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        for (unsigned int face_index = 0;
             face_index < GeometryInfo<dim>::faces_per_cell;
             ++face_index)
          {
            if (cell->face(face_index)->at_boundary())
              continue;

            fe_interface_values.reinit(
              cell,
              face_index,
              /* sub_face_no */ numbers::invalid_unsigned_int,
              /* cell_neighbor */ cell->neighbor(face_index),
              /* face_no_neighbor */ cell->neighbor_face_no(face_index),
              /* sub_face_no_neighbor */ numbers::invalid_unsigned_int);

            switch (fe_ptr->conforming_space)
              {
                case FiniteElementData<dim>::Conformity::H1:
                  {
                    get_function_jump(fe_interface_values,
                                      random_fe_function,
                                      interface_jumps);

                    deallog << "Function jumps (at quad points) in cell   "
                            << cell->id().to_string() << "   at face   "
                            << face_index << std::endl;

                    break;
                  } // H1

                case FiniteElementData<dim>::Conformity::Hdiv:
                  {
                    get_normal_jump(fe_interface_values,
                                    random_fe_function,
                                    interface_jumps);

                    deallog << "Normal jumps (at quad points) in cell   "
                            << cell->id().to_string() << "   at face   "
                            << face_index << std::endl;

                    break;
                  } // case H(div)

                case FiniteElementData<dim>::Conformity::Hcurl:
                  {
                    get_tangential_jump(fe_interface_values,
                                        random_fe_function,
                                        interface_jumps);

                    deallog << "Tangential jumps (at quad points) in cell   "
                            << cell->id().to_string() << "   at face   "
                            << face_index << std::endl;
                    break;
                  } // case H(curl)

                default:
                  break;
              } // switch (fe_ptr->conforming_space)

            // loop over quad_points and print jumps
            deallog << "   interface jumps:   ";
            for (unsigned int q = 0; q < n_interface_q_points; ++q)
              {
                deallog << interface_jumps[q] << "   ";
              }
            deallog << std::endl;
          } // face_index
      }     // cell
  }         // run_test()



  template <int dim>
  void
  FEConformityTest<dim>::run()
  {
    make_grid();
    make_dofs();
    run_test();
  }
} // namespace FEConforimityTest

#endif

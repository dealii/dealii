// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Output tensor-valued data in vtu

#include <deal.II/base/logstream.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iomanip>
#include <string>

#include "../tests.h"



template <int dim>
void
check()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);

  FESystem<dim> fe(FE_Q<dim>(1),
                   dim, // vector valued
                   FE_Q<dim>(1),
                   dim * dim); // tensor valued

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  Vector<double> v(dof_handler.n_dofs());


  // 1 cell is enough to test:
  const unsigned int dofs_per_cell = v.size();

  // take 8 tensors from tensors.vtk file inside
  // http://paraview.org/Wiki/images/8/81/TensorGlyph_v_1_0_3.zip or
  // http://paraview.org/Wiki/images/3/39/SuperquadricTensorGlyph_v3.zip
  const std::vector<std::vector<std::vector<double>>> tensors_3d = {
    {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
    {{-2, 0, 0}, {0, -2, 0}, {0, 0, 3}},
    {{1, 0, 0}, {0, 1, 0}, {0, 0, .001}},
    {{2, 1, 1}, {1, 2, 1}, {1, 1, 3}},
    {{0.247186415568, 0.490995206139, 0.131324884836},
     {0.490995206139, -0.371055707211, 0.719071682671},
     {0.131324884836, 0.719071682671, -0.156008182087}},
    {{0.280587657181, 0.467438945439, 0.934953136331},
     {0.467438945439, 0.0600321140579, -0.211376327727},
     {0.934953136331, -0.211376327727, 0.962975830149}},
    {{0.0628609549443, -0.0117908465212, -0.667617347633},
     {-0.0117908465212, -0.390601011028, -0.95780533241},
     {-0.667617347633, -0.95780533241, -0.193319773383}},
    {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};

  // coords:
  const std::vector<std::vector<double>> coords_3d = {{0, 0, 0},
                                                      {1, 0, 0},
                                                      {0, 1, 0},
                                                      {1, 1, 0},
                                                      {0, 0, 1},
                                                      {1, 0, 1},
                                                      {0, 1, 1},
                                                      {1, 1, 1}};

  const std::vector<std::vector<std::vector<double>>> tensors_2d = {
    {{1, 0}, {0, 1}}, {{3, 0}, {0, 1}}, {{1, 0}, {0, 3}}, {{3, 1}, {1, 2}}};

  // coords:
  const std::vector<std::vector<double>> coords_2d = {{0, 0},
                                                      {1, 0},
                                                      {0, 1},
                                                      {1, 1}};

  const std::vector<std::vector<std::vector<double>>> &tensors =
    (dim == 2 ? tensors_2d : tensors_3d);
  const std::vector<std::vector<double>> &coords =
    (dim == 2 ? coords_2d : coords_3d);

  const Quadrature<dim> support_quadrature(fe.get_unit_support_points());
  FEValues<dim> fe_values(fe, support_quadrature, update_quadrature_points);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      const auto &qp = fe_values.get_quadrature_points();

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          const auto        &this_qp = qp[i];
          const unsigned int i_group = fe.system_to_base_index(i).first.first;
          const unsigned int i_comp  = fe.system_to_component_index(i).first;
          if (i_group == 0)
            // vector
            {
              if (i_comp == 0)
                v(i) = 1;
              else if (i_comp == 1)
                v(i) = this_qp[1];
              else
                v(i) = 0.;
            }
          else
            // tensor
            {
              // find matching location
              unsigned int p = 0;
              for (const auto &c : coords)
                {
                  Point<dim> point;
                  for (unsigned int d = 0; d < dim; ++d)
                    point[d] = c[d];

                  const auto diff = point - this_qp;
                  if (diff.norm() < 1e-10)
                    break;

                  ++p;
                }

              const unsigned int tens_comp = i_comp - dim;
              const unsigned int ii =
                Tensor<2, dim>::unrolled_to_component_indices(tens_comp)[0];
              const unsigned int jj =
                Tensor<2, dim>::unrolled_to_component_indices(tens_comp)[1];
              v(i) = tensors[p][ii][jj];
            }
        }
    }

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      dim + dim * dim,
      DataComponentInterpretation::component_is_part_of_tensor);
  std::fill(data_component_interpretation.begin(),
            data_component_interpretation.begin() + dim,
            DataComponentInterpretation::component_is_part_of_vector);

  std::vector<std::string> component_name(dim + dim * dim, "tensor");
  std::fill(component_name.begin(), component_name.begin() + dim, "vector");

  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.compression_level = DataOutBase::CompressionLevel::best_compression;

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(v,
                           component_name,
                           DataOut<dim>::type_dof_data,
                           data_component_interpretation);
  data_out.build_patches();
  data_out.set_flags(vtk_flags);

  /*
  std::ofstream out("output_" + Utilities::int_to_string(dim) + "d.vtu");
  data_out.write_vtu(out);
  */

  data_out.write_vtu(deallog.get_file_stream());
}



int
main()
{
  initlog();
  check<2>();
  check<3>();
}

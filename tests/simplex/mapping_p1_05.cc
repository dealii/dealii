// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
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

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_p1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools_topology.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


// Test MappingP1::transform() by comparing with MappingFE.


// To call Mapping::transform() we have to get the mapping data, which requires
// accessing a protected member of FEValuesBase:
template <int dim, int spacedim>
class FEValues2 : public FEValues<dim, spacedim>
{
public:
  using FEValues<dim, spacedim>::FEValues;

  const typename Mapping<dim, spacedim>::InternalDataBase &
  get_mapping_data() const
  {
    Assert(this->mapping_data, ExcInternalError());
    return *this->mapping_data;
  }
};

using namespace dealii;

void
fill_tensor(double &t)
{
  t = random_value<double>();
}

template <int rank, int dim>
void
fill_tensor(Tensor<rank, dim> &t)
{
  for (unsigned int d = 0; d < dim; ++d)
    fill_tensor(t[d]);
}

template <int dim>
void
fill_tensor(Tensor<1, dim> &t)
{
  for (unsigned int d = 0; d < dim; ++d)
    t[d] = random_value<double>();
}

template <int dim, int spacedim = dim>
void
test()
{
  deallog << "MappingP1<" << dim << ", " << spacedim << ">" << std::endl;

  const auto reference_cell = ReferenceCells::get_simplex<dim>();

  Triangulation<dim, spacedim> tria;
  if constexpr (dim == 1)
    GridGenerator::subdivided_hyper_cube(tria, 1, 1.0, 2.0, true);
  else if constexpr (dim == 2 && spacedim == 3)
    {
      Triangulation<dim> tria_2d;
      GridGenerator::subdivided_hyper_cube_with_simplices(
        tria_2d, 1, 1.0, 2.0, true);
      std::vector<Point<dim>>    points_2d;
      std::vector<CellData<dim>> cell_data;
      SubCellData                subcell_data;
      std::tie(points_2d, cell_data, subcell_data) =
        GridTools::get_coarse_mesh_description(tria_2d);

      std::vector<Point<spacedim>> points_3d;
      for (const auto &p : points_2d)
        points_3d.emplace_back(p[0], p[1], 1.0 + p[0] * p[0] + p[1] * p[1]);

      tria.create_triangulation(points_3d, cell_data, subcell_data);
    }
  else
    GridGenerator::subdivided_hyper_cube_with_simplices(
      tria, 1, 1.0, 2.0, true);
  FE_SimplexP<dim, spacedim> fe(3);
  MappingP1<dim, spacedim>   mapping;
  MappingFE<dim, spacedim>   mapping_fe(fe);

  QGaussSimplex<dim> quadrature(2);

  UpdateFlags flags = update_values | update_gradients | update_hessians |
                      update_3rd_derivatives |
                      update_contravariant_transformation;

  FEValues2<dim, spacedim> fe_values(mapping, fe, quadrature, flags);
  FEValues2<dim, spacedim> fe_values_2(mapping_fe, fe, quadrature, flags);

  for (const auto &cell : tria.active_cell_iterators())
    {
      deallog << std::endl;

      // We need to update the internal states before we can
      // compute values with transform():
      fe_values.reinit(cell);
      fe_values_2.reinit(cell);

      deallog << "cell vertices" << std::endl;
      for (const unsigned int vertex_no : cell->vertex_indices())
        deallog << "  " << cell->vertex(vertex_no) << std::endl;
      deallog << std::endl;

      // First transform():
      {
        std::vector<Tensor<1, dim>> input(quadrature.size());
        for (auto &v : input)
          fill_tensor(v);
        auto input_2 = input;

        auto test_mapping = [&](const MappingKind &mapping_kind,
                                const std::string &mapping_name) {
          deallog << "transform Tensor<1, " << spacedim << "> with "
                  << mapping_name << ":" << std::endl;

          std::vector<Tensor<1, spacedim>> output(quadrature.size());
          auto                             output_2 = output;

          mapping.transform(make_array_view(input),
                            mapping_kind,
                            fe_values.get_mapping_data(),
                            make_array_view(output));
          mapping_fe.transform(make_array_view(input_2),
                               mapping_kind,
                               fe_values_2.get_mapping_data(),
                               make_array_view(output_2));

          for (unsigned int i = 0; i < input.size(); ++i)
            {
              Assert((output[i] - output_2[i]).norm() < 1e-12,
                     ExcInternalError());
              if (i == input.size() / 2)
                deallog << output[i] << std::endl;
            }
        };

        test_mapping(MappingKind::mapping_covariant, "mapping_covariant");
        test_mapping(MappingKind::mapping_contravariant,
                     "mapping_contravariant");
        test_mapping(MappingKind::mapping_piola, "mapping_piola");
      }
      deallog << std::endl;

      // Second transform():
      {
        std::vector<DerivativeForm<1, dim, spacedim>> input(quadrature.size());
        for (auto &v : input)
          for (unsigned int d = 0; d < spacedim; ++d)
            fill_tensor(v[d]);
        auto input_2 = input;

        auto test_mapping = [&](const MappingKind &mapping_kind,
                                const std::string &mapping_name) {
          deallog << "transform DerivativeForm<1, " << dim << ", " << spacedim
                  << "> with " << mapping_name << ":" << std::endl;

          std::vector<Tensor<2, spacedim>> output(quadrature.size());
          auto                             output_2 = output;

          mapping.transform(make_array_view(input),
                            mapping_kind,
                            fe_values.get_mapping_data(),
                            make_array_view(output));
          mapping_fe.transform(make_array_view(input_2),
                               mapping_kind,
                               fe_values_2.get_mapping_data(),
                               make_array_view(output_2));

          for (unsigned int i = 0; i < input.size(); ++i)
            {
              Assert((output[i] - output_2[i]).norm() < 1e-12,
                     ExcInternalError());
              if (i == input.size() / 2)
                deallog << output[i] << std::endl;
            }
        };

        test_mapping(MappingKind::mapping_covariant, "mapping_covariant");
        // Other mappings are not yet implemented for MappingFE so we do not
        // have them in MappingP1 (or tests here)
#if 0
        test_mapping(MappingKind::mapping_contravariant,
                     "mapping_contravariant");
        test_mapping(MappingKind::mapping_piola, "mapping_piola");
#endif
      }
      deallog << std::endl;

      // Third transform():
      //
      // these transformations are only implemented for codim 0
      if (dim == spacedim)
        {
          std::vector<Tensor<2, dim>> input(quadrature.size());
          for (auto &v : input)
            fill_tensor(v);
          auto input_2 = input;

          auto test_mapping = [&](const MappingKind &mapping_kind,
                                  const std::string &mapping_name) {
            deallog << "transform Tensor<2, " << dim << "> with "
                    << mapping_name << ":" << std::endl;

            std::vector<Tensor<2, spacedim>> output(quadrature.size());
            auto                             output_2 = output;

            mapping.transform(make_array_view(input),
                              mapping_kind,
                              fe_values.get_mapping_data(),
                              make_array_view(output));
            mapping_fe.transform(make_array_view(input_2),
                                 mapping_kind,
                                 fe_values_2.get_mapping_data(),
                                 make_array_view(output_2));

            for (unsigned int i = 0; i < input.size(); ++i)
              {
                Assert((output[i] - output_2[i]).norm() < 1e-12,
                       ExcInternalError());
                if (i == input.size() / 2)
                  deallog << output[i] << std::endl;
              }
          };

          // not yet implemented by MappingFE
#if 0
        test_mapping(MappingKind::mapping_covariant, "mapping_covariant");
#endif
          test_mapping(MappingKind::mapping_contravariant,
                       "mapping_contravariant");
          test_mapping(MappingKind::mapping_covariant_gradient,
                       "mapping_covariant_gradient");
          test_mapping(MappingKind::mapping_contravariant_gradient,
                       "mapping_contravariant_gradient");
          test_mapping(MappingKind::mapping_piola_gradient,
                       "mapping_piola_gradient");
        }
      deallog << std::endl;

      {
        std::vector<DerivativeForm<2, dim, spacedim>> input(quadrature.size());
        for (auto &v : input)
          for (unsigned int d = 0; d < spacedim; ++d)
            fill_tensor(v[d]);
        auto input_2 = input;

        auto test_mapping = [&](const MappingKind &mapping_kind,
                                const std::string &mapping_name) {
          deallog << "transform DerivativeForm<2, " << dim << ", " << spacedim
                  << "> with " << mapping_name << ":" << std::endl;

          std::vector<Tensor<3, spacedim>> output(quadrature.size());
          auto                             output_2 = output;

          mapping.transform(make_array_view(input),
                            mapping_kind,
                            fe_values.get_mapping_data(),
                            make_array_view(output));
          mapping_fe.transform(make_array_view(input_2),
                               mapping_kind,
                               fe_values_2.get_mapping_data(),
                               make_array_view(output_2));

          for (unsigned int i = 0; i < input.size(); ++i)
            {
              Assert((output[i] - output_2[i]).norm() < 1e-12,
                     ExcInternalError());
              if (i == input.size() / 2)
                deallog << output[i] << std::endl;
            }
        };

        // This is the only mapping implemented for MappingFE
        test_mapping(MappingKind::mapping_covariant_gradient,
                     "mapping_covariant_gradient");
      }
      deallog << std::endl;

      {
        std::vector<Tensor<3, dim>> input(quadrature.size());
        for (auto &v : input)
          fill_tensor(v);
        auto input_2 = input;

        auto test_mapping = [&](const MappingKind &mapping_kind,
                                const std::string &mapping_name) {
          deallog << "transform Tensor<3, " << dim << "> with " << mapping_name
                  << ":" << std::endl;

          std::vector<Tensor<3, spacedim>> output(quadrature.size());
          auto                             output_2 = output;

          mapping.transform(make_array_view(input),
                            mapping_kind,
                            fe_values.get_mapping_data(),
                            make_array_view(output));
          mapping_fe.transform(make_array_view(input_2),
                               mapping_kind,
                               fe_values_2.get_mapping_data(),
                               make_array_view(output_2));

          for (unsigned int i = 0; i < input.size(); ++i)
            {
              Assert((output[i] - output_2[i]).norm() < 1e-12,
                     ExcInternalError());
              if (i == input.size() / 2)
                deallog << output[i] << std::endl;
            }
        };

        // This is the only mapping implemented for MappingFE
        test_mapping(MappingKind::mapping_contravariant_hessian,
                     "mapping_contravariant_hessian");
        test_mapping(MappingKind::mapping_covariant_hessian,
                     "mapping_covariant_hessian");
        test_mapping(MappingKind::mapping_piola_hessian,
                     "mapping_piola_hessian");
      }
      deallog << std::endl;
    }

  deallog << std::endl;
}

int
main()
{
  initlog();

  test<1>();
  test<1, 2>();
  test<2>();
  test<2, 3>();
  test<3>();

  deallog << "OK" << std::endl;
}

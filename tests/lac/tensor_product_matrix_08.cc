// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test compression of TensorProductMatrixSymmetricSumCollection.

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q_cache.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/tensor_product_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/tensor_product_matrix_creator.h>

#include "../tests.h"

#include "./tensor_product_matrix.h"

template <int dim, typename Number>
void
do_test_mesh(const Mapping<dim> &mapping, const Triangulation<dim> &tria)
{
  using VectorizedArrayTrait = dealii::internal::VectorizedArrayTrait<Number>;
  using ScalarNumber         = typename VectorizedArrayTrait::value_type;
  static constexpr std::size_t width = VectorizedArrayTrait::width();

  using FDM = TensorProductMatrixSymmetricSumCollection<dim, Number>;

  const unsigned int fe_degree = 3;
  const unsigned int n_overlap = 1;

  FE_Q<dim>       fe(fe_degree);
  FE_Q<1>         fe_1D(fe_degree);
  QGauss<1>       quadrature_1D(fe_degree + 1);
  QGauss<dim - 1> quadrature_face(fe_degree + 1);

  const auto harmonic_patch_extent =
    GridTools::compute_harmonic_patch_extent(mapping, tria, quadrature_face);

  FDM collection_0(typename FDM::AdditionalData(true, false));
  FDM collection_1(typename FDM::AdditionalData(false, false));
  FDM collection_2(typename FDM::AdditionalData(true, true));
  FDM collection_3(typename FDM::AdditionalData(false, true));

  const auto n_cell_batches = (tria.n_active_cells() + width - 1) / width;

  collection_0.reserve(n_cell_batches);
  collection_1.reserve(n_cell_batches);
  collection_2.reserve(n_cell_batches);
  collection_3.reserve(n_cell_batches);

  auto cell = tria.begin_active();

  for (unsigned int counter = 0; counter < n_cell_batches; ++counter)
    {
      std::array<Table<2, Number>, dim> Ms;
      std::array<Table<2, Number>, dim> Ks;

      for (unsigned int v = 0; (v < width) && (cell != tria.end()); ++v, ++cell)
        {
          const auto &patch_extent =
            harmonic_patch_extent[cell->active_cell_index()];

          const auto M_and_K_scalar = TensorProductMatrixCreator::
            create_laplace_tensor_product_matrix<dim, ScalarNumber>(
              cell,
              std::set<unsigned int>{0},
              std::set<unsigned int>{},
              fe_1D,
              quadrature_1D,
              patch_extent,
              n_overlap);

          const auto Ms_scalar = M_and_K_scalar.first;
          const auto Ks_scalar = M_and_K_scalar.second;

          for (unsigned int d = 0; d < dim; ++d)
            {
              if (Ms[d].size(0) == 0 || Ms[d].size(1) == 0)
                {
                  Ms[d].reinit(Ms_scalar[d].size(0), Ms_scalar[d].size(1));
                  Ks[d].reinit(Ks_scalar[d].size(0), Ks_scalar[d].size(1));
                }

              for (unsigned int i = 0; i < Ms_scalar[d].size(0); ++i)
                for (unsigned int j = 0; j < Ms_scalar[d].size(0); ++j)
                  VectorizedArrayTrait::get(Ms[d][i][j], v) =
                    Ms_scalar[d][i][j];

              for (unsigned int i = 0; i < Ks_scalar[d].size(0); ++i)
                for (unsigned int j = 0; j < Ks_scalar[d].size(0); ++j)
                  VectorizedArrayTrait::get(Ks[d][i][j], v) =
                    Ks_scalar[d][i][j];
            }
        }

      collection_0.insert(counter, Ms, Ks);
      collection_1.insert(counter, Ms, Ks);
      collection_2.insert(counter, Ms, Ks);
      collection_3.insert(counter, Ms, Ks);
    }

  collection_0.finalize();
  collection_1.finalize();
  collection_2.finalize();
  collection_3.finalize();

  deallog << "Storage sizes: " << collection_0.storage_size() << " "
          << collection_1.storage_size() << std::endl;

  AlignedVector<Number> src(fe.n_dofs_per_cell());
  AlignedVector<Number> dst(fe.n_dofs_per_cell());
  Table<2, Number>      matrix_0(fe.n_dofs_per_cell(), fe.n_dofs_per_cell());
  Table<2, Number>      matrix_1(fe.n_dofs_per_cell(), fe.n_dofs_per_cell());
  Table<2, Number>      matrix_2(fe.n_dofs_per_cell(), fe.n_dofs_per_cell());
  Table<2, Number>      matrix_3(fe.n_dofs_per_cell(), fe.n_dofs_per_cell());

  for (unsigned int cell = 0; cell < n_cell_batches; ++cell)
    {
      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
        {
          for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
            src[j] = i == j;

          collection_0.apply_inverse(cell,
                                     make_array_view(dst.begin(), dst.end()),
                                     make_array_view(src.begin(), src.end()));
          for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
            matrix_0[j][i] = dst[j];

          collection_1.apply_inverse(cell,
                                     make_array_view(dst.begin(), dst.end()),
                                     make_array_view(src.begin(), src.end()));
          for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
            matrix_1[j][i] = dst[j];

          collection_2.apply_inverse(cell,
                                     make_array_view(dst.begin(), dst.end()),
                                     make_array_view(src.begin(), src.end()));
          for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
            matrix_2[j][i] = dst[j];

          collection_3.apply_inverse(cell,
                                     make_array_view(dst.begin(), dst.end()),
                                     make_array_view(src.begin(), src.end()));
          for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
            matrix_3[j][i] = dst[j];
        }


      FloatingPointComparator<Number> comp(1e-5, true);

      Assert((comp.compare(matrix_0, matrix_1) ==
              FloatingPointComparator<Number>::ComparisonResult::equal),
             ExcInternalError());
      Assert((comp.compare(matrix_0, matrix_2) ==
              FloatingPointComparator<Number>::ComparisonResult::equal),
             ExcInternalError());
      Assert((comp.compare(matrix_0, matrix_3) ==
              FloatingPointComparator<Number>::ComparisonResult::equal),
             ExcInternalError());
    }

  deallog << "OK!" << std::endl;
  deallog << std::endl;
}


template <typename Number>
void
do_test()
{
  const int dim = 2;

  MappingQ1<dim> mapping;

  const std::vector<std::tuple<unsigned int, bool, bool>> settings = {
    std::tuple<unsigned int, bool, bool>{0, false, false}, // hyper-cube
    std::tuple<unsigned int, bool, bool>{0, true, false}, // hyper-cube with PBC
    std::tuple<unsigned int, bool, bool>{1, false, false}, // hyper-rectangle
    std::tuple<unsigned int, bool, bool>{1,
                                         true,
                                         false}, // hyper-rectangle with PBC
    std::tuple<unsigned int, bool, bool>{0, false, true} // deformed hyper-cube
  };

  for (const auto &setting : settings)
    {
      Triangulation<dim> tria;

      const auto do_pbc = std::get<1>(setting);

      switch (std::get<0>(setting))
        {
          case 0:
            GridGenerator::hyper_cube(tria, 0.0, 1.0, do_pbc);
            break;
          case 1:
            GridGenerator::subdivided_hyper_rectangle(
              tria, {1, 2}, {0, 0}, {1.0, 1.0}, do_pbc);
            break;
          default:
            AssertThrow(false, ExcNotImplemented());
        }

      if (do_pbc)
        {
          std::vector<GridTools::PeriodicFacePair<
            typename Triangulation<dim>::cell_iterator>>
            periodic_faces;
          for (unsigned int d = 0; d < dim; ++d)
            GridTools::collect_periodic_faces(
              tria, 2 * d, 2 * d + 1, d, periodic_faces);
          tria.add_periodicity(periodic_faces);
        }

      tria.refine_global(3);

      MappingQCache<dim> mapping(3);

      mapping.initialize(
        MappingQ1<dim>(),
        tria,
        [&](const typename Triangulation<dim>::cell_iterator &,
            const Point<dim> &point) -> Point<dim> {
          Point<dim> result;

          if (std::get<2>(setting) == false) // Cartesian mesh
            return result;

          for (unsigned int d = 0; d < dim; ++d)
            result[d] = std::pow(point[d], d + 0.9) *
                        std::pow(point[(d + 1) % dim], d + 1.1);

          return result;
        },
        true);


      do_test_mesh<dim, Number>(mapping, tria);
    }
}


int
main()
{
  initlog();

  {
    deallog.push("v=0");
    do_test<double>();
    deallog << std::endl;
    do_test<float>();
    deallog.pop();
    deallog << std::endl;
  }

  {
    deallog.push("v=64");
    do_test<VectorizedArray<double, 1>>();
    deallog << std::endl;
    do_test<VectorizedArray<float, 1>>();
    deallog.pop();
    deallog << std::endl;
  }

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128
  {
    deallog.push("v=128");
    do_test<VectorizedArray<double, 2>>();
    deallog << std::endl;
    do_test<VectorizedArray<float, 4>>();
    deallog.pop();
    deallog << std::endl;
  }
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256
  {
    deallog.push("v=256");
    do_test<VectorizedArray<double, 4>>();
    deallog << std::endl;
    do_test<VectorizedArray<float, 8>>();
    deallog.pop();
    deallog << std::endl;
  }
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512
  {
    deallog.push("v=512");
    do_test<VectorizedArray<double, 8>>();
    deallog << std::endl;
    do_test<VectorizedArray<float, 16>>();
    deallog.pop();
    deallog << std::endl;
  }
#endif

  return 0;
}

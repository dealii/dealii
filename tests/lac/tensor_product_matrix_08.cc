// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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
  using FDM = TensorProductMatrixSymmetricSumCollection<dim, Number>;

  const unsigned int fe_degree = 3;
  const unsigned int n_overlap = 1;

  FE_Q<dim>       fe(fe_degree);
  FE_Q<1>         fe_1D(fe_degree);
  QGauss<1>       quadrature_1D(fe_degree + 1);
  QGauss<dim - 1> quadrature_face(fe_degree + 1);

  const auto harmonic_patch_extent =
    GridTools::compute_harmonic_patch_extent(mapping, tria, quadrature_face);

  FDM collection_0(typename FDM::AdditionalData(true));
  FDM collection_1(typename FDM::AdditionalData(false));

  collection_0.reserve(tria.n_active_cells());
  collection_1.reserve(tria.n_active_cells());

  for (const auto &cell : tria.active_cell_iterators())
    {
      const auto &patch_extent =
        harmonic_patch_extent[cell->active_cell_index()];

      const auto M_and_K = TensorProductMatrixCreator::
        create_laplace_tensor_product_matrix<dim, Number>(
          cell,
          std::set<unsigned int>{0},
          std::set<unsigned int>{},
          fe_1D,
          quadrature_1D,
          patch_extent,
          n_overlap);

      collection_0.insert(cell->active_cell_index(),
                          M_and_K.first,
                          M_and_K.second);
      collection_1.insert(cell->active_cell_index(),
                          M_and_K.first,
                          M_and_K.second);
    }

  collection_0.finalize();
  collection_1.finalize();

  deallog << "Storage sizes: " << collection_0.storage_size() << " "
          << collection_1.storage_size() << std::endl;

  Vector<Number>        src(fe.n_dofs_per_cell());
  Vector<Number>        dst(fe.n_dofs_per_cell());
  AlignedVector<Number> tmp;
  FullMatrix<Number>    matrix_0(fe.n_dofs_per_cell(), fe.n_dofs_per_cell());
  FullMatrix<Number>    matrix_1(fe.n_dofs_per_cell(), fe.n_dofs_per_cell());

  for (unsigned int cell = 0; cell < tria.n_active_cells(); ++cell)
    {
      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
        {
          for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
            src[j] = i == j;

          collection_0.apply_inverse(cell,
                                     make_array_view(dst),
                                     make_array_view(src),
                                     tmp);
          for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
            matrix_0[j][i] = dst[j];

          collection_1.apply_inverse(cell,
                                     make_array_view(dst),
                                     make_array_view(src),
                                     tmp);
          for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
            matrix_1[j][i] = dst[j];
        }


      FloatingPointComparator<Number> comp(1e-5, false);

      Assert((comp.compare(matrix_0, matrix_1) ==
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

  for (const auto setting : settings)
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

  do_test<double>();
  do_test<float>();

  return 0;
}

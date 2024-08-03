// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Compare the result of FEEvaluationImplHangingNodes with the CPU version
// of GPU implementation.

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/evaluation_kernels_hanging_nodes.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


namespace dealii
{
  namespace internal
  {
    /**
     * Evaluation kernels similar to them used in the GPU implementation
     * (see include/deal.II/matrix_free/hanging_nodes_internal.h).
     */
    template <int dim, typename Number, bool is_face>
    struct FEEvaluationImplHangingNodesReference
    {
      template <int fe_degree>
      static bool
      run(
        const FEEvaluationData<dim, Number, is_face> &fe_eval,
        const bool                                    transpose,
        const std::array<dealii::internal::MatrixFreeFunctions::ConstraintKinds,
                         Number::size()>             &c_mask,
        Number                                       *values)
      {
        Assert(is_face == false, ExcInternalError());

        if (dim == 2)
          {
            if (transpose)
              {
                run_2D<fe_degree, 0, true>(fe_eval, c_mask, values);
                run_2D<fe_degree, 1, true>(fe_eval, c_mask, values);
              }
            else
              {
                run_2D<fe_degree, 0, false>(fe_eval, c_mask, values);
                run_2D<fe_degree, 1, false>(fe_eval, c_mask, values);
              }
          }
        else if (dim == 3)
          {
            if (transpose)
              {
                run_3D<fe_degree, 0, true>(fe_eval, c_mask, values);
                run_3D<fe_degree, 1, true>(fe_eval, c_mask, values);
                run_3D<fe_degree, 2, true>(fe_eval, c_mask, values);
              }
            else
              {
                run_3D<fe_degree, 0, false>(fe_eval, c_mask, values);
                run_3D<fe_degree, 1, false>(fe_eval, c_mask, values);
                run_3D<fe_degree, 2, false>(fe_eval, c_mask, values);
              }
          }

        return false; // TODO
      }

    private:
      static unsigned int
      index2(unsigned int size, unsigned int i, unsigned int j)
      {
        return i + size * j;
      }

      static inline unsigned int
      index3(unsigned int size, unsigned int i, unsigned int j, unsigned int k)
      {
        return i + size * j + size * size * k;
      }

      template <int fe_degree_, unsigned int direction, bool transpose>
      static void
      run_2D(
        const FEEvaluationData<dim, Number, is_face> &fe_eval,
        const std::array<dealii::internal::MatrixFreeFunctions::ConstraintKinds,
                         Number::size()>             &constraint_mask,
        Number                                       *values)
      {
        const auto &constraint_weights = fe_eval.get_shape_info()
                                           .data.front()
                                           .subface_interpolation_matrices[0];

        const unsigned int fe_degree =
          fe_degree_ != -1 ? fe_degree_ :
                             fe_eval.get_shape_info().data.front().fe_degree;

        const unsigned int n_dofs =
          Utilities::pow<unsigned int>(fe_degree + 1, 2);

        AlignedVector<Number> values_temp(n_dofs);

        for (unsigned int i = 0; i < n_dofs; ++i)
          values_temp[i] = values[i];

        for (unsigned int v = 0; v < Number::size(); ++v)
          {
            if (constraint_mask[v] == dealii::internal::MatrixFreeFunctions::
                                        ConstraintKinds::unconstrained)
              continue;

            const auto this_type = (direction == 0) ?
                                     dealii::internal::MatrixFreeFunctions::
                                       ConstraintKinds::subcell_x :
                                     dealii::internal::MatrixFreeFunctions::
                                       ConstraintKinds::subcell_y;

            const bool constrained_face =
              (constraint_mask[v] &
               (((direction == 0) ? dealii::internal::MatrixFreeFunctions::
                                      ConstraintKinds::face_y :
                                    dealii::internal::MatrixFreeFunctions::
                                      ConstraintKinds::unconstrained) |
                ((direction == 1) ? dealii::internal::MatrixFreeFunctions::
                                      ConstraintKinds::face_x :
                                    dealii::internal::MatrixFreeFunctions::
                                      ConstraintKinds::unconstrained))) !=
              dealii::internal::MatrixFreeFunctions::ConstraintKinds::
                unconstrained;

            const bool type = (constraint_mask[v] & this_type) !=
                              dealii::internal::MatrixFreeFunctions::
                                ConstraintKinds::unconstrained;

            for (unsigned int x_idx = 0; x_idx < fe_degree + 1; ++x_idx)
              for (unsigned int y_idx = 0; y_idx < fe_degree + 1; ++y_idx)
                {
                  const unsigned int interp_idx =
                    (direction == 0) ? x_idx : y_idx;

                  typename Number::value_type t = 0.0;
                  // Flag is true if dof is constrained for the given direction
                  // and the given face.

                  // Flag is true if for the given direction, the dof is
                  // constrained with the right type and is on the correct side
                  // (left (= 0) or right (= fe_degree))
                  const bool constrained_dof =
                    ((direction == 0) &&
                     (((constraint_mask[v] &
                        dealii::internal::MatrixFreeFunctions::ConstraintKinds::
                          subcell_y) != dealii::internal::MatrixFreeFunctions::
                                          ConstraintKinds::unconstrained) ?
                        (y_idx == 0) :
                        (y_idx == fe_degree))) ||
                    ((direction == 1) &&
                     (((constraint_mask[v] &
                        dealii::internal::MatrixFreeFunctions::ConstraintKinds::
                          subcell_x) != dealii::internal::MatrixFreeFunctions::
                                          ConstraintKinds::unconstrained) ?
                        (x_idx == 0) :
                        (x_idx == fe_degree)));

                  if (constrained_face && constrained_dof)
                    {
                      const unsigned int real_idx =
                        index2(fe_degree + 1, x_idx, y_idx);

                      // deallog << "dir=" << direction << " real=" << real_idx
                      // << std::endl;


                      if (type)
                        {
                          for (unsigned int i = 0; i <= fe_degree; ++i)
                            {
                              const unsigned int real_idx =
                                (direction == 0) ?
                                  index2(fe_degree + 1, i, y_idx) :
                                  index2(fe_degree + 1, x_idx, i);

                              const auto w =
                                transpose ?
                                  constraint_weights[i * (fe_degree + 1) +
                                                     interp_idx] :
                                  constraint_weights[interp_idx *
                                                       (fe_degree + 1) +
                                                     i];
                              t += w * values_temp[real_idx][v];
                            }
                        }
                      else
                        {
                          for (unsigned int i = 0; i <= fe_degree; ++i)
                            {
                              const unsigned int real_idx =
                                (direction == 0) ?
                                  index2(fe_degree + 1, i, y_idx) :
                                  index2(fe_degree + 1, x_idx, i);

                              const auto w =
                                transpose ?
                                  constraint_weights[(fe_degree - i) *
                                                       (fe_degree + 1) +
                                                     fe_degree - interp_idx] :
                                  constraint_weights[(fe_degree - interp_idx) *
                                                       (fe_degree + 1) +
                                                     fe_degree - i];
                              t += w * values_temp[real_idx][v];
                            }
                        }

                      values[index2(fe_degree + 1, x_idx, y_idx)][v] = t;
                    }
                }
          }
      }

      template <int fe_degree_, unsigned int direction, bool transpose>
      static void
      run_3D(
        const FEEvaluationData<dim, Number, is_face> &fe_eval,
        const std::array<dealii::internal::MatrixFreeFunctions::ConstraintKinds,
                         Number::size()>             &constraint_mask,
        Number                                       *values)
      {
        const auto &constraint_weights = fe_eval.get_shape_info()
                                           .data.front()
                                           .subface_interpolation_matrices[0];

        const unsigned int fe_degree =
          fe_degree_ != -1 ? fe_degree_ :
                             fe_eval.get_shape_info().data.front().fe_degree;

        const unsigned int n_dofs =
          Utilities::pow<unsigned int>(fe_degree + 1, 3);

        AlignedVector<Number> values_temp(n_dofs);

        for (unsigned int i = 0; i < n_dofs; ++i)
          values_temp[i] = values[i];

        const auto this_type =
          (direction == 0) ?
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x :
          (direction == 1) ?
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y :
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_z;
        const auto face1_type =
          (direction == 0) ?
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y :
          (direction == 1) ?
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_z :
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x;
        const auto face2_type =
          (direction == 0) ?
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_z :
          (direction == 1) ?
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x :
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y;

        // If computing in x-direction, need to match against
        // dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_y or
        // dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_z
        const auto face1 =
          (direction == 0) ?
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_y :
          (direction == 1) ?
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_z :
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_x;
        const auto face2 =
          (direction == 0) ?
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_z :
          (direction == 1) ?
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_x :
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_y;
        const auto edge =
          (direction == 0) ?
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_x :
          (direction == 1) ?
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_y :
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_z;

        for (unsigned int v = 0; v < Number::size(); ++v)
          {
            if (constraint_mask[v] == dealii::internal::MatrixFreeFunctions::
                                        ConstraintKinds::unconstrained)
              continue;

            const auto constrained_face =
              constraint_mask[v] & (face1 | face2 | edge);

            const bool type = (constraint_mask[v] & this_type) !=
                              dealii::internal::MatrixFreeFunctions::
                                ConstraintKinds::unconstrained;

            for (unsigned int x_idx = 0; x_idx < fe_degree + 1; ++x_idx)
              for (unsigned int y_idx = 0; y_idx < fe_degree + 1; ++y_idx)
                for (unsigned int z_idx = 0; z_idx < fe_degree + 1; ++z_idx)
                  {
                    const unsigned int interp_idx = (direction == 0) ? x_idx :
                                                    (direction == 1) ? y_idx :
                                                                       z_idx;
                    const unsigned int face1_idx  = (direction == 0) ? y_idx :
                                                    (direction == 1) ? z_idx :
                                                                       x_idx;
                    const unsigned int face2_idx  = (direction == 0) ? z_idx :
                                                    (direction == 1) ? x_idx :
                                                                       y_idx;

                    typename Number::value_type t = 0;
                    const bool                  on_face1 =
                      ((constraint_mask[v] & face1_type) !=
                       dealii::internal::MatrixFreeFunctions::ConstraintKinds::
                         unconstrained) ?
                                         (face1_idx == 0) :
                                         (face1_idx == fe_degree);
                    const bool on_face2 =
                      ((constraint_mask[v] & face2_type) !=
                       dealii::internal::MatrixFreeFunctions::ConstraintKinds::
                         unconstrained) ?
                        (face2_idx == 0) :
                        (face2_idx == fe_degree);
                    const bool constrained_dof =
                      ((((constraint_mask[v] & face1) !=
                         dealii::internal::MatrixFreeFunctions::
                           ConstraintKinds::unconstrained) &&
                        on_face1) ||
                       (((constraint_mask[v] & face2) !=
                         dealii::internal::MatrixFreeFunctions::
                           ConstraintKinds::unconstrained) &&
                        on_face2) ||
                       (((constraint_mask[v] & edge) !=
                         dealii::internal::MatrixFreeFunctions::
                           ConstraintKinds::unconstrained) &&
                        on_face1 && on_face2));

                    if ((constrained_face !=
                         dealii::internal::MatrixFreeFunctions::
                           ConstraintKinds::unconstrained) &&
                        constrained_dof)
                      {
                        if (type)
                          {
                            for (unsigned int i = 0; i <= fe_degree; ++i)
                              {
                                const unsigned int real_idx =
                                  (direction == 0) ?
                                    index3(fe_degree + 1, i, y_idx, z_idx) :
                                  (direction == 1) ?
                                    index3(fe_degree + 1, x_idx, i, z_idx) :
                                    index3(fe_degree + 1, x_idx, y_idx, i);

                                const auto w =
                                  transpose ?
                                    constraint_weights[i * (fe_degree + 1) +
                                                       interp_idx] :
                                    constraint_weights[interp_idx *
                                                         (fe_degree + 1) +
                                                       i];
                                t += w * values_temp[real_idx][v];
                              }
                          }
                        else
                          {
                            for (unsigned int i = 0; i <= fe_degree; ++i)
                              {
                                const unsigned int real_idx =
                                  (direction == 0) ?
                                    index3(fe_degree + 1, i, y_idx, z_idx) :
                                  (direction == 1) ?
                                    index3(fe_degree + 1, x_idx, i, z_idx) :
                                    index3(fe_degree + 1, x_idx, y_idx, i);

                                const auto w =
                                  transpose ?
                                    constraint_weights[(fe_degree - i) *
                                                         (fe_degree + 1) +
                                                       fe_degree - interp_idx] :
                                    constraint_weights[(fe_degree -
                                                        interp_idx) *
                                                         (fe_degree + 1) +
                                                       fe_degree - i];
                                t += w * values_temp[real_idx][v];
                              }
                          }

                        values[index3(fe_degree + 1, x_idx, y_idx, z_idx)][v] =
                          t;
                      }
                  }
          }
      }
    };
  } // namespace internal
} // namespace dealii

template <int dim>
void
test(const unsigned int                                           degree,
     const dealii::internal::MatrixFreeFunctions::ConstraintKinds mask_value)
{
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube(tria, 2);
  tria.begin()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  QGauss<dim>    quad(degree + 1);
  FE_Q<dim>      fe(degree);
  MappingQ1<dim> mapping;

  DoFHandler<dim> dof_handler;
  dof_handler.reinit(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<double> constraints;

  typename MatrixFree<dim, double, VectorizedArray<double>>::AdditionalData
    additional_data;
  additional_data.mapping_update_flags = update_values | update_gradients |
                                         update_JxW_values |
                                         dealii::update_quadrature_points;

  MatrixFree<dim, double, VectorizedArray<double>> matrix_free;
  matrix_free.reinit(mapping, dof_handler, constraints, quad, additional_data);

  FEEvaluation<dim, -1, 0, 1, double> eval(matrix_free);
  eval.reinit(0);

  std::array<dealii::internal::MatrixFreeFunctions::compressed_constraint_kind,
             VectorizedArray<double>::size()>
    cmask;
  std::fill(cmask.begin(),
            cmask.end(),
            dealii::internal::MatrixFreeFunctions::
              unconstrained_compressed_constraint_kind);
  cmask[0] = dealii::internal::MatrixFreeFunctions::compress(mask_value, dim);

  std::array<dealii::internal::MatrixFreeFunctions::ConstraintKinds,
             VectorizedArray<double>::size()>
    cmask_;
  std::fill(
    cmask_.begin(),
    cmask_.end(),
    dealii::internal::MatrixFreeFunctions::ConstraintKinds::unconstrained);
  cmask_[0] = mask_value;

  for (unsigned int b = 0; b < 2; ++b)
    {
      AlignedVector<VectorizedArray<double>> values1(fe.n_dofs_per_cell());
      AlignedVector<VectorizedArray<double>> values2(fe.n_dofs_per_cell());

      for (unsigned int i = 0; i < values1.size(); ++i)
        {
          values1[i][0] = i;
          values2[i][0] = i;
        }

      for (const auto i : values1)
        deallog << i[0] << ' ';
      deallog << std::endl;

      internal::FEEvaluationImplHangingNodesReference<
        dim,
        VectorizedArray<double>,
        false>::template run<-1>(eval, b == 1, cmask_, values1.data());
      internal::FEEvaluationImplHangingNodes<dim, VectorizedArray<double>>::
        template run<-1>(
          1, eval.get_shape_info(), b == 1, cmask, values2.data());

      for (const auto i : values1)
        deallog << i[0] << ' ';
      deallog << std::endl;

      for (const auto i : values2)
        deallog << i[0] << ' ';
      deallog << std::endl;
      deallog << std::endl;

      for (unsigned int i = 0; i < values1.size(); ++i)
        Assert(std::abs(values1[i][0] - values2[i][0]) < 1e-5,
               ExcInternalError());
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  using namespace dealii::internal;

  for (unsigned int degree = 1; degree <= 3; ++degree)
    {
      test<2>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::unconstrained);
      deallog << std::endl;

      test<2>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_x |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::
            subcell_y); // face 0/0
      test<2>(degree,
              dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_x |
                dealii::internal::MatrixFreeFunctions::ConstraintKinds::
                  subcell_x); // face 0/1
      test<2>(degree,
              dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_x |
                dealii::internal::MatrixFreeFunctions::ConstraintKinds::
                  subcell_y); // face 1/0
      test<2>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_x); // face
                                                                         // 1/1
      deallog << std::endl;

      test<2>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_y |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::
            subcell_x); // face 2/0
      test<2>(degree,
              dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_y |
                dealii::internal::MatrixFreeFunctions::ConstraintKinds::
                  subcell_y); // face 2/1
      test<2>(degree,
              dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_y |
                dealii::internal::MatrixFreeFunctions::ConstraintKinds::
                  subcell_x); // face 3/0
      test<2>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_y); // face
                                                                         // 3/1
      deallog << std::endl;
    }

  for (unsigned int degree = 1; degree <= 3; ++degree)
    {
      // edge 2
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_x |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_z);
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_x |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_z |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x);

      // edge 3
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_x |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_z);
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_x |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_z |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x);

      // edge 6
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_x |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y);
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_x |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x);

      // edge 7
      test<3>(degree,
              dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_x);
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_x |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x);


      // edge 0
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_y |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_z);
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_y |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_z |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y);

      // edge 1
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_y |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_z);
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_y |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_z |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y);

      // edge 4
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_y |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x);
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_y |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y);

      // edge 5
      test<3>(degree,
              dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_y);
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_y |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y);


      // edge 8
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_z |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y);
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_z |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_z);

      // edge 9
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_z |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y);
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_z |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_z);

      // edge 10
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_z |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x);
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_z |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_z);

      // edge 11
      test<3>(degree,
              dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_z);
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_z |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_z);


      // face 0
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_x |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x);

      // face 1
      test<3>(degree,
              dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_x);

      // face 2
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_y |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y);

      // face 3
      test<3>(degree,
              dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_y);

      // face 4
      test<3>(
        degree,
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_z |
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_z);

      // face 5
      test<3>(degree,
              dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_z);
    }
}

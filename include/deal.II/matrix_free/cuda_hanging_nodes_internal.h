// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii_cuda_hanging_nodes_internal_h
#define dealii_cuda_hanging_nodes_internal_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_COMPILER_CUDA_AWARE

#  include <deal.II/base/cuda_size.h>

#  include <deal.II/matrix_free/hanging_nodes_internal.h>

DEAL_II_NAMESPACE_OPEN
namespace CUDAWrappers
{
  namespace internal
  {
    __constant__ double
      constraint_weights[(CUDAWrappers::mf_max_elem_degree + 1) *
                         (CUDAWrappers::mf_max_elem_degree + 1)];

    //------------------------------------------------------------------------//
    // Functions for resolving the hanging node constraints on the GPU        //
    //------------------------------------------------------------------------//
    template <unsigned int size>
    __device__ inline unsigned int
    index2(unsigned int i, unsigned int j)
    {
      return i + size * j;
    }



    template <unsigned int size>
    __device__ inline unsigned int
    index3(unsigned int i, unsigned int j, unsigned int k)
    {
      return i + size * j + size * size * k;
    }



    template <unsigned int fe_degree,
              unsigned int direction,
              bool         transpose,
              typename Number>
    __device__ inline void
    interpolate_boundary_2d(
      const dealii::internal::MatrixFreeFunctions::ConstraintKinds
              constraint_mask,
      Number *values)
    {
      const unsigned int x_idx = threadIdx.x % (fe_degree + 1);
      const unsigned int y_idx = threadIdx.y;

      const auto this_type =
        (direction == 0) ?
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x :
          dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y;

      const unsigned int interp_idx = (direction == 0) ? x_idx : y_idx;

      Number t = 0;
      // Flag is true if dof is constrained for the given direction and the
      // given face.
      const bool constrained_face =
        (constraint_mask &
         (((direction == 0) ?
             dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_y :
             dealii::internal::MatrixFreeFunctions::ConstraintKinds::
               unconstrained) |
          ((direction == 1) ?
             dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_x :
             dealii::internal::MatrixFreeFunctions::ConstraintKinds::
               unconstrained))) !=
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::unconstrained;

      // Flag is true if for the given direction, the dof is constrained with
      // the right type and is on the correct side (left (= 0) or right (=
      // fe_degree))
      const bool constrained_dof =
        ((direction == 0) &&
         (((constraint_mask & dealii::internal::MatrixFreeFunctions::
                                ConstraintKinds::subcell_y) !=
           dealii::internal::MatrixFreeFunctions::ConstraintKinds::
             unconstrained) ?
            (y_idx == 0) :
            (y_idx == fe_degree))) ||
        ((direction == 1) &&
         (((constraint_mask & dealii::internal::MatrixFreeFunctions::
                                ConstraintKinds::subcell_x) !=
           dealii::internal::MatrixFreeFunctions::ConstraintKinds::
             unconstrained) ?
            (x_idx == 0) :
            (x_idx == fe_degree)));

      if (constrained_face && constrained_dof)
        {
          const bool type = (constraint_mask & this_type) !=
                            dealii::internal::MatrixFreeFunctions::
                              ConstraintKinds::unconstrained;

          if (type)
            {
              for (unsigned int i = 0; i <= fe_degree; ++i)
                {
                  const unsigned int real_idx =
                    (direction == 0) ? index2<fe_degree + 1>(i, y_idx) :
                                       index2<fe_degree + 1>(x_idx, i);

                  const Number w =
                    transpose ?
                      constraint_weights[i * (fe_degree + 1) + interp_idx] :
                      constraint_weights[interp_idx * (fe_degree + 1) + i];
                  t += w * values[real_idx];
                }
            }
          else
            {
              for (unsigned int i = 0; i <= fe_degree; ++i)
                {
                  const unsigned int real_idx =
                    (direction == 0) ? index2<fe_degree + 1>(i, y_idx) :
                                       index2<fe_degree + 1>(x_idx, i);

                  const Number w =
                    transpose ?
                      constraint_weights[(fe_degree - i) * (fe_degree + 1) +
                                         fe_degree - interp_idx] :
                      constraint_weights[(fe_degree - interp_idx) *
                                           (fe_degree + 1) +
                                         fe_degree - i];
                  t += w * values[real_idx];
                }
            }
        }

      // The synchronization is done for all the threads in one block with
      // each block being assigned to one element.
      __syncthreads();
      if (constrained_face && constrained_dof)
        values[index2<fe_degree + 1>(x_idx, y_idx)] = t;

      __syncthreads();
    }



    template <unsigned int fe_degree,
              unsigned int direction,
              bool         transpose,
              typename Number>
    __device__ inline void
    interpolate_boundary_3d(
      const dealii::internal::MatrixFreeFunctions::ConstraintKinds
              constraint_mask,
      Number *values)
    {
      const unsigned int x_idx = threadIdx.x % (fe_degree + 1);
      const unsigned int y_idx = threadIdx.y;
      const unsigned int z_idx = threadIdx.z;

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

      // If computing in x-direction, need to match against face_y or
      // face_z
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
      const auto constrained_face = constraint_mask & (face1 | face2 | edge);

      const unsigned int interp_idx = (direction == 0) ? x_idx :
                                      (direction == 1) ? y_idx :
                                                         z_idx;
      const unsigned int face1_idx  = (direction == 0) ? y_idx :
                                      (direction == 1) ? z_idx :
                                                         x_idx;
      const unsigned int face2_idx  = (direction == 0) ? z_idx :
                                      (direction == 1) ? x_idx :
                                                         y_idx;

      Number     t        = 0;
      const bool on_face1 = ((constraint_mask & face1_type) !=
                             dealii::internal::MatrixFreeFunctions::
                               ConstraintKinds::unconstrained) ?
                              (face1_idx == 0) :
                              (face1_idx == fe_degree);
      const bool on_face2 = ((constraint_mask & face2_type) !=
                             dealii::internal::MatrixFreeFunctions::
                               ConstraintKinds::unconstrained) ?
                              (face2_idx == 0) :
                              (face2_idx == fe_degree);
      const bool constrained_dof =
        ((((constraint_mask & face1) != dealii::internal::MatrixFreeFunctions::
                                          ConstraintKinds::unconstrained) &&
          on_face1) ||
         (((constraint_mask & face2) != dealii::internal::MatrixFreeFunctions::
                                          ConstraintKinds::unconstrained) &&
          on_face2) ||
         (((constraint_mask & edge) != dealii::internal::MatrixFreeFunctions::
                                         ConstraintKinds::unconstrained) &&
          on_face1 && on_face2));

      if ((constrained_face != dealii::internal::MatrixFreeFunctions::
                                 ConstraintKinds::unconstrained) &&
          constrained_dof)
        {
          const bool type = (constraint_mask & this_type) !=
                            dealii::internal::MatrixFreeFunctions::
                              ConstraintKinds::unconstrained;
          if (type)
            {
              for (unsigned int i = 0; i <= fe_degree; ++i)
                {
                  const unsigned int real_idx =
                    (direction == 0) ? index3<fe_degree + 1>(i, y_idx, z_idx) :
                    (direction == 1) ? index3<fe_degree + 1>(x_idx, i, z_idx) :
                                       index3<fe_degree + 1>(x_idx, y_idx, i);

                  const Number w =
                    transpose ?
                      constraint_weights[i * (fe_degree + 1) + interp_idx] :
                      constraint_weights[interp_idx * (fe_degree + 1) + i];
                  t += w * values[real_idx];
                }
            }
          else
            {
              for (unsigned int i = 0; i <= fe_degree; ++i)
                {
                  const unsigned int real_idx =
                    (direction == 0) ? index3<fe_degree + 1>(i, y_idx, z_idx) :
                    (direction == 1) ? index3<fe_degree + 1>(x_idx, i, z_idx) :
                                       index3<fe_degree + 1>(x_idx, y_idx, i);

                  const Number w =
                    transpose ?
                      constraint_weights[(fe_degree - i) * (fe_degree + 1) +
                                         fe_degree - interp_idx] :
                      constraint_weights[(fe_degree - interp_idx) *
                                           (fe_degree + 1) +
                                         fe_degree - i];
                  t += w * values[real_idx];
                }
            }
        }

      // The synchronization is done for all the threads in one block with
      // each block being assigned to one element.
      __syncthreads();

      if ((constrained_face != dealii::internal::MatrixFreeFunctions::
                                 ConstraintKinds::unconstrained) &&
          constrained_dof)
        values[index3<fe_degree + 1>(x_idx, y_idx, z_idx)] = t;

      __syncthreads();
    }



    /**
     * This function resolves the hanging nodes using tensor product.
     *
     * The implementation of this class is explained in Section 3 of
     * @cite ljungkvist2017matrix and in Section 3.4 of
     * @cite kronbichler2019multigrid.
     */
    template <int dim, int fe_degree, bool transpose, typename Number>
    __device__ void
    resolve_hanging_nodes(
      const dealii::internal::MatrixFreeFunctions::ConstraintKinds
              constraint_mask,
      Number *values)
    {
      if (dim == 2)
        {
          interpolate_boundary_2d<fe_degree, 0, transpose>(constraint_mask,
                                                           values);

          interpolate_boundary_2d<fe_degree, 1, transpose>(constraint_mask,
                                                           values);
        }
      else if (dim == 3)
        {
          // Interpolate y and z faces (x-direction)
          interpolate_boundary_3d<fe_degree, 0, transpose>(constraint_mask,
                                                           values);
          // Interpolate x and z faces (y-direction)
          interpolate_boundary_3d<fe_degree, 1, transpose>(constraint_mask,
                                                           values);
          // Interpolate x and y faces (z-direction)
          interpolate_boundary_3d<fe_degree, 2, transpose>(constraint_mask,
                                                           values);
        }
    }
  } // namespace internal
} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE
#endif

#endif

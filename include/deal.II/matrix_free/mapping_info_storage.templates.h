// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2021 by the deal.II authors
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

#ifndef dealii_matrix_free_mapping_info_storage_templates_h
#define dealii_matrix_free_mapping_info_storage_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/matrix_free/evaluation_template_factory.h>
#include <deal.II/matrix_free/mapping_info_storage.h>
#include <deal.II/matrix_free/task_info.h>
#include <deal.II/matrix_free/util.h>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace MatrixFreeFunctions
  {
    template <int structdim,
              int spacedim,
              typename Number,
              typename VectorizedArrayType>
    MappingInfoStorage<structdim, spacedim, Number, VectorizedArrayType>::
      QuadratureDescriptor::QuadratureDescriptor()
      : n_q_points(numbers::invalid_unsigned_int)
    {}



    template <int structdim,
              int spacedim,
              typename Number,
              typename VectorizedArrayType>
    template <int dim_q>
    void
    MappingInfoStorage<structdim, spacedim, Number, VectorizedArrayType>::
      QuadratureDescriptor::initialize(
        const Quadrature<dim_q> &quadrature,
        const UpdateFlags        update_flags_inner_faces)
    {
      Assert(structdim + 1 <= spacedim ||
               update_flags_inner_faces == update_default,
             ExcMessage("Volume cells do not allow for setting inner faces"));
      this->quadrature = quadrature;
      n_q_points       = quadrature.size();
      quadrature_weights.resize(n_q_points);
      for (unsigned int i = 0; i < n_q_points; ++i)
        quadrature_weights[i] = quadrature.weight(i);

      // note: quadrature_1d and tensor_quadrature_weights are not set up

      // TODO: set up face_orientations
      (void)update_flags_inner_faces;
    }



    template <int structdim,
              int spacedim,
              typename Number,
              typename VectorizedArrayType>
    void
    MappingInfoStorage<structdim, spacedim, Number, VectorizedArrayType>::
      QuadratureDescriptor::initialize(
        const Quadrature<1> &quadrature_1d,
        const UpdateFlags    update_flags_inner_faces)
    {
      Assert(structdim + 1 <= spacedim ||
               update_flags_inner_faces == update_default,
             ExcMessage("Volume cells do not allow for setting inner faces"));
      this->quadrature_1d = quadrature_1d;
      quadrature          = Quadrature<structdim>(quadrature_1d);
      n_q_points          = quadrature.size();
      quadrature_weights.resize(n_q_points);
      for (unsigned int i = 0; i < n_q_points; ++i)
        quadrature_weights[i] = quadrature.weight(i);

      for (int d = 0; d < structdim; ++d)
        {
          tensor_quadrature_weights[d].resize(quadrature_1d.size());
          for (unsigned int i = 0; i < quadrature_1d.size(); ++i)
            tensor_quadrature_weights[d][i] = quadrature_1d.weight(i);
        }

      // face orientation for faces in 3D
      if (structdim == spacedim - 1 && spacedim == 3 &&
          update_flags_inner_faces != update_default)
        {
          const unsigned int n = quadrature_1d.size();
          face_orientations.reinit(8, n * n);
          for (unsigned int j = 0, i = 0; j < n; ++j)
            for (unsigned int k = 0; k < n; ++k, ++i)
              {
                // face_orientation=true,  face_flip=false, face_rotation=false
                face_orientations[0][i] = i;
                // face_orientation=false, face_flip=false, face_rotation=false
                face_orientations[1][i] = j + k * n;
                // face_orientation=true,  face_flip=true,  face_rotation=false
                face_orientations[2][i] = (n - 1 - k) + (n - 1 - j) * n;
                // face_orientation=false, face_flip=true,  face_rotation=false
                face_orientations[3][i] = (n - 1 - j) + (n - 1 - k) * n;
                // face_orientation=true,  face_flip=false, face_rotation=true
                face_orientations[4][i] = j + (n - 1 - k) * n;
                // face_orientation=false, face_flip=false, face_rotation=true
                face_orientations[5][i] = k + (n - 1 - j) * n;
                // face_orientation=true,  face_flip=true,  face_rotation=true
                face_orientations[6][i] = (n - 1 - j) + k * n;
                // face_orientation=false, face_flip=true,  face_rotation=true
                face_orientations[7][i] = (n - 1 - k) + j * n;
              }
        }
    }



    template <int structdim,
              int spacedim,
              typename Number,
              typename VectorizedArrayType>
    std::size_t
    MappingInfoStorage<structdim, spacedim, Number, VectorizedArrayType>::
      QuadratureDescriptor::memory_consumption() const
    {
      std::size_t memory = sizeof(this) + quadrature.memory_consumption() +
                           quadrature_weights.memory_consumption() +
                           face_orientations.memory_consumption();
      for (int d = 0; d < structdim; ++d)
        memory += tensor_quadrature_weights[d].memory_consumption();
      return memory;
    }



    template <int structdim,
              int spacedim,
              typename Number,
              typename VectorizedArrayType>
    void
    MappingInfoStorage<structdim, spacedim, Number, VectorizedArrayType>::
      clear_data_fields()
    {
      data_index_offsets.clear();
      JxW_values.clear();
      normal_vectors.clear();
      for (unsigned int i = 0; i < 2; ++i)
        {
          jacobians[i].clear();
          jacobian_gradients[i].clear();
          normals_times_jacobians[i].clear();
        }
      quadrature_point_offsets.clear();
      quadrature_points.clear();
    }



    template <int structdim,
              int spacedim,
              typename Number,
              typename VectorizedArrayType>
    std::size_t
    MappingInfoStorage<structdim, spacedim, Number, VectorizedArrayType>::
      memory_consumption() const
    {
      return MemoryConsumption::memory_consumption(descriptor) +
             MemoryConsumption::memory_consumption(data_index_offsets) +
             MemoryConsumption::memory_consumption(JxW_values) +
             MemoryConsumption::memory_consumption(normal_vectors) +
             MemoryConsumption::memory_consumption(jacobians[0]) +
             MemoryConsumption::memory_consumption(jacobians[1]) +
             MemoryConsumption::memory_consumption(jacobian_gradients[0]) +
             MemoryConsumption::memory_consumption(jacobian_gradients[1]) +
             MemoryConsumption::memory_consumption(normals_times_jacobians[0]) +
             MemoryConsumption::memory_consumption(normals_times_jacobians[1]) +
             MemoryConsumption::memory_consumption(quadrature_point_offsets) +
             MemoryConsumption::memory_consumption(quadrature_points);
    }



    template <int structdim,
              int spacedim,
              typename Number,
              typename VectorizedArrayType>
    template <typename StreamType>
    void
    MappingInfoStorage<structdim, spacedim, Number, VectorizedArrayType>::
      print_memory_consumption(StreamType &out, const TaskInfo &task_info) const
    {
      // print_memory_statistics involves global communication, so we can
      // disable the check here only if no processor has any such data
      const std::size_t size =
        Utilities::MPI::sum(jacobians[0].size(), task_info.communicator);
      if (size > 0)
        {
          out << "      Memory JxW data:               ";
          task_info.print_memory_statistics(
            out,
            MemoryConsumption::memory_consumption(data_index_offsets) +
              MemoryConsumption::memory_consumption(JxW_values));
          out << "      Memory Jacobian data:          ";
          task_info.print_memory_statistics(
            out,
            MemoryConsumption::memory_consumption(jacobians[0]) +
              MemoryConsumption::memory_consumption(jacobians[1]));
          out << "      Memory second derivative data: ";
          task_info.print_memory_statistics(
            out,
            MemoryConsumption::memory_consumption(jacobian_gradients[0]) +
              MemoryConsumption::memory_consumption(jacobian_gradients[1]));
        }
      const std::size_t normal_size =
        Utilities::MPI::sum(normal_vectors.size(), task_info.communicator);
      if (normal_size > 0)
        {
          out << "      Memory normal vectors data:    ";
          task_info.print_memory_statistics(
            out,
            MemoryConsumption::memory_consumption(normal_vectors) +
              MemoryConsumption::memory_consumption(
                normals_times_jacobians[0]) +
              MemoryConsumption::memory_consumption(
                normals_times_jacobians[1]));
        }

      const std::size_t quad_size =
        Utilities::MPI::sum(quadrature_points.size(), task_info.communicator);
      if (quad_size > 0)
        {
          out << "      Memory quadrature points:      ";
          task_info.print_memory_statistics(
            out,
            MemoryConsumption::memory_consumption(quadrature_point_offsets) +
              MemoryConsumption::memory_consumption(quadrature_points));
        }
    }

  } // namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif

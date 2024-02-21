// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_matrix_free_mapping_info_storage_templates_h
#define dealii_matrix_free_mapping_info_storage_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/matrix_free/evaluation_template_factory.h>
#include <deal.II/matrix_free/mapping_info_storage.h>
#include <deal.II/matrix_free/task_info.h>
#include <deal.II/matrix_free/util.h>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace MatrixFreeFunctions
  {
    template <int structdim, int spacedim, typename Number>
    MappingInfoStorage<structdim, spacedim, Number>::QuadratureDescriptor::
      QuadratureDescriptor()
      : n_q_points(numbers::invalid_unsigned_int)
    {}



    template <int structdim, int spacedim, typename Number>
    template <int dim_q>
    void
    MappingInfoStorage<structdim, spacedim, Number>::QuadratureDescriptor::
      initialize(const Quadrature<dim_q> &quadrature)
    {
      this->quadrature = quadrature;
      n_q_points       = quadrature.size();
      quadrature_weights.resize(n_q_points);
      for (unsigned int i = 0; i < n_q_points; ++i)
        quadrature_weights[i] = quadrature.weight(i);

      // note: quadrature_1d and tensor_quadrature_weights are not set up
    }



    template <int structdim, int spacedim, typename Number>
    void
    MappingInfoStorage<structdim, spacedim, Number>::QuadratureDescriptor::
      initialize(const Quadrature<1> &quadrature_1d)
    {
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
    }



    template <int structdim, int spacedim, typename Number>
    std::size_t
    MappingInfoStorage<structdim, spacedim, Number>::QuadratureDescriptor::
      memory_consumption() const
    {
      std::size_t memory = sizeof(this) + quadrature.memory_consumption() +
                           quadrature_weights.memory_consumption();
      for (int d = 0; d < structdim; ++d)
        memory += tensor_quadrature_weights[d].memory_consumption();
      return memory;
    }



    template <int structdim, int spacedim, typename Number>
    void
    MappingInfoStorage<structdim, spacedim, Number>::clear_data_fields()
    {
      data_index_offsets.clear();
      JxW_values.clear();
      normal_vectors.clear();
      for (unsigned int i = 0; i < 2; ++i)
        {
          jacobians[i].clear();
          jacobian_gradients[i].clear();
          jacobian_gradients_non_inverse[i].clear();
          normals_times_jacobians[i].clear();
        }
      quadrature_point_offsets.clear();
      quadrature_points.clear();
    }



    template <int structdim, int spacedim, typename Number>
    UpdateFlags
    MappingInfoStorage<structdim, spacedim, Number>::compute_update_flags(
      const UpdateFlags                                     update_flags,
      const std::vector<dealii::hp::QCollection<spacedim>> &quads,
      const bool                                            piola_transform)
    {
      // this class is build around the evaluation of jacobians, so compute
      // them in any case. The Jacobians will be inverted manually. Since we
      // always do support integration, we also include the JxW values
      UpdateFlags new_flags = update_jacobians | update_JxW_values;

      // for Hessian information, need inverse Jacobians and the derivative of
      // Jacobians (these two together will give use the gradients of the
      // inverse Jacobians, which is what we need)
      if ((update_flags & update_hessians) != 0u ||
          (update_flags & update_jacobian_grads) != 0u ||
          (piola_transform &&
           ((update_flags &
             (update_gradients | update_contravariant_transformation)) != 0u)))
        new_flags |= update_jacobian_grads;

      if ((update_flags & update_quadrature_points) != 0u)
        new_flags |= update_quadrature_points;

      // there is one more thing: if we have a quadrature formula with only
      // one quadrature point on the first component, but more points on later
      // components, we need to have Jacobian gradients anyway in order to
      // determine whether the Jacobian is constant throughout a cell
      if (quads.empty() == false)
        {
          bool formula_with_one_point = false;
          for (unsigned int i = 0; i < quads[0].size(); ++i)
            if (quads[0][i].size() == 1)
              {
                formula_with_one_point = true;
                break;
              }
          if (formula_with_one_point == true)
            for (unsigned int comp = 1; comp < quads.size(); ++comp)
              for (unsigned int i = 0; i < quads[comp].size(); ++i)
                if (quads[comp][i].size() > 1)
                  {
                    new_flags |= update_jacobian_grads;
                  }
        }
      return new_flags;
    }



    template <int structdim, int spacedim, typename Number>
    std::size_t
    MappingInfoStorage<structdim, spacedim, Number>::memory_consumption() const
    {
      return MemoryConsumption::memory_consumption(descriptor) +
             MemoryConsumption::memory_consumption(data_index_offsets) +
             MemoryConsumption::memory_consumption(JxW_values) +
             MemoryConsumption::memory_consumption(normal_vectors) +
             MemoryConsumption::memory_consumption(jacobians[0]) +
             MemoryConsumption::memory_consumption(jacobians[1]) +
             MemoryConsumption::memory_consumption(jacobian_gradients[0]) +
             MemoryConsumption::memory_consumption(jacobian_gradients[1]) +
             MemoryConsumption::memory_consumption(
               jacobian_gradients_non_inverse[0]) +
             MemoryConsumption::memory_consumption(
               jacobian_gradients_non_inverse[1]) +
             MemoryConsumption::memory_consumption(normals_times_jacobians[0]) +
             MemoryConsumption::memory_consumption(normals_times_jacobians[1]) +
             MemoryConsumption::memory_consumption(quadrature_point_offsets) +
             MemoryConsumption::memory_consumption(quadrature_points);
    }



    template <int structdim, int spacedim, typename Number>
    template <typename StreamType>
    void
    MappingInfoStorage<structdim, spacedim, Number>::print_memory_consumption(
      StreamType     &out,
      const TaskInfo &task_info) const
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
              MemoryConsumption::memory_consumption(jacobian_gradients[1]) +
              MemoryConsumption::memory_consumption(
                jacobian_gradients_non_inverse[0]) +
              MemoryConsumption::memory_consumption(
                jacobian_gradients_non_inverse[1]));
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

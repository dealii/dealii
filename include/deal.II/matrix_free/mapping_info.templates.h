// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2018 by the deal.II authors
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

#ifndef dealii_matrix_free_mapping_info_templates_h
#define dealii_matrix_free_mapping_info_templates_h

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/matrix_free/mapping_info.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MatrixFreeFunctions
  {
    /* ------------------------ MappingInfoStorage implementation ---------- */

    template <int structdim, int spacedim, typename Number>
    MappingInfoStorage<structdim, spacedim, Number>::QuadratureDescriptor ::
      QuadratureDescriptor()
      : n_q_points(numbers::invalid_unsigned_int)
    {}

    template <int structdim, int spacedim, typename Number>
    void
    MappingInfoStorage<structdim, spacedim, Number>::QuadratureDescriptor ::
      initialize(const Quadrature<1>& quadrature_1d,
                 const UpdateFlags    update_flags_inner_faces)
    {
      Assert(structdim + 1 <= spacedim
               || update_flags_inner_faces == update_default,
             ExcMessage("Volume cells do not allow for setting inner faces"));
      quadrature = Quadrature<structdim>(quadrature_1d);
      n_q_points = quadrature.size();
      quadrature_weights.resize(n_q_points);
      for(unsigned int i = 0; i < n_q_points; ++i)
        quadrature_weights[i] = quadrature.weight(i);

      for(unsigned int d = 0; d < structdim; ++d)
        {
          tensor_quadrature_weights[d].resize(quadrature_1d.size());
          for(unsigned int i = 0; i < quadrature_1d.size(); ++i)
            tensor_quadrature_weights[d][i] = quadrature_1d.weight(i);
        }

      // face orientation for faces in 3D
      if(structdim == spacedim - 1 && spacedim == 3
         && update_flags_inner_faces != update_default)
        {
          const unsigned int n = quadrature_1d.size();
          face_orientations.reinit(8, n * n);
          for(unsigned int j = 0, i = 0; j < n; ++j)
            for(unsigned int k = 0; k < n; ++k, ++i)
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

    template <int structdim, int spacedim, typename Number>
    std::size_t
    MappingInfoStorage<structdim, spacedim, Number>::QuadratureDescriptor ::
      memory_consumption() const
    {
      std::size_t memory = sizeof(this) + quadrature.memory_consumption()
                           + quadrature_weights.memory_consumption()
                           + face_orientations.memory_consumption();
      for(unsigned int d = 0; d < structdim; ++d)
        memory += tensor_quadrature_weights[d].memory_consumption();
      return memory;
    }

    template <int structdim, int spacedim, typename Number>
    std::size_t
    MappingInfoStorage<structdim, spacedim, Number>::memory_consumption() const
    {
      return MemoryConsumption::memory_consumption(descriptor)
             + MemoryConsumption::memory_consumption(data_index_offsets)
             + MemoryConsumption::memory_consumption(JxW_values)
             + MemoryConsumption::memory_consumption(normal_vectors)
             + MemoryConsumption::memory_consumption(jacobians[0])
             + MemoryConsumption::memory_consumption(jacobians[1])
             + MemoryConsumption::memory_consumption(jacobian_gradients[0])
             + MemoryConsumption::memory_consumption(jacobian_gradients[1])
             + MemoryConsumption::memory_consumption(normals_times_jacobians[0])
             + MemoryConsumption::memory_consumption(normals_times_jacobians[1])
             + MemoryConsumption::memory_consumption(quadrature_point_offsets)
             + MemoryConsumption::memory_consumption(quadrature_points);
    }

    template <int structdim, int spacedim, typename Number>
    template <typename StreamType>
    void
    MappingInfoStorage<structdim, spacedim, Number>::print_memory_consumption(
      StreamType&     out,
      const SizeInfo& task_info) const
    {
      // print_memory_statistics involves global communication, so we can
      // disable the check here only if no processor has any such data
      const std::size_t size
        = Utilities::MPI::sum(jacobians[0].size(), task_info.communicator);
      if(size > 0)
        {
          out << "      Memory JxW data:               ";
          task_info.print_memory_statistics(
            out,
            MemoryConsumption::memory_consumption(data_index_offsets)
              + MemoryConsumption::memory_consumption(JxW_values));
          out << "      Memory Jacobian data:          ";
          task_info.print_memory_statistics(
            out,
            MemoryConsumption::memory_consumption(jacobians[0])
              + MemoryConsumption::memory_consumption(jacobians[1]));
          out << "      Memory second derivative data: ";
          task_info.print_memory_statistics(
            out,
            MemoryConsumption::memory_consumption(jacobian_gradients[0])
              + MemoryConsumption::memory_consumption(jacobian_gradients[1]));
        }
      const std::size_t normal_size
        = Utilities::MPI::sum(normal_vectors.size(), task_info.communicator);
      if(normal_size > 0)
        {
          out << "      Memory normal vectors data:    ";
          task_info.print_memory_statistics(
            out,
            MemoryConsumption::memory_consumption(normal_vectors)
              + MemoryConsumption::memory_consumption(
                  normals_times_jacobians[0])
              + MemoryConsumption::memory_consumption(
                  normals_times_jacobians[1]));
        }

      const std::size_t quad_size
        = Utilities::MPI::sum(quadrature_points.size(), task_info.communicator);
      if(quad_size > 0)
        {
          out << "      Memory quadrature points:      ";
          task_info.print_memory_statistics(
            out,
            MemoryConsumption::memory_consumption(quadrature_point_offsets)
              + MemoryConsumption::memory_consumption(quadrature_points));
        }
    }

    /* ------------------------ MappingInfo implementation ----------------- */

    template <int dim, typename Number>
    MappingInfo<dim, Number>::MappingInfo()
    {}

    template <int dim, typename Number>
    void
    MappingInfo<dim, Number>::clear()
    {
      cell_data.clear();
      face_data.clear();
      face_data_by_cells.clear();
      cell_type.clear();
      face_type.clear();
    }

    template <int dim, typename Number>
    UpdateFlags
    MappingInfo<dim, Number>::compute_update_flags(
      const UpdateFlags                              update_flags,
      const std::vector<dealii::hp::QCollection<1>>& quad)
    {
      // this class is build around the evaluation of jacobians, so compute
      // them in any case. The Jacobians will be inverted manually. Since we
      // always do support integration, we also include the JxW values
      UpdateFlags new_flags = update_jacobians | update_JxW_values;

      // for Hessian information, need inverse Jacobians and the derivative of
      // Jacobians (these two together will give use the gradients of the
      // inverse Jacobians, which is what we need)
      if(update_flags & update_hessians || update_flags & update_jacobian_grads)
        new_flags |= update_jacobian_grads;

      if(update_flags & update_quadrature_points)
        new_flags |= update_quadrature_points;

      // there is one more thing: if we have a quadrature formula with only
      // one quadrature point on the first component, but more points on later
      // components, we need to have Jacobian gradients anyway in order to
      // determine whether the Jacobian is constant throughout a cell
      if(quad.empty() == false)
        {
          bool formula_with_one_point = false;
          for(unsigned int i = 0; i < quad[0].size(); ++i)
            if(quad[0][i].size() == 1)
              {
                formula_with_one_point = true;
                break;
              }
          if(formula_with_one_point == true)
            for(unsigned int comp = 1; comp < quad.size(); ++comp)
              for(unsigned int i = 0; i < quad[comp].size(); ++i)
                if(quad[comp][i].size() > 1)
                  {
                    new_flags |= update_jacobian_grads;
                  }
        }
      return new_flags;
    }

    template <int dim, typename Number>
    void
    MappingInfo<dim, Number>::initialize(
      const dealii::Triangulation<dim>&                          tria,
      const std::vector<std::pair<unsigned int, unsigned int>>&  cells,
      const FaceInfo<VectorizedArray<Number>::n_array_elements>& face_info,
      const std::vector<unsigned int>&               active_fe_index,
      const Mapping<dim>&                            mapping,
      const std::vector<dealii::hp::QCollection<1>>& quad,
      const UpdateFlags                              update_flags_cells,
      const UpdateFlags update_flags_boundary_faces,
      const UpdateFlags update_flags_inner_faces,
      const UpdateFlags update_flags_faces_by_cells)
    {
      clear();

      // Could call these functions in parallel, but not useful because the
      // work inside is nicely split up already
      initialize_cells(
        tria, cells, active_fe_index, mapping, quad, update_flags_cells);
      initialize_faces(tria,
                       cells,
                       face_info.faces,
                       mapping,
                       quad,
                       update_flags_boundary_faces,
                       update_flags_inner_faces);
      initialize_faces_by_cells(
        tria, cells, mapping, quad, update_flags_faces_by_cells);
    }

    /* ------------------------- initialization of cells ------------------- */

    // Anonymous namespace with implementation of extraction of values on cell
    // range
    namespace
    {
      template <int dim>
      double
      get_jacobian_size(const dealii::Triangulation<dim>& tria)
      {
        if(tria.n_cells() == 0)
          return 1;
        else
          return tria.begin()->diameter();
      }

      template <int dim, typename Number>
      struct CompressedCellData
      {
        CompressedCellData(const double expected_size)
          : data(FPArrayComparator<Number>(expected_size))
        {}

        std::map<
          Tensor<2,
                 dim,
                 Tensor<1, VectorizedArray<Number>::n_array_elements, Number>>,
          unsigned int,
          FPArrayComparator<Number>>
          data;
      };

      /**
       * Internal temporary data used for the initialization.
       */
      template <int dim, typename Number>
      struct LocalData
      {
        LocalData(const double jac_size);
        void
        resize(const unsigned int size);

        AlignedVector<Point<dim, VectorizedArray<Number>>> quadrature_points;
        AlignedVector<Tensor<2, dim, VectorizedArray<Number>>> general_jac;
        AlignedVector<VectorizedArray<Number>>                 JxW_values;
        AlignedVector<Tensor<3, dim, VectorizedArray<Number>>> general_jac_grad;
        AlignedVector<Tensor<1, dim, VectorizedArray<Number>>> normal_vectors;
        Tensor<2, dim, VectorizedArray<Number>>                const_jac;
        const double                                           jac_size;
      };

      template <int dim, typename Number>
      LocalData<dim, Number>::LocalData(const double jac_size_in)
        : jac_size(jac_size_in)
      {}

      template <int dim, typename Number>
      void
      LocalData<dim, Number>::resize(const unsigned int size)
      {
        if(JxW_values.size() != size)
          {
            quadrature_points.resize_fast(size);
            general_jac.resize_fast(size * 2);
            JxW_values.resize_fast(size);
            general_jac_grad.resize_fast(size * 2);
            normal_vectors.resize_fast(size);
          }
      }

      /**
       * Helper function called internally during the initialize function.
       */
      template <int dim, typename Number>
      void
      evaluate_on_cell(
        const dealii::Triangulation<dim>&            tria,
        const std::pair<unsigned int, unsigned int>* cells,
        const unsigned int                           my_q,
        GeometryType&                                cell_t_prev,
        GeometryType (&cell_t)[VectorizedArray<Number>::n_array_elements],
        dealii::FEValues<dim, dim>& fe_val,
        LocalData<dim, Number>&     cell_data)
      {
        const unsigned int n_q_points   = fe_val.n_quadrature_points;
        const UpdateFlags  update_flags = fe_val.get_update_flags();

        cell_data.const_jac = Tensor<2, dim, VectorizedArray<Number>>();

        // this should be the same value as used in HashValue::scaling (but we
        // not have that field here)
        const double zero_tolerance_double
          = cell_data.jac_size * std::numeric_limits<double>::epsilon() * 1024.;
        for(unsigned int j = 0; j < VectorizedArray<Number>::n_array_elements;
            ++j)
          {
            typename dealii::Triangulation<dim>::cell_iterator cell_it(
              &tria, cells[j].first, cells[j].second);
            fe_val.reinit(cell_it);
            cell_t[j] = general;

            // extract quadrature points and store them temporarily. if we have
            // Cartesian cells, we can compress the indices
            if(update_flags & update_quadrature_points)
              for(unsigned int q = 0; q < n_q_points; ++q)
                {
                  const Point<dim>& point = fe_val.quadrature_point(q);
                  for(unsigned int d = 0; d < dim; ++d)
                    cell_data.quadrature_points[q][d][j] = point[d];
                }

            // if this is not the first quadrature formula and we already have
            // determined that this cell is either Cartesian or with constant
            // Jacobian, we have nothing more to do.
            if(my_q > 0 && cell_t_prev <= affine)
              continue;

            // first round: if the transformation is detected to be the same as
            // on the old cell, we only need to copy over the data.
            if(fe_val.get_cell_similarity() == CellSimilarity::translation
               && my_q == 0)
              {
                if(j == 0)
                  cell_t[j] = cell_t_prev;
                else
                  cell_t[j] = cell_t[j - 1];
              }

            const DerivativeForm<1, dim, dim>& jac_0 = fe_val.jacobian(0);

            if(my_q == 0)
              {
                // check whether the Jacobian is constant on this cell the first
                // time we come around here
                if(cell_t[j] == general)
                  {
                    bool jacobian_constant = true;
                    for(unsigned int q = 1; q < n_q_points; ++q)
                      {
                        const DerivativeForm<1, dim, dim>& jac
                          = fe_val.jacobian(q);
                        for(unsigned int d = 0; d < dim; ++d)
                          for(unsigned int e = 0; e < dim; ++e)
                            if(std::fabs(jac_0[d][e] - jac[d][e])
                               > zero_tolerance_double)
                              jacobian_constant = false;
                        if(jacobian_constant == false)
                          break;
                      }

                    // check whether the Jacobian is diagonal to machine
                    // accuracy
                    bool cell_cartesian = jacobian_constant;
                    for(unsigned int d = 0; d < dim; ++d)
                      for(unsigned int e = 0; e < dim; ++e)
                        if(d != e)
                          if(std::fabs(jac_0[d][e]) > zero_tolerance_double)
                            {
                              cell_cartesian = false;
                              break;
                            }

                    // in case we have only one quadrature point, we can have
                    // non-constant Jacobians, but we cannot detect it by
                    // comparison from one quadrature point to the next: in that
                    // case, need to look at second derivatives and see whether
                    // there are some non-zero entries (this is necessary since
                    // we determine the constness of the Jacobian for the first
                    // quadrature formula and might not look at them any more
                    // for the second, third quadrature formula). in any case,
                    // the flag update_jacobian_grads will be set in that case
                    if(cell_cartesian == false && n_q_points == 1
                       && update_flags & update_jacobian_grads)
                      {
                        const DerivativeForm<1, dim, dim>& jac
                          = fe_val.jacobian(0);
                        const DerivativeForm<2, dim, dim>& jacobian_grad
                          = fe_val.jacobian_grad(0);
                        for(unsigned int d = 0; d < dim; ++d)
                          for(unsigned int e = 0; e < dim; ++e)
                            for(unsigned int f = 0; f < dim; ++f)
                              {
                                double jac_grad_comp
                                  = (jac[f][0] * jacobian_grad[d][e][0]);
                                for(unsigned int g = 1; g < dim; ++g)
                                  jac_grad_comp
                                    += (jac[f][g] * jacobian_grad[d][e][g]);
                                if(std::fabs(jac_grad_comp)
                                   > zero_tolerance_double)
                                  jacobian_constant = false;
                              }
                      }
                    // set cell type
                    if(cell_cartesian == true)
                      cell_t[j] = cartesian;
                    else if(jacobian_constant == true)
                      cell_t[j] = affine;
                    else
                      cell_t[j] = general;
                  }

                // Cartesian cell
                if(cell_t[j] == cartesian)
                  {
                    // set Jacobian into diagonal (off-diagonal part is already
                    // zeroed out)
                    for(unsigned int d = 0; d < dim; ++d)
                      cell_data.const_jac[d][d][j] = jac_0[d][d];
                    continue;
                  }

                // cell with affine mapping
                else if(cell_t[j] == affine)
                  {
                    // compress out very small values
                    for(unsigned int d = 0; d < dim; ++d)
                      for(unsigned int e = 0; e < dim; ++e)
                        if(std::fabs(jac_0[d][e]))
                          cell_data.const_jac[d][e][j] = jac_0[d][e];
                    continue;
                  }
              }

            // general cell case

            // go through all quadrature points and fill in the data into the
            // temporary data structures with slots for the vectorized data
            // types
            for(unsigned int q = 0; q < n_q_points; ++q)
              {
                // compress out very small numbers which are only noise. Then it
                // is cleaner to use zero straight away (though it does not save
                // any memory)
                const DerivativeForm<1, dim, dim>& jac = fe_val.jacobian(q);
                for(unsigned int d = 0; d < dim; ++d)
                  for(unsigned int e = 0; e < dim; ++e)
                    cell_data.general_jac[q][d][e][j]
                      = std::fabs(jac[d][e]) < zero_tolerance_double ?
                          0. :
                          jac[d][e];

                // need to do some calculus based on the gradient of the
                // Jacobian, in order to find the gradient of the inverse
                // Jacobian which is needed in user code. however, we would like
                // to perform that on vectorized data types instead of doubles
                // or floats. to this end, copy the gradients first
                if(update_flags & update_jacobian_grads)
                  {
                    const DerivativeForm<2, dim, dim>& jacobian_grad
                      = fe_val.jacobian_grad(q);
                    for(unsigned int d = 0; d < dim; ++d)
                      for(unsigned int e = 0; e < dim; ++e)
                        for(unsigned int f = 0; f < dim; ++f)
                          cell_data.general_jac_grad[q][d][e][f][j]
                            = jacobian_grad[d][e][f];
                  }
              }
          } // end loop over entries of vectorization (n_array_elements cells)

        // set information for next cell
        cell_t_prev = cell_t[VectorizedArray<Number>::n_array_elements - 1];
      }

      template <int dim, typename Number>
      void
      initialize_cell_range(
        const std::pair<unsigned int, unsigned int>               cell_range,
        const dealii::Triangulation<dim>&                         tria,
        const std::vector<std::pair<unsigned int, unsigned int>>& cells,
        const std::vector<unsigned int>&               active_fe_index,
        const Mapping<dim>&                            mapping,
        const std::vector<dealii::hp::QCollection<1>>& quad,
        const UpdateFlags                              update_flags,
        MappingInfo<dim, Number>&                      mapping_info,
        std::pair<std::vector<MappingInfoStorage<dim, dim, Number>>,
                  CompressedCellData<dim, Number>>&    data)
      {
        FE_Nothing<dim> dummy_fe;

        Tensor<3, dim, VectorizedArray<Number>> jac_grad, grad_jac_inv;
        Tensor<1, dim, VectorizedArray<Number>> tmp;

        // when we make comparisons about the size of Jacobians we need to
        // know the approximate size of typical entries in Jacobians. We need
        // to fix the Jacobian size once and for all. We choose the diameter
        // of the first cell (on level zero, which is the best accuracy we can
        // hope for, since diameters on finer levels are computed by
        // differences of nearby cells) as the order of magnitude by which we
        // make comparisons "relative."
        const double jacobian_size = get_jacobian_size(tria);

        // objects that hold the data for up to vectorization_width cells while
        // we fill them up. Only after all vectorization_width cells have been
        // processed, we can insert the data into the data structures of this
        // class
        LocalData<dim, Number> cell_data(jacobian_size);

        // encodes the cell types of the current cell. Since several cells
        // must be considered together, this variable holds the individual
        // info of the last chunk of cells
        GeometryType cell_t[VectorizedArray<Number>::n_array_elements];
        GeometryType cell_t_prev = general;

        // fe_values object that is used to compute the mapping data. for
        // the hp case there might be more than one finite element. since we
        // manually select the active FE index and not via a
        // hp::DoFHandler<dim>::active_cell_iterator, we need to manually
        // select the correct finite element, so just hold a vector of
        // FEValues
        std::vector<std::vector<std::shared_ptr<dealii::FEValues<dim>>>>
          fe_values(mapping_info.cell_data.size());
        for(unsigned int i = 0; i < fe_values.size(); ++i)
          fe_values[i].resize(mapping_info.cell_data[i].descriptor.size());
        UpdateFlags update_flags_feval
          = (update_flags & update_jacobians ? update_jacobians :
                                               update_default)
            | (update_flags & update_jacobian_grads ? update_jacobian_grads :
                                                      update_default)
            | (update_flags & update_quadrature_points ?
                 update_quadrature_points :
                 update_default);

        std::vector<std::vector<unsigned int>> n_q_points_1d(quad.size()),
          step_size_cartesian(quad.size());
        for(unsigned int my_q = 0; my_q < quad.size(); ++my_q)
          {
            n_q_points_1d[my_q].resize(quad[my_q].size());
            step_size_cartesian[my_q].resize(quad[my_q].size());
            for(unsigned int hpq = 0; hpq < quad[my_q].size(); ++hpq)
              {
                n_q_points_1d[my_q][hpq] = quad[my_q][hpq].size();

                // To walk on the diagonal for lexicographic ordering, we have
                // to jump one index ahead in each direction. For direction 0,
                // this is just the next point, for direction 1, it means adding
                // n_q_points_1d, and so on.
                step_size_cartesian[my_q][hpq] = 0;
                unsigned int factor            = 1;
                for(unsigned int d = 0; d < dim; ++d)
                  {
                    step_size_cartesian[my_q][hpq] += factor;
                    factor *= n_q_points_1d[my_q][hpq];
                  }
              }
          }

        const unsigned int end_cell = std::min(mapping_info.cell_type.size(),
                                               std::size_t(cell_range.second));
        // loop over given cells
        for(unsigned int cell = cell_range.first; cell < end_cell; ++cell)
          for(unsigned int my_q = 0; my_q < mapping_info.cell_data.size();
              ++my_q)
            {
              // GENERAL OUTLINE: First generate the data in format "number"
              // for vectorization_width cells, and then find the most
              // general type of cell for appropriate vectorized formats. then
              // fill this data in
              const unsigned int fe_index
                = active_fe_index.size() > 0 ? active_fe_index[cell] : 0;
              const unsigned int n_q_points
                = mapping_info.cell_data[my_q].descriptor[fe_index].n_q_points;
              if(fe_values[my_q][fe_index].get() == nullptr)
                fe_values[my_q][fe_index].reset(new dealii::FEValues<dim>(
                  mapping,
                  dummy_fe,
                  mapping_info.cell_data[my_q].descriptor[fe_index].quadrature,
                  update_flags_feval));
              dealii::FEValues<dim>& fe_val = *fe_values[my_q][fe_index];
              cell_data.resize(n_q_points);

              // if the fe index has changed from the previous cell, set the
              // old cell type to invalid (otherwise, we might detect
              // similarity due to some cells further ahead)
              if(my_q > 0)
                cell_t_prev = GeometryType(mapping_info.cell_type[cell]);
              else if(cell > cell_range.first && active_fe_index.size() > 0
                      && active_fe_index[cell] != active_fe_index[cell - 1])
                cell_t_prev = general;

              evaluate_on_cell(
                tria,
                &cells[cell * VectorizedArray<Number>::n_array_elements],
                my_q,
                cell_t_prev,
                cell_t,
                fe_val,
                cell_data);

              // now reorder the data into vectorized types. if we are here
              // for the first time, we need to find out whether the Jacobian
              // allows for some simplification (Cartesian, affine) taking
              // vectorization_width cell together

              if(my_q == 0)
                {
                  // find the most general cell type (most general type is 3
                  // (general cell))
                  GeometryType most_general_type = cartesian;
                  for(unsigned int j = 0;
                      j < VectorizedArray<Number>::n_array_elements;
                      ++j)
                    if(cell_t[j] > most_general_type)
                      most_general_type = cell_t[j];
                  AssertIndexRange((unsigned int) most_general_type, 4U);
                  mapping_info.cell_type[cell] = most_general_type;
                }

              AssertThrow(
                data.first[my_q].JxW_values.size()
                  < static_cast<std::size_t>(
                      std::numeric_limits<unsigned int>::max()),
                ExcMessage(
                  "Index overflow. Cannot fit data in 32 bit integers"));

              unsigned int insert_position = data.first[my_q].JxW_values.size();
              // Cartesian/affine cell with constant Jacobians throughout the
              // cell. We need to store the data in another data field because
              // std::map cannot store data based on VectorizedArray directly
              // (alignment issue).
              if(mapping_info.cell_type[cell] <= affine)
                {
                  if(my_q == 0)
                    {
                      std::pair<
                        Tensor<2,
                               dim,
                               Tensor<1,
                                      VectorizedArray<Number>::n_array_elements,
                                      Number>>,
                        unsigned int>
                        new_entry;
                      // This number overlaps with the general data but we
                      // take care of that when we merge data from different
                      // threads
                      new_entry.second = data.second.data.size();
                      for(unsigned int d = 0; d < dim; ++d)
                        for(unsigned int e = 0; e < dim; ++e)
                          for(unsigned int v = 0;
                              v < VectorizedArray<Number>::n_array_elements;
                              ++v)
                            new_entry.first[d][e][v]
                              = cell_data.const_jac[d][e][v];

                      insert_position
                        = data.second.data.insert(new_entry).first->second;
                    }
                  else
                    insert_position
                      = data.first[0]
                          .data_index_offsets[cell - cell_range.first];
                }

              // general cell case: now go through all quadrature points and
              // collect the data. done for all different quadrature formulas,
              // so do it outside the above loop.
              data.first[my_q].data_index_offsets.push_back(insert_position);
              if(mapping_info.get_cell_type(cell) == general)
                {
                  for(unsigned int q = 0; q < n_q_points; ++q)
                    {
                      Tensor<2, dim, VectorizedArray<Number>>& jac
                        = cell_data.general_jac[q];
                      Tensor<3, dim, VectorizedArray<Number>>& jacobian_grad
                        = cell_data.general_jac_grad[q];
                      for(unsigned int j = 0;
                          j < VectorizedArray<Number>::n_array_elements;
                          ++j)
                        if(cell_t[j] < general)
                          {
                            for(unsigned int d = 0; d < dim; ++d)
                              for(unsigned int e = 0; e < dim; ++e)
                                {
                                  jac[d][e][j] = cell_data.const_jac[d][e][j];
                                  for(unsigned int f = 0; f < dim; ++f)
                                    jacobian_grad[d][e][f][j] = 0.;
                                }
                          }

                      data.first[my_q].JxW_values.push_back(
                        determinant(jac) * fe_val.get_quadrature().weight(q));
                      Tensor<2, dim, VectorizedArray<Number>> inv_jac
                        = transpose(invert(jac));
                      data.first[my_q].jacobians[0].push_back(inv_jac);

                      if(update_flags & update_jacobian_grads)
                        {
                          // for second derivatives on the real cell, need
                          // also the gradient of the inverse Jacobian J. This
                          // involves some calculus and is done
                          // vectorized. This is very cheap compared to what
                          // fe_values does (in early 2011). If L is the
                          // gradient of the jacobian on the unit cell, the
                          // gradient of the inverse is given by
                          // (multidimensional calculus) - J * (J * L) * J
                          // (the third J is because we need to transform the
                          // gradient L from the unit to the real cell, and
                          // then apply the inverse Jacobian). Compare this
                          // with 1D with j(x) = 1/k(phi(x)), where j = phi'
                          // is the inverse of the jacobian and k is the
                          // derivative of the jacobian on the unit cell. Then
                          // j' = phi' k'/k^2 = j k' j^2.

                          // compute: jac_grad = J*grad_unit(J^-1)
                          for(unsigned int d = 0; d < dim; ++d)
                            for(unsigned int e = 0; e < dim; ++e)
                              for(unsigned int f = 0; f < dim; ++f)
                                {
                                  jac_grad[f][e][d]
                                    = (inv_jac[f][0] * jacobian_grad[d][e][0]);
                                  for(unsigned int g = 1; g < dim; ++g)
                                    jac_grad[f][e][d]
                                      += (inv_jac[f][g]
                                          * jacobian_grad[d][e][g]);
                                }

                          // compute: transpose (-jac * jac_grad[d] * jac)
                          for(unsigned int d = 0; d < dim; ++d)
                            for(unsigned int e = 0; e < dim; ++e)
                              {
                                for(unsigned int f = 0; f < dim; ++f)
                                  {
                                    tmp[f] = VectorizedArray<Number>();
                                    for(unsigned int g = 0; g < dim; ++g)
                                      tmp[f]
                                        -= jac_grad[d][f][g] * inv_jac[g][e];
                                  }

                                // needed for non-diagonal part of Jacobian
                                // grad
                                for(unsigned int f = 0; f < dim; ++f)
                                  {
                                    grad_jac_inv[f][d][e]
                                      = inv_jac[f][0] * tmp[0];
                                    for(unsigned int g = 1; g < dim; ++g)
                                      grad_jac_inv[f][d][e]
                                        += inv_jac[f][g] * tmp[g];
                                  }
                              }

                          // the diagonal part of Jacobian gradient comes first
                          Tensor<1,
                                 dim*(dim + 1) / 2,
                                 Tensor<1, dim, VectorizedArray<Number>>>
                            final_grad;
                          for(unsigned int d = 0; d < dim; ++d)
                            for(unsigned int e = 0; e < dim; ++e)
                              final_grad[d][e] = grad_jac_inv[d][d][e];

                          // then the upper-diagonal part
                          for(unsigned int d = 0, count = 0; d < dim; ++d)
                            for(unsigned int e = d + 1; e < dim; ++e, ++count)
                              for(unsigned int f = 0; f < dim; ++f)
                                final_grad[dim + count][f]
                                  = grad_jac_inv[d][e][f];
                          data.first[my_q].jacobian_gradients[0].push_back(
                            final_grad);
                        }
                    }
                }

              if(update_flags & update_quadrature_points)
                {
                  // eventually we turn to the quadrature points that we can
                  // compress in case we have Cartesian cells. we also need to
                  // reorder them into arrays of vectorized data types.  first
                  // go through the cells and find out how much memory we need
                  // to allocate for the quadrature points. We store
                  // n_q_points for all cells but Cartesian cells. For
                  // Cartesian cells, only need to store the values on a
                  // diagonal through the cell (n_q_points_1d). This will give
                  // (with some little indexing) the location of all
                  // quadrature points.
                  const unsigned int old_size
                    = data.first[my_q].quadrature_points.size();
                  data.first[my_q].quadrature_point_offsets.push_back(old_size);

                  if(mapping_info.get_cell_type(cell) == cartesian)
                    {
                      for(unsigned int q = 0; q < n_q_points_1d[my_q][fe_index];
                          ++q)
                        {
                          Point<dim, VectorizedArray<Number>> quad_point;
                          for(unsigned int d = 0; d < dim; ++d)
                            quad_point[d]
                              = cell_data.quadrature_points
                                  [q * step_size_cartesian[my_q][fe_index]][d];
                          data.first[my_q].quadrature_points.push_back(
                            quad_point);
                        }
                    }
                  else
                    {
                      for(unsigned int q = 0; q < n_q_points; ++q)
                        data.first[my_q].quadrature_points.push_back(
                          cell_data.quadrature_points[q]);
                    }
                }
            } // end for ( cell < end_cells )
      }

      template <typename CONTAINER>
      void
      merge_compressed_data(const CONTAINER&           source,
                            CONTAINER&                 destination,
                            std::vector<unsigned int>& indices)
      {
        indices.resize(source.size());
        typename CONTAINER::iterator lookup = destination.begin();
        for(typename CONTAINER::const_iterator it = source.begin();
            it != source.end();
            ++it)
          {
            typename CONTAINER::value_type entry = *it;
            entry.second                         = destination.size();
            lookup = destination.insert(lookup, entry);
            AssertIndexRange(it->second, indices.size());
            indices[it->second] = lookup->second;
            // best guess for insert position of next item
            ++lookup;
          }
      }

      template <int structdim, int dim, typename Number>
      void
      copy_data(const unsigned int                          first_cell,
                const std::array<std::size_t, 2>&           data_shift,
                const std::vector<unsigned int>&            indices_compressed,
                const std::vector<GeometryType>&            cell_type,
                MappingInfoStorage<structdim, dim, Number>& data_cells_local,
                MappingInfoStorage<structdim, dim, Number>& data_cells)
      {
        // Copy the index offsets and shift by the appropriate value
        for(unsigned int lcell = 0;
            lcell < data_cells_local.data_index_offsets.size();
            ++lcell)
          {
            const unsigned int cell = lcell + first_cell;
            data_cells.data_index_offsets[cell]
              = cell_type[cell] <= static_cast<unsigned int>(affine) ?
                  indices_compressed[data_cells_local
                                       .data_index_offsets[lcell]] :
                  data_cells_local.data_index_offsets[lcell] + data_shift[0];
            if(data_cells_local.quadrature_point_offsets.size() > lcell)
              data_cells.quadrature_point_offsets[cell]
                = data_cells_local.quadrature_point_offsets[lcell]
                  + data_shift[1];
          }

        // Copy quadrature points
        if(data_cells.quadrature_point_offsets.empty() == false)
          {
            Point<dim, VectorizedArray<Number>>* out_point
              = &data_cells.quadrature_points[data_shift[1]];
            for(const Point<dim, VectorizedArray<Number>>* point
                = data_cells_local.quadrature_points.begin();
                point != data_cells_local.quadrature_points.end();
                ++point, ++out_point)
              *out_point = *point;
            data_cells_local.quadrature_points.clear();
          }

        // If we have collected Jacobian data, copy Jacobians, JxW values,
        // Jacobian gradients
        if(data_cells_local.JxW_values.empty())
          return;

        std::copy(data_cells_local.JxW_values.begin(),
                  data_cells_local.JxW_values.end(),
                  data_cells.JxW_values.begin() + data_shift[0]);
        data_cells_local.JxW_values.clear();
        std::copy(data_cells_local.normal_vectors.begin(),
                  data_cells_local.normal_vectors.end(),
                  data_cells.normal_vectors.begin() + data_shift[0]);
        data_cells_local.normal_vectors.clear();
        for(unsigned int i = 0; i < 2; ++i)
          {
            std::copy(data_cells_local.jacobians[i].begin(),
                      data_cells_local.jacobians[i].end(),
                      data_cells.jacobians[i].begin() + data_shift[0]);
            data_cells_local.jacobians[i].clear();
            std::copy(data_cells_local.jacobian_gradients[i].begin(),
                      data_cells_local.jacobian_gradients[i].end(),
                      data_cells.jacobian_gradients[i].begin() + data_shift[0]);
            data_cells_local.jacobian_gradients[i].clear();
            std::copy(data_cells_local.normals_times_jacobians[i].begin(),
                      data_cells_local.normals_times_jacobians[i].end(),
                      data_cells.normals_times_jacobians[i].begin()
                        + data_shift[0]);
            data_cells_local.normals_times_jacobians[i].clear();
          }
      }

    } // end of anonymous namespace

    template <int dim, typename Number>
    void
    MappingInfo<dim, Number>::initialize_cells(
      const dealii::Triangulation<dim>&                         tria,
      const std::vector<std::pair<unsigned int, unsigned int>>& cells,
      const std::vector<unsigned int>&                          active_fe_index,
      const Mapping<dim>&                                       mapping,
      const std::vector<dealii::hp::QCollection<1>>&            quad,
      const UpdateFlags update_flags_input)
    {
      const unsigned int n_quads = quad.size();
      const unsigned int n_cells = cells.size();
      const unsigned int vectorization_width
        = VectorizedArray<Number>::n_array_elements;
      Assert(n_cells % vectorization_width == 0, ExcInternalError());
      const unsigned int n_macro_cells = n_cells / vectorization_width;
      cell_data.resize(n_quads);
      cell_type.resize(n_macro_cells);

      // dummy FE that is used to set up an FEValues object. Do not need the
      // actual finite element because we will only evaluate quantities for
      // the mapping that are independent of the FE
      UpdateFlags update_flags = compute_update_flags(update_flags_input, quad);

      for(unsigned int my_q = 0; my_q < n_quads; ++my_q)
        {
          const unsigned int n_hp_quads = quad[my_q].size();
          AssertIndexRange(0, n_hp_quads);
          cell_data[my_q].descriptor.resize(n_hp_quads);
          for(unsigned int q = 0; q < n_hp_quads; ++q)
            cell_data[my_q].descriptor[q].initialize(quad[my_q][q],
                                                     update_default);
        }

      if(n_macro_cells == 0)
        return;

      // Create as many chunks of cells as we have threads and spawn the work
      unsigned int work_per_chunk
        = std::max(8U,
                   (n_macro_cells + MultithreadInfo::n_threads() - 1)
                     / MultithreadInfo::n_threads());

      std::vector<std::pair<std::vector<MappingInfoStorage<dim, dim, Number>>,
                            CompressedCellData<dim, Number>>>
        data_cells_local;
      // Reserve enough space to avoid re-allocation (which would break the
      // references to the data fields passed to the tasks!)
      data_cells_local.reserve(MultithreadInfo::n_threads());

      {
        Threads::TaskGroup<>                  tasks;
        std::pair<unsigned int, unsigned int> cell_range(0U, work_per_chunk);
        while(cell_range.first < n_macro_cells)
          {
            data_cells_local.push_back(std::make_pair(
              std::vector<MappingInfoStorage<dim, dim, Number>>(n_quads),
              CompressedCellData<dim, Number>(get_jacobian_size(tria))));
            tasks += Threads::new_task(&initialize_cell_range<dim, Number>,
                                       cell_range,
                                       tria,
                                       cells,
                                       active_fe_index,
                                       mapping,
                                       quad,
                                       update_flags,
                                       *this,
                                       data_cells_local.back());
            cell_range.first = cell_range.second;
            cell_range.second += work_per_chunk;
          }
        tasks.join_all();
      }

      // Fill in each thread's constant Jacobians into the data of the zeroth
      // chunk in serial
      std::vector<std::vector<unsigned int>> indices_compressed(
        data_cells_local.size());
      for(unsigned int i = 0; i < data_cells_local.size(); ++i)
        merge_compressed_data(data_cells_local[i].second.data,
                              data_cells_local[0].second.data,
                              indices_compressed[i]);

      // Collect all data in the final data fields.
      // First allocate the memory
      const unsigned int n_constant_jacobians
        = data_cells_local[0].second.data.size();
      for(unsigned int my_q = 0; my_q < cell_data.size(); ++my_q)
        {
          cell_data[my_q].data_index_offsets.resize(cell_type.size());
          std::vector<std::array<std::size_t, 2>> shift(
            data_cells_local.size());
          shift[0][0] = n_constant_jacobians;
          shift[0][1] = 0;
          for(unsigned int i = 1; i < data_cells_local.size(); ++i)
            {
              shift[i][0]
                = shift[i - 1][0]
                  + data_cells_local[i - 1].first[my_q].JxW_values.size();
              shift[i][1] = shift[i - 1][1]
                            + data_cells_local[i - 1]
                                .first[my_q]
                                .quadrature_points.size();
            }
          cell_data[my_q].JxW_values.resize_fast(
            shift.back()[0]
            + data_cells_local.back().first[my_q].JxW_values.size());
          cell_data[my_q].jacobians[0].resize_fast(
            cell_data[my_q].JxW_values.size());
          if(update_flags & update_jacobian_grads)
            cell_data[my_q].jacobian_gradients[0].resize_fast(
              cell_data[my_q].JxW_values.size());
          if(update_flags & update_quadrature_points)
            {
              cell_data[my_q].quadrature_point_offsets.resize(cell_type.size());
              cell_data[my_q].quadrature_points.resize_fast(
                shift.back()[1]
                + data_cells_local.back().first[my_q].quadrature_points.size());
            }

          // Start tasks that copy the local data
          Threads::TaskGroup<> tasks;
          for(unsigned int i = 0; i < data_cells_local.size(); ++i)
            tasks += Threads::new_task(&copy_data<dim, dim, Number>,
                                       work_per_chunk * i,
                                       shift[i],
                                       indices_compressed[i],
                                       cell_type,
                                       data_cells_local[i].first[my_q],
                                       cell_data[my_q]);

          // finally, insert the constant cell data at the beginning (the
          // other tasks can already start copying the non-constant data)
          if(my_q == 0)
            {
              for(auto& it : data_cells_local[0].second.data)
                {
                  Tensor<2, dim, VectorizedArray<Number>> jac;
                  for(unsigned int d = 0; d < dim; ++d)
                    for(unsigned int e = 0; e < dim; ++e)
                      for(unsigned int v = 0;
                          v < VectorizedArray<Number>::n_array_elements;
                          ++v)
                        jac[d][e][v] = it.first[d][e][v];
                  AssertIndexRange(it.second, n_constant_jacobians);
                  const std::size_t index           = it.second;
                  cell_data[my_q].JxW_values[index] = determinant(jac);
                  // invert and transpose jac
                  cell_data[my_q].jacobians[0][index] = transpose(invert(jac));
                  // second derivative of transformation is zero on affine cells
                }
            }
          else
            {
              for(unsigned int i = 0; i < n_constant_jacobians; ++i)
                {
                  cell_data[my_q].JxW_values[i] = cell_data[0].JxW_values[i];
                  cell_data[my_q].jacobians[0][i]
                    = cell_data[0].jacobians[0][i];
                }
            }

          // ... wait for the parallel work to finish
          tasks.join_all();
        }
    }

    /* ------------------------- initialization of faces ------------------- */

    // Anonymous namespace with implementation of extraction of values on cell
    // range
    namespace
    {
      template <int dim, typename Number>
      struct CompressedFaceData
      {
        // Constructor. As a scaling factor for the FPArrayComparator, we
        // select the inverse of the Jacobian (not the Jacobian as in the
        // CompressedCellData) and add another factor of 512 to account for
        // some roundoff effects.
        CompressedFaceData(const Number jacobian_size)
          : data(FPArrayComparator<Number>(512. / jacobian_size)),
            jacobian_size(jacobian_size)
        {}

        // Store the Jacobians on both sides of a face (2*(dim*dim) entries),
        // the normal vector (dim entries), and the Jacobian determinant on
        // the face (first tensor) for each entry in the vectorized array
        // (inner tensor). We cannot choose a VectorizedArray type directly
        // because std::map does not provide the necessary alignment upon
        // memory allocation.
        std::map<
          Tensor<1,
                 2 * dim * dim + dim + 1,
                 Tensor<1, VectorizedArray<Number>::n_array_elements, Number>>,
          unsigned int,
          FPArrayComparator<Number>>
          data;

        // Store the scaling factor
        const Number jacobian_size;
      };

      // We always put the derivative normal to the face in the last slot for
      // simpler unit cell gradient computations. This function reorders the
      // indices of the Jacobian appropriately.
      template <int dim>
      unsigned int
      reorder_face_derivative_indices(const unsigned int face_no,
                                      const unsigned int index)
      {
        Assert(index < dim, ExcInternalError());
        if(dim == 3)
          {
            unsigned int table[3][3] = {{1, 2, 0}, {2, 0, 1}, {0, 1, 2}};
            return table[face_no / 2][index];
          }
        else if(dim == 2)
          {
            unsigned int table[2][2] = {{1, 0}, {0, 1}};
            return table[face_no / 2][index];
          }
        else if(dim == 1)
          return 0;
        else
          Assert(
            false,
            ExcNotImplemented("Not possible in dim=" + std::to_string(dim)));

        return numbers::invalid_unsigned_int;
      }

      template <int dim, typename Number>
      void
      initialize_face_range(
        const std::pair<unsigned int, unsigned int>               face_range,
        const dealii::Triangulation<dim>&                         tria,
        const std::vector<std::pair<unsigned int, unsigned int>>& cells,
        const std::vector<
          FaceToCellTopology<VectorizedArray<Number>::n_array_elements>>& faces,
        const Mapping<dim>&                         mapping,
        const UpdateFlags                           update_flags_boundary,
        const UpdateFlags                           update_flags_inner,
        MappingInfo<dim, Number>&                   mapping_info,
        std::pair<std::vector<MappingInfoStorage<dim - 1, dim, Number>>,
                  CompressedFaceData<dim, Number>>& data)
      {
        FE_Nothing<dim> dummy_fe;

        std::vector<std::vector<std::shared_ptr<FEFaceValues<dim>>>>
          fe_face_values_container(mapping_info.face_data.size());
        for(unsigned int my_q = 0; my_q < mapping_info.face_data.size(); ++my_q)
          fe_face_values_container[my_q].resize(
            mapping_info.face_data[my_q].descriptor.size());

        std::vector<std::vector<std::shared_ptr<FEFaceValues<dim>>>>
          fe_boundary_face_values_container(mapping_info.face_data.size());
        for(unsigned int my_q = 0; my_q < mapping_info.face_data.size(); ++my_q)
          fe_boundary_face_values_container[my_q].resize(
            mapping_info.face_data[my_q].descriptor.size());

        std::vector<std::vector<std::shared_ptr<FESubfaceValues<dim>>>>
          fe_subface_values_container(mapping_info.face_data.size());
        for(unsigned int my_q = 0; my_q < mapping_info.face_data.size(); ++my_q)
          fe_subface_values_container[my_q].resize(
            mapping_info.face_data[my_q].descriptor.size());

        LocalData<dim, Number> face_data(get_jacobian_size(tria));

        const unsigned int end_face
          = std::min(std::size_t(face_range.second), faces.size());
        for(unsigned int face = face_range.first; face < end_face; ++face)
          for(unsigned int my_q = 0; my_q < mapping_info.face_data.size();
              ++my_q)
            {
              // currently only non-hp case...
              Assert(mapping_info.face_data[my_q].descriptor.size() == 1,
                     ExcNotImplemented());
              const Quadrature<dim - 1>& quadrature
                = mapping_info.face_data[my_q].descriptor[0].quadrature;

              const bool is_boundary_face = faces[face].cells_exterior[0]
                                            == numbers::invalid_unsigned_int;

              if(is_boundary_face
                 && fe_boundary_face_values_container[my_q][0] == nullptr)
                fe_boundary_face_values_container[my_q][0].reset(
                  new FEFaceValues<dim>(
                    mapping, dummy_fe, quadrature, update_flags_boundary));
              else if(fe_face_values_container[my_q][0] == nullptr)
                fe_face_values_container[my_q][0].reset(new FEFaceValues<dim>(
                  mapping, dummy_fe, quadrature, update_flags_inner));

              FEFaceValues<dim>& fe_face_values
                = is_boundary_face ?
                    *fe_boundary_face_values_container[my_q][0] :
                    *fe_face_values_container[my_q][0];
              const unsigned int n_q_points
                = fe_face_values.n_quadrature_points;
              face_data.resize(n_q_points);

              bool normal_is_similar = true;
              bool JxW_is_similar    = true;
              bool cell_is_cartesian = true;
              for(unsigned int v = 0;
                  v < VectorizedArray<Number>::n_array_elements;
                  ++v)
                {
                  Tensor<2, dim> jacobian_0;
                  double         compare_norm_jac = 1.;
                  if(faces[face].cells_interior[v]
                     != numbers::invalid_unsigned_int)
                    {
                      typename dealii::Triangulation<dim>::cell_iterator
                        cell_it(&tria,
                                cells[faces[face].cells_interior[v]].first,
                                cells[faces[face].cells_interior[v]].second);

                      fe_face_values.reinit(cell_it,
                                            faces[face].interior_face_no);

                      for(unsigned int q = 0; q < n_q_points; ++q)
                        {
                          if(std::abs(
                               fe_face_values.JxW(q) * quadrature.weight(0)
                               - fe_face_values.JxW(0) * quadrature.weight(q))
                             > 2048. * std::numeric_limits<double>::epsilon()
                                 * fe_face_values.JxW(0) * quadrature.weight(q))
                            JxW_is_similar = false;
                          face_data.JxW_values[q][v] = fe_face_values.JxW(q);

                          DerivativeForm<1, dim, dim> inv_jac
                            = fe_face_values.jacobian(q).covariant_form();
                          Tensor<1, dim> normal
                            = fe_face_values.normal_vector(q);

                          // Filter out very small values in normal. No need
                          // to re-normalize because these values cannot enter
                          // the norm significantly: Total size is 1 but 1e-13
                          // squared is 1e-26.
                          for(unsigned int d = 0; d < dim; ++d)
                            if(std::abs(normal[d])
                               < 1024. * std::numeric_limits<double>::epsilon())
                              normal[d] = 0.;

                          if(q == 0)
                            {
                              jacobian_0       = inv_jac;
                              compare_norm_jac = jacobian_0.norm();
                            }
                          for(unsigned int d = 0; d < dim; ++d)
                            for(unsigned int e = 0; e < dim; ++e)
                              {
                                if(std::abs(inv_jac[d][e] - jacobian_0[d][e])
                                   > 2048.
                                       * std::numeric_limits<double>::epsilon()
                                       * compare_norm_jac)
                                  JxW_is_similar = false;
                                const unsigned int ee
                                  = reorder_face_derivative_indices<dim>(
                                    faces[face].interior_face_no, e);
                                face_data.general_jac[q][d][e][v]
                                  = inv_jac[d][ee];
                              }

                          for(unsigned int d = 0; d < dim; ++d)
                            {
                              face_data.normal_vectors[q][d][v] = normal[d];
                              if(std::abs(normal[d]
                                          - fe_face_values.normal_vector(0)[d])
                                 > 1024.
                                     * std::numeric_limits<double>::epsilon())
                                normal_is_similar = false;
                            }
                          if(std::abs(
                               std::abs(
                                 normal[faces[face].interior_face_no / 2])
                               - 1.)
                             > 1024. * std::numeric_limits<double>::epsilon())
                            cell_is_cartesian = false;
                          for(unsigned int d = 0; d < dim; ++d)
                            for(unsigned int e = 0; e < dim; ++e)
                              if(e != d
                                 && std::abs(inv_jac[d][e])
                                      > 2048
                                          * std::numeric_limits<
                                              double>::epsilon()
                                          * compare_norm_jac)
                                cell_is_cartesian = false;

                          if(fe_face_values.get_update_flags()
                             & update_quadrature_points)
                            for(unsigned int d = 0; d < dim; ++d)
                              face_data.quadrature_points[q][d][v]
                                = fe_face_values.quadrature_point(q)[d];
                        }
                    }
                  // Fill up with data of the zeroth component to avoid
                  // false negatives when checking for similarities
                  else
                    {
                      for(unsigned int q = 0; q < n_q_points; ++q)
                        {
                          face_data.JxW_values[q][v]
                            = face_data.JxW_values[q][0];
                          for(unsigned int d = 0; d < dim; ++d)
                            for(unsigned int e = 0; e < dim; ++e)
                              {
                                face_data.general_jac[q][d][e][v]
                                  = face_data.general_jac[q][d][e][0];
                                face_data.normal_vectors[q][d][v]
                                  = face_data.normal_vectors[q][d][0];
                              }
                        }
                    }
                  if(is_boundary_face == false
                     && faces[face].cells_exterior[v]
                          != numbers::invalid_unsigned_int)
                    {
                      typename dealii::Triangulation<dim>::cell_iterator
                        cell_it(&tria,
                                cells[faces[face].cells_exterior[v]].first,
                                cells[faces[face].cells_exterior[v]].second);

                      const FEValuesBase<dim>* actual_fe_face_values = nullptr;
                      if(faces[face].subface_index
                         >= GeometryInfo<dim>::max_children_per_cell)
                        {
                          fe_face_values.reinit(cell_it,
                                                faces[face].exterior_face_no);
                          actual_fe_face_values = &fe_face_values;
                        }
                      else
                        {
                          if(fe_subface_values_container[my_q][0] == nullptr)
                            fe_subface_values_container[my_q][0].reset(
                              new FESubfaceValues<dim>(mapping,
                                                       dummy_fe,
                                                       quadrature,
                                                       update_flags_inner));
                          fe_subface_values_container[my_q][0]->reinit(
                            cell_it,
                            faces[face].exterior_face_no,
                            faces[face].subface_index);
                          actual_fe_face_values
                            = fe_subface_values_container[my_q][0].get();
                        }
                      for(unsigned int q = 0; q < n_q_points; ++q)
                        {
                          DerivativeForm<1, dim, dim> inv_jac
                            = actual_fe_face_values->jacobian(q)
                                .covariant_form();
                          if(q == 0)
                            {
                              jacobian_0       = inv_jac;
                              compare_norm_jac = jacobian_0.norm();
                            }
                          for(unsigned int d = 0; d < dim; ++d)
                            for(unsigned int e = 0; e < dim; ++e)
                              {
                                if(std::abs(inv_jac[d][e] - jacobian_0[d][e])
                                   > 2048.
                                       * std::numeric_limits<double>::epsilon()
                                       * compare_norm_jac)
                                  JxW_is_similar = false;
                                const unsigned int ee
                                  = reorder_face_derivative_indices<dim>(
                                    faces[face].exterior_face_no, e);
                                face_data.general_jac[n_q_points + q][d][e][v]
                                  = inv_jac[d][ee];
                              }
                          for(unsigned int d = 0; d < dim; ++d)
                            for(unsigned int e = 0; e < dim; ++e)
                              if(e != d
                                 && std::abs(inv_jac[d][e])
                                      > 2048
                                          * std::numeric_limits<
                                              double>::epsilon()
                                          * compare_norm_jac)
                                cell_is_cartesian = false;
                        }
                    }
                  // Fill up with 'known' values
                  else if(is_boundary_face == false)
                    {
                      Assert(faces[face].cells_exterior[0]
                               != numbers::invalid_unsigned_int,
                             ExcInternalError());
                      for(unsigned int q = 0; q < n_q_points; ++q)
                        for(unsigned int d = 0; d < dim; ++d)
                          for(unsigned int e = 0; e < dim; ++e)
                            face_data.general_jac[n_q_points + q][d][e][v]
                              = face_data.general_jac[n_q_points + q][d][e][0];
                    }
                  // If boundary face, simply set the data to zero (will not
                  // be used). Note that faces over periodic boundary
                  // conditions will be treated as interior ones in this setup.
                  else
                    for(unsigned int q = 0; q < n_q_points; ++q)
                      for(unsigned int d = 0; d < dim; ++d)
                        for(unsigned int e = 0; e < dim; ++e)
                          face_data.general_jac[n_q_points + q][d][e][v] = 0.;
                }

              // check if face is affine or at least if it is flat
              // (i.e., all normal vectors are the same)
              if(my_q == 0)
                {
                  GeometryType face_type = affine;
                  if(JxW_is_similar == false)
                    face_type = flat_faces;
                  if(normal_is_similar == false)
                    face_type = general;
                  if(face_type == affine && cell_is_cartesian)
                    face_type = cartesian;
                  mapping_info.face_type[face] = face_type;
                }

              // Fill in quadrature points
              if(fe_face_values.get_update_flags() & update_quadrature_points)
                {
                  data.first[my_q].quadrature_point_offsets.push_back(
                    data.first[my_q].quadrature_points.size());
                  if(fe_face_values.get_update_flags()
                     & update_quadrature_points)
                    for(unsigned int q = 0; q < n_q_points; ++q)
                      data.first[my_q].quadrature_points.push_back(
                        face_data.quadrature_points[q]);
                }

              typedef Tensor<1,
                             VectorizedArray<Number>::n_array_elements,
                             Number>
                           VEC_ARRAY;
              unsigned int insert_position = data.first[my_q].JxW_values.size();

              // Fill in JxW values, apply compression
              if(mapping_info.face_type[face] <= affine)
                {
                  if(my_q == 0)
                    {
                      // find out if we already had the same JxW values before
                      std::pair<Tensor<1, 2 * dim * dim + dim + 1, VEC_ARRAY>,
                                unsigned int>
                        new_entry;
                      new_entry.second = data.second.data.size();
                      for(unsigned int v = 0;
                          v < VectorizedArray<Number>::n_array_elements;
                          ++v)
                        new_entry.first[2 * dim * dim + dim][v]
                          = face_data.JxW_values[0][v] / quadrature.weight(0)
                            / Utilities::fixed_power<dim>(face_data.jac_size);

                      new_entry.second = data.second.data.size();
                      for(unsigned int d = 0; d < dim; ++d)
                        for(unsigned int e = 0; e < dim; ++e)
                          for(unsigned int v = 0;
                              v < VectorizedArray<Number>::n_array_elements;
                              ++v)
                            new_entry.first[d * dim + e][v]
                              = face_data.general_jac[0][d][e][v];
                      if(is_boundary_face == false)
                        for(unsigned int d = 0; d < dim; ++d)
                          for(unsigned int e = 0; e < dim; ++e)
                            for(unsigned int v = 0;
                                v < VectorizedArray<Number>::n_array_elements;
                                ++v)
                              new_entry.first[dim * dim + d * dim + e][v]
                                = face_data.general_jac[n_q_points][d][e][v];
                      // we need to add the normal vector here because we
                      // store both the inverse jacobian and the normal vector
                      // times the jacobian; of course, there will be
                      // different values in their product for normal vectors
                      // oriented in different ways (the memory saving is
                      // still significant); we need to divide by the jacobian
                      // size to get the right scaling
                      for(unsigned int d = 0; d < dim; ++d)
                        for(unsigned int v = 0;
                            v < VectorizedArray<Number>::n_array_elements;
                            ++v)
                          new_entry.first[2 * dim * dim + d][v]
                            = face_data.normal_vectors[0][d][v]
                              / face_data.jac_size;

                      insert_position
                        = data.second.data.insert(new_entry).first->second;
                    }
                  else
                    insert_position = data.first[0].data_index_offsets[face];
                }
              data.first[my_q].data_index_offsets.push_back(insert_position);
              if(mapping_info.face_type[face] > affine)
                {
                  for(unsigned int q = 0; q < n_q_points; ++q)
                    data.first[my_q].JxW_values.push_back(
                      face_data.JxW_values[q]);
                  for(unsigned int q = 0; q < n_q_points; ++q)
                    data.first[my_q].normal_vectors.push_back(
                      face_data.normal_vectors[q]);
                  for(unsigned int q = 0; q < n_q_points; ++q)
                    data.first[my_q].jacobians[0].push_back(
                      face_data.general_jac[q]);
                  if(is_boundary_face == false)
                    for(unsigned int q = 0; q < n_q_points; ++q)
                      data.first[my_q].jacobians[1].push_back(
                        face_data.general_jac[n_q_points + q]);
                }
            }
      }

      template <int dim, typename Number>
      void
      compute_normal_times_jacobian(
        const unsigned int               first_face,
        const unsigned int               last_face,
        const std::vector<GeometryType>& face_type,
        const std::vector<
          FaceToCellTopology<VectorizedArray<Number>::n_array_elements>>& faces,
        MappingInfoStorage<dim - 1, dim, Number>& data_faces)
      {
        for(unsigned int face = first_face; face < last_face; ++face)
          {
            const bool is_boundary_face
              = faces[face].cells_exterior[0] == numbers::invalid_unsigned_int;
            const unsigned int n_q_points_work
              = face_type[face] > affine ? data_faces.descriptor[0].n_q_points :
                                           1;
            const unsigned int offset = data_faces.data_index_offsets[face];

            for(unsigned int q = 0; q < n_q_points_work; ++q)
              {
                data_faces.normals_times_jacobians[0][offset + q]
                  = data_faces.normal_vectors[offset + q]
                    * data_faces.jacobians[0][offset + q];
                if(is_boundary_face == false)
                  data_faces.normals_times_jacobians[1][offset + q]
                    = data_faces.normal_vectors[offset + q]
                      * data_faces.jacobians[1][offset + q];
              }
          }
      }

    } // end of anonymous namespace

    template <int dim, typename Number>
    void
    MappingInfo<dim, Number>::initialize_faces(
      const dealii::Triangulation<dim>&                         tria,
      const std::vector<std::pair<unsigned int, unsigned int>>& cells,
      const std::vector<
        FaceToCellTopology<VectorizedArray<Number>::n_array_elements>>& faces,
      const Mapping<dim>&                                               mapping,
      const std::vector<dealii::hp::QCollection<1>>&                    quad,
      const UpdateFlags update_flags_boundary_faces,
      const UpdateFlags update_flags_inner_faces)
    {
      face_type.resize(faces.size(), general);
      face_data.resize(quad.size());

      // We currently always set the same flags on both inner and boundary
      // faces to the same value. At some point, we might want to separate the
      // two.
      UpdateFlags update_flags_compute_boundary
        = ((update_flags_inner_faces | update_flags_boundary_faces)
               & update_quadrature_points ?
             update_quadrature_points :
             update_default)
          | update_normal_vectors | update_JxW_values | update_jacobians;
      UpdateFlags update_flags_compute_inner
        = ((update_flags_inner_faces | update_flags_boundary_faces)
               & update_quadrature_points ?
             update_quadrature_points :
             update_default)
          | update_normal_vectors | update_JxW_values | update_jacobians;
      UpdateFlags update_flags_common
        = update_flags_inner_faces | update_flags_boundary_faces;

      for(unsigned int my_q = 0; my_q < quad.size(); ++my_q)
        {
          face_data[my_q].descriptor.resize(quad[my_q].size());
          for(unsigned int hpq = 0; hpq < quad[my_q].size(); ++hpq)
            face_data[my_q].descriptor[hpq].initialize(
              quad[my_q][hpq], update_flags_compute_inner);
        }

      if(faces.size() == 0)
        return;

      // Create as many chunks of cells as we have threads and spawn the work
      unsigned int work_per_chunk
        = std::max(std::size_t(8),
                   (faces.size() + MultithreadInfo::n_threads() - 1)
                     / MultithreadInfo::n_threads());

      std::vector<
        std::pair<std::vector<MappingInfoStorage<dim - 1, dim, Number>>,
                  CompressedFaceData<dim, Number>>>
        data_faces_local;
      // Reserve enough space to avoid re-allocation (which would destroy the
      // references passed to the tasks!)
      data_faces_local.reserve(MultithreadInfo::n_threads());

      {
        Threads::TaskGroup<>                  tasks;
        std::pair<unsigned int, unsigned int> face_range(0U, work_per_chunk);
        while(face_range.first < faces.size())
          {
            data_faces_local.push_back(std::make_pair(
              std::vector<MappingInfoStorage<dim - 1, dim, Number>>(
                quad.size()),
              CompressedFaceData<dim, Number>(get_jacobian_size(tria))));
            tasks += Threads::new_task(&initialize_face_range<dim, Number>,
                                       face_range,
                                       tria,
                                       cells,
                                       faces,
                                       mapping,
                                       update_flags_compute_boundary,
                                       update_flags_compute_inner,
                                       *this,
                                       data_faces_local.back());
            face_range.first = face_range.second;
            face_range.second += work_per_chunk;
          }

        tasks.join_all();
      }

      // Fill in each thread's constant data (normals, JxW, normals x
      // Jacobian) into the data of the zeroth chunk in serial
      std::vector<std::vector<unsigned int>> indices_compressed(
        data_faces_local.size());
      for(unsigned int i = 0; i < data_faces_local.size(); ++i)
        merge_compressed_data(data_faces_local[i].second.data,
                              data_faces_local[0].second.data,
                              indices_compressed[i]);

      // Collect all data in the final data fields.
      // First allocate the memory
      const unsigned int n_constant_data
        = data_faces_local[0].second.data.size();
      for(unsigned int my_q = 0; my_q < face_data.size(); ++my_q)
        {
          face_data[my_q].data_index_offsets.resize(face_type.size());
          std::vector<std::array<std::size_t, 2>> shift(
            data_faces_local.size());
          shift[0][0] = n_constant_data;
          shift[0][1] = 0;
          for(unsigned int i = 1; i < data_faces_local.size(); ++i)
            {
              shift[i][0]
                = shift[i - 1][0]
                  + data_faces_local[i - 1].first[my_q].JxW_values.size();
              shift[i][1] = shift[i - 1][1]
                            + data_faces_local[i - 1]
                                .first[my_q]
                                .quadrature_points.size();
            }
          face_data[my_q].JxW_values.resize_fast(
            shift.back()[0]
            + data_faces_local.back().first[my_q].JxW_values.size());
          face_data[my_q].normal_vectors.resize_fast(
            face_data[my_q].JxW_values.size());
          face_data[my_q].jacobians[0].resize_fast(
            face_data[my_q].JxW_values.size());
          face_data[my_q].jacobians[1].resize_fast(
            face_data[my_q].JxW_values.size());
          if(update_flags_common & update_jacobian_grads)
            {
              face_data[my_q].jacobian_gradients[0].resize_fast(
                face_data[my_q].JxW_values.size());
              face_data[my_q].jacobian_gradients[1].resize_fast(
                face_data[my_q].JxW_values.size());
            }
          face_data[my_q].normals_times_jacobians[0].resize_fast(
            face_data[my_q].JxW_values.size());
          face_data[my_q].normals_times_jacobians[1].resize_fast(
            face_data[my_q].JxW_values.size());
          if(update_flags_common & update_quadrature_points)
            {
              face_data[my_q].quadrature_point_offsets.resize(face_type.size());
              face_data[my_q].quadrature_points.resize_fast(
                shift.back()[1]
                + data_faces_local.back().first[my_q].quadrature_points.size());
            }

          // start the tasks to gather the data in parallel
          Threads::TaskGroup<> tasks;
          for(unsigned int i = 0; i < data_faces_local.size(); ++i)
            tasks += Threads::new_task(&copy_data<dim - 1, dim, Number>,
                                       work_per_chunk * i,
                                       shift[i],
                                       indices_compressed[i],
                                       face_type,
                                       data_faces_local[i].first[my_q],
                                       face_data[my_q]);

          // fill the constant data fields (in parallel to the loop above)
          if(my_q == 0)
            {
              const Number jac_size = get_jacobian_size(tria);
              for(auto& it : data_faces_local[0].second.data)
                {
                  // JxW values; invert previously applied scaling
                  for(unsigned int v = 0;
                      v < VectorizedArray<Number>::n_array_elements;
                      ++v)
                    face_data[my_q].JxW_values[it.second][v]
                      = it.first[2 * dim * dim + dim][v]
                        * Utilities::fixed_power<dim>(jac_size);

                  // inverse Jacobians
                  for(unsigned int i = 0; i < 2; ++i)
                    for(unsigned int d = 0; d < dim; ++d)
                      for(unsigned int e = 0; e < dim; ++e)
                        for(unsigned int v = 0;
                            v < VectorizedArray<Number>::n_array_elements;
                            ++v)
                          face_data[my_q].jacobians[i][it.second][d][e][v]
                            = it.first[i * dim * dim + d * dim + e][v];

                  // normal vectors; invert previously applied scaling
                  for(unsigned int d = 0; d < dim; ++d)
                    for(unsigned int v = 0;
                        v < VectorizedArray<Number>::n_array_elements;
                        ++v)
                      face_data[my_q].normal_vectors[it.second][d][v]
                        = it.first[2 * dim * dim + d][v] * jac_size;
                }
            }
          else
            {
              for(unsigned int i = 0; i < n_constant_data; ++i)
                {
                  face_data[my_q].JxW_values[i] = face_data[0].JxW_values[i];
                  face_data[my_q].normal_vectors[i]
                    = face_data[0].normal_vectors[i];
                  for(unsigned k = 0; k < 2; ++k)
                    face_data[my_q].jacobians[k][i]
                      = face_data[0].jacobians[k][i];
                }
            }

          // ... wait for the parallel work to finish
          tasks.join_all();

          // finally compute the normal times the jacobian
          for(unsigned int i = 0; i < data_faces_local.size(); ++i)
            tasks += Threads::new_task(
              &compute_normal_times_jacobian<dim, Number>,
              work_per_chunk * i,
              std::min(work_per_chunk * (i + 1), (unsigned int) faces.size()),
              face_type,
              faces,
              face_data[my_q]);
          tasks.join_all();
        }
    }

    template <int dim, typename Number>
    void
    MappingInfo<dim, Number>::initialize_faces_by_cells(
      const dealii::Triangulation<dim>&                         tria,
      const std::vector<std::pair<unsigned int, unsigned int>>& cells,
      const Mapping<dim>&                                       mapping,
      const std::vector<dealii::hp::QCollection<1>>&            quad,
      const UpdateFlags update_flags_faces_by_cells)
    {
      if(update_flags_faces_by_cells == update_default)
        return;

      face_data_by_cells.resize(quad.size());
      const unsigned int n_quads = quad.size();
      const unsigned int vectorization_width
        = VectorizedArray<Number>::n_array_elements;
      UpdateFlags update_flags
        = (update_flags_faces_by_cells & update_quadrature_points ?
             update_quadrature_points :
             update_default)
          | update_normal_vectors | update_JxW_values | update_jacobians;

      for(unsigned int my_q = 0; my_q < n_quads; ++my_q)
        {
          const unsigned int n_hp_quads = quad[my_q].size();
          AssertIndexRange(0, n_hp_quads);
          face_data_by_cells[my_q].descriptor.resize(n_hp_quads);
          for(unsigned int q = 0; q < n_hp_quads; ++q)
            face_data_by_cells[my_q].descriptor[q].initialize(quad[my_q][q],
                                                              update_default);

          // since we already know the cell type, we can pre-allocate the right
          // amount of data straight away and we just need to do some basic counting
          AssertDimension(cell_type.size(), cells.size() / vectorization_width);
          face_data_by_cells[my_q].data_index_offsets.resize(
            cell_type.size() * GeometryInfo<dim>::faces_per_cell);
          if(update_flags & update_quadrature_points)
            face_data_by_cells[my_q].quadrature_point_offsets.resize(
              cell_type.size() * GeometryInfo<dim>::faces_per_cell);
          std::size_t storage_length = 0;
          for(unsigned int i = 0; i < cell_type.size(); ++i)
            for(unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
                ++face)
              {
                if(cell_type[i] <= affine)
                  {
                    face_data_by_cells[my_q]
                      .data_index_offsets[i * GeometryInfo<dim>::faces_per_cell
                                          + face]
                      = storage_length;
                    ++storage_length;
                  }
                else
                  {
                    face_data_by_cells[my_q]
                      .data_index_offsets[i * GeometryInfo<dim>::faces_per_cell
                                          + face]
                      = storage_length;
                    storage_length
                      += face_data_by_cells[my_q].descriptor[0].n_q_points;
                  }
                if(update_flags & update_quadrature_points)
                  face_data_by_cells[my_q].quadrature_point_offsets
                    [i * GeometryInfo<dim>::faces_per_cell + face]
                    = (i * GeometryInfo<dim>::faces_per_cell + face)
                      * face_data_by_cells[my_q].descriptor[0].n_q_points;
              }
          face_data_by_cells[my_q].JxW_values.resize_fast(
            storage_length * GeometryInfo<dim>::faces_per_cell);
          face_data_by_cells[my_q].jacobians[0].resize_fast(
            storage_length * GeometryInfo<dim>::faces_per_cell);
          if(update_flags & update_normal_vectors)
            face_data_by_cells[my_q].normal_vectors.resize_fast(
              storage_length * GeometryInfo<dim>::faces_per_cell);
          if(update_flags & update_normal_vectors
             && update_flags & update_jacobians)
            face_data_by_cells[my_q].normals_times_jacobians[0].resize_fast(
              storage_length * GeometryInfo<dim>::faces_per_cell);
          if(update_flags & update_jacobian_grads)
            face_data_by_cells[my_q].jacobian_gradients[0].resize_fast(
              storage_length * GeometryInfo<dim>::faces_per_cell);

          if(update_flags & update_quadrature_points)
            face_data_by_cells[my_q].quadrature_points.resize_fast(
              cell_type.size() * GeometryInfo<dim>::faces_per_cell
              * face_data_by_cells[my_q].descriptor[0].n_q_points);
        }

      FE_Nothing<dim> dummy_fe;
      // currently no hp-indices implemented
      const unsigned int fe_index = 0;
      std::vector<std::vector<std::shared_ptr<dealii::FEFaceValues<dim>>>>
        fe_face_values(face_data_by_cells.size());
      for(unsigned int i = 0; i < fe_face_values.size(); ++i)
        fe_face_values[i].resize(face_data_by_cells[i].descriptor.size());
      for(unsigned int cell = 0; cell < cell_type.size(); ++cell)
        for(unsigned int my_q = 0; my_q < face_data_by_cells.size(); ++my_q)
          for(unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
              ++face)
            {
              if(fe_face_values[my_q][fe_index].get() == nullptr)
                fe_face_values[my_q][fe_index].reset(
                  new dealii::FEFaceValues<dim>(
                    mapping,
                    dummy_fe,
                    face_data_by_cells[my_q].descriptor[fe_index].quadrature,
                    update_flags));
              dealii::FEFaceValues<dim>& fe_val
                = *fe_face_values[my_q][fe_index];
              const unsigned int offset
                = face_data_by_cells[my_q]
                    .data_index_offsets[cell * GeometryInfo<dim>::faces_per_cell
                                        + face];

              for(unsigned int v = 0; v < vectorization_width; ++v)
                {
                  typename dealii::Triangulation<dim>::cell_iterator cell_it(
                    &tria,
                    cells[cell * vectorization_width + v].first,
                    cells[cell * vectorization_width + v].second);
                  fe_val.reinit(cell_it, face);

                  // copy data for affine data type
                  if(cell_type[cell] <= affine)
                    {
                      if(update_flags & update_JxW_values)
                        face_data_by_cells[my_q].JxW_values[offset][v]
                          = fe_val.JxW(0)
                            / face_data_by_cells[my_q]
                                .descriptor[fe_index]
                                .quadrature.weight(0);
                      if(update_flags & update_jacobians)
                        {
                          DerivativeForm<1, dim, dim> inv_jac
                            = fe_val.jacobian(0).covariant_form();
                          for(unsigned int d = 0; d < dim; ++d)
                            for(unsigned int e = 0; e < dim; ++e)
                              {
                                const unsigned int ee
                                  = reorder_face_derivative_indices<dim>(face,
                                                                         e);
                                face_data_by_cells[my_q]
                                  .jacobians[0][offset][d][e][v]
                                  = inv_jac[d][ee];
                              }
                        }
                      if(update_flags & update_jacobian_grads)
                        {
                          Assert(false, ExcNotImplemented());
                        }
                      if(update_flags & update_normal_vectors)
                        for(unsigned int d = 0; d < dim; ++d)
                          face_data_by_cells[my_q].normal_vectors[offset][d][v]
                            = fe_val.normal_vector(0)[d];
                    }
                  // copy data for general data type
                  else
                    {
                      if(update_flags & update_JxW_values)
                        for(unsigned int q = 0; q < fe_val.n_quadrature_points;
                            ++q)
                          face_data_by_cells[my_q].JxW_values[offset + q][v]
                            = fe_val.JxW(q);
                      if(update_flags & update_jacobians)
                        for(unsigned int q = 0; q < fe_val.n_quadrature_points;
                            ++q)
                          {
                            DerivativeForm<1, dim, dim> inv_jac
                              = fe_val.jacobian(q).covariant_form();
                            for(unsigned int d = 0; d < dim; ++d)
                              for(unsigned int e = 0; e < dim; ++e)
                                {
                                  const unsigned int ee
                                    = reorder_face_derivative_indices<dim>(face,
                                                                           e);
                                  face_data_by_cells[my_q]
                                    .jacobians[0][offset + q][d][e][v]
                                    = inv_jac[d][ee];
                                }
                          }
                      if(update_flags & update_jacobian_grads)
                        {
                          Assert(false, ExcNotImplemented());
                        }
                      if(update_flags & update_normal_vectors)
                        for(unsigned int q = 0; q < fe_val.n_quadrature_points;
                            ++q)
                          for(unsigned int d = 0; d < dim; ++d)
                            face_data_by_cells[my_q]
                              .normal_vectors[offset + q][d][v]
                              = fe_val.normal_vector(q)[d];
                    }
                  if(update_flags & update_quadrature_points)
                    for(unsigned int q = 0; q < fe_val.n_quadrature_points; ++q)
                      for(unsigned int d = 0; d < dim; ++d)
                        face_data_by_cells[my_q].quadrature_points
                          [face_data_by_cells[my_q].quadrature_point_offsets
                             [cell * GeometryInfo<dim>::faces_per_cell + face]
                           + q][d][v]
                          = fe_val.quadrature_point(q)[d];
                }
              if(update_flags & update_normal_vectors
                 && update_flags & update_jacobians)
                for(unsigned int q = 0; q < (cell_type[cell] <= affine ?
                                               1 :
                                               fe_val.n_quadrature_points);
                    ++q)
                  face_data_by_cells[my_q]
                    .normals_times_jacobians[0][offset + q]
                    = face_data_by_cells[my_q].normal_vectors[offset + q]
                      * face_data_by_cells[my_q].jacobians[0][offset + q];
            }
    }

    template <int dim, typename Number>
    std::size_t
    MappingInfo<dim, Number>::memory_consumption() const
    {
      std::size_t memory = MemoryConsumption::memory_consumption(cell_data);
      memory += MemoryConsumption::memory_consumption(face_data);
      memory += cell_type.capacity() * sizeof(GeometryType);
      memory += face_type.capacity() * sizeof(GeometryType);
      memory += sizeof(*this);
      return memory;
    }

    template <int dim, typename Number>
    template <typename StreamType>
    void
    MappingInfo<dim, Number>::print_memory_consumption(
      StreamType&     out,
      const SizeInfo& task_info) const
    {
      out << "    Cell types:                      ";
      task_info.print_memory_statistics(
        out, cell_type.capacity() * sizeof(GeometryType));
      out << "    Face types:                      ";
      task_info.print_memory_statistics(
        out, face_type.capacity() * sizeof(GeometryType));
      for(unsigned int j = 0; j < cell_data.size(); ++j)
        {
          out << "    Data component " << j << std::endl;
          cell_data[j].print_memory_consumption(out, task_info);
          face_data[j].print_memory_consumption(out, task_info);
        }
    }

    /* ------------------------------------------------------------------ */

    template <typename Number>
    FPArrayComparator<Number>::FPArrayComparator(const Number scaling)
      : tolerance(scaling * std::numeric_limits<double>::epsilon() * 1024.)
    {}

    template <typename Number>
    bool
    FPArrayComparator<Number>::operator()(const std::vector<Number>& v1,
                                          const std::vector<Number>& v2) const
    {
      const unsigned int s1 = v1.size(), s2 = v2.size();
      if(s1 < s2)
        return true;
      else if(s1 > s2)
        return false;
      else
        for(unsigned int i = 0; i < s1; ++i)
          if(v1[i] < v2[i] - tolerance)
            return true;
          else if(v1[i] > v2[i] + tolerance)
            return false;
      return false;
    }

    template <typename Number>
    bool
    FPArrayComparator<Number>::operator()(
      const Tensor<1, VectorizedArray<Number>::n_array_elements, Number>& t1,
      const Tensor<1, VectorizedArray<Number>::n_array_elements, Number>& t2)
      const
    {
      for(unsigned int k = 0; k < VectorizedArray<Number>::n_array_elements;
          ++k)
        if(t1[k] < t2[k] - tolerance)
          return true;
        else if(t1[k] > t2[k] + tolerance)
          return false;
      return false;
    }

    template <typename Number>
    template <int dim>
    bool
    FPArrayComparator<Number>::operator()(
      const Tensor<
        1,
        dim,
        Tensor<1, VectorizedArray<Number>::n_array_elements, Number>>& t1,
      const Tensor<
        1,
        dim,
        Tensor<1, VectorizedArray<Number>::n_array_elements, Number>>& t2) const
    {
      for(unsigned int d = 0; d < dim; ++d)
        for(unsigned int k = 0; k < VectorizedArray<Number>::n_array_elements;
            ++k)
          if(t1[d][k] < t2[d][k] - tolerance)
            return true;
          else if(t1[d][k] > t2[d][k] + tolerance)
            return false;
      return false;
    }

    template <typename Number>
    template <int dim>
    bool
    FPArrayComparator<Number>::operator()(
      const Tensor<
        2,
        dim,
        Tensor<1, VectorizedArray<Number>::n_array_elements, Number>>& t1,
      const Tensor<
        2,
        dim,
        Tensor<1, VectorizedArray<Number>::n_array_elements, Number>>& t2) const
    {
      for(unsigned int d = 0; d < dim; ++d)
        for(unsigned int e = 0; e < dim; ++e)
          for(unsigned int k = 0; k < VectorizedArray<Number>::n_array_elements;
              ++k)
            if(t1[d][e][k] < t2[d][e][k] - tolerance)
              return true;
            else if(t1[d][e][k] > t2[d][e][k] + tolerance)
              return false;
      return false;
    }

  } // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif

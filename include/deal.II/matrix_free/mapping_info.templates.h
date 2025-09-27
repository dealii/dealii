// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_matrix_free_mapping_info_templates_h
#define dealii_matrix_free_mapping_info_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/floating_point_comparator.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/matrix_free/evaluation_template_factory.h>
#include <deal.II/matrix_free/fe_evaluation_data.h>
#include <deal.II/matrix_free/mapping_info.h>
#include <deal.II/matrix_free/mapping_info_storage.templates.h>
#include <deal.II/matrix_free/util.h>

#include <limits>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace MatrixFreeFunctions
  {
    template <int dim, typename Number, typename VectorizedArrayType>
    void
    MappingInfo<dim, Number, VectorizedArrayType>::clear()
    {
      cell_data.clear();
      face_data.clear();
      face_data_by_cells.clear();
      cell_type.clear();
      face_type.clear();
      mapping_collection = nullptr;
      mapping            = nullptr;
    }



    template <int dim, typename Number, typename VectorizedArrayType>
    void
    MappingInfo<dim, Number, VectorizedArrayType>::initialize(
      const dealii::Triangulation<dim>                         &tria,
      const std::vector<std::pair<unsigned int, unsigned int>> &cells,
      const FaceInfo<VectorizedArrayType::size()>              &face_info,
      const std::vector<unsigned int>                          &active_fe_index,
      const std::shared_ptr<dealii::hp::MappingCollection<dim>> &mapping,
      const std::vector<dealii::hp::QCollection<dim>>           &quad,
      const UpdateFlags update_flags_cells,
      const UpdateFlags update_flags_boundary_faces,
      const UpdateFlags update_flags_inner_faces,
      const UpdateFlags update_flags_faces_by_cells,
      const bool        piola_transform)
    {
      clear();
      this->mapping_collection = mapping;
      this->mapping            = &mapping->operator[](0);

      cell_data.resize(quad.size());
      face_data.resize(quad.size());
      face_data_by_cells.resize(quad.size());

      const bool is_mixed_mesh = tria.is_mixed_mesh();

      // dummy FE that is used to set up an FEValues object. Do not need the
      // actual finite element because we will only evaluate quantities for
      // the mapping that are independent of the FE
      this->update_flags_cells =
        MappingInfoStorage<dim, dim, VectorizedArrayType>::compute_update_flags(
          update_flags_cells, quad, piola_transform);

      this->update_flags_boundary_faces =
        ((update_flags_inner_faces | update_flags_boundary_faces) &
             update_quadrature_points ?
           update_quadrature_points :
           update_default) |
        ((((update_flags_inner_faces | update_flags_boundary_faces) &
           (update_jacobian_grads | update_hessians)) != 0u ||
          (piola_transform &&
           ((update_flags_inner_faces | update_flags_boundary_faces) &
            (update_gradients | update_contravariant_transformation)) != 0u)) ?
           update_jacobian_grads :
           update_default) |
        update_normal_vectors | update_JxW_values | update_jacobians;
      this->update_flags_inner_faces    = this->update_flags_boundary_faces;
      this->update_flags_faces_by_cells = update_flags_faces_by_cells;

      reference_cell_types.resize(quad.size());

      for (unsigned int my_q = 0; my_q < quad.size(); ++my_q)
        {
          const unsigned int n_hp_quads = quad[my_q].size();
          AssertIndexRange(0, n_hp_quads);

          const unsigned int scale = std::max<unsigned int>(1, dim - 1);

          cell_data[my_q].descriptor.resize(n_hp_quads);
          face_data[my_q].descriptor.resize(n_hp_quads * scale);
          face_data[my_q].q_collection.resize(n_hp_quads);
          face_data_by_cells[my_q].descriptor.resize(n_hp_quads * scale);
          reference_cell_types[my_q].resize(n_hp_quads);

          for (unsigned int hpq = 0; hpq < n_hp_quads; ++hpq)
            {
              bool flag = quad[my_q][hpq].is_tensor_product();

              if (flag)
                for (unsigned int i = 1; i < dim; ++i)
                  flag &= quad[my_q][hpq].get_tensor_basis()[0] ==
                          quad[my_q][hpq].get_tensor_basis()[i];

              if (flag == false)
                {
                  cell_data[my_q].descriptor[hpq].initialize(quad[my_q][hpq]);
                  const auto quad_face =
                    get_unique_face_quadratures(quad[my_q][hpq]);

                  if (quad_face.first.size() > 0) // line, quad
                    {
                      face_data[my_q].descriptor[hpq * scale].initialize(
                        quad_face.first);
                      face_data_by_cells[my_q]
                        .descriptor[hpq * scale]
                        .initialize(quad_face.first);
                    }

                  if (quad_face.second.size() > 0) // triangle
                    {
                      AssertDimension(dim, 3);
                      face_data[my_q]
                        .descriptor[hpq * scale + (is_mixed_mesh ? 1 : 0)]
                        .initialize(quad_face.second);
                      face_data_by_cells[my_q]
                        .descriptor[hpq * scale + (is_mixed_mesh ? 1 : 0)]
                        .initialize(quad_face.second);
                    }

                  const auto face_quadrature_collection =
                    get_face_quadrature_collection(quad[my_q][hpq]);

                  face_data[my_q].q_collection[hpq] =
                    face_quadrature_collection.second;

                  reference_cell_types[my_q][hpq] =
                    face_quadrature_collection.first;
                }
              else
                {
                  cell_data[my_q].descriptor[hpq].initialize(
                    quad[my_q][hpq].get_tensor_basis()[0]);

                  const auto quad_face = quad[my_q][hpq].get_tensor_basis()[0];
                  face_data[my_q].descriptor[hpq * scale].initialize(quad_face);
                  face_data[my_q].q_collection[hpq] =
                    dealii::hp::QCollection<dim - 1>(quad_face);
                  face_data_by_cells[my_q].descriptor[hpq * scale].initialize(
                    quad_face);
                  reference_cell_types[my_q][hpq] =
                    ReferenceCells::get_hypercube<dim>();
                }
            }
        }

      // In case we have no hp-adaptivity (active_fe_index is empty), we have
      // cells, and the mapping is MappingQ or a derived class, we can
      // use the fast method.
      if (active_fe_index.empty() && !cells.empty() && mapping->size() == 1 &&
          dynamic_cast<const MappingQ<dim> *>(&mapping->operator[](0)))
        compute_mapping_q(tria, cells, face_info);
      else
        {
          // Could call these functions in parallel, but not useful because
          // the work inside is nicely split up already
          initialize_cells(tria, cells, active_fe_index, *mapping);
          initialize_faces(
            tria, cells, face_info.faces, active_fe_index, *mapping);
          initialize_faces_by_cells(tria, cells, face_info, *mapping);
        }
    }



    template <int dim, typename Number, typename VectorizedArrayType>
    void
    MappingInfo<dim, Number, VectorizedArrayType>::update_mapping(
      const dealii::Triangulation<dim>                         &tria,
      const std::vector<std::pair<unsigned int, unsigned int>> &cells,
      const FaceInfo<VectorizedArrayType::size()>              &face_info,
      const std::vector<unsigned int>                          &active_fe_index,
      const std::shared_ptr<dealii::hp::MappingCollection<dim>> &mapping)
    {
      AssertDimension(cells.size() / VectorizedArrayType::size(),
                      cell_type.size());

      for (auto &data : cell_data)
        data.clear_data_fields();
      for (auto &data : face_data)
        data.clear_data_fields();
      for (auto &data : face_data_by_cells)
        data.clear_data_fields();

      this->mapping_collection = mapping;
      this->mapping            = &mapping->operator[](0);

      if (active_fe_index.empty() && !cells.empty() && mapping->size() == 1 &&
          dynamic_cast<const MappingQ<dim> *>(&mapping->operator[](0)))
        compute_mapping_q(tria, cells, face_info);
      else
        {
          // Could call these functions in parallel, but not useful because
          // the work inside is nicely split up already
          initialize_cells(tria, cells, active_fe_index, *mapping);
          initialize_faces(
            tria, cells, face_info.faces, active_fe_index, *mapping);
          initialize_faces_by_cells(tria, cells, face_info, *mapping);
        }
    }



    // Copy a vectorized array of one type to another type
    template <typename VectorizedArrayType1, typename VectorizedArrayType2>
    inline DEAL_II_ALWAYS_INLINE void
    store_vectorized_array(const VectorizedArrayType1 value,
                           const unsigned int         offset,
                           VectorizedArrayType2      &result)
    {
      static_assert(VectorizedArrayType2::size() >=
                      VectorizedArrayType1::size(),
                    "Cannot convert to vectorized array of wider number type");

      DEAL_II_OPENMP_SIMD_PRAGMA
      for (unsigned int v = 0; v < VectorizedArrayType1::size(); ++v)
        result[offset + v] = value[v];
    }



    // For second derivatives on the real cell, we need the gradient of the
    // inverse Jacobian J. This involves some calculus and is done
    // vectorized. If L is the gradient of the Jacobian on the unit cell,
    // the gradient of the inverse is given by (multidimensional calculus) -
    // J * (J * L) * J (the third J is because we need to transform the
    // gradient L from the unit to the real cell, and then apply the inverse
    // Jacobian). Compare this with 1d with j(x) = 1/k(phi(x)), where j =
    // phi' is the inverse of the Jacobian and k is the derivative of the
    // Jacobian on the unit cell. Then j' = phi' k'/k^2 = j k' j^2.
    template <int dim, typename Number>
    Tensor<1, dim *(dim + 1) / 2, Tensor<1, dim, Number>>
    process_jacobian_gradient(const Tensor<2, dim, Number> &inv_jac_permutation,
                              const Tensor<2, dim, Number> &inv_jac,
                              const Tensor<3, dim, Number> &jac_grad)
    {
      Number inv_jac_grad[dim][dim][dim];

      // compute: inv_jac_grad = inv_jac_permutation * grad_unit(jac)
      for (unsigned int d = 0; d < dim; ++d)
        for (unsigned int e = 0; e < dim; ++e)
          for (unsigned int f = 0; f < dim; ++f)
            {
              inv_jac_grad[f][e][d] =
                (inv_jac_permutation[f][0] * jac_grad[d][e][0]);
              for (unsigned int g = 1; g < dim; ++g)
                inv_jac_grad[f][e][d] +=
                  (inv_jac_permutation[f][g] * jac_grad[d][e][g]);
            }

      // compute: transpose (-inv_jac_permutation * inv_jac_grad[d] * inv_jac)
      Number tmp[dim];
      Number grad_jac_inv[dim][dim][dim];
      for (unsigned int d = 0; d < dim; ++d)
        for (unsigned int e = 0; e < dim; ++e)
          {
            for (unsigned int f = 0; f < dim; ++f)
              {
                tmp[f] = -inv_jac_grad[d][f][0] * inv_jac[0][e];
                for (unsigned int g = 1; g < dim; ++g)
                  tmp[f] -= inv_jac_grad[d][f][g] * inv_jac[g][e];
              }

            // needed for non-diagonal part of Jacobian grad
            for (unsigned int f = 0; f < dim; ++f)
              {
                grad_jac_inv[f][d][e] = inv_jac_permutation[f][0] * tmp[0];
                for (unsigned int g = 1; g < dim; ++g)
                  grad_jac_inv[f][d][e] += inv_jac_permutation[f][g] * tmp[g];
              }
          }

      Tensor<1, dim *(dim + 1) / 2, Tensor<1, dim, Number>> result;

      // the diagonal part of Jacobian gradient comes first
      for (unsigned int d = 0; d < dim; ++d)
        for (unsigned int f = 0; f < dim; ++f)
          result[d][f] = grad_jac_inv[d][d][f];

      // then the upper-diagonal part
      for (unsigned int d = 0, count = 0; d < dim; ++d)
        for (unsigned int e = d + 1; e < dim; ++e, ++count)
          for (unsigned int f = 0; f < dim; ++f)
            result[dim + count][f] = grad_jac_inv[d][e][f];
      return result;
    }



    /* ------------------------- initialization of cells ------------------- */

    // Namespace with implementation of extraction of values on cell
    // range
    namespace ExtractCellHelper
    {
      template <int dim>
      double
      get_jacobian_size(const dealii::Triangulation<dim> &tria)
      {
        if (tria.n_cells() == 0)
          return 1;
        else
          return tria.begin()->diameter();
      }



      template <int dim, typename Number, typename VectorizedArrayType>
      struct CompressedCellData
      {
        CompressedCellData(const double expected_size)
          : data(FloatingPointComparator<VectorizedArrayType>(
              expected_size * std::numeric_limits<double>::epsilon() * 1024.))
        {}

        std::map<Tensor<2, dim, Tensor<1, VectorizedArrayType::size(), Number>>,
                 unsigned int,
                 FloatingPointComparator<VectorizedArrayType>>
          data;
      };

      /**
       * Internal temporary data used for the initialization.
       */
      template <int dim, typename Number, typename VectorizedArrayType>
      struct LocalData
      {
        LocalData(const double jac_size);
        void
        resize(const unsigned int size);

        AlignedVector<Point<dim, VectorizedArrayType>>     quadrature_points;
        AlignedVector<Tensor<2, dim, VectorizedArrayType>> general_jac;
        AlignedVector<VectorizedArrayType>                 JxW_values;
        AlignedVector<Tensor<3, dim, VectorizedArrayType>> general_jac_grad;
        AlignedVector<Tensor<1, dim, VectorizedArrayType>> normal_vectors;
        Tensor<2, dim, VectorizedArrayType>                const_jac;
        const double                                       jac_size;
      };



      template <int dim, typename Number, typename VectorizedArrayType>
      LocalData<dim, Number, VectorizedArrayType>::LocalData(
        const double jac_size_in)
        : jac_size(jac_size_in)
      {}



      template <int dim, typename Number, typename VectorizedArrayType>
      void
      LocalData<dim, Number, VectorizedArrayType>::resize(
        const unsigned int size)
      {
        if (JxW_values.size() != size)
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
      template <int dim, typename VectorizedArrayType>
      void
      evaluate_on_cell(const dealii::Triangulation<dim>            &tria,
                       const std::pair<unsigned int, unsigned int> *cells,
                       const unsigned int                           my_q,
                       GeometryType                                &cell_t_prev,
                       GeometryType                                *cell_t,
                       dealii::FEValues<dim, dim>                  &fe_val,
                       LocalData<dim,
                                 typename VectorizedArrayType::value_type,
                                 VectorizedArrayType>              &cell_data)
      {
        const unsigned int n_q_points   = fe_val.n_quadrature_points;
        const UpdateFlags  update_flags = fe_val.get_update_flags();

        cell_data.const_jac = Tensor<2, dim, VectorizedArrayType>();

        // this should be the same value as used in HashValue::scaling (but we
        // not have that field here)
        const double zero_tolerance_double =
          cell_data.jac_size * std::numeric_limits<double>::epsilon() * 1024.;
        for (unsigned int j = 0; j < VectorizedArrayType::size(); ++j)
          {
            typename dealii::Triangulation<dim>::cell_iterator cell_it(
              &tria, cells[j].first, cells[j].second);
            fe_val.reinit(cell_it);
            cell_t[j] = general;

            // extract quadrature points and store them temporarily. if we have
            // Cartesian cells, we can compress the indices
            if (update_flags & update_quadrature_points)
              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  const Point<dim> &point = fe_val.quadrature_point(q);
                  for (unsigned int d = 0; d < dim; ++d)
                    cell_data.quadrature_points[q][d][j] = point[d];
                }

            // if this is not the first quadrature formula and we already have
            // determined that this cell is either Cartesian or with constant
            // Jacobian, we have nothing more to do.
            if (my_q > 0 && cell_t_prev <= affine)
              continue;

            // first round: if the transformation is detected to be the same as
            // on the old cell, we only need to copy over the data.
            if (fe_val.get_cell_similarity() == CellSimilarity::translation &&
                my_q == 0)
              {
                if (j == 0)
                  cell_t[j] = cell_t_prev;
                else
                  cell_t[j] = cell_t[j - 1];
              }

            const DerivativeForm<1, dim, dim> &jac_0 = fe_val.jacobian(0);

            if (my_q == 0)
              {
                // check whether the Jacobian is constant on this cell the first
                // time we come around here
                if (cell_t[j] == general)
                  {
                    bool jacobian_constant = true;
                    for (unsigned int q = 1; q < n_q_points; ++q)
                      {
                        const DerivativeForm<1, dim, dim> &jac =
                          fe_val.jacobian(q);
                        for (unsigned int d = 0; d < dim; ++d)
                          for (unsigned int e = 0; e < dim; ++e)
                            if (std::fabs(jac_0[d][e] - jac[d][e]) >
                                zero_tolerance_double)
                              jacobian_constant = false;
                        if (jacobian_constant == false)
                          break;
                      }

                    // check whether the Jacobian is diagonal to machine
                    // accuracy
                    bool cell_cartesian = jacobian_constant;
                    for (unsigned int d = 0; d < dim; ++d)
                      for (unsigned int e = 0; e < dim; ++e)
                        if (d != e)
                          if (std::fabs(jac_0[d][e]) > zero_tolerance_double)
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
                    if (cell_cartesian == false && n_q_points == 1 &&
                        update_flags & update_jacobian_grads)
                      {
                        const DerivativeForm<1, dim, dim> &jac =
                          fe_val.jacobian(0);
                        const DerivativeForm<2, dim, dim> &jacobian_grad =
                          fe_val.jacobian_grad(0);
                        for (unsigned int d = 0; d < dim; ++d)
                          for (unsigned int e = 0; e < dim; ++e)
                            for (unsigned int f = 0; f < dim; ++f)
                              {
                                double jac_grad_comp =
                                  (jac[f][0] * jacobian_grad[d][e][0]);
                                for (unsigned int g = 1; g < dim; ++g)
                                  jac_grad_comp +=
                                    (jac[f][g] * jacobian_grad[d][e][g]);
                                if (std::fabs(jac_grad_comp) >
                                    zero_tolerance_double)
                                  jacobian_constant = false;
                              }
                      }
                    // set cell type
                    if (cell_cartesian == true)
                      cell_t[j] = cartesian;
                    else if (jacobian_constant == true)
                      cell_t[j] = affine;
                    else
                      cell_t[j] = general;
                  }

                // Cartesian cell
                if (cell_t[j] == cartesian)
                  {
                    // set Jacobian into diagonal (off-diagonal part is already
                    // zeroed out)
                    for (unsigned int d = 0; d < dim; ++d)
                      cell_data.const_jac[d][d][j] = jac_0[d][d];
                    continue;
                  }

                // cell with affine mapping
                else if (cell_t[j] == affine)
                  {
                    // compress out very small values
                    for (unsigned int d = 0; d < dim; ++d)
                      for (unsigned int e = 0; e < dim; ++e)
                        if (std::fabs(jac_0[d][e]) != 0.)
                          cell_data.const_jac[d][e][j] = jac_0[d][e];
                    continue;
                  }
              }

            // general cell case

            // go through all quadrature points and fill in the data into the
            // temporary data structures with slots for the vectorized data
            // types
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                // compress out very small numbers which are only noise. Then it
                // is cleaner to use zero straight away (though it does not save
                // any memory)
                const DerivativeForm<1, dim, dim> &jac = fe_val.jacobian(q);
                for (unsigned int d = 0; d < dim; ++d)
                  for (unsigned int e = 0; e < dim; ++e)
                    cell_data.general_jac[q][d][e][j] =
                      std::fabs(jac[d][e]) < zero_tolerance_double ? 0. :
                                                                     jac[d][e];

                // need to do some calculus based on the gradient of the
                // Jacobian, in order to find the gradient of the inverse
                // Jacobian which is needed in user code. however, we would like
                // to perform that on vectorized data types instead of doubles
                // or floats. to this end, copy the gradients first
                if (update_flags & update_jacobian_grads)
                  {
                    const DerivativeForm<2, dim, dim> &jacobian_grad =
                      fe_val.jacobian_grad(q);
                    for (unsigned int d = 0; d < dim; ++d)
                      for (unsigned int e = 0; e < dim; ++e)
                        for (unsigned int f = 0; f < dim; ++f)
                          cell_data.general_jac_grad[q][d][e][f][j] =
                            jacobian_grad[d][e][f];
                  }
              }
          } // end loop over entries of vectorization (size() cells)

        // set information for next cell
        cell_t_prev = cell_t[VectorizedArrayType::size() - 1];
      }



      /**
       * This function computes and tabulates the mapping information for
       * generic mappings on a range of cells by a call to FEValues with the
       * respective UpdateFlags set. There is a specialized function
       * mapping_q_compute_range function that provides a faster
       * initialization for mappings derived from MappingQ and suitable
       * quadrature formulas.
       */
      template <int dim, typename Number, typename VectorizedArrayType>
      void
      initialize_cell_range(
        const std::pair<unsigned int, unsigned int>               cell_range,
        const dealii::Triangulation<dim>                         &tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cells,
        const std::vector<unsigned int>               &active_fe_index,
        const dealii::hp::MappingCollection<dim>      &mapping,
        MappingInfo<dim, Number, VectorizedArrayType> &mapping_info,
        std::pair<
          std::vector<MappingInfoStorage<dim, dim, VectorizedArrayType>>,
          CompressedCellData<dim, Number, VectorizedArrayType>> &data)
      {
        FE_Nothing<dim> dummy_fe;

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
        LocalData<dim, Number, VectorizedArrayType> cell_data(jacobian_size);

        // encodes the cell types of the current cell. Since several cells
        // must be considered together, this variable holds the individual
        // info of the last chunk of cells
        GeometryType cell_t[VectorizedArrayType::size()];
        GeometryType cell_t_prev = general;

        // fe_values object that is used to compute the mapping data. for
        // the hp-case there might be more than one finite element. since we
        // manually select the active FE index and not via a
        // DoFHandler<dim>::active_cell_iterator, we need to manually
        // select the correct finite element, so just hold a vector of
        // FEValues
        std::vector<std::vector<std::shared_ptr<dealii::FEValues<dim>>>>
          fe_values(mapping_info.cell_data.size());

        const unsigned int max_active_fe_index =
          active_fe_index.size() > 0 ?
            *std::max_element(active_fe_index.begin(), active_fe_index.end()) :
            0;

        for (unsigned int i = 0; i < fe_values.size(); ++i)
          fe_values[i].resize(max_active_fe_index + 1);
        const UpdateFlags update_flags = mapping_info.update_flags_cells;
        const UpdateFlags update_flags_feval =
          (update_flags & update_jacobians ? update_jacobians :
                                             update_default) |
          (update_flags & update_jacobian_grads ? update_jacobian_grads :
                                                  update_default) |
          (update_flags & update_quadrature_points ? update_quadrature_points :
                                                     update_default);

        const unsigned int end_cell = std::min(mapping_info.cell_type.size(),
                                               std::size_t(cell_range.second));
        // loop over given cells
        for (unsigned int cell = cell_range.first; cell < end_cell; ++cell)
          for (unsigned int my_q = 0; my_q < mapping_info.cell_data.size();
               ++my_q)
            {
              // GENERAL OUTLINE: First generate the data in format "number"
              // for vectorization_width cells, and then find the most
              // general type of cell for appropriate vectorized formats. then
              // fill this data in
              const unsigned int fe_index =
                active_fe_index.size() > 0 ? active_fe_index[cell] : 0;
              const unsigned int hp_quad_index =
                mapping_info.cell_data[my_q].descriptor.size() == 1 ? 0 :
                                                                      fe_index;
              const unsigned int hp_mapping_index =
                mapping.size() == 1 ? 0 : fe_index;
              const unsigned int n_q_points = mapping_info.cell_data[my_q]
                                                .descriptor[hp_quad_index]
                                                .n_q_points;
              if (fe_values[my_q][fe_index].get() == nullptr)
                fe_values[my_q][fe_index] =
                  std::make_shared<dealii::FEValues<dim>>(
                    mapping[hp_mapping_index],
                    dummy_fe,
                    mapping_info.cell_data[my_q]
                      .descriptor[hp_quad_index]
                      .quadrature,
                    update_flags_feval);
              dealii::FEValues<dim> &fe_val = *fe_values[my_q][fe_index];
              cell_data.resize(n_q_points);

              // if the FE index has changed from the previous cell, set the
              // old cell type to invalid (otherwise, we might detect
              // similarity due to some cells further ahead)
              if (my_q > 0)
                cell_t_prev = GeometryType(mapping_info.cell_type[cell]);
              else if (cell > cell_range.first && active_fe_index.size() > 0 &&
                       active_fe_index[cell] != active_fe_index[cell - 1])
                cell_t_prev = general;

              evaluate_on_cell(tria,
                               &cells[cell * VectorizedArrayType::size()],
                               my_q,
                               cell_t_prev,
                               cell_t,
                               fe_val,
                               cell_data);

              // now reorder the data into vectorized types. if we are here
              // for the first time, we need to find out whether the Jacobian
              // allows for some simplification (Cartesian, affine) taking
              // vectorization_width cell together

              if (my_q == 0)
                {
                  // find the most general cell type (most general type is 3
                  // (general cell))
                  GeometryType most_general_type = cartesian;
                  for (unsigned int j = 0; j < VectorizedArrayType::size(); ++j)
                    if (cell_t[j] > most_general_type)
                      most_general_type = cell_t[j];
                  AssertIndexRange(most_general_type, 4U);
                  mapping_info.cell_type[cell] = most_general_type;
                }

              AssertThrow(
                data.first[my_q].JxW_values.size() <
                  static_cast<std::size_t>(
                    std::numeric_limits<unsigned int>::max()),
                ExcMessage(
                  "Index overflow. Cannot fit data in 32 bit integers"));

              unsigned int insert_position = data.first[my_q].JxW_values.size();
              // Cartesian/affine cell with constant Jacobians throughout the
              // cell. We need to store the data in another data field because
              // std::map cannot store data based on VectorizedArray directly
              // (alignment issue).
              if (mapping_info.cell_type[cell] <= affine)
                {
                  if (my_q == 0)
                    {
                      std::pair<
                        Tensor<2,
                               dim,
                               Tensor<1, VectorizedArrayType::size(), Number>>,
                        unsigned int>
                        new_entry;
                      // This number overlaps with the general data but we
                      // take care of that when we merge data from different
                      // threads
                      new_entry.second = data.second.data.size();
                      for (unsigned int d = 0; d < dim; ++d)
                        for (unsigned int e = 0; e < dim; ++e)
                          for (unsigned int v = 0;
                               v < VectorizedArrayType::size();
                               ++v)
                            new_entry.first[d][e][v] =
                              cell_data.const_jac[d][e][v];

                      insert_position =
                        data.second.data.insert(new_entry).first->second;
                    }
                  else
                    insert_position =
                      data.first[0].data_index_offsets[cell - cell_range.first];
                }

              // general cell case: now go through all quadrature points and
              // collect the data. done for all different quadrature formulas,
              // so do it outside the above loop.
              data.first[my_q].data_index_offsets.push_back(insert_position);
              if (mapping_info.get_cell_type(cell) == general)
                {
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    {
                      Tensor<2, dim, VectorizedArrayType> &jac =
                        cell_data.general_jac[q];
                      Tensor<3, dim, VectorizedArrayType> &jacobian_grad =
                        cell_data.general_jac_grad[q];
                      for (unsigned int j = 0; j < VectorizedArrayType::size();
                           ++j)
                        if (cell_t[j] < general)
                          {
                            for (unsigned int d = 0; d < dim; ++d)
                              for (unsigned int e = 0; e < dim; ++e)
                                {
                                  jac[d][e][j] = cell_data.const_jac[d][e][j];
                                  for (unsigned int f = 0; f < dim; ++f)
                                    jacobian_grad[d][e][f][j] = 0.;
                                }
                          }

                      data.first[my_q].JxW_values.push_back(
                        determinant(jac) * fe_val.get_quadrature().weight(q));
                      Tensor<2, dim, VectorizedArrayType> inv_jac =
                        transpose(invert(jac));
                      data.first[my_q].jacobians[0].push_back(inv_jac);

                      if (update_flags & update_jacobian_grads)
                        {
                          data.first[my_q].jacobian_gradients[0].push_back(
                            process_jacobian_gradient(inv_jac,
                                                      inv_jac,
                                                      jacobian_grad));
                          Tensor<1,
                                 dim *(dim + 1) / 2,
                                 Tensor<1, dim, VectorizedArrayType>>
                            jac_grad_sym;
                          // the diagonal part of Jacobian gradient comes
                          // first
                          for (unsigned int d = 0; d < dim; ++d)
                            for (unsigned int f = 0; f < dim; ++f)
                              jac_grad_sym[d][f] = jacobian_grad[f][d][d];

                          // then the upper-diagonal part
                          for (unsigned int d = 0, count = dim; d < dim; ++d)
                            for (unsigned int e = d + 1; e < dim; ++e, ++count)
                              for (unsigned int f = 0; f < dim; ++f)
                                jac_grad_sym[count][f] = jacobian_grad[f][d][e];

                          data.first[my_q]
                            .jacobian_gradients_non_inverse[0]
                            .push_back(jac_grad_sym);
                        }
                    }
                }

              if (update_flags & update_quadrature_points)
                {
                  // eventually we turn to the quadrature points that we can
                  // compress in case we have affine cells. we also need to
                  // reorder them into arrays of vectorized data types.  first
                  // go through the cells and find out how much memory we need
                  // to allocate for the quadrature points. We store
                  // n_q_points for general cells and a single value for
                  // Cartesian and affine cells (the position of the (0,0)
                  // point from the reference coordinates)
                  const unsigned int old_size =
                    data.first[my_q].quadrature_points.size();
                  data.first[my_q].quadrature_point_offsets.push_back(old_size);

                  if (mapping_info.get_cell_type(cell) < general)
                    {
                      Point<dim, VectorizedArrayType> quad_point;
                      for (unsigned int v = 0; v < VectorizedArrayType::size();
                           ++v)
                        {
                          typename dealii::Triangulation<dim>::cell_iterator
                            cell_it(
                              &tria,
                              cells[cell * VectorizedArrayType::size() + v]
                                .first,
                              cells[cell * VectorizedArrayType::size() + v]
                                .second);
                          const Point<dim> p =
                            mapping[hp_mapping_index]
                              .transform_unit_to_real_cell(cell_it,
                                                           Point<dim>());
                          for (unsigned int d = 0; d < dim; ++d)
                            quad_point[d][v] = p[d];
                        }
                      data.first[my_q].quadrature_points.push_back(quad_point);
                    }
                  else
                    {
                      for (unsigned int q = 0; q < n_q_points; ++q)
                        data.first[my_q].quadrature_points.push_back(
                          cell_data.quadrature_points[q]);
                    }
                }
            } // end for ( cell < end_cells )
      }



      template <typename CONTAINER>
      void
      merge_compressed_data(const CONTAINER           &source,
                            CONTAINER                 &destination,
                            std::vector<unsigned int> &indices)
      {
        indices.resize(source.size());
        typename CONTAINER::iterator lookup = destination.begin();
        for (typename CONTAINER::const_iterator it = source.begin();
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



      template <int structdim,
                int dim,
                typename Number,
                typename VectorizedArrayType>
      void
      copy_data(
        const unsigned int                first_cell,
        const std::array<std::size_t, 2> &data_shift,
        const std::vector<unsigned int>  &indices_compressed,
        const std::vector<GeometryType>  &cell_type,
        MappingInfoStorage<structdim, dim, VectorizedArrayType>
          &data_cells_local,
        MappingInfoStorage<structdim, dim, VectorizedArrayType> &data_cells)
      {
        // Copy the index offsets and shift by the appropriate value
        for (unsigned int lcell = 0;
             lcell < data_cells_local.data_index_offsets.size();
             ++lcell)
          {
            const unsigned int cell = lcell + first_cell;
            data_cells.data_index_offsets[cell] =
              cell_type[cell] <= affine ?
                (dim == structdim ? 2 : 1) *
                  indices_compressed[data_cells_local
                                       .data_index_offsets[lcell]] :
                data_cells_local.data_index_offsets[lcell] + data_shift[0];
            if (data_cells_local.quadrature_point_offsets.size() > lcell)
              data_cells.quadrature_point_offsets[cell] =
                data_cells_local.quadrature_point_offsets[lcell] +
                data_shift[1];
          }

        // Copy quadrature points
        if (data_cells.quadrature_point_offsets.empty() == false)
          {
            Point<dim, VectorizedArrayType> *out_point =
              &data_cells.quadrature_points[data_shift[1]];
            for (const Point<dim, VectorizedArrayType> *point =
                   data_cells_local.quadrature_points.begin();
                 point != data_cells_local.quadrature_points.end();
                 ++point, ++out_point)
              *out_point = *point;
            data_cells_local.quadrature_points.clear();
          }

        // If we have collected Jacobian data, copy Jacobians, JxW values,
        // Jacobian gradients
        if (data_cells_local.JxW_values.empty())
          return;

        std::copy(data_cells_local.JxW_values.begin(),
                  data_cells_local.JxW_values.end(),
                  data_cells.JxW_values.begin() + data_shift[0]);
        data_cells_local.JxW_values.clear();
        std::copy(data_cells_local.normal_vectors.begin(),
                  data_cells_local.normal_vectors.end(),
                  data_cells.normal_vectors.begin() + data_shift[0]);
        data_cells_local.normal_vectors.clear();
        for (unsigned int i = 0; i < 2; ++i)
          {
            std::copy(data_cells_local.jacobians[i].begin(),
                      data_cells_local.jacobians[i].end(),
                      data_cells.jacobians[i].begin() + data_shift[0]);
            data_cells_local.jacobians[i].clear();
            std::copy(data_cells_local.jacobian_gradients[i].begin(),
                      data_cells_local.jacobian_gradients[i].end(),
                      data_cells.jacobian_gradients[i].begin() + data_shift[0]);
            data_cells_local.jacobian_gradients[i].clear();
            std::copy(
              data_cells_local.jacobian_gradients_non_inverse[i].begin(),
              data_cells_local.jacobian_gradients_non_inverse[i].end(),
              data_cells.jacobian_gradients_non_inverse[i].begin() +
                data_shift[0]);
            data_cells_local.jacobian_gradients_non_inverse[i].clear();
            std::copy(data_cells_local.normals_times_jacobians[i].begin(),
                      data_cells_local.normals_times_jacobians[i].end(),
                      data_cells.normals_times_jacobians[i].begin() +
                        data_shift[0]);
            data_cells_local.normals_times_jacobians[i].clear();
          }
      }



      /**
       * This function is the preparatory step for the faster MappingQ-based
       * setup of the data structures in Mapping Info, using FEValues as an
       * initial query mechanism into MappingQ, storing the resulting
       * quadrature points and an initial representation of the Jacobians in
       * two arrays.
       */
      template <int dim>
      void
      mapping_q_query_fe_values(
        const unsigned int                                        begin_cell,
        const unsigned int                                        end_cell,
        const MappingQ<dim>                                      &mapping_q,
        const dealii::Triangulation<dim>                         &tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cell_array,
        const double                                              jacobian_size,
        std::vector<GeometryType> &preliminary_cell_type,
        AlignedVector<double>     &plain_quadrature_points,
        AlignedVector<std::array<Tensor<2, dim>, dim + 1>>
          &jacobians_on_stencil)
      {
        if (begin_cell == end_cell)
          return;

        const unsigned int mapping_degree = mapping_q.get_degree();
        FE_Nothing<dim>    dummy_fe;
        QGaussLobatto<dim> quadrature(mapping_degree + 1);
        const unsigned int n_mapping_points =
          Utilities::pow(mapping_degree + 1, dim);

        FEValues<dim> fe_values(mapping_q,
                                dummy_fe,
                                quadrature,
                                update_quadrature_points | update_jacobians);

        for (unsigned int cell = begin_cell; cell < end_cell; ++cell)
          {
            typename dealii::Triangulation<dim>::cell_iterator cell_it(
              &tria, cell_array[cell].first, cell_array[cell].second);
            fe_values.reinit(cell_it);
            for (unsigned int d = 0; d < dim; ++d)
              for (unsigned int q = 0; q < n_mapping_points; ++q)
                plain_quadrature_points[(cell * dim + d) * n_mapping_points +
                                        q] = fe_values.quadrature_point(q)[d];

            // store the first, second, n-th and n^2-th one along a
            // stencil-like pattern
            std::array<Tensor<2, dim, double>, dim + 1> &my_jacobians =
              jacobians_on_stencil[cell];
            my_jacobians[0] = Tensor<2, dim, double>(fe_values.jacobian(0));
            for (unsigned int d = 0, skip = 1; d < dim;
                 ++d, skip *= (mapping_degree + 1))
              my_jacobians[1 + d] =
                Tensor<2, dim, double>(fe_values.jacobian(skip));

            // check whether cell is Cartesian/affine/general
            GeometryType type = cartesian;
            for (unsigned int d = 0; d < dim; ++d)
              for (unsigned int e = 0; e < dim; ++e)
                if (d != e)
                  if (std::abs(my_jacobians[0][d][e]) > 1e-12 * jacobian_size)
                    type = affine;

            for (unsigned int q = 1; q < n_mapping_points; ++q)
              for (unsigned int d = 0; d < dim; ++d)
                for (unsigned int e = 0; e < dim; ++e)
                  if (std::abs(fe_values.jacobian(q)[d][e] -
                               fe_values.jacobian(0)[d][e]) >
                      1e-12 * jacobian_size)
                    {
                      type = general;
                      goto endloop;
                    }
          endloop:
            preliminary_cell_type[cell] = type;
          }
      }



      template <int dim>
      std::vector<unsigned int>
      mapping_q_find_compression(
        const double jacobian_size,
        const AlignedVector<std::array<Tensor<2, dim>, dim + 1>>
                                    &jacobians_on_stencil,
        const unsigned int           n_mapping_points,
        const AlignedVector<double> &plain_quadrature_points,
        std::vector<GeometryType>   &preliminary_cell_type)
      {
        std::vector<unsigned int> cell_data_index(jacobians_on_stencil.size());

        // we include a map to store some compressed information about the
        // Jacobians which we collect by a stencil-like pattern around the
        // first quadrature point on the cell - we use a relatively coarse
        // tolerance to account for some inaccuracies in the manifold
        // evaluation
        std::map<std::array<Tensor<2, dim>, dim + 1>,
                 unsigned int,
                 FloatingPointComparator<double>>
          compressed_jacobians(FloatingPointComparator<double>(
            1e4 * jacobian_size * std::numeric_limits<double>::epsilon() *
            1024.));

        unsigned int n_data_buckets = 0;
        for (unsigned int cell = 0; cell < jacobians_on_stencil.size(); ++cell)
          {
            // check in the map for the index of this cell
            const auto position =
              compressed_jacobians.find(jacobians_on_stencil[cell]);
            bool add_this_cell = position == compressed_jacobians.end();
            if (!add_this_cell)
              {
                // check if the found duplicate really is a translation and
                // the similarity identified by the map is not by accident
                double        max_distance = 0;
                const double *ptr_origin =
                  plain_quadrature_points.data() +
                  position->second * dim * n_mapping_points;
                const double *ptr_mine = plain_quadrature_points.data() +
                                         cell * dim * n_mapping_points;
                for (unsigned int d = 0; d < dim; ++d)
                  {
                    const double translate_d =
                      ptr_origin[d * n_mapping_points] -
                      ptr_mine[d * n_mapping_points];
                    for (unsigned int q = 1; q < n_mapping_points; ++q)
                      max_distance =
                        std::max(std::abs(ptr_origin[d * n_mapping_points + q] -
                                          ptr_mine[d * n_mapping_points + q] -
                                          translate_d),
                                 max_distance);
                  }

                // this is not a duplicate, must add it again
                if (max_distance > 1e-10 * jacobian_size)
                  add_this_cell = true;
              }
            else
              compressed_jacobians.insert(
                std::make_pair(jacobians_on_stencil[cell], cell));

            if (add_this_cell)
              cell_data_index[cell] = n_data_buckets++;
            else
              {
                cell_data_index[cell] = cell_data_index[position->second];
                // make sure that the cell type is the same as in the original
                // field, despite possibly small differences due to roundoff
                // and the tolerances we use
                preliminary_cell_type[cell] =
                  preliminary_cell_type[position->second];
              }
          }
        return cell_data_index;
      }


      /**
       * This function computes and tabulates the mapping information for
       * MappingQ-derived mappings on a range of cells calling into the tensor
       * product evaluators of the matrix-free framework, using a
       * polynomial expansion of the cell geometry underlying the MappingQ
       * class.
       */
      template <int dim,
                typename Number,
                typename VectorizedArrayType,
                typename VectorizedDouble>
      void
      mapping_q_compute_range(
        const unsigned int                                        begin_cell,
        const unsigned int                                        end_cell,
        const dealii::Triangulation<dim>                         &tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cell_array,
        const std::vector<GeometryType>                          &cell_type,
        const std::vector<bool>                                  &process_cell,
        const UpdateFlags            update_flags_cells,
        const AlignedVector<double> &plain_quadrature_points,
        const ShapeInfo<double>     &shape_info,
        MappingInfoStorage<dim, dim, VectorizedArrayType> &my_data)
      {
        constexpr unsigned int n_lanes   = VectorizedArrayType::size();
        constexpr unsigned int n_lanes_d = VectorizedDouble::size();

        const unsigned int n_q_points = my_data.descriptor[0].n_q_points;
        const unsigned int n_mapping_points =
          shape_info.dofs_per_component_on_cell;
        constexpr unsigned int hess_dim = dim * (dim + 1) / 2;

        FEEvaluationData<dim, VectorizedDouble, false> eval(shape_info);

        AlignedVector<VectorizedDouble> evaluation_data;
        eval.set_data_pointers(&evaluation_data, dim);

        for (unsigned int cell = begin_cell; cell < end_cell; ++cell)
          for (unsigned vv = 0; vv < n_lanes; vv += n_lanes_d)
            {
              if (cell_type[cell] > affine || process_cell[cell])
                {
                  unsigned int start_indices[n_lanes_d];
                  for (unsigned int v = 0; v < n_lanes_d; ++v)
                    start_indices[v] =
                      (cell * n_lanes + vv + v) * n_mapping_points * dim;
                  vectorized_load_and_transpose(n_mapping_points * dim,
                                                plain_quadrature_points.data(),
                                                start_indices,
                                                eval.begin_dof_values());

                  FEEvaluationFactory<dim, VectorizedDouble>::evaluate(
                    dim,
                    EvaluationFlags::values | EvaluationFlags::gradients |
                      (update_flags_cells & update_jacobian_grads ?
                         EvaluationFlags::hessians :
                         EvaluationFlags::nothing),
                    eval.begin_dof_values(),
                    eval);
                }
              if (update_flags_cells & update_quadrature_points)
                {
                  Point<dim, VectorizedArrayType> *quadrature_points =
                    my_data.quadrature_points.data() +
                    my_data.quadrature_point_offsets[cell];
                  if (cell_type[cell] <= affine)
                    for (unsigned int d = 0; d < dim; ++d)
                      for (unsigned int v = 0; v < n_lanes_d; ++v)
                        quadrature_points[0][d][vv + v] =
                          plain_quadrature_points
                            [(dim * (cell * n_lanes + vv + v) + d) *
                             n_mapping_points];
                  else
                    for (unsigned int d = 0; d < dim; ++d)
                      for (unsigned int q = 0; q < n_q_points; ++q)
                        store_vectorized_array(
                          eval.begin_values()[q + d * n_q_points],
                          vv,
                          quadrature_points[q][d]);
                }

              const unsigned int n_points =
                cell_type[cell] <= affine ? 1 : n_q_points;
              if (process_cell[cell])
                for (unsigned int q = 0; q < n_points; ++q)
                  {
                    const unsigned int idx =
                      my_data.data_index_offsets[cell] + q;
                    Tensor<2, dim, VectorizedDouble> jac;
                    for (unsigned int d = 0; d < dim; ++d)
                      for (unsigned int e = 0; e < dim; ++e)
                        jac[d][e] =
                          eval
                            .begin_gradients()[e + (d * n_q_points + q) * dim];

                    // eliminate roundoff errors
                    if (cell_type[cell] == cartesian)
                      for (unsigned int d = 0; d < dim; ++d)
                        for (unsigned int e = 0; e < dim; ++e)
                          if (d != e)
                            jac[d][e] = 0.;

                    const VectorizedDouble jac_det = determinant(jac);

                    if constexpr (running_in_debug_mode())
                      {
                        for (unsigned int v = 0; v < n_lanes_d; ++v)
                          {
                            const typename Triangulation<dim>::cell_iterator
                              cell_iterator(
                                &tria,
                                cell_array[cell * n_lanes + vv + v].first,
                                cell_array[cell * n_lanes + vv + v].second);

                            Assert(
                              jac_det[v] > 1e-12 * Utilities::fixed_power<dim>(
                                                     cell_iterator->diameter() /
                                                     std::sqrt(double(dim))),
                              (typename Mapping<dim>::ExcDistortedMappedCell(
                                cell_iterator->center(), jac_det[v], q)));
                          }
                      }
                    else
                      {
                        (void)tria;
                        (void)cell_array;
                      }

                    const Tensor<2, dim, VectorizedDouble> inv_jac =
                      transpose(invert(jac));

                    if (cell_type[cell] <= affine)
                      {
                        store_vectorized_array(jac_det,
                                               vv,
                                               my_data.JxW_values[idx]);

                        for (unsigned int d = 0; d < dim; ++d)
                          for (unsigned int e = 0; e < dim; ++e)
                            store_vectorized_array(
                              jac[d][e],
                              vv,
                              my_data.jacobians[0][idx + 1][d][e]);
                      }
                    else
                      {
                        const double weight =
                          my_data.descriptor[0].quadrature.weight(q);
                        store_vectorized_array(jac_det * weight,
                                               vv,
                                               my_data.JxW_values[idx]);
                      }
                    for (unsigned int d = 0; d < dim; ++d)
                      for (unsigned int e = 0; e < dim; ++e)
                        store_vectorized_array(inv_jac[d][e],
                                               vv,
                                               my_data.jacobians[0][idx][d][e]);

                    if (update_flags_cells & update_jacobian_grads &&
                        cell_type[cell] > affine)
                      {
                        Tensor<3, dim, VectorizedDouble> jac_grad;
                        for (unsigned int d = 0; d < dim; ++d)
                          {
                            for (unsigned int e = 0; e < dim; ++e)
                              jac_grad[d][e][e] =
                                eval.begin_hessians()[q + (d * hess_dim + e) *
                                                            n_q_points];
                            for (unsigned int c = dim, e = 0; e < dim; ++e)
                              for (unsigned int f = e + 1; f < dim; ++f, ++c)
                                jac_grad[d][e][f] = jac_grad[d][f][e] =
                                  eval.begin_hessians()[q + (d * hess_dim + c) *
                                                              n_q_points];
                            const auto inv_jac_grad =
                              process_jacobian_gradient(inv_jac,
                                                        inv_jac,
                                                        jac_grad);
                            for (unsigned int d = 0; d < hess_dim; ++d)
                              for (unsigned int e = 0; e < dim; ++e)
                                store_vectorized_array(
                                  inv_jac_grad[d][e],
                                  vv,
                                  my_data.jacobian_gradients[0][idx][d][e]);

                            // Also store the non-inverse Jacobian gradient.
                            // the diagonal part of Jacobian gradient comes
                            // first
                            for (unsigned int d = 0; d < dim; ++d)
                              for (unsigned int f = 0; f < dim; ++f)
                                store_vectorized_array(
                                  jac_grad[f][d][d],
                                  vv,
                                  my_data.jacobian_gradients_non_inverse[0][idx]
                                                                        [d][f]);

                            // then the upper-diagonal part
                            for (unsigned int d = 0, count = dim; d < dim; ++d)
                              for (unsigned int e = d + 1; e < dim;
                                   ++e, ++count)
                                for (unsigned int f = 0; f < dim; ++f)
                                  store_vectorized_array(
                                    jac_grad[f][d][e],
                                    vv,
                                    my_data.jacobian_gradients_non_inverse
                                      [0][idx][count][f]);
                          }
                      }
                  }
            }
      }

    } // namespace ExtractCellHelper



    template <int dim, typename Number, typename VectorizedArrayType>
    void
    MappingInfo<dim, Number, VectorizedArrayType>::initialize_cells(
      const dealii::Triangulation<dim>                         &tria,
      const std::vector<std::pair<unsigned int, unsigned int>> &cells,
      const std::vector<unsigned int>                          &active_fe_index,
      const dealii::hp::MappingCollection<dim>                 &mapping)
    {
      const unsigned int n_cells = cells.size();
      const unsigned int n_lanes = VectorizedArrayType::size();
      Assert(n_cells % n_lanes == 0, ExcInternalError());
      const unsigned int n_cell_batches = n_cells / n_lanes;
      cell_type.resize(n_cell_batches);

      if (n_cell_batches == 0)
        return;

      // Create as many chunks of cells as we have threads and spawn the work
      unsigned int work_per_chunk =
        std::max(8U,
                 (n_cell_batches + MultithreadInfo::n_threads() - 1) /
                   MultithreadInfo::n_threads());

      std::vector<std::pair<
        std::vector<MappingInfoStorage<dim, dim, VectorizedArrayType>>,
        ExtractCellHelper::
          CompressedCellData<dim, Number, VectorizedArrayType>>>
        data_cells_local;
      // Reserve enough space to avoid re-allocation (which would break the
      // references to the data fields passed to the tasks!)
      data_cells_local.reserve(MultithreadInfo::n_threads());

      {
        Threads::TaskGroup<>                  tasks;
        std::pair<unsigned int, unsigned int> cell_range(0U, work_per_chunk);
        while (cell_range.first < n_cell_batches)
          {
            data_cells_local.push_back(std::make_pair(
              std::vector<MappingInfoStorage<dim, dim, VectorizedArrayType>>(
                cell_data.size()),
              ExtractCellHelper::
                CompressedCellData<dim, Number, VectorizedArrayType>(
                  ExtractCellHelper::get_jacobian_size(tria))));
            tasks += Threads::new_task(
              &ExtractCellHelper::
                initialize_cell_range<dim, Number, VectorizedArrayType>,
              cell_range,
              tria,
              cells,
              active_fe_index,
              mapping,
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
      for (unsigned int i = 0; i < data_cells_local.size(); ++i)
        ExtractCellHelper::merge_compressed_data(
          data_cells_local[i].second.data,
          data_cells_local[0].second.data,
          indices_compressed[i]);

      // Collect all data in the final data fields.
      // First allocate the memory
      const unsigned int n_constant_jacobians =
        data_cells_local[0].second.data.size();
      for (unsigned int my_q = 0; my_q < cell_data.size(); ++my_q)
        {
          cell_data[my_q].data_index_offsets.resize(cell_type.size());
          std::vector<std::array<std::size_t, 2>> shift(
            data_cells_local.size());
          shift[0][0] = 2 * n_constant_jacobians;
          shift[0][1] = 0;
          for (unsigned int i = 1; i < data_cells_local.size(); ++i)
            {
              shift[i][0] =
                shift[i - 1][0] +
                data_cells_local[i - 1].first[my_q].JxW_values.size();
              shift[i][1] =
                shift[i - 1][1] +
                data_cells_local[i - 1].first[my_q].quadrature_points.size();
            }
          cell_data[my_q].JxW_values.resize_fast(
            shift.back()[0] +
            data_cells_local.back().first[my_q].JxW_values.size());
          cell_data[my_q].jacobians[0].resize_fast(
            cell_data[my_q].JxW_values.size());
          if (update_flags_cells & update_jacobian_grads)
            {
              cell_data[my_q].jacobian_gradients[0].resize_fast(
                cell_data[my_q].JxW_values.size());
              cell_data[my_q].jacobian_gradients_non_inverse[0].resize_fast(
                cell_data[my_q].JxW_values.size());
            }
          if (update_flags_cells & update_quadrature_points)
            {
              cell_data[my_q].quadrature_point_offsets.resize(cell_type.size());
              cell_data[my_q].quadrature_points.resize_fast(
                shift.back()[1] +
                data_cells_local.back().first[my_q].quadrature_points.size());
            }

          // Start tasks that copy the local data
          Threads::TaskGroup<> tasks;
          for (unsigned int i = 0; i < data_cells_local.size(); ++i)
            tasks += Threads::new_task(
              &ExtractCellHelper::
                copy_data<dim, dim, Number, VectorizedArrayType>,
              work_per_chunk * i,
              shift[i],
              indices_compressed[i],
              cell_type,
              data_cells_local[i].first[my_q],
              cell_data[my_q]);

          // finally, insert the constant cell data at the beginning (the
          // other tasks can already start copying the non-constant
          // data). Note that we use two slots for the constant data to
          // accommodate for both the inverse transposed Jacobian (that we
          // need for derivatives) and the Jacobian (that we need for
          // quadrature points)
          if (my_q == 0)
            {
              for (const auto &it : data_cells_local[0].second.data)
                {
                  Tensor<2, dim, VectorizedArrayType> jac;
                  for (unsigned int d = 0; d < dim; ++d)
                    for (unsigned int e = 0; e < dim; ++e)
                      for (unsigned int v = 0; v < VectorizedArrayType::size();
                           ++v)
                        jac[d][e][v] = it.first[d][e][v];
                  AssertIndexRange(it.second, n_constant_jacobians);
                  const std::size_t index               = it.second;
                  cell_data[my_q].JxW_values[2 * index] = determinant(jac);
                  // invert and transpose jac
                  cell_data[my_q].jacobians[0][2 * index] =
                    transpose(invert(jac));
                  cell_data[my_q].jacobians[0][2 * index + 1] = jac;
                  // second derivative of transformation is zero on affine cells
                }
            }
          else
            {
              for (unsigned int i = 0; i < 2 * n_constant_jacobians; ++i)
                {
                  cell_data[my_q].JxW_values[i] = cell_data[0].JxW_values[i];
                  cell_data[my_q].jacobians[0][i] =
                    cell_data[0].jacobians[0][i];
                }
            }

          // ... wait for the parallel work to finish
          tasks.join_all();
        }
    }



    /* ------------------------- initialization of faces ------------------- */

    // Namespace with implementation of extraction of values on face
    // range
    namespace ExtractFaceHelper
    {
      template <int dim, typename Number, typename VectorizedArrayType>
      struct CompressedFaceData
      {
        // Constructor. As a scaling factor for the FloatingPointComparator, we
        // select the inverse of the Jacobian (not the Jacobian as in the
        // CompressedCellData) and add another factor of 512 to account for
        // some roundoff effects.
        CompressedFaceData(const Number jacobian_size)
          : data(FloatingPointComparator<VectorizedArrayType>(
              512. / jacobian_size * std::numeric_limits<double>::epsilon() *
              1024.))
          , jacobian_size(jacobian_size)
        {}

        // Store the Jacobians on both sides of a face (2*(dim*dim) entries),
        // the normal vector (dim entries), and the Jacobian determinant on
        // the face (first tensor) for each entry in the vectorized array
        // (inner tensor). We cannot choose a VectorizedArray type directly
        // because std::map does not provide the necessary alignment upon
        // memory allocation.
        std::map<Tensor<1,
                        2 * dim * dim + dim + 1,
                        Tensor<1, VectorizedArrayType::size(), Number>>,
                 unsigned int,
                 FloatingPointComparator<VectorizedArrayType>>
          data;

        // Store the scaling factor
        const Number jacobian_size;
      };



      // We always put the derivative normal to the face in the last slot for
      // simpler unit cell gradient computations. This function reorders the
      // indices of the Jacobian appropriately.
      template <int dim>
      unsigned int
      reorder_face_derivative_indices(
        const unsigned int  face_no,
        const unsigned int  index,
        const ReferenceCell reference_cell = ReferenceCells::Invalid)
      {
        Assert(index < dim, ExcInternalError());

        if ((reference_cell == ReferenceCells::Invalid ||
             reference_cell == ReferenceCells::get_hypercube<dim>()) == false)
          {
            return index;
          }

        if (dim == 3)
          {
            unsigned int table[3][3] = {{1, 2, 0}, {2, 0, 1}, {0, 1, 2}};
            return table[face_no / 2][index];
          }
        else if (dim == 2)
          {
            unsigned int table[2][2] = {{1, 0}, {0, 1}};
            return table[face_no / 2][index];
          }
        else if (dim == 1)
          return 0;
        else
          Assert(false,
                 ExcNotImplemented("Not possible in dim=" +
                                   std::to_string(dim)));

        return numbers::invalid_unsigned_int;
      }



      /**
       * This function computes and tabulates the mapping information for
       * generic mappings on a range of faces by a call to FEFaceValues with
       * the respective UpdateFlags set. There is a specialized function
       * mapping_q_compute_range function that provides a faster
       * initialization for mappings derived from MappingQ and suitable
       * quadrature formulas.
       */
      template <int dim, typename Number, typename VectorizedArrayType>
      void
      initialize_face_range(
        const std::pair<unsigned int, unsigned int>               face_range,
        const dealii::Triangulation<dim>                         &tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cells,
        const std::vector<FaceToCellTopology<VectorizedArrayType::size()>>
                                                      &faces,
        const std::vector<unsigned int>               &active_fe_index,
        const dealii::hp::MappingCollection<dim>      &mapping_in,
        MappingInfo<dim, Number, VectorizedArrayType> &mapping_info,
        std::pair<
          std::vector<MappingInfoStorage<dim - 1, dim, VectorizedArrayType>>,
          CompressedFaceData<dim, Number, VectorizedArrayType>> &data)
      {
        std::vector<std::vector<std::shared_ptr<FE_Nothing<dim>>>> dummy_fe(
          mapping_info.reference_cell_types.size());
        for (unsigned int my_q = 0;
             my_q < mapping_info.reference_cell_types.size();
             ++my_q)
          {
            const unsigned int n_hp_quads =
              mapping_info.reference_cell_types[my_q].size();
            dummy_fe[my_q].resize(n_hp_quads);

            for (unsigned int hpq = 0; hpq < n_hp_quads; ++hpq)
              dummy_fe[my_q][hpq] = std::make_shared<FE_Nothing<dim>>(
                mapping_info.reference_cell_types[my_q][hpq]);
          }

        const unsigned int max_active_fe_index =
          active_fe_index.size() > 0 ?
            *std::max_element(active_fe_index.begin(), active_fe_index.end()) :
            0;

        std::vector<std::vector<std::shared_ptr<FEFaceValues<dim>>>>
          fe_face_values_container(mapping_info.face_data.size());
        for (unsigned int my_q = 0; my_q < mapping_info.face_data.size();
             ++my_q)
          fe_face_values_container[my_q].resize(max_active_fe_index + 1);

        std::vector<std::vector<std::shared_ptr<FEFaceValues<dim>>>>
          fe_boundary_face_values_container(mapping_info.face_data.size());
        for (unsigned int my_q = 0; my_q < mapping_info.face_data.size();
             ++my_q)
          fe_boundary_face_values_container[my_q].resize(max_active_fe_index +
                                                         1);

        std::vector<std::vector<std::shared_ptr<FESubfaceValues<dim>>>>
          fe_subface_values_container(mapping_info.face_data.size());
        for (unsigned int my_q = 0; my_q < mapping_info.face_data.size();
             ++my_q)
          fe_subface_values_container[my_q].resize(
            mapping_info.face_data[my_q].descriptor.size());

        ExtractCellHelper::LocalData<dim, Number, VectorizedArrayType>
          face_data(ExtractCellHelper::get_jacobian_size(tria));

        const unsigned int end_face =
          std::min(std::size_t(face_range.second), faces.size());
        for (unsigned int face = face_range.first; face < end_face; ++face)
          for (unsigned int my_q = 0; my_q < mapping_info.face_data.size();
               ++my_q)
            {
              // We assume that we have the faces sorted by the active FE
              // indices so that the active FE index of the interior side of the
              // face batch is the same as the FE index of the interior side of
              // its first entry.
              const unsigned int fe_index =
                active_fe_index.size() > 0 ?
                  active_fe_index[faces[face].cells_interior[0] /
                                  VectorizedArrayType::size()] :
                  0;
              const unsigned int hp_quad_index =
                mapping_info.cell_data[my_q].descriptor.size() == 1 ? 0 :
                                                                      fe_index;
              const unsigned int hp_mapping_index =
                mapping_in.size() == 1 ? 0 : fe_index;

              const auto &mapping = mapping_in[hp_mapping_index];
              const auto &quadrature =
                mapping_info.face_data[my_q].q_collection[hp_quad_index];

              const bool is_boundary_face =
                faces[face].cells_exterior[0] == numbers::invalid_unsigned_int;

              if (is_boundary_face &&
                  fe_boundary_face_values_container[my_q][fe_index] == nullptr)
                fe_boundary_face_values_container[my_q][fe_index] =
                  std::make_shared<FEFaceValues<dim>>(
                    mapping,
                    *dummy_fe[my_q][hp_quad_index],
                    quadrature,
                    mapping_info.update_flags_boundary_faces);
              else if (fe_face_values_container[my_q][fe_index] == nullptr)
                fe_face_values_container[my_q][fe_index] =
                  std::make_shared<FEFaceValues<dim>>(
                    mapping,
                    *dummy_fe[my_q][hp_quad_index],
                    quadrature,
                    mapping_info.update_flags_inner_faces);

              FEFaceValues<dim> &fe_face_values =
                is_boundary_face ?
                  *fe_boundary_face_values_container[my_q][fe_index] :
                  *fe_face_values_container[my_q][fe_index];

              unsigned int n_q_points = 0; // will be override once FEFaceValues
                                           // is set up

              bool normal_is_similar = true;
              bool JxW_is_similar    = true;
              bool cell_is_cartesian = true;
              for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
                {
                  Tensor<2, dim> jacobian_0;
                  double         compare_norm_jac = 1.;
                  if (faces[face].cells_interior[v] !=
                      numbers::invalid_unsigned_int)
                    {
                      typename dealii::Triangulation<dim>::cell_iterator
                        cell_it(&tria,
                                cells[faces[face].cells_interior[v]].first,
                                cells[faces[face].cells_interior[v]].second);

                      fe_face_values.reinit(cell_it,
                                            faces[face].interior_face_no);

                      if (v == 0)
                        {
                          n_q_points = fe_face_values.n_quadrature_points;
                          face_data.resize(n_q_points);
                        }

                      for (unsigned int q = 0; q < n_q_points; ++q)
                        {
                          if (std::abs(
                                fe_face_values.JxW(q) *
                                  fe_face_values.get_quadrature().weight(0) -
                                fe_face_values.JxW(0) *
                                  fe_face_values.get_quadrature().weight(q)) >
                              2048. * std::numeric_limits<double>::epsilon() *
                                fe_face_values.JxW(0) *
                                fe_face_values.get_quadrature().weight(q))
                            JxW_is_similar = false;
                          face_data.JxW_values[q][v] = fe_face_values.JxW(q);

                          DerivativeForm<1, dim, dim> inv_jac =
                            fe_face_values.jacobian(q).covariant_form();
                          Tensor<1, dim> normal =
                            fe_face_values.normal_vector(q);

                          // Filter out very small values in normal. No need
                          // to re-normalize because these values cannot enter
                          // the norm significantly: Total size is 1 but 1e-13
                          // squared is 1e-26.
                          for (unsigned int d = 0; d < dim; ++d)
                            if (std::abs(normal[d]) <
                                1024. * std::numeric_limits<double>::epsilon())
                              normal[d] = 0.;

                          if (q == 0)
                            {
                              jacobian_0       = inv_jac;
                              compare_norm_jac = jacobian_0.norm();
                            }
                          for (unsigned int d = 0; d < dim; ++d)
                            for (unsigned int e = 0; e < dim; ++e)
                              {
                                if (std::abs(inv_jac[d][e] - jacobian_0[d][e]) >
                                    2048. *
                                      std::numeric_limits<double>::epsilon() *
                                      compare_norm_jac)
                                  JxW_is_similar = false;
                                const unsigned int ee =
                                  reorder_face_derivative_indices<dim>(
                                    faces[face].interior_face_no,
                                    e,
                                    cell_it->reference_cell());
                                face_data.general_jac[q][d][e][v] =
                                  inv_jac[d][ee];
                              }

                          for (unsigned int d = 0; d < dim; ++d)
                            {
                              face_data.normal_vectors[q][d][v] = normal[d];
                              if (std::abs(normal[d] -
                                           fe_face_values.normal_vector(0)[d]) >
                                  1024. *
                                    std::numeric_limits<double>::epsilon())
                                normal_is_similar = false;
                            }
                          if (std::abs(
                                std::abs(
                                  normal[faces[face].interior_face_no / 2]) -
                                1.) >
                              1024. * std::numeric_limits<double>::epsilon())
                            cell_is_cartesian = false;
                          for (unsigned int d = 0; d < dim; ++d)
                            for (unsigned int e = 0; e < dim; ++e)
                              if (e != d &&
                                  std::abs(inv_jac[d][e]) >
                                    2048 *
                                      std::numeric_limits<double>::epsilon() *
                                      compare_norm_jac)
                                cell_is_cartesian = false;

                          if (fe_face_values.get_update_flags() &
                              update_quadrature_points)
                            for (unsigned int d = 0; d < dim; ++d)
                              face_data.quadrature_points[q][d][v] =
                                fe_face_values.quadrature_point(q)[d];
                        }
                    }
                  // Fill up with data of the zeroth component to avoid
                  // false negatives when checking for similarities
                  else
                    {
                      for (unsigned int q = 0; q < n_q_points; ++q)
                        {
                          face_data.JxW_values[q][v] =
                            face_data.JxW_values[q][0];
                          for (unsigned int d = 0; d < dim; ++d)
                            for (unsigned int e = 0; e < dim; ++e)
                              {
                                face_data.general_jac[q][d][e][v] =
                                  face_data.general_jac[q][d][e][0];
                                face_data.normal_vectors[q][d][v] =
                                  face_data.normal_vectors[q][d][0];
                              }
                          if (fe_face_values.get_update_flags() &
                              update_quadrature_points)
                            for (unsigned int d = 0; d < dim; ++d)
                              face_data.quadrature_points[q][d][v] =
                                face_data.quadrature_points[q][d][0];
                        }
                    }
                  if (is_boundary_face == false &&
                      faces[face].cells_exterior[v] !=
                        numbers::invalid_unsigned_int)
                    {
                      typename dealii::Triangulation<dim>::cell_iterator
                        cell_it(&tria,
                                cells[faces[face].cells_exterior[v]].first,
                                cells[faces[face].cells_exterior[v]].second);

                      const FEValuesBase<dim> *actual_fe_face_values = nullptr;
                      if (faces[face].subface_index >=
                          GeometryInfo<dim>::max_children_per_cell)
                        {
                          // We assume that we have the faces sorted by the
                          // active FE indices so that the active FE index of
                          // the exterior side of the face batch is the same as
                          // the FE index of the exterior side of its first
                          // entry.
                          const unsigned int fe_index =
                            active_fe_index.size() > 0 ?
                              active_fe_index[faces[face].cells_exterior[0] /
                                              VectorizedArrayType::size()] :
                              0;
                          const unsigned int hp_quad_index =
                            mapping_info.cell_data[my_q].descriptor.size() ==
                                1 ?
                              0 :
                              fe_index;
                          const unsigned int hp_mapping_index =
                            mapping_in.size() == 1 ? 0 : fe_index;

                          const auto &mapping = mapping_in[hp_mapping_index];
                          const auto &quadrature =
                            mapping_info.face_data[my_q]
                              .q_collection[hp_quad_index];

                          if (fe_face_values_container[my_q][fe_index] ==
                              nullptr)
                            fe_face_values_container[my_q][fe_index] =
                              std::make_shared<FEFaceValues<dim>>(
                                mapping,
                                *dummy_fe[my_q][hp_quad_index],
                                quadrature,
                                mapping_info.update_flags_boundary_faces);

                          fe_face_values_container[my_q][fe_index]->reinit(
                            cell_it, faces[face].exterior_face_no);

                          actual_fe_face_values =
                            fe_face_values_container[my_q][fe_index].get();
                        }
                      else
                        {
                          if (fe_subface_values_container[my_q][0] == nullptr)
                            fe_subface_values_container[my_q][0] =
                              std::make_shared<FESubfaceValues<dim>>(
                                mapping,
                                *dummy_fe[my_q][hp_quad_index],
                                quadrature,
                                mapping_info.update_flags_inner_faces);
                          fe_subface_values_container[my_q][0]->reinit(
                            cell_it,
                            faces[face].exterior_face_no,
                            faces[face].subface_index);
                          actual_fe_face_values =
                            fe_subface_values_container[my_q][0].get();
                        }
                      for (unsigned int q = 0; q < n_q_points; ++q)
                        {
                          DerivativeForm<1, dim, dim> inv_jac =
                            actual_fe_face_values->jacobian(q).covariant_form();
                          if (q == 0)
                            {
                              jacobian_0       = inv_jac;
                              compare_norm_jac = jacobian_0.norm();
                            }
                          for (unsigned int d = 0; d < dim; ++d)
                            for (unsigned int e = 0; e < dim; ++e)
                              {
                                if (std::abs(inv_jac[d][e] - jacobian_0[d][e]) >
                                    2048. *
                                      std::numeric_limits<double>::epsilon() *
                                      compare_norm_jac)
                                  JxW_is_similar = false;
                                const unsigned int ee =
                                  reorder_face_derivative_indices<dim>(
                                    faces[face].exterior_face_no,
                                    e,
                                    cell_it->reference_cell());
                                face_data.general_jac[n_q_points + q][d][e][v] =
                                  inv_jac[d][ee];
                              }
                          for (unsigned int d = 0; d < dim; ++d)
                            for (unsigned int e = 0; e < dim; ++e)
                              if (e != d &&
                                  std::abs(inv_jac[d][e]) >
                                    2048 *
                                      std::numeric_limits<double>::epsilon() *
                                      compare_norm_jac)
                                cell_is_cartesian = false;
                        }
                    }
                  // Fill up with 'known' values
                  else if (is_boundary_face == false)
                    {
                      Assert(faces[face].cells_exterior[0] !=
                               numbers::invalid_unsigned_int,
                             ExcInternalError());
                      for (unsigned int q = 0; q < n_q_points; ++q)
                        for (unsigned int d = 0; d < dim; ++d)
                          for (unsigned int e = 0; e < dim; ++e)
                            face_data.general_jac[n_q_points + q][d][e][v] =
                              face_data.general_jac[n_q_points + q][d][e][0];
                    }
                  // If boundary face, simply set the data to zero (will not
                  // be used). Note that faces over periodic boundary
                  // conditions will be treated as interior ones in this setup.
                  else
                    for (unsigned int q = 0; q < n_q_points; ++q)
                      for (unsigned int d = 0; d < dim; ++d)
                        for (unsigned int e = 0; e < dim; ++e)
                          face_data.general_jac[n_q_points + q][d][e][v] = 0.;
                }

              // check if face is affine or at least if it is flat
              // (i.e., all normal vectors are the same)
              if (my_q == 0)
                {
                  GeometryType face_type = affine;
                  if (JxW_is_similar == false)
                    face_type = flat_faces;
                  if (normal_is_similar == false)
                    face_type = general;
                  if (face_type == affine && cell_is_cartesian)
                    face_type = cartesian;
                  mapping_info.face_type[face] = face_type;
                }

              // Fill in quadrature points
              if (fe_face_values.get_update_flags() & update_quadrature_points)
                {
                  data.first[my_q].quadrature_point_offsets.push_back(
                    data.first[my_q].quadrature_points.size());
                  if (fe_face_values.get_update_flags() &
                      update_quadrature_points)
                    for (unsigned int q = 0; q < n_q_points; ++q)
                      data.first[my_q].quadrature_points.push_back(
                        face_data.quadrature_points[q]);
                }

              using VEC_ARRAY = Tensor<1, VectorizedArrayType::size(), Number>;
              unsigned int insert_position = data.first[my_q].JxW_values.size();

              // Fill in JxW values, apply compression
              if (mapping_info.face_type[face] <= affine)
                {
                  if (my_q == 0)
                    {
                      // find out if we already had the same JxW values before
                      std::pair<Tensor<1, 2 * dim * dim + dim + 1, VEC_ARRAY>,
                                unsigned int>
                        new_entry;
                      new_entry.second = data.second.data.size();
                      for (unsigned int v = 0; v < VectorizedArrayType::size();
                           ++v)
                        new_entry.first[2 * dim * dim + dim][v] =
                          face_data.JxW_values[0][v] /
                          fe_face_values.get_quadrature().weight(0) /
                          Utilities::fixed_power<dim>(face_data.jac_size);

                      new_entry.second = data.second.data.size();
                      for (unsigned int d = 0; d < dim; ++d)
                        for (unsigned int e = 0; e < dim; ++e)
                          for (unsigned int v = 0;
                               v < VectorizedArrayType::size();
                               ++v)
                            new_entry.first[d * dim + e][v] =
                              face_data.general_jac[0][d][e][v];
                      if (is_boundary_face == false)
                        for (unsigned int d = 0; d < dim; ++d)
                          for (unsigned int e = 0; e < dim; ++e)
                            for (unsigned int v = 0;
                                 v < VectorizedArrayType::size();
                                 ++v)
                              new_entry.first[dim * dim + d * dim + e][v] =
                                face_data.general_jac[n_q_points][d][e][v];
                      // we need to add the normal vector here because we
                      // store both the inverse jacobian and the normal vector
                      // times the Jacobian; of course, there will be
                      // different values in their product for normal vectors
                      // oriented in different ways (the memory saving is
                      // still significant); we need to divide by the Jacobian
                      // size to get the right scaling
                      for (unsigned int d = 0; d < dim; ++d)
                        for (unsigned int v = 0;
                             v < VectorizedArrayType::size();
                             ++v)
                          new_entry.first[2 * dim * dim + d][v] =
                            face_data.normal_vectors[0][d][v] /
                            face_data.jac_size;

                      insert_position =
                        data.second.data.insert(new_entry).first->second;
                    }
                  else
                    insert_position = data.first[0].data_index_offsets[face];
                }
              data.first[my_q].data_index_offsets.push_back(insert_position);
              if (mapping_info.face_type[face] > affine)
                {
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    data.first[my_q].JxW_values.push_back(
                      face_data.JxW_values[q]);
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    data.first[my_q].normal_vectors.push_back(
                      face_data.normal_vectors[q]);
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    data.first[my_q].jacobians[0].push_back(
                      face_data.general_jac[q]);
                  if (is_boundary_face == false)
                    for (unsigned int q = 0; q < n_q_points; ++q)
                      data.first[my_q].jacobians[1].push_back(
                        face_data.general_jac[n_q_points + q]);
                }
            }
      }

      template <int dim, typename Number, typename VectorizedArrayType>
      void
      compute_normal_times_jacobian(
        const unsigned int               first_face,
        const unsigned int               last_face,
        const std::vector<GeometryType> &face_type,
        const std::vector<FaceToCellTopology<VectorizedArrayType::size()>>
                                                              &faces,
        MappingInfoStorage<dim - 1, dim, VectorizedArrayType> &data_faces)
      {
        for (unsigned int face = first_face; face < last_face; ++face)
          {
            const bool is_boundary_face =
              faces[face].cells_exterior[0] == numbers::invalid_unsigned_int;
            const unsigned int n_q_points_work =
              face_type[face] > affine ? data_faces.descriptor[0].n_q_points :
                                         1;
            const unsigned int offset = data_faces.data_index_offsets[face];

            for (unsigned int q = 0; q < n_q_points_work; ++q)
              {
                data_faces.normals_times_jacobians[0][offset + q] =
                  data_faces.normal_vectors[offset + q] *
                  data_faces.jacobians[0][offset + q];
                if (is_boundary_face == false)
                  data_faces.normals_times_jacobians[1][offset + q] =
                    data_faces.normal_vectors[offset + q] *
                    data_faces.jacobians[1][offset + q];
              }
          }
      }



      /**
       * This function computes and tabulates the mapping information for
       * MappingQ-derived mappings on a range of faces calling into the tensor
       * product face evaluators interpolators of the matrix-free framework,
       * using a polynomial expansion of the cell geometry underlying the
       * MappingQ class.
       */
      template <int dim,
                typename Number,
                typename VectorizedArrayType,
                typename VectorizedDouble>
      void
      mapping_q_compute_range(
        const unsigned int begin_face,
        const unsigned int end_face,
        const std::vector<FaceToCellTopology<VectorizedArrayType::size()>>
                                        &faces,
        const std::vector<GeometryType> &face_type,
        const std::vector<bool>         &process_face,
        const UpdateFlags                update_flags_faces,
        const AlignedVector<double>     &plain_quadrature_points,
        const ShapeInfo<double>         &shape_info,
        MappingInfoStorage<dim - 1, dim, VectorizedArrayType> &my_data)
      {
        constexpr unsigned int n_lanes   = VectorizedArrayType::size();
        constexpr unsigned int n_lanes_d = VectorizedDouble::size();

        const unsigned int n_q_points = my_data.descriptor[0].n_q_points;
        const unsigned int n_mapping_points =
          shape_info.dofs_per_component_on_cell;
        constexpr unsigned int hess_dim = dim * (dim + 1) / 2;

        FEEvaluationData<dim, VectorizedDouble, true> eval_int(shape_info,
                                                               true);
        FEEvaluationData<dim, VectorizedDouble, true> eval_ext(shape_info,
                                                               false);

        AlignedVector<VectorizedDouble> evaluation_data_int,
          evaluation_data_ext;
        eval_int.set_data_pointers(&evaluation_data_int, dim);
        eval_ext.set_data_pointers(&evaluation_data_ext, dim);

        for (unsigned int face = begin_face; face < end_face; ++face)
          for (unsigned vv = 0; vv < n_lanes; vv += n_lanes_d)
            {
              FaceToCellTopology<VectorizedDouble::size()> face_double = {};
              face_double.interior_face_no = faces[face].interior_face_no;
              face_double.exterior_face_no = faces[face].exterior_face_no;
              face_double.face_orientation = faces[face].face_orientation;
              face_double.subface_index    = faces[face].subface_index;

              // load the geometry field for all SIMD lanes
              unsigned int start_indices[n_lanes_d];
              for (unsigned int v = 0; v < n_lanes_d; ++v)
                if (faces[face].cells_interior[vv + v] !=
                    numbers::invalid_unsigned_int)
                  start_indices[v] =
                    faces[face].cells_interior[vv + v] * n_mapping_points * dim;
                else
                  start_indices[v] =
                    faces[face].cells_interior[0] * n_mapping_points * dim;

              eval_int.reinit_face(face_double);
              vectorized_load_and_transpose(n_mapping_points * dim,
                                            plain_quadrature_points.data(),
                                            start_indices,
                                            eval_int.begin_dof_values());

              // now let the matrix-free evaluators provide us with the
              // data on faces
              FEFaceEvaluationFactory<dim, VectorizedDouble>::evaluate(
                dim,
                EvaluationFlags::values | EvaluationFlags::gradients |
                  (update_flags_faces & update_jacobian_grads ?
                     EvaluationFlags::hessians :
                     EvaluationFlags::nothing),
                eval_int.begin_dof_values(),
                eval_int);

              if (update_flags_faces & update_quadrature_points)
                for (unsigned int q = 0; q < n_q_points; ++q)
                  for (unsigned int d = 0; d < dim; ++d)
                    store_vectorized_array(
                      eval_int.begin_values()[d * n_q_points + q],
                      vv,
                      my_data.quadrature_points
                        [my_data.quadrature_point_offsets[face] + q][d]);

              if (process_face[face] == false)
                continue;

              const unsigned int offset = my_data.data_index_offsets[face];

              const auto compute_jacobian_grad =
                [&](unsigned int                                   face_no,
                    int                                            is_exterior,
                    unsigned int                                   q,
                    Tensor<2, dim, VectorizedDouble>               inv_jac,
                    FEEvaluationData<dim, VectorizedDouble, true> &eval) {
                  Tensor<2, dim, VectorizedDouble> inv_transp_jac_permutation;
                  for (unsigned int d = 0; d < dim; ++d)
                    for (unsigned int e = 0; e < dim; ++e)
                      {
                        const unsigned int ee =
                          ExtractFaceHelper::reorder_face_derivative_indices<
                            dim>(face_no, e);
                        inv_transp_jac_permutation[d][e] = inv_jac[ee][d];
                      }
                  Tensor<2, dim, VectorizedDouble> jacobi;
                  for (unsigned int e = 0; e < dim; ++e)
                    for (unsigned int d = 0; d < dim; ++d)
                      jacobi[d][e] =
                        eval.begin_gradients()[(d * n_q_points + q) * dim + e];
                  Tensor<2, dim, VectorizedDouble> inv_transp_jac =
                    transpose(invert(jacobi));
                  Tensor<3, dim, VectorizedDouble> jac_grad;
                  for (unsigned int d = 0; d < dim; ++d)
                    {
                      for (unsigned int e = 0; e < dim; ++e)
                        jac_grad[d][e][e] =
                          eval.begin_hessians()[q + (d * hess_dim + e) *
                                                      n_q_points];
                      for (unsigned int c = dim, e = 0; e < dim; ++e)
                        for (unsigned int f = e + 1; f < dim; ++f, ++c)
                          jac_grad[d][e][f] = jac_grad[d][f][e] =
                            eval.begin_hessians()[q + (d * hess_dim + c) *
                                                        n_q_points];
                      const auto inv_jac_grad =
                        process_jacobian_gradient(inv_transp_jac_permutation,
                                                  inv_transp_jac,
                                                  jac_grad);
                      for (unsigned int e = 0; e < dim; ++e)
                        {
                          for (unsigned int d = 0; d < hess_dim; ++d)
                            store_vectorized_array(
                              inv_jac_grad[d][e],
                              vv,
                              my_data.jacobian_gradients[is_exterior]
                                                        [offset + q][d][e]);
                        }

                      // Also store the non-inverse Jacobian gradient.
                      // the diagonal part of Jacobian gradient comes first.
                      // jac_grad already has its derivatives reordered,
                      // so no need to compensate for this here
                      for (unsigned int d = 0; d < dim; ++d)
                        for (unsigned int f = 0; f < dim; ++f)
                          store_vectorized_array(
                            jac_grad[f][d][d],
                            vv,
                            my_data.jacobian_gradients_non_inverse[is_exterior]
                                                                  [offset + q]
                                                                  [d][f]);

                      // then the upper-diagonal part
                      for (unsigned int d = 0, count = dim; d < dim; ++d)
                        for (unsigned int e = d + 1; e < dim; ++e, ++count)
                          for (unsigned int f = 0; f < dim; ++f)
                            store_vectorized_array(
                              jac_grad[f][d][e],
                              vv,
                              my_data.jacobian_gradients_non_inverse
                                [is_exterior][offset + q][count][f]);
                    }
                };

              // go through the faces and fill the result
              const unsigned int n_points_compute =
                face_type[face] <= affine ? 1 : n_q_points;
              for (unsigned int q = 0; q < n_points_compute; ++q)
                {
                  const unsigned int interior_face_no =
                    faces[face].interior_face_no;
                  Tensor<2, dim, VectorizedDouble> jac;
                  for (unsigned int e = 0; e < dim; ++e)
                    {
                      const unsigned int ee =
                        ExtractFaceHelper::reorder_face_derivative_indices<dim>(
                          interior_face_no, e);
                      for (unsigned int d = 0; d < dim; ++d)
                        jac[d][ee] =
                          eval_int
                            .begin_gradients()[(d * n_q_points + q) * dim + e];
                    }
                  Tensor<2, dim, VectorizedDouble> inv_jac = invert(jac);
                  for (unsigned int e = 0; e < dim; ++e)
                    {
                      const unsigned int ee =
                        ExtractFaceHelper::reorder_face_derivative_indices<dim>(
                          interior_face_no, e);
                      for (unsigned int d = 0; d < dim; ++d)
                        store_vectorized_array(
                          inv_jac[ee][d],
                          vv,
                          my_data.jacobians[0][offset + q][d][e]);
                    }
                  if (face_type[face] <= affine)
                    for (unsigned int e = 0; e < dim; ++e)
                      {
                        const unsigned int ee =
                          ExtractFaceHelper::reorder_face_derivative_indices<
                            dim>(interior_face_no, e);
                        for (unsigned int d = 0; d < dim; ++d)
                          store_vectorized_array(
                            jac[d][ee],
                            vv,
                            my_data.jacobians[0][offset + q + 1][d][e]);
                      }

                  if (update_flags_faces & update_jacobian_grads)
                    {
                      compute_jacobian_grad(
                        interior_face_no, 0, q, inv_jac, eval_int);
                    }

                  std::array<Tensor<1, dim, VectorizedDouble>, dim - 1>
                    tangential_vectors;
                  for (unsigned int d = 0; d != dim - 1; ++d)
                    for (unsigned int e = 0; e < dim; ++e)
                      for (unsigned int f = 0; f < dim; ++f)
                        tangential_vectors[d][e] +=
                          jac[e][f] * GeometryInfo<dim>::unit_tangential_vectors
                                        [interior_face_no][d][f];

                  Tensor<1, dim, VectorizedDouble> boundary_form;
                  if (dim == 1)
                    boundary_form[0] = interior_face_no == 0 ? -1. : 1.;
                  else if (dim == 2)
                    boundary_form = cross_product_2d(tangential_vectors[0]);
                  else if (dim == 3)
                    boundary_form = cross_product_3d(tangential_vectors[0],
                                                     tangential_vectors[1]);
                  else
                    DEAL_II_NOT_IMPLEMENTED();

                  const VectorizedDouble JxW =
                    boundary_form.norm() *
                    (face_type[face] <= affine ?
                       1. :
                       my_data.descriptor[0].quadrature.weight(q));

                  if constexpr (running_in_debug_mode())
                    {
                      for (unsigned int v = 0; v < n_lanes_d; ++v)
                        Assert(JxW[v] > 0.0, ExcInternalError());
                    }

                  store_vectorized_array(JxW,
                                         vv,
                                         my_data.JxW_values[offset + q]);

                  const Tensor<1, dim, VectorizedDouble> normal =
                    boundary_form / boundary_form.norm();

                  for (unsigned int d = 0; d < dim; ++d)
                    store_vectorized_array(
                      normal[d], vv, my_data.normal_vectors[offset + q][d]);

                  my_data.normals_times_jacobians[0][offset + q] =
                    my_data.normal_vectors[offset + q] *
                    my_data.jacobians[0][offset + q];
                }

              if (faces[face].cells_exterior[0] !=
                  numbers::invalid_unsigned_int)
                {
                  for (unsigned int v = 0; v < n_lanes_d; ++v)
                    if (faces[face].cells_exterior[vv + v] !=
                        numbers::invalid_unsigned_int)
                      start_indices[v] = faces[face].cells_exterior[vv + v] *
                                         n_mapping_points * dim;
                    else
                      start_indices[v] =
                        faces[face].cells_exterior[0] * n_mapping_points * dim;

                  eval_ext.reinit_face(face_double);
                  vectorized_load_and_transpose(n_mapping_points * dim,
                                                plain_quadrature_points.data(),
                                                start_indices,
                                                eval_ext.begin_dof_values());

                  FEFaceEvaluationFactory<dim, VectorizedDouble>::evaluate(
                    dim,
                    EvaluationFlags::values | EvaluationFlags::gradients |
                      (update_flags_faces & update_jacobian_grads ?
                         EvaluationFlags::hessians :
                         EvaluationFlags::nothing),
                    eval_ext.begin_dof_values(),
                    eval_ext);

                  for (unsigned int q = 0; q < n_points_compute; ++q)
                    {
                      const unsigned int exterior_face_no =
                        faces[face].exterior_face_no;
                      Tensor<2, dim, VectorizedDouble> jac;
                      for (unsigned int e = 0; e < dim; ++e)
                        {
                          const unsigned int ee =
                            ExtractFaceHelper::reorder_face_derivative_indices<
                              dim>(exterior_face_no, e);
                          for (unsigned int d = 0; d < dim; ++d)
                            jac[d][ee] =
                              eval_ext
                                .begin_gradients()[(d * n_q_points + q) * dim +
                                                   e];
                        }
                      Tensor<2, dim, VectorizedDouble> inv_jac = invert(jac);
                      for (unsigned int e = 0; e < dim; ++e)
                        {
                          const unsigned int ee =
                            ExtractFaceHelper::reorder_face_derivative_indices<
                              dim>(exterior_face_no, e);
                          for (unsigned int d = 0; d < dim; ++d)
                            store_vectorized_array(
                              inv_jac[ee][d],
                              vv,
                              my_data.jacobians[1][offset + q][d][e]);
                        }
                      if (face_type[face] <= affine)
                        for (unsigned int e = 0; e < dim; ++e)
                          {
                            const unsigned int ee = ExtractFaceHelper::
                              reorder_face_derivative_indices<dim>(
                                exterior_face_no, e);
                            for (unsigned int d = 0; d < dim; ++d)
                              store_vectorized_array(
                                jac[d][ee],
                                vv,
                                my_data.jacobians[1][offset + q + 1][d][e]);
                          }

                      if (update_flags_faces & update_jacobian_grads)
                        {
                          compute_jacobian_grad(
                            exterior_face_no, 1, q, inv_jac, eval_ext);
                        }

                      my_data.normals_times_jacobians[1][offset + q] =
                        my_data.normal_vectors[offset + q] *
                        my_data.jacobians[1][offset + q];
                    }
                }
            }
      }

    } // namespace ExtractFaceHelper



    template <int dim, typename Number, typename VectorizedArrayType>
    void
    MappingInfo<dim, Number, VectorizedArrayType>::initialize_faces(
      const dealii::Triangulation<dim>                                   &tria,
      const std::vector<std::pair<unsigned int, unsigned int>>           &cells,
      const std::vector<FaceToCellTopology<VectorizedArrayType::size()>> &faces,
      const std::vector<unsigned int>          &active_fe_index,
      const dealii::hp::MappingCollection<dim> &mapping)
    {
      face_type.resize(faces.size(), general);

      if (faces.empty())
        return;

      // Create as many chunks of cells as we have threads and spawn the
      // work
      unsigned int work_per_chunk =
        std::max(std::size_t(8),
                 (faces.size() + MultithreadInfo::n_threads() - 1) /
                   MultithreadInfo::n_threads());

      std::vector<std::pair<
        std::vector<MappingInfoStorage<dim - 1, dim, VectorizedArrayType>>,
        ExtractFaceHelper::
          CompressedFaceData<dim, Number, VectorizedArrayType>>>
        data_faces_local;
      // Reserve enough space to avoid re-allocation (which would destroy
      // the references passed to the tasks!)
      data_faces_local.reserve(MultithreadInfo::n_threads());

      {
        Threads::TaskGroup<>                  tasks;
        std::pair<unsigned int, unsigned int> face_range(0U, work_per_chunk);
        while (face_range.first < faces.size())
          {
            data_faces_local.push_back(std::make_pair(
              std::vector<
                MappingInfoStorage<dim - 1, dim, VectorizedArrayType>>(
                face_data.size()),
              ExtractFaceHelper::
                CompressedFaceData<dim, Number, VectorizedArrayType>(
                  ExtractCellHelper::get_jacobian_size(tria))));
            tasks += Threads::new_task(
              &ExtractFaceHelper::
                initialize_face_range<dim, Number, VectorizedArrayType>,
              face_range,
              tria,
              cells,
              faces,
              active_fe_index,
              mapping,
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
      for (unsigned int i = 0; i < data_faces_local.size(); ++i)
        ExtractCellHelper::merge_compressed_data(
          data_faces_local[i].second.data,
          data_faces_local[0].second.data,
          indices_compressed[i]);

      const UpdateFlags update_flags_common =
        update_flags_boundary_faces | update_flags_inner_faces;

      // Collect all data in the final data fields.
      // First allocate the memory
      const unsigned int n_constant_data =
        data_faces_local[0].second.data.size();
      for (unsigned int my_q = 0; my_q < face_data.size(); ++my_q)
        {
          face_data[my_q].data_index_offsets.resize(face_type.size());
          std::vector<std::array<std::size_t, 2>> shift(
            data_faces_local.size());
          shift[0][0] = n_constant_data;
          shift[0][1] = 0;
          for (unsigned int i = 1; i < data_faces_local.size(); ++i)
            {
              shift[i][0] =
                shift[i - 1][0] +
                data_faces_local[i - 1].first[my_q].JxW_values.size();
              shift[i][1] =
                shift[i - 1][1] +
                data_faces_local[i - 1].first[my_q].quadrature_points.size();
            }
          face_data[my_q].JxW_values.resize_fast(
            shift.back()[0] +
            data_faces_local.back().first[my_q].JxW_values.size());
          face_data[my_q].normal_vectors.resize_fast(
            face_data[my_q].JxW_values.size());
          face_data[my_q].jacobians[0].resize_fast(
            face_data[my_q].JxW_values.size());
          face_data[my_q].jacobians[1].resize_fast(
            face_data[my_q].JxW_values.size());
          if (update_flags_common & update_jacobian_grads)
            {
              face_data[my_q].jacobian_gradients[0].resize_fast(
                face_data[my_q].JxW_values.size());
              face_data[my_q].jacobian_gradients[1].resize_fast(
                face_data[my_q].JxW_values.size());
              face_data[my_q].jacobian_gradients_non_inverse[0].resize_fast(
                face_data[my_q].JxW_values.size());
              face_data[my_q].jacobian_gradients_non_inverse[1].resize_fast(
                face_data[my_q].JxW_values.size());
            }
          face_data[my_q].normals_times_jacobians[0].resize_fast(
            face_data[my_q].JxW_values.size());
          face_data[my_q].normals_times_jacobians[1].resize_fast(
            face_data[my_q].JxW_values.size());
          if (update_flags_common & update_quadrature_points)
            {
              face_data[my_q].quadrature_point_offsets.resize(face_type.size());
              face_data[my_q].quadrature_points.resize_fast(
                shift.back()[1] +
                data_faces_local.back().first[my_q].quadrature_points.size());
            }

          // start the tasks to gather the data in parallel
          Threads::TaskGroup<> tasks;
          for (unsigned int i = 0; i < data_faces_local.size(); ++i)
            tasks += Threads::new_task(
              &ExtractCellHelper::
                copy_data<dim - 1, dim, Number, VectorizedArrayType>,
              work_per_chunk * i,
              shift[i],
              indices_compressed[i],
              face_type,
              data_faces_local[i].first[my_q],
              face_data[my_q]);

          // fill the constant data fields (in parallel to the loop above)
          if (my_q == 0)
            {
              const Number jac_size =
                ExtractCellHelper::get_jacobian_size(tria);
              for (const auto &it : data_faces_local[0].second.data)
                {
                  // JxW values; invert previously applied scaling
                  for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
                    face_data[my_q].JxW_values[it.second][v] =
                      it.first[2 * dim * dim + dim][v] *
                      Utilities::fixed_power<dim>(jac_size);

                  // inverse Jacobians
                  for (unsigned int i = 0; i < 2; ++i)
                    for (unsigned int d = 0; d < dim; ++d)
                      for (unsigned int e = 0; e < dim; ++e)
                        for (unsigned int v = 0;
                             v < VectorizedArrayType::size();
                             ++v)
                          face_data[my_q].jacobians[i][it.second][d][e][v] =
                            it.first[i * dim * dim + d * dim + e][v];

                  // normal vectors; invert previously applied scaling
                  for (unsigned int d = 0; d < dim; ++d)
                    for (unsigned int v = 0; v < VectorizedArrayType::size();
                         ++v)
                      face_data[my_q].normal_vectors[it.second][d][v] =
                        it.first[2 * dim * dim + d][v] * jac_size;
                }
            }
          else
            {
              for (unsigned int i = 0; i < n_constant_data; ++i)
                {
                  face_data[my_q].JxW_values[i] = face_data[0].JxW_values[i];
                  face_data[my_q].normal_vectors[i] =
                    face_data[0].normal_vectors[i];
                  for (unsigned k = 0; k < 2; ++k)
                    face_data[my_q].jacobians[k][i] =
                      face_data[0].jacobians[k][i];
                }
            }

          // ... wait for the parallel work to finish
          tasks.join_all();

          // finally compute the normal times the jacobian
          for (unsigned int i = 0; i < data_faces_local.size(); ++i)
            tasks += Threads::new_task(
              &ExtractFaceHelper::
                compute_normal_times_jacobian<dim, Number, VectorizedArrayType>,
              work_per_chunk * i,
              std::min<unsigned int>(work_per_chunk * (i + 1), faces.size()),
              face_type,
              faces,
              face_data[my_q]);
          tasks.join_all();
        }
    }



    template <int dim, typename Number, typename VectorizedArrayType>
    void
    MappingInfo<dim, Number, VectorizedArrayType>::compute_mapping_q(
      const dealii::Triangulation<dim>                         &tria,
      const std::vector<std::pair<unsigned int, unsigned int>> &cell_array,
      const FaceInfo<VectorizedArrayType::size()>              &face_info)
    {
      // step 1: extract quadrature point data with the data appropriate for
      // MappingQ
      AssertDimension(this->mapping_collection->size(), 1);

      const MappingQ<dim> *mapping_q = dynamic_cast<const MappingQ<dim> *>(
        &this->mapping_collection->operator[](0));
      Assert(mapping_q != nullptr, ExcInternalError());

      const unsigned int mapping_degree = mapping_q->get_degree();
      const unsigned int n_mapping_points =
        Utilities::pow(mapping_degree + 1, dim);
      AlignedVector<double> plain_quadrature_points(cell_array.size() *
                                                    n_mapping_points * dim);

      const double jacobian_size = ExtractCellHelper::get_jacobian_size(tria);

      std::vector<unsigned int> cell_data_index;
      std::vector<GeometryType> preliminary_cell_type(cell_array.size());
      {
        AlignedVector<std::array<Tensor<2, dim>, dim + 1>> jacobians_on_stencil(
          cell_array.size());

        // Create as many chunks of cells as we have threads and spawn the
        // work
        unsigned int work_per_chunk =
          std::max(std::size_t(1),
                   (cell_array.size() + MultithreadInfo::n_threads() - 1) /
                     MultithreadInfo::n_threads());

        // we manually use tasks here rather than parallel::apply_to_subranges
        // because we want exactly as many loops as we have threads - the
        // initialization of the loops with FEValues is expensive
        std::size_t          offset = 0;
        Threads::TaskGroup<> tasks;
        for (unsigned int t = 0; t < MultithreadInfo::n_threads();
             ++t, offset += work_per_chunk)
          tasks += Threads::new_task(
            &ExtractCellHelper::mapping_q_query_fe_values<dim>,
            offset,
            std::min(cell_array.size(), offset + work_per_chunk),
            *mapping_q,
            tria,
            cell_array,
            jacobian_size,
            preliminary_cell_type,
            plain_quadrature_points,
            jacobians_on_stencil);
        tasks.join_all();
        cell_data_index =
          ExtractCellHelper::mapping_q_find_compression(jacobian_size,
                                                        jacobians_on_stencil,
                                                        n_mapping_points,
                                                        plain_quadrature_points,
                                                        preliminary_cell_type);
      }

      // step 2: compute the appropriate evaluation matrices for cells and
      // faces

      // We want to use vectorization for computing the quantities, but must
      // evaluate the geometry in double precision; thus, for floats we need
      // to do things in two sweeps and convert the final result.
      constexpr unsigned int n_lanes = VectorizedArrayType::size();

      using VectorizedDouble =
        VectorizedArray<double,
                        ((std::is_same_v<Number, float> &&
                          VectorizedArrayType::size() > 1) ?
                           VectorizedArrayType::size() / 2 :
                           VectorizedArrayType::size())>;

      // Create a ShapeInfo object to provide the necessary interpolators to
      // the various quadrature points. Note that it is initialized with the
      // finite element fe_geometry using the degree of the mapping, which is
      // not the same as the degree of the underlying finite element shape
      // functions or the quadrature points; shape info is merely a vehicle to
      // return us the right interpolation matrices from the cell support
      // points to the cell and face quadrature points.
      std::vector<ShapeInfo<double>> shape_infos(cell_data.size());
      {
        FE_DGQ<dim> fe_geometry(mapping_degree);
        for (unsigned int my_q = 0; my_q < cell_data.size(); ++my_q)
          shape_infos[my_q].reinit(cell_data[my_q].descriptor[0].quadrature,
                                   fe_geometry);
      }

      // step 3: find compression of cells with vectorization
      std::map<std::array<unsigned int, n_lanes>, unsigned int> compressed_data;

      cell_type.resize(cell_array.size() / n_lanes);
      std::vector<bool>         process_cell(cell_type.size());
      std::vector<unsigned int> cell_data_index_vect(cell_type.size());

      for (unsigned int cell = 0; cell < cell_array.size(); cell += n_lanes)
        {
          std::pair<std::array<unsigned int, n_lanes>, unsigned int>
            data_indices;
          for (unsigned int i = 0; i < n_lanes; ++i)
            data_indices.first[i] = cell_data_index[cell + i];
          data_indices.second = cell / n_lanes;

          auto inserted = compressed_data.insert(data_indices);

          process_cell[cell / n_lanes] = inserted.second;
          if (inserted.second == true)
            cell_data_index_vect[cell / n_lanes] = data_indices.second;
          else
            cell_data_index_vect[cell / n_lanes] = inserted.first->second;

          cell_type[cell / n_lanes] =
            *std::max_element(preliminary_cell_type.data() + cell,
                              preliminary_cell_type.data() + cell + n_lanes);
        }

      // step 4: compute the data on cells from the cached quadrature
      // points, filling up all SIMD lanes as appropriate
      for (unsigned int my_q = 0; my_q < cell_data.size(); ++my_q)
        {
          MappingInfoStorage<dim, dim, VectorizedArrayType> &my_data =
            cell_data[my_q];

          // step 4a: set the index offsets, find out how much to allocate,
          // and allocate the memory
          const unsigned int n_q_points = my_data.descriptor[0].n_q_points;
          unsigned int       max_size   = 0;
          my_data.data_index_offsets.resize(cell_type.size());
          for (unsigned int cell = 0; cell < cell_type.size(); ++cell)
            {
              if (process_cell[cell] == false)
                my_data.data_index_offsets[cell] =
                  my_data.data_index_offsets[cell_data_index_vect[cell]];
              else
                my_data.data_index_offsets[cell] = max_size;
              max_size =
                std::max(max_size,
                         my_data.data_index_offsets[cell] +
                           (cell_type[cell] <= affine ? 2 : n_q_points));
            }

          my_data.JxW_values.resize_fast(max_size);
          my_data.jacobians[0].resize_fast(max_size);
          if (update_flags_cells & update_jacobian_grads)
            {
              my_data.jacobian_gradients[0].resize_fast(max_size);
              my_data.jacobian_gradients_non_inverse[0].resize_fast(max_size);
            }

          if (update_flags_cells & update_quadrature_points)
            {
              my_data.quadrature_point_offsets.resize(cell_type.size());
              for (unsigned int cell = 1; cell < cell_type.size(); ++cell)
                if (cell_type[cell - 1] <= affine)
                  my_data.quadrature_point_offsets[cell] =
                    my_data.quadrature_point_offsets[cell - 1] + 1;
                else
                  my_data.quadrature_point_offsets[cell] =
                    my_data.quadrature_point_offsets[cell - 1] + n_q_points;
              my_data.quadrature_points.resize_fast(
                my_data.quadrature_point_offsets.back() +
                (cell_type.back() <= affine ? 1 : n_q_points));
            }

          // step 4b: go through the cells and compute the information using
          // similar evaluators as for the matrix-free integrals
          dealii::parallel::apply_to_subranges(
            0U,
            cell_type.size(),
            [&](const unsigned int begin, const unsigned int end) {
              ExtractCellHelper::mapping_q_compute_range<dim,
                                                         Number,
                                                         VectorizedArrayType,
                                                         VectorizedDouble>(
                begin,
                end,
                tria,
                cell_array,
                cell_type,
                process_cell,
                update_flags_cells,
                plain_quadrature_points,
                shape_infos[my_q],
                my_data);
            },
            std::max(cell_type.size() / MultithreadInfo::n_threads() / 2,
                     std::size_t(2U)));
        }

      const std::vector<FaceToCellTopology<VectorizedArrayType::size()>>
        &faces = face_info.faces;
      if (faces.empty())
        return;

      // step 5: find compression of faces with vectorization
      std::map<std::array<unsigned int, 2 * n_lanes + 3>, unsigned int>
        compressed_faces;

      face_type.resize(faces.size());
      std::vector<bool>         process_face(face_type.size());
      std::vector<unsigned int> face_data_index_vect(face_type.size());

      for (unsigned int face = 0; face < faces.size(); ++face)
        {
          std::pair<std::array<unsigned int, 2 * n_lanes + 3>, unsigned int>
            data_indices;
          for (unsigned int i = 0; i < n_lanes; ++i)
            if (faces[face].cells_interior[i] != numbers::invalid_unsigned_int)
              data_indices.first[i] =
                cell_data_index[faces[face].cells_interior[i]];
            else
              data_indices.first[i] = data_indices.first[0];
          for (unsigned int i = 0; i < n_lanes; ++i)
            data_indices.first[n_lanes + i] = data_indices.first[i];
          for (unsigned int i = 0; i < n_lanes; ++i)
            if (faces[face].cells_exterior[i] != numbers::invalid_unsigned_int)
              data_indices.first[n_lanes + i] =
                cell_data_index[faces[face].cells_exterior[i]];
          data_indices.first[2 * n_lanes]     = faces[face].interior_face_no;
          data_indices.first[2 * n_lanes + 1] = faces[face].exterior_face_no;
          data_indices.first[2 * n_lanes + 2] = faces[face].subface_index;

          data_indices.second = face;

          auto inserted = compressed_faces.insert(data_indices);

          process_face[face] = inserted.second;
          if (inserted.second == true)
            face_data_index_vect[face] = face;
          else
            face_data_index_vect[face] = inserted.first->second;

          face_type[face] = cartesian;
          for (unsigned int i = 0; i < n_lanes; ++i)
            if (faces[face].cells_interior[i] != numbers::invalid_unsigned_int)
              face_type[face] =
                std::max(face_type[face],
                         preliminary_cell_type[faces[face].cells_interior[i]]);
          for (unsigned int i = 0; i < n_lanes; ++i)
            if (faces[face].cells_exterior[i] != numbers::invalid_unsigned_int)
              face_type[face] =
                std::max(face_type[face],
                         preliminary_cell_type[faces[face].cells_exterior[i]]);
        }

      // step 6: compute the data on faces from the cached cell quadrature
      // points, filling up all SIMD lanes as appropriate
      for (unsigned int my_q = 0; my_q < face_data.size(); ++my_q)
        {
          MappingInfoStorage<dim - 1, dim, VectorizedArrayType> &my_data =
            face_data[my_q];

          // step 6a: set the index offsets, find out how much to allocate,
          // and allocate the memory
          const unsigned int n_q_points = my_data.descriptor[0].n_q_points;
          unsigned int       max_size   = 0;
          my_data.data_index_offsets.resize(face_type.size());
          for (unsigned int face = 0; face < face_type.size(); ++face)
            {
              if (process_face[face] == false)
                my_data.data_index_offsets[face] =
                  my_data.data_index_offsets[face_data_index_vect[face]];
              else
                my_data.data_index_offsets[face] = max_size;
              max_size =
                std::max(max_size,
                         my_data.data_index_offsets[face] +
                           (face_type[face] <= affine ? 2 : n_q_points));
            }

          const UpdateFlags update_flags_common =
            update_flags_boundary_faces | update_flags_inner_faces;

          my_data.JxW_values.resize_fast(max_size);
          my_data.normal_vectors.resize_fast(max_size);
          my_data.jacobians[0].resize_fast(max_size);
          my_data.jacobians[1].resize_fast(max_size);
          if (update_flags_common & update_jacobian_grads)
            {
              my_data.jacobian_gradients[0].resize_fast(max_size);
              my_data.jacobian_gradients[1].resize_fast(max_size);
              my_data.jacobian_gradients_non_inverse[0].resize_fast(max_size);
              my_data.jacobian_gradients_non_inverse[1].resize_fast(max_size);
            }
          my_data.normals_times_jacobians[0].resize_fast(max_size);
          my_data.normals_times_jacobians[1].resize_fast(max_size);

          if (update_flags_common & update_quadrature_points)
            {
              my_data.quadrature_point_offsets.resize(face_type.size());
              my_data.quadrature_point_offsets[0] = 0;
              for (unsigned int face = 1; face < faces.size(); ++face)
                my_data.quadrature_point_offsets[face] =
                  n_q_points + my_data.quadrature_point_offsets[face - 1];
              my_data.quadrature_points.resize_fast(face_type.size() *
                                                    n_q_points);
            }

          // step 6b: go through the faces and compute the information using
          // similar evaluators as for the matrix-free face integrals
          dealii::parallel::apply_to_subranges(
            0U,
            face_type.size(),
            [&](const unsigned int begin, const unsigned int end) {
              ExtractFaceHelper::mapping_q_compute_range<dim,
                                                         Number,
                                                         VectorizedArrayType,
                                                         VectorizedDouble>(
                begin,
                end,
                faces,
                face_type,
                process_face,
                update_flags_common,
                plain_quadrature_points,
                shape_infos[my_q],
                my_data);
            },
            std::max(face_type.size() / MultithreadInfo::n_threads() / 2,
                     std::size_t(2U)));
        }

      // step 6c: figure out if normal vectors are the same on some of the
      // faces which allows us to set the flat_faces face type
      unsigned int quad_with_most_points = 0;
      for (unsigned int my_q = 1; my_q < face_data.size(); ++my_q)
        if (face_data[my_q].descriptor[0].n_q_points >
            face_data[quad_with_most_points].descriptor[0].n_q_points)
          quad_with_most_points = my_q;
      dealii::parallel::apply_to_subranges(
        0U,
        face_type.size(),
        [&](const unsigned int begin, const unsigned int end) {
          for (unsigned int face = begin; face < end; ++face)
            if (face_type[face] == general)
              {
                const unsigned int n_q_points =
                  face_data[quad_with_most_points].descriptor[0].n_q_points;
                const Tensor<1, dim, VectorizedArrayType> *normals =
                  face_data[quad_with_most_points].normal_vectors.data() +
                  face_data[quad_with_most_points].data_index_offsets[face];
                VectorizedArrayType distance = 0.;
                for (unsigned int q = 1; q < n_q_points; ++q)
                  distance += (normals[q] - normals[0]).norm_square();
                bool all_small = true;
                for (unsigned int v = 0; v < n_lanes; ++v)
                  if (distance[v] >
                      50. * std::numeric_limits<Number>::epsilon() *
                        std::numeric_limits<Number>::epsilon() * n_q_points)
                    all_small = false;
                if (all_small)
                  face_type[face] = flat_faces;
              }
        },
        std::max(face_type.size() / MultithreadInfo::n_threads() / 2,
                 std::size_t(2U)));

      // step 7: compute the face data by cells. This still needs to be
      // transitioned to extracting the information from cell quadrature
      // points but we need to figure out the correct indices of neighbors
      // within the list of arrays still
      initialize_faces_by_cells(tria,
                                cell_array,
                                face_info,
                                *this->mapping_collection);
    }



    template <int dim, typename Number, typename VectorizedArrayType>
    void
    MappingInfo<dim, Number, VectorizedArrayType>::initialize_faces_by_cells(
      const dealii::Triangulation<dim>                         &tria,
      const std::vector<std::pair<unsigned int, unsigned int>> &cells,
      const FaceInfo<VectorizedArrayType::size()>              &face_info,
      const dealii::hp::MappingCollection<dim>                 &mapping_in)
    {
      if (update_flags_faces_by_cells == update_default)
        return;


      AssertDimension(mapping_in.size(), 1);
      const auto &mapping = mapping_in[0];

      const unsigned int     n_quads = face_data_by_cells.size();
      constexpr unsigned int n_lanes = VectorizedArrayType::size();
      UpdateFlags            update_flags =
        (update_flags_faces_by_cells & update_quadrature_points ?
           update_quadrature_points :
           update_default) |
        update_normal_vectors | update_JxW_values | update_jacobians;

      const auto compute_neighbor_index =
        [&face_info](const unsigned int cell,
                     const unsigned int face,
                     const unsigned int lane) {
          constexpr unsigned int n_lanes   = VectorizedArrayType::size();
          const unsigned int     cell_this = cell * n_lanes + lane;
          unsigned int           face_index =
            face_info.cell_and_face_to_plain_faces(cell, face, lane);
          if (face_index == numbers::invalid_unsigned_int)
            return numbers::invalid_unsigned_int;

          unsigned int cell_neighbor = face_info.faces[face_index / n_lanes]
                                         .cells_interior[face_index % n_lanes];
          if (cell_neighbor == cell_this)
            cell_neighbor = face_info.faces[face_index / n_lanes]
                              .cells_exterior[face_index % n_lanes];
          return cell_neighbor;
        };

      faces_by_cells_type.resize(cell_type.size());
      for (unsigned int cell = 0; cell < cell_type.size(); ++cell)
        for (const unsigned int face : GeometryInfo<dim>::face_indices())
          {
            faces_by_cells_type[cell][face] = cell_type[cell];
            for (unsigned int v = 0; v < n_lanes; ++v)
              {
                const unsigned int cell_neighbor =
                  compute_neighbor_index(cell, face, v);
                if (cell_neighbor != numbers::invalid_unsigned_int)
                  faces_by_cells_type[cell][face] =
                    std::max(faces_by_cells_type[cell][face],
                             cell_type[cell_neighbor / n_lanes]);
              }
          }

      for (unsigned int my_q = 0; my_q < n_quads; ++my_q)
        {
          // since we already know the cell type, we can pre-allocate the right
          // amount of data straight away and we just need to do some basic
          // counting
          AssertDimension(cell_type.size(), cells.size() / n_lanes);
          face_data_by_cells[my_q].data_index_offsets.resize(
            cell_type.size() * ReferenceCells::max_n_faces<dim>());
          if (update_flags & update_quadrature_points)
            face_data_by_cells[my_q].quadrature_point_offsets.resize(
              cell_type.size() * ReferenceCells::max_n_faces<dim>());
          std::size_t storage_length = 0;
          for (unsigned int i = 0; i < cell_type.size(); ++i)
            for (const unsigned int face : GeometryInfo<dim>::face_indices())
              {
                if (faces_by_cells_type[i][face] <= affine)
                  {
                    face_data_by_cells[my_q].data_index_offsets
                      [i * ReferenceCells::max_n_faces<dim>() + face] =
                      storage_length;
                    ++storage_length;
                  }
                else
                  {
                    face_data_by_cells[my_q].data_index_offsets
                      [i * ReferenceCells::max_n_faces<dim>() + face] =
                      storage_length;
                    storage_length +=
                      face_data_by_cells[my_q].descriptor[0].n_q_points;
                  }
                if (update_flags & update_quadrature_points)
                  face_data_by_cells[my_q].quadrature_point_offsets
                    [i * ReferenceCells::max_n_faces<dim>() + face] =
                    (i * ReferenceCells::max_n_faces<dim>() + face) *
                    face_data_by_cells[my_q].descriptor[0].n_q_points;
              }
          face_data_by_cells[my_q].JxW_values.resize_fast(
            storage_length * ReferenceCells::max_n_faces<dim>());
          face_data_by_cells[my_q].jacobians[0].resize_fast(
            storage_length * ReferenceCells::max_n_faces<dim>());
          face_data_by_cells[my_q].jacobians[1].resize_fast(
            storage_length * ReferenceCells::max_n_faces<dim>());
          if (update_flags & update_normal_vectors)
            face_data_by_cells[my_q].normal_vectors.resize_fast(
              storage_length * ReferenceCells::max_n_faces<dim>());
          if (update_flags & update_normal_vectors &&
              update_flags & update_jacobians)
            face_data_by_cells[my_q].normals_times_jacobians[0].resize_fast(
              storage_length * ReferenceCells::max_n_faces<dim>());
          if (update_flags & update_normal_vectors &&
              update_flags & update_jacobians)
            face_data_by_cells[my_q].normals_times_jacobians[1].resize_fast(
              storage_length * ReferenceCells::max_n_faces<dim>());
          if (update_flags & update_jacobian_grads)
            {
              face_data_by_cells[my_q].jacobian_gradients[0].resize_fast(
                storage_length * ReferenceCells::max_n_faces<dim>());
              face_data_by_cells[my_q]
                .jacobian_gradients_non_inverse[0]
                .resize_fast(storage_length *
                             ReferenceCells::max_n_faces<dim>());
            }

          if (update_flags & update_quadrature_points)
            face_data_by_cells[my_q].quadrature_points.resize_fast(
              cell_type.size() * ReferenceCells::max_n_faces<dim>() *
              face_data_by_cells[my_q].descriptor[0].n_q_points);
        }

      FE_Nothing<dim> dummy_fe;
      // currently no hp-indices implemented
      const unsigned int fe_index = 0;
      std::vector<std::vector<std::shared_ptr<dealii::FEFaceValues<dim>>>>
        fe_face_values(face_data_by_cells.size());
      for (unsigned int i = 0; i < fe_face_values.size(); ++i)
        fe_face_values[i].resize(face_data_by_cells[i].descriptor.size());
      std::vector<std::vector<std::shared_ptr<dealii::FEFaceValues<dim>>>>
        fe_face_values_neigh(face_data_by_cells.size());
      for (unsigned int i = 0; i < fe_face_values_neigh.size(); ++i)
        fe_face_values_neigh[i].resize(face_data_by_cells[i].descriptor.size());
      for (unsigned int cell = 0; cell < cell_type.size(); ++cell)
        for (unsigned int my_q = 0; my_q < face_data_by_cells.size(); ++my_q)
          for (const unsigned int face : GeometryInfo<dim>::face_indices())
            {
              if (fe_face_values[my_q][fe_index].get() == nullptr)
                fe_face_values[my_q][fe_index] =
                  std::make_shared<dealii::FEFaceValues<dim>>(
                    mapping,
                    dummy_fe,
                    face_data[my_q].q_collection[fe_index],
                    update_flags);
              if (fe_face_values_neigh[my_q][fe_index].get() == nullptr)
                fe_face_values_neigh[my_q][fe_index] =
                  std::make_shared<dealii::FEFaceValues<dim>>(
                    mapping,
                    dummy_fe,
                    face_data[my_q].q_collection[fe_index],
                    update_flags);
              dealii::FEFaceValues<dim> &fe_val =
                *fe_face_values[my_q][fe_index];
              dealii::FEFaceValues<dim> &fe_val_neigh =
                *fe_face_values_neigh[my_q][fe_index];
              const unsigned int offset =
                face_data_by_cells[my_q].data_index_offsets
                  [cell * ReferenceCells::max_n_faces<dim>() + face];

              const GeometryType my_cell_type = faces_by_cells_type[cell][face];

              for (unsigned int v = 0; v < n_lanes; ++v)
                {
                  typename dealii::Triangulation<dim>::cell_iterator cell_it(
                    &tria,
                    cells[cell * n_lanes + v].first,
                    cells[cell * n_lanes + v].second);
                  fe_val.reinit(cell_it, face);

                  const unsigned int cell_neighbor =
                    compute_neighbor_index(cell, face, v);

                  if (cell_neighbor != numbers::invalid_unsigned_int)
                    {
                      typename dealii::Triangulation<dim>::cell_iterator
                        cell_it_neigh(&tria,
                                      cells[cell_neighbor].first,
                                      cells[cell_neighbor].second);
                      fe_val_neigh.reinit(cell_it_neigh,
                                          cell_it->at_boundary(face) ?
                                            cell_it->periodic_neighbor_face_no(
                                              face) :
                                            cell_it->neighbor_face_no(face));
                    }

                  // copy data for affine data type
                  if (my_cell_type <= affine)
                    {
                      if (update_flags & update_JxW_values)
                        face_data_by_cells[my_q].JxW_values[offset][v] =
                          fe_val.JxW(0) / fe_val.get_quadrature().weight(0);
                      if (update_flags & update_jacobians)
                        {
                          DerivativeForm<1, dim, dim> inv_jac =
                            fe_val.jacobian(0).covariant_form();
                          for (unsigned int d = 0; d < dim; ++d)
                            for (unsigned int e = 0; e < dim; ++e)
                              {
                                const unsigned int ee = ExtractFaceHelper::
                                  reorder_face_derivative_indices<dim>(face, e);
                                face_data_by_cells[my_q]
                                  .jacobians[0][offset][d][e][v] =
                                  inv_jac[d][ee];
                              }
                        }
                      if (cell_neighbor != numbers::invalid_unsigned_int &&
                          (update_flags & update_jacobians))
                        for (unsigned int q = 0; q < fe_val.n_quadrature_points;
                             ++q)
                          {
                            DerivativeForm<1, dim, dim> inv_jac =
                              fe_val_neigh.jacobian(q).covariant_form();
                            for (unsigned int d = 0; d < dim; ++d)
                              for (unsigned int e = 0; e < dim; ++e)
                                {
                                  const unsigned int ee = ExtractFaceHelper::
                                    reorder_face_derivative_indices<dim>(
                                      fe_val_neigh.get_face_number(), e);
                                  face_data_by_cells[my_q]
                                    .jacobians[1][offset][d][e][v] =
                                    inv_jac[d][ee];
                                }
                          }
                      if (update_flags & update_jacobian_grads)
                        {
                          DEAL_II_NOT_IMPLEMENTED();
                        }
                      if (update_flags & update_normal_vectors)
                        for (unsigned int d = 0; d < dim; ++d)
                          face_data_by_cells[my_q]
                            .normal_vectors[offset][d][v] =
                            fe_val.normal_vector(0)[d];
                    }
                  // copy data for general data type
                  else
                    {
                      if (update_flags & update_JxW_values)
                        for (unsigned int q = 0; q < fe_val.n_quadrature_points;
                             ++q)
                          face_data_by_cells[my_q].JxW_values[offset + q][v] =
                            fe_val.JxW(q);
                      if (update_flags & update_jacobians)
                        for (unsigned int q = 0; q < fe_val.n_quadrature_points;
                             ++q)
                          {
                            DerivativeForm<1, dim, dim> inv_jac =
                              fe_val.jacobian(q).covariant_form();
                            for (unsigned int d = 0; d < dim; ++d)
                              for (unsigned int e = 0; e < dim; ++e)
                                {
                                  const unsigned int ee = ExtractFaceHelper::
                                    reorder_face_derivative_indices<dim>(face,
                                                                         e);
                                  face_data_by_cells[my_q]
                                    .jacobians[0][offset + q][d][e][v] =
                                    inv_jac[d][ee];
                                }
                          }
                      if (cell_neighbor != numbers::invalid_unsigned_int &&
                          (update_flags & update_jacobians))
                        for (unsigned int q = 0; q < fe_val.n_quadrature_points;
                             ++q)
                          {
                            DerivativeForm<1, dim, dim> inv_jac =
                              fe_val_neigh.jacobian(q).covariant_form();
                            for (unsigned int d = 0; d < dim; ++d)
                              for (unsigned int e = 0; e < dim; ++e)
                                {
                                  const unsigned int ee = ExtractFaceHelper::
                                    reorder_face_derivative_indices<dim>(
                                      fe_val_neigh.get_face_number(), e);
                                  face_data_by_cells[my_q]
                                    .jacobians[1][offset + q][d][e][v] =
                                    inv_jac[d][ee];
                                }
                          }
                      if (update_flags & update_jacobian_grads)
                        {
                          DEAL_II_NOT_IMPLEMENTED();
                        }
                      if (update_flags & update_normal_vectors)
                        for (unsigned int q = 0; q < fe_val.n_quadrature_points;
                             ++q)
                          for (unsigned int d = 0; d < dim; ++d)
                            face_data_by_cells[my_q]
                              .normal_vectors[offset + q][d][v] =
                              fe_val.normal_vector(q)[d];
                    }
                  if (update_flags & update_quadrature_points)
                    for (unsigned int q = 0; q < fe_val.n_quadrature_points;
                         ++q)
                      for (unsigned int d = 0; d < dim; ++d)
                        face_data_by_cells[my_q].quadrature_points
                          [face_data_by_cells[my_q].quadrature_point_offsets
                             [cell * ReferenceCells::max_n_faces<dim>() +
                              face] +
                           q][d][v] = fe_val.quadrature_point(q)[d];
                }
              if (update_flags & update_normal_vectors &&
                  update_flags & update_jacobians)
                for (unsigned int q = 0;
                     q <
                     (my_cell_type <= affine ? 1 : fe_val.n_quadrature_points);
                     ++q)
                  face_data_by_cells[my_q]
                    .normals_times_jacobians[0][offset + q] =
                    face_data_by_cells[my_q].normal_vectors[offset + q] *
                    face_data_by_cells[my_q].jacobians[0][offset + q];
              if (update_flags & update_normal_vectors &&
                  update_flags & update_jacobians)
                for (unsigned int q = 0;
                     q <
                     (my_cell_type <= affine ? 1 : fe_val.n_quadrature_points);
                     ++q)
                  face_data_by_cells[my_q]
                    .normals_times_jacobians[1][offset + q] =
                    face_data_by_cells[my_q].normal_vectors[offset + q] *
                    face_data_by_cells[my_q].jacobians[1][offset + q];
            }
    }



    template <int dim, typename Number, typename VectorizedArrayType>
    std::size_t
    MappingInfo<dim, Number, VectorizedArrayType>::memory_consumption() const
    {
      std::size_t memory = MemoryConsumption::memory_consumption(cell_data);
      memory += MemoryConsumption::memory_consumption(face_data);
      memory += MemoryConsumption::memory_consumption(face_data_by_cells);
      memory += cell_type.capacity() * sizeof(GeometryType);
      memory += face_type.capacity() * sizeof(GeometryType);
      memory += faces_by_cells_type.capacity() *
                ReferenceCells::max_n_faces<dim>() * sizeof(GeometryType);
      memory += sizeof(*this);
      return memory;
    }



    template <int dim, typename Number, typename VectorizedArrayType>
    template <typename StreamType>
    void
    MappingInfo<dim, Number, VectorizedArrayType>::print_memory_consumption(
      StreamType     &out,
      const TaskInfo &task_info) const
    {
      out << "    Cell types:                      ";
      task_info.print_memory_statistics(out,
                                        cell_type.capacity() *
                                          sizeof(GeometryType));
      out << "    Face types:                      ";
      task_info.print_memory_statistics(out,
                                        face_type.capacity() *
                                          sizeof(GeometryType));

      out << "    Faces by cells types:            ";
      task_info.print_memory_statistics(out,
                                        faces_by_cells_type.capacity() *
                                          ReferenceCells::max_n_faces<dim>() *
                                          sizeof(GeometryType));

      for (unsigned int j = 0; j < cell_data.size(); ++j)
        {
          out << "    Data component " << j << std::endl;
          cell_data[j].print_memory_consumption(out, task_info);
          face_data[j].print_memory_consumption(out, task_info);
          face_data_by_cells[j].print_memory_consumption(out, task_info);
        }
    }

  } // namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif

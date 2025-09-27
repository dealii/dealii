// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_vector_tools_interpolate_templates_h
#define dealii_vector_tools_interpolate_templates_h


#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/intergrid_map.h>

#include <deal.II/hp/fe_values.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include <deal.II/numerics/vector_tools_interpolate.h>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace EpetraWrappers
  {
    class Vector;
  }
} // namespace LinearAlgebra


namespace VectorTools
{
  // This namespace contains the actual implementation called
  // by VectorTools::interpolate and variants (such as
  // VectorTools::interpolate_by_material_id).
  namespace internal
  {
    // A small helper function to transform a component range starting
    // at offset from the real to the unit cell according to the
    // supplied conformity. The function_values vector is transformed
    // in place.
    //
    // FIXME: This should be refactored into the mapping (i.e.
    // implement the inverse function of Mapping<dim, spacedim>::transform).
    // Further, the finite element should make the information about
    // the correct mapping directly accessible (i.e. which MappingKind
    // should be used). Using fe.conforming_space might be a bit of a
    // problem because we only support doing nothing, Hcurl, and Hdiv
    // conforming mappings.
    //
    // Input:
    //  conformity: conformity of the finite element, used to select
    //              appropriate type of transformation
    //  fe_values_jacobians: used for jacobians (and inverses of
    //                        jacobians). the object is supposed to be
    //                        reinit()'d for the current cell
    //  function_values, offset: function_values is manipulated in place
    //                           starting at position offset
    template <int dim, int spacedim, typename FEValuesType, typename T3>
    void
    transform(const typename FiniteElementData<dim>::Conformity conformity,
              const unsigned int                                offset,
              const FEValuesType &fe_values_jacobians,
              T3                 &function_values)
    {
      switch (conformity)
        {
          case FiniteElementData<dim>::Hcurl:
            // See Monk, Finite Element Methods for Maxwell's Equations,
            // p. 77ff, formula (3.76) and Corollary 3.58.
            // For given mapping F_K: \hat K \to K, we have to transform
            //  \hat u = (dF_K)^T u\circ F_K

            for (unsigned int i = 0; i < function_values.size(); ++i)
              {
                const auto &jacobians =
                  fe_values_jacobians.get_present_fe_values().get_jacobians();

                const ArrayView<typename T3::value_type::value_type> source(
                  &function_values[i][0] + offset, dim);

                Tensor<1,
                       dim,
                       typename ProductType<typename T3::value_type::value_type,
                                            double>::type>
                  destination;

                // value[m] <- sum jacobian_transpose[m][n] * old_value[n]:
                TensorAccessors::contract<1, 2, 1, dim>(
                  destination, jacobians[i].transpose(), source);

                // now copy things back into the input=output vector
                for (unsigned int d = 0; d < dim; ++d)
                  source[d] = destination[d];
              }
            break;

          case FiniteElementData<dim>::Hdiv:
            // See Monk, Finite Element Methods for Maxwell's Equations,
            // p. 79ff, formula (3.77) and Lemma 3.59.
            // For given mapping F_K: \hat K \to K, we have to transform
            //  \hat w = det(dF_K) (dF_K)^{-1} w\circ F_K

            for (unsigned int i = 0; i < function_values.size(); ++i)
              {
                const auto &jacobians =
                  fe_values_jacobians.get_present_fe_values().get_jacobians();
                const auto &inverse_jacobians =
                  fe_values_jacobians.get_present_fe_values()
                    .get_inverse_jacobians();

                const ArrayView<typename T3::value_type::value_type> source(
                  &function_values[i][0] + offset, dim);

                Tensor<1,
                       dim,
                       typename ProductType<typename T3::value_type::value_type,
                                            double>::type>
                  destination;

                // value[m] <- sum inverse_jacobians[m][n] * old_value[n]:
                TensorAccessors::contract<1, 2, 1, dim>(destination,
                                                        inverse_jacobians[i],
                                                        source);
                destination *= jacobians[i].determinant();

                // now copy things back into the input=output vector
                for (unsigned int d = 0; d < dim; ++d)
                  source[d] = destination[d];
              }
            break;

          case FiniteElementData<dim>::H1:
            DEAL_II_FALLTHROUGH;
          case FiniteElementData<dim>::L2:
            // See Monk, Finite Element Methods for Maxwell's Equations,
            // p. 77ff, formula (3.74).
            // For given mapping F_K: \hat K \to K, we have to transform
            //  \hat p = p\circ F_K
            //  i.e., do nothing.
            break;

          default:
            // In case we deal with an unknown conformity, just assume we
            // deal with a Lagrange element and do nothing.
            break;

        } /*switch*/
    }


    // A small helper function that iteratively applies above transform
    // function to a vector function_values recursing over a given finite
    // element decomposing it into base elements:
    //
    // Input
    //   fe: the full finite element corresponding to function_values
    //   [ rest see above]
    // Output: the offset after we have handled the element at
    //   a given offset
    template <int dim, int spacedim, typename FEValuesType, typename T3>
    unsigned int
    apply_transform(const FiniteElement<dim, spacedim> &fe,
                    const unsigned int                  offset,
                    const FEValuesType                 &fe_values_jacobians,
                    T3                                 &function_values)
    {
      if (fe.n_base_elements() > 1 || fe.element_multiplicity(0) > 1)
        {
          // In case of an FESystem transform every (vector) component
          // separately:
          unsigned current_offset = offset;
          for (unsigned int i = 0; i < fe.n_base_elements(); ++i)
            {
              const auto &base_fe      = fe.base_element(i);
              const auto  multiplicity = fe.element_multiplicity(i);
              for (unsigned int m = 0; m < multiplicity; ++m)
                {
                  // recursively call apply_transform to make sure to
                  // correctly handle nested FE systems.
                  current_offset = apply_transform(base_fe,
                                                   current_offset,
                                                   fe_values_jacobians,
                                                   function_values);
                }
            }
          return current_offset;
        }
      else
        {
          transform<dim, spacedim>(fe.conforming_space,
                                   offset,
                                   fe_values_jacobians,
                                   function_values);
          return (offset + fe.n_components());
        }
    }


    // Internal implementation of interpolate that takes a generic functor
    // function such that function(cell) is of type
    // Function<spacedim, typename VectorType::value_type>*
    //
    // A given cell is skipped if function(cell) == nullptr
    template <int dim, int spacedim, typename VectorType, typename T>
    DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
    void interpolate(
      const hp::MappingCollection<dim, spacedim> &mapping_collection,
      const DoFHandler<dim, spacedim>            &dof_handler,
      T                                          &function,
      VectorType                                 &vec,
      const ComponentMask                        &component_mask,
      const unsigned int level = numbers::invalid_unsigned_int)
    {
      Assert(component_mask.represents_n_components(
               dof_handler.get_fe_collection().n_components()),
             ExcMessage(
               "The number of components in the mask has to be either "
               "zero or equal to the number of components in the finite "
               "element."));

      if (level == numbers::invalid_unsigned_int)
        AssertDimension(vec.size(), dof_handler.n_dofs());
      else
        AssertDimension(vec.size(), dof_handler.n_dofs(level));

      Assert(component_mask.n_selected_components(
               dof_handler.get_fe_collection().n_components()) > 0,
             ComponentMask::ExcNoComponentSelected());

      //
      // Computing the generalized interpolant isn't quite as straightforward
      // as for classical Lagrange elements. A major complication is the fact
      // it generally doesn't hold true that a function evaluates to the same
      // dof coefficient on different cells. This means *setting* the value
      // of a (global) degree of freedom computed on one cell doesn't
      // necessarily lead to the same result when computed on a neighboring
      // cell (that shares the same global degree of freedom).
      //
      // We thus, do the following operation:
      //
      // On each cell:
      //
      //  - We first determine all function values u(x_i) in generalized
      //    support points
      //
      //  - We transform these function values back to the unit cell
      //    according to the conformity of the component (scalar, Hdiv, or
      //    Hcurl conforming); see [Monk, Finite Element Methods for Maxwell's
      //    Equations, p.77ff Section 3.9] for details. This results in
      //    \hat u(\hat x_i)
      //
      //  - We convert these generalized support point values to nodal values
      //
      //  - For every global dof we take the average 1 / n_K \sum_{K} dof_K
      //    where n_K is the number of cells sharing the global dof and dof_K
      //    is the computed value on the cell K.
      //
      // For every degree of freedom that is shared by k cells, we compute
      // its value on all k cells and take the weighted average with respect
      // to the JxW values.
      //

      using number = typename VectorType::value_type;

      const hp::FECollection<dim, spacedim> &fe(
        dof_handler.get_fe_collection());

      std::vector<types::global_dof_index> dofs_on_cell(fe.max_dofs_per_cell());

      // Temporary storage for cell-wise interpolation operation. We store a
      // variant for every FE we encounter to speed up resizing operations.
      // The first vector is used for local function evaluation. The vector
      // dof_values is used to store intermediate cell-wise interpolation
      // results (see the detailed explanation in the for loop further down
      // below).

      std::vector<std::vector<Vector<number>>> fe_function_values(fe.size());
      std::vector<std::vector<number>>         fe_dof_values(fe.size());

      // We will need two temporary global vectors that store the new values
      // and weights.
      VectorType interpolation;
      VectorType weights;
      interpolation.reinit(vec);
      weights.reinit(vec);

      // Store locally owned dofs, so that we can skip all non-local dofs,
      // if they do not need to be interpolated.
      const IndexSet locally_owned_dofs = vec.locally_owned_elements();

      // We use an FEValues object to transform all generalized support
      // points from the unit cell to the real cell coordinates. Thus,
      // initialize a quadrature with all generalized support points and
      // create an FEValues object with it.

      std::vector<bool>    needs_expensive_algorithm(fe.size(), true);
      hp::QCollection<dim> support_quadrature;
      for (unsigned int fe_index = 0; fe_index < fe.size(); ++fe_index)
        {
          const auto &fe_i = fe[fe_index];
          // If the finite element has no dofs, we can skip it
          if (fe_i.dofs_per_cell == 0)
            {
              support_quadrature.push_back(Quadrature<dim>());
              continue;
            }
          Assert(fe_i.has_generalized_support_points(),
                 ExcMessage(
                   "The finite element does not have generalized support "
                   "points. This is required for interpolation."));
          const auto &points = fe_i.get_generalized_support_points();
          support_quadrature.push_back(Quadrature<dim>(points));
          if (fe_i.n_base_elements() == 1 &&
              fe_i.element_multiplicity(0) == fe.n_components() &&
              fe_i.has_support_points())
            {
              const auto &fe_base          = fe_i.base_element(0);
              bool        all_points_equal = true;
              // Check points for exact equality - they are either copied
              // inside an FESystem or genuinely different, so no need for a
              // tolerance
              for (unsigned int i = 0; i < fe_base.n_dofs_per_cell(); ++i)
                if (fe_base.get_unit_support_points()[i].distance(points[i]) >
                    0.)
                  {
                    all_points_equal = false;
                    break;
                  }
              if (all_points_equal)
                needs_expensive_algorithm[fe_index] = false;
            }
        }

      // An FEValues object to evaluate (generalized) support point
      // locations as well as Jacobians and their inverses.
      // The latter are only needed for Hcurl or Hdiv conforming elements,
      // but we'll just always include them.
      hp::FEValues<dim, spacedim> fe_values(mapping_collection,
                                            fe,
                                            support_quadrature,
                                            update_quadrature_points |
                                              update_jacobians |
                                              update_inverse_jacobians);

      //
      // Now loop over all locally owned, active cells.
      //
      const auto runner = [&](const auto &cell) {
        const unsigned int fe_index = cell->active_fe_index();

        // Do nothing if there are no local degrees of freedom.
        if (fe[fe_index].n_dofs_per_cell() == 0)
          return;

        // Skip processing of the current cell if the function object is
        // invalid. This is used by interpolate_by_material_id to skip
        // interpolating over cells with unknown material id.
        if (!function(cell))
          return;

        // Get transformed, generalized support points
        fe_values.reinit(cell);
        const std::vector<Point<spacedim>> &generalized_support_points =
          fe_values.get_present_fe_values().get_quadrature_points();

        // Get indices of the dofs on this cell
        const auto n_dofs = fe[fe_index].n_dofs_per_cell();
        dofs_on_cell.resize(n_dofs);
        cell->get_active_or_mg_dof_indices(dofs_on_cell);

        // Prepare temporary storage
        auto &function_values = fe_function_values[fe_index];
        auto &dof_values      = fe_dof_values[fe_index];

        const auto n_components = fe[fe_index].n_components();
        // Only resize (and create sample entry) if sizes do not match
        if (function_values.size() != generalized_support_points.size())
          function_values.resize(generalized_support_points.size(),
                                 Vector<number>(n_components));

        // Get all function values:
        AssertDimension(n_components, function(cell)->n_components);
        function(cell)->vector_value_list(generalized_support_points,
                                          function_values);

        // For the simple case with elements with support points, we will
        // simply use the interpolated DoF values in the access loop further
        // down. Otherwise, we have to transform all function values from
        // the real cell back to the unit cell. We query the finite element
        // for the correct transformation. Matters get a bit more
        // complicated because we have to apply said transformation for
        // every base element.
        if (needs_expensive_algorithm[fe_index])
          {
            dof_values.resize(n_dofs);
            const unsigned int offset =
              apply_transform(fe[fe_index],
                              /* starting_offset = */ 0,
                              fe_values,
                              function_values);
            Assert(offset == n_components, ExcInternalError());

            FETools::convert_generalized_support_point_values_to_dof_values(
              fe[fe_index], function_values, dof_values);
          }

        for (unsigned int i = 0; i < n_dofs; ++i)
          {
            const auto &nonzero_components =
              fe[fe_index].get_nonzero_components(i);

            // Figure out whether the component mask applies. We assume
            // that we are allowed to set degrees of freedom if at least
            // one of the components (of the dof) is selected.
            bool selected = false;
            for (unsigned int c = 0; c < nonzero_components.size(); ++c)
              selected =
                selected || (nonzero_components[c] && component_mask[c]);

            if (selected)
              {
                if constexpr (running_in_debug_mode())
                  {
                    // make sure that all selected base elements are indeed
                    // interpolatory

                    if (const auto fe_system =
                          dynamic_cast<const FESystem<dim> *>(&fe[fe_index]))
                      {
                        const auto index =
                          fe_system->system_to_base_index(i).first.first;
                        Assert(fe_system->base_element(index)
                                 .has_generalized_support_points(),
                               ExcMessage("The component mask supplied to "
                                          "VectorTools::interpolate selects a "
                                          "non-interpolatory element."));
                      }
                  }

                // Add local values to the global vectors
                if (needs_expensive_algorithm[fe_index])
                  ::dealii::internal::ElementAccess<VectorType>::add(
                    dof_values[i], dofs_on_cell[i], interpolation);
                else
                  {
                    const auto base_index =
                      fe[fe_index].system_to_base_index(i);
                    ::dealii::internal::ElementAccess<VectorType>::add(
                      function_values[base_index.second]
                                     [base_index.first.second],
                      dofs_on_cell[i],
                      interpolation);
                  }
                ::dealii::internal::ElementAccess<VectorType>::add(
                  typename VectorType::value_type(1.0),
                  dofs_on_cell[i],
                  weights);
              }
            else
              {
                // If a component is ignored, copy the dof values
                // from the vector "vec", but only if they are locally
                // available
                if (locally_owned_dofs.is_element(dofs_on_cell[i]))
                  {
                    const auto value =
                      ::dealii::internal::ElementAccess<VectorType>::get(
                        vec, dofs_on_cell[i]);
                    ::dealii::internal::ElementAccess<VectorType>::add(
                      value, dofs_on_cell[i], interpolation);
                    ::dealii::internal::ElementAccess<VectorType>::add(
                      typename VectorType::value_type(1.0),
                      dofs_on_cell[i],
                      weights);
                  }
              }
          }
      };

      if (level == numbers::invalid_unsigned_int)
        {
          for (const auto &cell : dof_handler.active_cell_iterators())
            {
              if (cell->is_locally_owned())
                runner(cell);
            }
        }
      else
        {
          for (const auto &cell : dof_handler.mg_cell_iterators_on_level(level))
            if (cell->is_locally_owned_on_level())
              runner(cell);
        }

      interpolation.compress(VectorOperation::add);
      weights.compress(VectorOperation::add);

      for (const auto i : interpolation.locally_owned_elements())
        {
          const auto weight =
            ::dealii::internal::ElementAccess<VectorType>::get(weights, i);

          // See if we touched this DoF at all. If so, set the average
          // of the value we computed in the output vector. Otherwise,
          // don't touch the value at all.
          if (weight != number(0))
            {
              const auto value =
                ::dealii::internal::ElementAccess<VectorType>::get(
                  interpolation, i);
              ::dealii::internal::ElementAccess<VectorType>::set(value / weight,
                                                                 i,
                                                                 vec);
            }
        }
      vec.compress(VectorOperation::insert);
    }

  } // namespace internal



  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void interpolate(
    const hp::MappingCollection<dim, spacedim>                &mapping,
    const DoFHandler<dim, spacedim>                           &dof_handler,
    const Function<spacedim, typename VectorType::value_type> &function,
    VectorType                                                &vec,
    const ComponentMask                                       &component_mask,
    const unsigned int                                         level)
  {
    AssertDimension(dof_handler.get_fe_collection().n_components(),
                    function.n_components);

    // Create a small lambda capture wrapping function and call the
    // internal implementation
    const auto function_map = [&function](const auto &)
      -> const Function<spacedim, typename VectorType::value_type> * {
      return &function;
    };

    internal::interpolate(
      mapping, dof_handler, function_map, vec, component_mask, level);
  }



  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void interpolate(
    const Mapping<dim, spacedim>                              &mapping,
    const DoFHandler<dim, spacedim>                           &dof_handler,
    const Function<spacedim, typename VectorType::value_type> &function,
    VectorType                                                &vec,
    const ComponentMask                                       &component_mask,
    const unsigned int                                         level)
  {
    interpolate(hp::MappingCollection<dim, spacedim>(mapping),
                dof_handler,
                function,
                vec,
                component_mask,
                level);
  }



  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void interpolate(
    const DoFHandler<dim, spacedim>                           &dof,
    const Function<spacedim, typename VectorType::value_type> &function,
    VectorType                                                &vec,
    const ComponentMask                                       &component_mask,
    const unsigned int                                         level)
  {
    AssertDimension(dof.get_fe_collection().n_components(),
                    function.n_components);
    interpolate(get_default_linear_mapping(dof.get_triangulation()),
                dof,
                function,
                vec,
                component_mask,
                level);
  }



  template <int dim, class InVector, class OutVector, int spacedim>
  DEAL_II_CXX20_REQUIRES(concepts::is_dealii_vector_type<InVector> &&
                           concepts::is_writable_dealii_vector_type<OutVector>)
  void interpolate(const DoFHandler<dim, spacedim> &dof_1,
                   const DoFHandler<dim, spacedim> &dof_2,
                   const FullMatrix<double>        &transfer,
                   const InVector                  &data_1,
                   OutVector                       &data_2)
  {
    using number = typename OutVector::value_type;
    Vector<number> cell_data_1(dof_1.get_fe().n_dofs_per_cell());
    Vector<number> cell_data_2(dof_2.get_fe().n_dofs_per_cell());

    // Reset output vector.
    data_2 = static_cast<number>(0);

    // Store how many cells share each dof (unghosted).
    OutVector touch_count;
    touch_count.reinit(data_2);

    std::vector<types::global_dof_index> local_dof_indices(
      dof_2.get_fe().n_dofs_per_cell());

    typename DoFHandler<dim, spacedim>::active_cell_iterator cell_1 =
      dof_1.begin_active();
    typename DoFHandler<dim, spacedim>::active_cell_iterator cell_2 =
      dof_2.begin_active();
    const typename DoFHandler<dim, spacedim>::cell_iterator end_1 = dof_1.end();

    for (; cell_1 != end_1; ++cell_1, ++cell_2)
      {
        if (cell_1->is_locally_owned())
          {
            Assert(cell_2->is_locally_owned(), ExcInternalError());

            // Perform dof interpolation.
            cell_1->get_dof_values(data_1, cell_data_1);
            transfer.vmult(cell_data_2, cell_data_1);

            cell_2->get_dof_indices(local_dof_indices);

            // Distribute cell vector.
            for (unsigned int j = 0; j < dof_2.get_fe().n_dofs_per_cell(); ++j)
              {
                ::dealii::internal::ElementAccess<OutVector>::add(
                  cell_data_2(j), local_dof_indices[j], data_2);

                // Count cells that share each dof.
                ::dealii::internal::ElementAccess<OutVector>::add(
                  static_cast<number>(1), local_dof_indices[j], touch_count);
              }
          }
      }

    // Collect information over all the parallel processes.
    data_2.compress(VectorOperation::add);
    touch_count.compress(VectorOperation::add);

    // Compute the mean value of the sum which has been placed in
    // each entry of the output vector only at locally owned elements.
    for (const auto i : data_2.locally_owned_elements())
      {
        const number touch_count_i =
          ::dealii::internal::ElementAccess<OutVector>::get(touch_count, i);

        Assert(touch_count_i != static_cast<number>(0), ExcInternalError());

        const number value =
          ::dealii::internal::ElementAccess<OutVector>::get(data_2, i) /
          touch_count_i;

        ::dealii::internal::ElementAccess<OutVector>::set(value, i, data_2);
      }

    // Compress data_2 to set the proper values on all the parallel processes.
    data_2.compress(VectorOperation::insert);
  }



  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void get_position_vector(const DoFHandler<dim, spacedim> &dh,
                           VectorType                      &vector,
                           const ComponentMask             &mask)
  {
    const FiniteElement<dim, spacedim> &fe = dh.get_fe();
    get_position_vector(
      *fe.reference_cell().template get_default_mapping<dim, spacedim>(
        fe.degree),
      dh,
      vector,
      mask);
  }



  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void get_position_vector(const Mapping<dim, spacedim>    &map_q,
                           const DoFHandler<dim, spacedim> &dh,
                           VectorType                      &vector,
                           const ComponentMask             &mask)
  {
    AssertDimension(vector.size(), dh.n_dofs());
    const FiniteElement<dim, spacedim> &fe = dh.get_fe();

    // Construct default fe_mask;
    const ComponentMask fe_mask(
      mask.size() ? mask :
                    ComponentMask(fe.get_nonzero_components(0).size(), true));

    AssertDimension(fe_mask.size(), fe.get_nonzero_components(0).size());

    std::vector<unsigned int> fe_to_real(fe_mask.size(),
                                         numbers::invalid_unsigned_int);
    unsigned int              size = 0;
    for (unsigned int i = 0; i < fe_mask.size(); ++i)
      {
        if (fe_mask[i])
          fe_to_real[i] = size++;
      }
    Assert(
      size == spacedim,
      ExcMessage(
        "The Component Mask you provided is invalid. It has to select exactly spacedim entries."));


    if (fe.has_support_points())
      {
        const Quadrature<dim> quad(fe.get_unit_support_points());

        FEValues<dim, spacedim> fe_v(map_q, fe, quad, update_quadrature_points);
        std::vector<types::global_dof_index> dofs(fe.n_dofs_per_cell());

        AssertDimension(fe.n_dofs_per_cell(),
                        fe.get_unit_support_points().size());
        Assert(fe.is_primitive(),
               ExcMessage("FE is not Primitive! This won't work."));

        for (const auto &cell :
             dh.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
          {
            fe_v.reinit(cell);
            cell->get_dof_indices(dofs);
            const std::vector<Point<spacedim>> &points =
              fe_v.get_quadrature_points();
            for (unsigned int q = 0; q < points.size(); ++q)
              {
                const unsigned int comp = fe.system_to_component_index(q).first;
                if (fe_mask[comp])
                  ::dealii::internal::ElementAccess<VectorType>::set(
                    points[q][fe_to_real[comp]], dofs[q], vector);
              }
          }
      }
    else
      {
        // Construct a FiniteElement with FE_Q^spacedim, and call this
        // function again.
        //
        // Once we have this, interpolate with the given finite element
        // to get a Mapping which is interpolatory at the support points
        // of FE_Q(fe.degree())
        const FESystem<dim, spacedim> *fe_system =
          dynamic_cast<const FESystem<dim, spacedim> *>(&fe);
        Assert(fe_system, ExcNotImplemented());
        unsigned int degree = numbers::invalid_unsigned_int;

        // Get information about the blocks
        for (unsigned int i = 0; i < fe_mask.size(); ++i)
          if (fe_mask[i])
            {
              const unsigned int base_i =
                fe_system->component_to_base_index(i).first;
              Assert(degree == numbers::invalid_unsigned_int ||
                       degree == fe_system->base_element(base_i).degree,
                     ExcNotImplemented());
              Assert(fe_system->base_element(base_i).is_primitive(),
                     ExcNotImplemented());
              degree = fe_system->base_element(base_i).degree;
            }

        // We create an intermediate FE_Q vector space, and then
        // interpolate from that vector space to this one, by
        // carefully selecting the right components.

        FESystem<dim, spacedim>   feq(FE_Q<dim, spacedim>(degree), spacedim);
        DoFHandler<dim, spacedim> dhq(dh.get_triangulation());
        dhq.distribute_dofs(feq);
        Vector<double>      eulerq(dhq.n_dofs());
        const ComponentMask maskq(spacedim, true);
        get_position_vector(map_q, dhq, eulerq);

        FullMatrix<double>             transfer(fe.n_dofs_per_cell(),
                                    feq.n_dofs_per_cell());
        FullMatrix<double>             local_transfer(feq.n_dofs_per_cell());
        const std::vector<Point<dim>> &points = feq.get_unit_support_points();

        // Here we construct the interpolation matrix from
        // FE_Q^spacedim to the FiniteElement used by
        // euler_dof_handler.
        //
        // In order to construct such interpolation matrix, we have to
        // solve the following system:
        //
        // v_j phi_j(q_i) = w_k psi_k(q_i) = w_k delta_ki = w_i
        //
        // where psi_k are the basis functions for fe_q, and phi_i are
        // the basis functions of the target space while q_i are the
        // support points for the fe_q space. With this choice of
        // interpolation points, on the matrices is the identity
        // matrix, and we have to invert only one matrix. The
        // resulting vector will be interpolatory at the support
        // points of fe_q, even if the finite element does not have
        // support points.
        //
        // Morally, we should invert the matrix T_ij = phi_i(q_j),
        // however in general this matrix is not invertible, since
        // there may be components which do not contribute to the
        // displacement vector. Since we are not interested in those
        // components, we construct a square matrix with the same
        // number of components of the FE_Q system. The FE_Q system
        // was constructed above in such a way that the polynomial
        // degree of the FE_Q system and that of the given FE are the
        // same on the cell, which should guarantee that, for the
        // displacement components only, the interpolation matrix is
        // invertible. We construct a mapping between indices first,
        // and check that this is the case. If not, we bail out, not
        // knowing what to do in this case.

        std::vector<unsigned int> fe_to_feq(fe.n_dofs_per_cell(),
                                            numbers::invalid_unsigned_int);
        unsigned int              index = 0;
        for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
          if (fe_mask[fe.system_to_component_index(i).first])
            fe_to_feq[i] = index++;

        // If index is not the same as feq.n_dofs_per_cell(), we won't
        // know how to invert the resulting matrix. Bail out.
        Assert(index == feq.n_dofs_per_cell(), ExcNotImplemented());

        for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
          {
            const unsigned int comp_j = fe.system_to_component_index(j).first;
            if (fe_mask[comp_j])
              for (unsigned int i = 0; i < points.size(); ++i)
                {
                  if (fe_to_real[comp_j] ==
                      feq.system_to_component_index(i).first)
                    local_transfer(i, fe_to_feq[j]) =
                      fe.shape_value(j, points[i]);
                }
          }

        // Now we construct the rectangular interpolation matrix. This
        // one is filled only with the information from the components
        // of the displacement. The rest is set to zero.
        local_transfer.invert(local_transfer);
        for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
          if (fe_to_feq[i] != numbers::invalid_unsigned_int)
            for (unsigned int j = 0; j < feq.n_dofs_per_cell(); ++j)
              transfer(i, j) = local_transfer(fe_to_feq[i], j);

        // The interpolation matrix is then passed to the
        // VectorTools::interpolate() function to generate the correct
        // interpolation.
        interpolate(dhq, dh, transfer, eulerq, vector);
      }
  }

  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void interpolate_based_on_material_id(
    const Mapping<dim, spacedim>    &mapping,
    const DoFHandler<dim, spacedim> &dof_handler,
    const std::map<types::material_id,
                   const Function<spacedim, typename VectorType::value_type> *>
                        &functions,
    VectorType          &vec,
    const ComponentMask &component_mask)
  {
    // Create a small lambda capture wrapping the function map and call the
    // internal implementation
    const auto function_map =
      [&functions](
        const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell)
      -> const Function<spacedim, typename VectorType::value_type> * {
      const auto function = functions.find(cell->material_id());
      if (function != functions.end())
        return function->second;
      else
        return nullptr;
    };

    internal::interpolate(hp::MappingCollection<dim, spacedim>(mapping),
                          dof_handler,
                          function_map,
                          vec,
                          component_mask);
  }

  namespace internal
  {
    /**
     * Return whether the cell and all of its descendants are locally owned.
     */
    template <typename cell_iterator>
    bool
    is_locally_owned(const cell_iterator &cell)
    {
      if (cell->is_active())
        return cell->is_locally_owned();

      for (unsigned int c = 0; c < cell->n_children(); ++c)
        if (!is_locally_owned(cell->child(c)))
          return false;

      return true;
    }
  } // namespace internal



  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void interpolate_to_different_mesh(
    const DoFHandler<dim, spacedim> &dof_handler_1,
    const VectorType                &u1,
    const DoFHandler<dim, spacedim> &dof_handler_2,
    VectorType                      &u2)
  {
    Assert(GridTools::have_same_coarse_mesh(dof_handler_1, dof_handler_2),
           ExcMessage("The two DoF handlers must represent triangulations that "
                      "have the same coarse meshes"));

    InterGridMap<DoFHandler<dim, spacedim>> intergridmap;
    intergridmap.make_mapping(dof_handler_1, dof_handler_2);

    AffineConstraints<typename VectorType::value_type> dummy;
    dummy.close();

    interpolate_to_different_mesh(intergridmap, u1, dummy, u2);
  }



  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void interpolate_to_different_mesh(
    const DoFHandler<dim, spacedim>                          &dof_handler_1,
    const VectorType                                         &u1,
    const DoFHandler<dim, spacedim>                          &dof_handler_2,
    const AffineConstraints<typename VectorType::value_type> &constraints,
    VectorType                                               &u2)
  {
    Assert(GridTools::have_same_coarse_mesh(dof_handler_1, dof_handler_2),
           ExcMessage("The two DoF handlers must represent triangulations that "
                      "have the same coarse meshes"));

    InterGridMap<DoFHandler<dim, spacedim>> intergridmap;
    intergridmap.make_mapping(dof_handler_1, dof_handler_2);

    interpolate_to_different_mesh(intergridmap, u1, constraints, u2);
  }



  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void interpolate_to_different_mesh(
    const InterGridMap<DoFHandler<dim, spacedim>>            &intergridmap,
    const VectorType                                         &u1,
    const AffineConstraints<typename VectorType::value_type> &constraints,
    VectorType                                               &u2)
  {
    const DoFHandler<dim, spacedim> &dof_handler_1 =
      intergridmap.get_source_grid();
    const DoFHandler<dim, spacedim> &dof_handler_2 =
      intergridmap.get_destination_grid();

    Assert(
      dof_handler_1.get_fe_collection() == dof_handler_2.get_fe_collection(),
      ExcMessage("The FECollections of both DoFHandler objects must match"));
    AssertDimension(u1.size(), dof_handler_1.n_dofs());
    AssertDimension(u2.size(), dof_handler_2.n_dofs());

    Vector<typename VectorType::value_type> cache;

    // Looping over the finest common
    // mesh, this means that source and
    // destination cells have to be on the
    // same level and at least one has to
    // be active.
    //
    // Therefore, loop over all cells
    // (active and inactive) of the source
    // grid ..
    for (const auto &cell1 : dof_handler_1.cell_iterators())
      {
        const typename DoFHandler<dim, spacedim>::cell_iterator cell2 =
          intergridmap[cell1];

        // .. and skip if source and destination
        // cells are not on the same level ..
        if (cell1->level() != cell2->level())
          continue;
        // .. or none of them is active.
        if (!cell1->is_active() && !cell2->is_active())
          continue;

        Assert(
          internal::is_locally_owned(cell1) ==
            internal::is_locally_owned(cell2),
          ExcMessage(
            "The two Triangulations are required to have the same parallel partitioning."));

        // Skip foreign cells.
        if (cell1->is_active() && !cell1->is_locally_owned())
          continue;
        if (cell2->is_active() && !cell2->is_locally_owned())
          continue;

        // Get and set the corresponding
        // dof_values by interpolation.
        if (cell1->is_active())
          {
            cache.reinit(cell1->get_fe().n_dofs_per_cell());
            cell1->get_interpolated_dof_values(u1,
                                               cache,
                                               cell1->active_fe_index());
            cell2->set_dof_values_by_interpolation(cache,
                                                   u2,
                                                   cell1->active_fe_index());
          }
        else
          {
            cache.reinit(cell2->get_fe().n_dofs_per_cell());
            cell1->get_interpolated_dof_values(u1,
                                               cache,
                                               cell2->active_fe_index());
            cell2->set_dof_values_by_interpolation(cache,
                                                   u2,
                                                   cell2->active_fe_index());
          }
      }

    // finish the work on parallel vectors
    u2.compress(VectorOperation::insert);
    // Apply hanging node constraints.
    constraints.distribute(u2);
  }


  // Some helper functions. Much of the functions here is taken from a similar
  // implementation in data_out_dof_data.templates.h.
  namespace InterpolateBetweenMeshes
  {
    template <int dim, int spacedim, typename Number>
    void
    create_vector(const DoFHandler<dim, spacedim>            &dof_handler,
                  LinearAlgebra::distributed::Vector<Number> &u)
    {
      const IndexSet &locally_owned_dofs = dof_handler.locally_owned_dofs();

      const IndexSet locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(dof_handler);

      Assert(locally_owned_dofs.is_contiguous(),
             ExcMessage("You are trying to use vectors with non-contiguous "
                        "locally-owned index sets. This is not possible."));

      u.reinit(locally_owned_dofs,
               locally_relevant_dofs,
               dof_handler.get_mpi_communicator());
    }


    /**
     * Copy the data from an arbitrary non-block vector to a
     * LinearAlgebra::distributed::Vector.
     */
    template <typename VectorType, typename Number>
    void
    copy_locally_owned_data_from(
      const VectorType                           &src,
      LinearAlgebra::distributed::Vector<Number> &dst)
    {
      if constexpr (!IsBlockVector<VectorType>::value)
        {
          // If source and destination vector have the same underlying scalar,
          // we can directly import elements by using only one temporary vector:
          if constexpr (std::is_same_v<typename VectorType::value_type, Number>)
            {
              LinearAlgebra::ReadWriteVector<typename VectorType::value_type>
                temp;
              temp.reinit(src.locally_owned_elements());
              temp.import_elements(src, VectorOperation::insert);

              dst.import_elements(temp, VectorOperation::insert);
            }
          else
            // The source and destination vector have different scalar types. We
            // need to split the parallel import and local copy operations into
            // two phases
            {
              LinearAlgebra::ReadWriteVector<typename VectorType::value_type>
                temp;
              temp.reinit(src.locally_owned_elements());
              temp.import_elements(src, VectorOperation::insert);

              LinearAlgebra::ReadWriteVector<Number> temp2;
              temp2.reinit(temp, true);
              temp2 = temp;

              dst.import_elements(temp2, VectorOperation::insert);
            }
        }
      else
        DEAL_II_NOT_IMPLEMENTED();
    }

#ifdef DEAL_II_WITH_TRILINOS
    template <typename Number>
    void
    copy_locally_owned_data_from(
      const TrilinosWrappers::MPI::Vector        &src,
      LinearAlgebra::distributed::Vector<Number> &dst)
    {
      // ReadWriteVector does not work for ghosted
      // TrilinosWrappers::MPI::Vector objects. Fall back to copy the
      // entries manually.
      for (const auto i : dst.locally_owned_elements())
        dst[i] = src[i];
    }
#endif

#ifdef DEAL_II_TRILINOS_WITH_TPETRA
    template <typename Number, typename MemorySpace>
    void
    copy_locally_owned_data_from(
      const LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace> &src,
      LinearAlgebra::distributed::Vector<Number>                       &dst)
    {
      // ReadWriteVector does not work for ghosted
      // TrilinosWrappers::MPI::Vector objects. Fall back to copy the
      // entries manually.
      for (const auto i : dst.locally_owned_elements())
        dst[i] = src[i];
    }
#endif
  } // namespace InterpolateBetweenMeshes


  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void interpolate_to_coarser_mesh(
    const DoFHandler<dim, spacedim> &dof_handler_fine,
    const VectorType                &u_fine,
    const DoFHandler<dim, spacedim> &dof_handler_coarse,
    const AffineConstraints<typename VectorType::value_type>
               &constraints_coarse,
    VectorType &u_coarse)
  {
    Assert(GridTools::have_same_coarse_mesh(dof_handler_coarse,
                                            dof_handler_fine),
           ExcMessage("The two DoF handlers must represent triangulations that "
                      "have the same coarse meshes"));
    AssertDimension(dof_handler_fine.n_dofs(), u_fine.size());
    AssertDimension(dof_handler_coarse.n_dofs(), u_coarse.size());

    using LAVector =
      LinearAlgebra::distributed::Vector<typename VectorType::value_type>;

    // Create and initialize a version of the coarse vector using the
    // p::d::Vector type:
    LAVector my_u_coarse;
    InterpolateBetweenMeshes::create_vector(dof_handler_coarse, my_u_coarse);

    // Then do the same for the fine vector, and initialize it with a copy
    // of the source vector
    LAVector my_u_fine;
    InterpolateBetweenMeshes::create_vector(dof_handler_fine, my_u_fine);
    InterpolateBetweenMeshes::copy_locally_owned_data_from(u_fine, my_u_fine);
    my_u_fine.update_ghost_values();


    // Set up transfer operator. The transfer object also wants a constraints
    // object for the fine level, because it implements both the up- and
    // down-transfers. But we here only ever need the down-transfer, which
    // never touches the fine constraints. So we can use a dummy object for
    // that.
    AffineConstraints<typename VectorType::value_type> constraints_fine(
      dof_handler_fine.locally_owned_dofs(),
      DoFTools::extract_locally_relevant_dofs(dof_handler_fine));
    constraints_fine.close();

    MGTwoLevelTransfer<dim, LAVector> transfer;
    transfer.reinit(dof_handler_fine,
                    dof_handler_coarse,
                    constraints_fine,
                    constraints_coarse);

    // Then perform the interpolation from fine to coarse mesh:
    transfer.interpolate(my_u_coarse, my_u_fine);

    // Finally, apply the constraints to make sure that the resulting
    // object is conforming. Then copy it back into the user vector.
    // (At the time of writing, LinearAlgebra::EpetraWrappers::Vector
    // does not allow individual element access via operator() or
    // operator[]. So disallow this vector type for the moment here.)
    constraints_coarse.distribute(my_u_coarse);
    if constexpr (!std::is_same_v<VectorType,
                                  LinearAlgebra::EpetraWrappers::Vector>)
      {
        for (const auto i : u_coarse.locally_owned_elements())
          u_coarse(i) = my_u_coarse(i);
      }
    else
      DEAL_II_NOT_IMPLEMENTED();
  }



  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void interpolate_to_finer_mesh(
    const DoFHandler<dim, spacedim> &dof_handler_coarse,
    const VectorType                &u_coarse,
    const DoFHandler<dim, spacedim> &dof_handler_fine,
    const AffineConstraints<typename VectorType::value_type> &constraints_fine,
    VectorType                                               &u_fine)
  {
    Assert(GridTools::have_same_coarse_mesh(dof_handler_fine,
                                            dof_handler_coarse),
           ExcMessage("The two DoF handlers must represent triangulations that "
                      "have the same coarse meshes"));
    AssertDimension(dof_handler_fine.n_dofs(), u_fine.size());
    AssertDimension(dof_handler_coarse.n_dofs(), u_coarse.size());

    using LAVector =
      LinearAlgebra::distributed::Vector<typename VectorType::value_type>;

    // Create and initialize a version of the coarse vector using the
    // p::d::Vector type. Copy the source vector
    LAVector my_u_coarse;
    InterpolateBetweenMeshes::create_vector(dof_handler_coarse, my_u_coarse);
    InterpolateBetweenMeshes::copy_locally_owned_data_from(u_coarse,
                                                           my_u_coarse);
    my_u_coarse.update_ghost_values();

    // Then do the same for the fine vector, and initialize it with a copy
    // of the source vector
    LAVector my_u_fine;
    InterpolateBetweenMeshes::create_vector(dof_handler_fine, my_u_fine);

    // Set up transfer operator. The transfer object also wants a constraints
    // object for the coarse level, because it implements both the up- and
    // down-transfers. But we here only ever need the up-transfer, which
    // never touches the coarse constraints. So we can use a dummy object for
    // that.
    AffineConstraints<typename VectorType::value_type> constraints_coarse(
      dof_handler_coarse.locally_owned_dofs(),
      DoFTools::extract_locally_relevant_dofs(dof_handler_coarse));
    constraints_coarse.close();

    MGTwoLevelTransfer<dim, LAVector> transfer;
    transfer.reinit(dof_handler_fine,
                    dof_handler_coarse,
                    constraints_fine,
                    constraints_coarse);

    // Then perform the interpolation from coarse to fine mesh. (Note that
    // we add to my_u_fine, but that vector is initially a zero vector,
    // so we are really just writing into it.)
    transfer.prolongate_and_add(my_u_fine, my_u_coarse);

    // Finally, apply the constraints to make sure that the resulting
    // object is conforming. Then copy it back into the user vector.
    // (At the time of writing, LinearAlgebra::EpetraWrappers::Vector
    // does not allow individual element access via operator() or
    // operator[]. So disallow this vector type for the moment here.)
    constraints_fine.distribute(my_u_fine);
    if constexpr (!std::is_same_v<VectorType,
                                  LinearAlgebra::EpetraWrappers::Vector>)
      {
        for (const auto i : u_fine.locally_owned_elements())
          u_fine(i) = my_u_fine(i);
      }
    else
      DEAL_II_NOT_IMPLEMENTED();
  }


} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_interpolate_templates_h

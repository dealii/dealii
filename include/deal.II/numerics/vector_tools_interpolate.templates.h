// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_vector_tools_interpolate_templates_h
#define dealii_vector_tools_interpolate_templates_h


#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/intergrid_map.h>

#include <deal.II/hp/fe_values.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/vector_tools_interpolate.h>

DEAL_II_NAMESPACE_OPEN

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
              T3 &                function_values)
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
                    const FEValuesType &                fe_values_jacobians,
                    T3 &                                function_values)
    {
      if (const auto *system =
            dynamic_cast<const FESystem<dim, spacedim> *>(&fe))
        {
          // In case of an FESystem transform every (vector) component
          // separately:
          unsigned current_offset = offset;
          for (unsigned int i = 0; i < system->n_base_elements(); ++i)
            {
              const auto &base_fe      = system->base_element(i);
              const auto  multiplicity = system->element_multiplicity(i);
              for (unsigned int m = 0; m < multiplicity; ++m)
                {
                  // recursively call apply_transform to make sure to
                  // correctly handle nested fe systems.
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
    template <int dim,
              int spacedim,
              typename VectorType,
              template <int, int> class DoFHandlerType,
              typename T>
    void
    interpolate(const Mapping<dim, spacedim> &       mapping,
                const DoFHandlerType<dim, spacedim> &dof_handler,
                T &                                  function,
                VectorType &                         vec,
                const ComponentMask &                component_mask)
    {
      Assert(component_mask.represents_n_components(
               dof_handler.get_fe_collection().n_components()),
             ExcMessage(
               "The number of components in the mask has to be either "
               "zero or equal to the number of components in the finite "
               "element."));

      Assert(vec.size() == dof_handler.n_dofs(),
             ExcDimensionMismatch(vec.size(), dof_handler.n_dofs()));

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
      // variant for every fe we encounter to speed up resizing operations.
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

      hp::QCollection<dim> support_quadrature;
      for (unsigned int fe_index = 0; fe_index < fe.size(); ++fe_index)
        {
          const auto &points = fe[fe_index].get_generalized_support_points();
          support_quadrature.push_back(Quadrature<dim>(points));
        }

      const hp::MappingCollection<dim, spacedim> mapping_collection(mapping);

      // An FEValues object to evaluate (generalized) support point
      // locations as well as Jacobians and their inverses.
      // the latter are only needed for Hcurl or Hdiv conforming elements,
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

      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          // If this cell is not locally owned, do nothing.
          if (!cell->is_locally_owned())
            continue;

          const unsigned int fe_index = cell->active_fe_index();

          // Do nothing if there are no local degrees of freedom.
          if (fe[fe_index].dofs_per_cell == 0)
            continue;

          // Skip processing of the current cell if the function object is
          // invalid. This is used by interpolate_by_material_id to skip
          // interpolating over cells with unknown material id.
          if (!function(cell))
            continue;

          // Get transformed, generalized support points
          fe_values.reinit(cell);
          const std::vector<Point<spacedim>> &generalized_support_points =
            fe_values.get_present_fe_values().get_quadrature_points();

          // Get indices of the dofs on this cell
          const auto n_dofs = fe[fe_index].dofs_per_cell;
          dofs_on_cell.resize(n_dofs);
          cell->get_dof_indices(dofs_on_cell);

          // Prepare temporary storage
          auto &function_values = fe_function_values[fe_index];
          auto &dof_values      = fe_dof_values[fe_index];

          const auto n_components = fe[fe_index].n_components();
          function_values.resize(generalized_support_points.size(),
                                 Vector<number>(n_components));
          dof_values.resize(n_dofs);

          // Get all function values:
          Assert(
            n_components == function(cell)->n_components,
            ExcDimensionMismatch(dof_handler.get_fe_collection().n_components(),
                                 function(cell)->n_components));
          function(cell)->vector_value_list(generalized_support_points,
                                            function_values);

          {
            // Before we can average, we have to transform all function values
            // from the real cell back to the unit cell. We query the finite
            // element for the correct transformation. Matters get a bit more
            // complicated because we have to apply said transformation for
            // every base element.

            const unsigned int offset =
              apply_transform(fe[fe_index],
                              /* starting_offset = */ 0,
                              fe_values,
                              function_values);
            (void)offset;
            Assert(offset == n_components, ExcInternalError());
          }

          FETools::convert_generalized_support_point_values_to_dof_values(
            fe[fe_index], function_values, dof_values);

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
#ifdef DEBUG
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
#endif

                  // Add local values to the global vectors
                  ::dealii::internal::ElementAccess<VectorType>::add(
                    dof_values[i], dofs_on_cell[i], interpolation);
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
        } /* loop over dof_handler.active_cell_iterators() */

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



  template <int dim,
            int spacedim,
            typename VectorType,
            template <int, int> class DoFHandlerType>
  void
  interpolate(
    const Mapping<dim, spacedim> &                             mapping,
    const DoFHandlerType<dim, spacedim> &                      dof_handler,
    const Function<spacedim, typename VectorType::value_type> &function,
    VectorType &                                               vec,
    const ComponentMask &                                      component_mask)
  {
    Assert(dof_handler.get_fe_collection().n_components() ==
             function.n_components,
           ExcDimensionMismatch(dof_handler.get_fe_collection().n_components(),
                                function.n_components));

    // Create a small lambda capture wrapping function and call the
    // internal implementation
    const auto function_map = [&function](
      const typename DoFHandlerType<dim, spacedim>::active_cell_iterator &)
      -> const Function<spacedim, typename VectorType::value_type> *
    {
      return &function;
    };

    internal::interpolate(
      mapping, dof_handler, function_map, vec, component_mask);
  }



  template <int dim,
            int spacedim,
            typename VectorType,
            template <int, int> class DoFHandlerType>
  void
  interpolate(
    const DoFHandlerType<dim, spacedim> &                      dof,
    const Function<spacedim, typename VectorType::value_type> &function,
    VectorType &                                               vec,
    const ComponentMask &                                      component_mask)
  {
    interpolate(StaticMappingQ1<dim, spacedim>::mapping,
                dof,
                function,
                vec,
                component_mask);
  }



  template <int dim, class InVector, class OutVector, int spacedim>
  void
  interpolate(const DoFHandler<dim, spacedim> &dof_1,
              const DoFHandler<dim, spacedim> &dof_2,
              const FullMatrix<double> &       transfer,
              const InVector &                 data_1,
              OutVector &                      data_2)
  {
    using number = typename OutVector::value_type;
    Vector<number> cell_data_1(dof_1.get_fe().dofs_per_cell);
    Vector<number> cell_data_2(dof_2.get_fe().dofs_per_cell);

    // Reset output vector.
    data_2 = static_cast<number>(0);

    // Store how many cells share each dof (unghosted).
    OutVector touch_count;
    touch_count.reinit(data_2);

    std::vector<types::global_dof_index> local_dof_indices(
      dof_2.get_fe().dofs_per_cell);

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
            for (unsigned int j = 0; j < dof_2.get_fe().dofs_per_cell; ++j)
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
    for (const auto &i : data_2.locally_owned_elements())
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


  template <int dim,
            int spacedim,
            template <int, int> class DoFHandlerType,
            typename VectorType>
  void
  get_position_vector(const DoFHandlerType<dim, spacedim> &dh,
                      VectorType &                         vector,
                      const ComponentMask &                mask)
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

        MappingQGeneric<dim, spacedim> map_q(fe.degree);
        FEValues<dim, spacedim> fe_v(map_q, fe, quad, update_quadrature_points);
        std::vector<types::global_dof_index> dofs(fe.dofs_per_cell);

        AssertDimension(fe.dofs_per_cell, fe.get_unit_support_points().size());
        Assert(fe.is_primitive(),
               ExcMessage("FE is not Primitive! This won't work."));

        for (const auto &cell : dh.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              fe_v.reinit(cell);
              cell->get_dof_indices(dofs);
              const std::vector<Point<spacedim>> &points =
                fe_v.get_quadrature_points();
              for (unsigned int q = 0; q < points.size(); ++q)
                {
                  const unsigned int comp =
                    fe.system_to_component_index(q).first;
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

        FESystem<dim, spacedim> feq(FE_Q<dim, spacedim>(degree), spacedim);
        DoFHandlerType<dim, spacedim> dhq(dh.get_triangulation());
        dhq.distribute_dofs(feq);
        Vector<double>      eulerq(dhq.n_dofs());
        const ComponentMask maskq(spacedim, true);
        get_position_vector(dhq, eulerq);

        FullMatrix<double> transfer(fe.dofs_per_cell, feq.dofs_per_cell);
        FullMatrix<double> local_transfer(feq.dofs_per_cell);
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

        std::vector<unsigned int> fe_to_feq(fe.dofs_per_cell,
                                            numbers::invalid_unsigned_int);
        unsigned int              index = 0;
        for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
          if (fe_mask[fe.system_to_component_index(i).first])
            fe_to_feq[i] = index++;

        // If index is not the same as feq.dofs_per_cell, we won't
        // know how to invert the resulting matrix. Bail out.
        Assert(index == feq.dofs_per_cell, ExcNotImplemented());

        for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
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
        for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
          if (fe_to_feq[i] != numbers::invalid_unsigned_int)
            for (unsigned int j = 0; j < feq.dofs_per_cell; ++j)
              transfer(i, j) = local_transfer(fe_to_feq[i], j);

        // The interpolation matrix is then passed to the
        // VectorTools::interpolate() function to generate the correct
        // interpolation.
        interpolate(dhq, dh, transfer, eulerq, vector);
      }
  }

  template <int dim,
            int spacedim,
            typename VectorType,
            template <int, int> class DoFHandlerType>
  void
  interpolate_based_on_material_id(
    const Mapping<dim, spacedim> &       mapping,
    const DoFHandlerType<dim, spacedim> &dof_handler,
    const std::map<types::material_id,
                   const Function<spacedim, typename VectorType::value_type> *>
      &                  functions,
    VectorType &         vec,
    const ComponentMask &component_mask)
  {
    // Create a small lambda capture wrapping the function map and call the
    // internal implementation
    const auto function_map = [&functions](
      const typename DoFHandlerType<dim, spacedim>::active_cell_iterator &cell)
      -> const Function<spacedim, typename VectorType::value_type> *
    {
      const auto function = functions.find(cell->material_id());
      if (function != functions.end())
        return function->second;
      else
        return nullptr;
    };

    internal::interpolate(
      mapping, dof_handler, function_map, vec, component_mask);
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

  template <int dim,
            int spacedim,
            typename VectorType,
            template <int, int> class DoFHandlerType>
  void
  interpolate_to_different_mesh(const DoFHandlerType<dim, spacedim> &dof1,
                                const VectorType &                   u1,
                                const DoFHandlerType<dim, spacedim> &dof2,
                                VectorType &                         u2)
  {
    Assert(GridTools::have_same_coarse_mesh(dof1, dof2),
           ExcMessage("The two DoF handlers must represent triangulations that "
                      "have the same coarse meshes"));

    InterGridMap<DoFHandlerType<dim, spacedim>> intergridmap;
    intergridmap.make_mapping(dof1, dof2);

    AffineConstraints<typename VectorType::value_type> dummy;
    dummy.close();

    interpolate_to_different_mesh(intergridmap, u1, dummy, u2);
  }



  template <int dim,
            int spacedim,
            typename VectorType,
            template <int, int> class DoFHandlerType>
  void
  interpolate_to_different_mesh(
    const DoFHandlerType<dim, spacedim> &                     dof1,
    const VectorType &                                        u1,
    const DoFHandlerType<dim, spacedim> &                     dof2,
    const AffineConstraints<typename VectorType::value_type> &constraints,
    VectorType &                                              u2)
  {
    Assert(GridTools::have_same_coarse_mesh(dof1, dof2),
           ExcMessage("The two DoF handlers must represent triangulations that "
                      "have the same coarse meshes"));

    InterGridMap<DoFHandlerType<dim, spacedim>> intergridmap;
    intergridmap.make_mapping(dof1, dof2);

    interpolate_to_different_mesh(intergridmap, u1, constraints, u2);
  }

  template <int dim,
            int spacedim,
            typename VectorType,
            template <int, int> class DoFHandlerType>
  void
  interpolate_to_different_mesh(
    const InterGridMap<DoFHandlerType<dim, spacedim>> &       intergridmap,
    const VectorType &                                        u1,
    const AffineConstraints<typename VectorType::value_type> &constraints,
    VectorType &                                              u2)
  {
    const DoFHandlerType<dim, spacedim> &dof1 = intergridmap.get_source_grid();
    const DoFHandlerType<dim, spacedim> &dof2 =
      intergridmap.get_destination_grid();
    (void)dof2;

    Assert(dof1.get_fe_collection() == dof2.get_fe_collection(),
           ExcMessage(
             "The FECollections of both DoFHandler objects must match"));
    Assert(u1.size() == dof1.n_dofs(),
           ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
    Assert(u2.size() == dof2.n_dofs(),
           ExcDimensionMismatch(u2.size(), dof2.n_dofs()));

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
    typename DoFHandlerType<dim, spacedim>::cell_iterator cell1 = dof1.begin();
    const typename DoFHandlerType<dim, spacedim>::cell_iterator endc1 =
      dof1.end();

    for (; cell1 != endc1; ++cell1)
      {
        const typename DoFHandlerType<dim, spacedim>::cell_iterator cell2 =
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
            cache.reinit(cell1->get_fe().dofs_per_cell);
            cell1->get_interpolated_dof_values(u1,
                                               cache,
                                               cell1->active_fe_index());
            cell2->set_dof_values_by_interpolation(cache,
                                                   u2,
                                                   cell1->active_fe_index());
          }
        else
          {
            cache.reinit(cell2->get_fe().dofs_per_cell);
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
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_interpolate_templates_h

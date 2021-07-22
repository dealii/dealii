/*
 * copy_to_different_dof_handler.h
 *
 *  Created on: Apr 29, 2020
 *      Author: mwichro
 */

#ifndef INCLUDE_INTERPOLATE_TO_DIFFERENT_DOF_HANDLER_H_
#define INCLUDE_INTERPOLATE_TO_DIFFERENT_DOF_HANDLER_H_

#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/la_parallel_vector.h>

// apply_transform is defined here
// I am probably including wrong file, so please check this
#include <deal.II/numerics/vector_tools_interpolate.templates.h>

namespace dealii
{
  namespace VectorTools
  {
    // Modyfied internal interpolate
    //
    template <int dim,
              int spacedim,
              typename VectorType,
              template <int, int> class DoFHandlerType>
    void
    interpolate_to_different_dof_handler(
      VectorType &                                      dst,
      const DoFHandlerType<dim, spacedim> &             dof_handler,
      const unsigned int &                              dst_component,
      const LinearAlgebra::distributed::Vector<double> &src,
      const DoFHandlerType<dim, spacedim> &             src_dofs,
      const unsigned int &                              src_component,
      const types::material_id on_material_id = numbers::invalid_material_id,
      const Mapping<dim, spacedim> &mapping   = StaticMappingQ1<dim>::mapping)
    {
      src.update_ghost_values();

      FEValuesExtractors::Scalar src_extractor(src_component);
      ComponentMask              component_mask(
        dof_handler.get_fe_collection().n_components(), false);
      component_mask.set(dst_component, true);

      Assert(component_mask.n_selected_components(
               dof_handler.get_fe_collection().n_components()),
             ComponentMask::ExcNoComponentSelected());

      Assert(dst.size() == dof_handler.n_dofs(),
             ExcDimensionMismatch(dst.size(), dof_handler.n_dofs()));

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
      //  and weights.
      VectorType interpolation;
      VectorType weights;
      interpolation.reinit(dst);
      weights.reinit(dst);


      // Store locally owned dofs, so that we can skip all non-local dofs,
      // if they do not need to be interpolated.
      const IndexSet locally_owned_dofs = dst.locally_owned_elements();

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
                                              update_values | update_jacobians |
                                              update_inverse_jacobians);

      const hp::FECollection<dim, spacedim> &fe_src(
        src_dofs.get_fe_collection());

      hp::FEValues<dim, spacedim> fe_values_source(mapping_collection,
                                                   fe_src,
                                                   support_quadrature,
                                                   update_values);

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
          if (cell->material_id() != on_material_id &&
              on_material_id != numbers::invalid_material_id)
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



          // Interpolate extisting values on cell with other material ID
          //          if (cell->material_id() != on_material_id
          //        		  && on_material_id != numbers::invalid_material_id){
          //        	  fe_values.get_present_fe_values()
          //        			  .get_function_values(dst,function_values  );
          //
          //          }else

          {
            // Get cell iterator with src DoFHandler

            typename DoFHandlerType<dim, spacedim>::active_cell_iterator
              src_cell(&(src_dofs.get_triangulation()),
                       cell->level(),
                       cell->index(),
                       &src_dofs);
            // reinit fe_values_source
            fe_values_source.reinit(src_cell);

            // And get all function values:
            std::vector<double> source_values(
              generalized_support_points.size());
            fe_values_source.get_present_fe_values()[src_extractor]
              .get_function_values(src, source_values);
            for (unsigned int i = 0; i < source_values.size(); ++i)
              function_values[i](dst_component) = source_values[i];
          }

          //          function(cell)->vector_value_list(generalized_support_points,
          //                                            function_values);

          {
            // Before we can average, we have to transform all function values
            // from the real cell back to the unit cell. We query the finite
            // element for the correct transformation. Matters get a bit more
            // complicated because we have to apply said transformation for
            // every base element.

            const unsigned int offset =
              internal::apply_transform(fe[fe_index],
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
                  // from the vector "dst", but only if they are locally
                  // available
                  if (locally_owned_dofs.is_element(dofs_on_cell[i]))
                    {
                      const auto value =
                        ::dealii::internal::ElementAccess<VectorType>::get(
                          dst, dofs_on_cell[i]);
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
                                                                 dst);
            }
        }
      //      dst.compress(VectorOperation::insert);
      dst.update_ghost_values();
    }

  } // namespace VectorTools
} // namespace dealii
#endif /* INCLUDE_INTERPOLATE_TO_DIFFERENT_DOF_HANDLER_H_ */

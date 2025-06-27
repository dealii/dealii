// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/distributed/field_transfer.h>

#ifdef DEAL_II_WITH_P4EST

#  include <deal.II/lac/lapack_full_matrix.h>

#  include <limits>

DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace distributed
  {
    namespace experimental
    {
      template <int dim, typename VectorType, int spacedim>
      FieldTransfer<dim, VectorType, spacedim>::FieldTransfer(
        const DoFHandler<dim, spacedim> &dof)
        : dof_handler(dof)
      {
        // When coarsening, we want to mimic the behavior of SolutionTransfer
        // and interpolate from child cells to parent. Define this strategy here
        // since it is not readily available
        const auto coarsening_strategy =
          [this](
            const typename dealii::Triangulation<dim, spacedim>::cell_iterator
                                              &parent,
            const std::vector<Vector<Number>> &children_values) {
            // get the equivalent DoFCellAccessor
            typename DoFHandler<dim, spacedim>::cell_iterator dof_cell_iterator(
              &dof_handler.get_triangulation(),
              parent->level(),
              parent->index(),
              &dof_handler);

            int fe_index = 0;
            if (dof_handler.has_hp_capabilities())
              fe_index = dealii::internal::hp::DoFHandlerImplementation::
                dominated_future_fe_on_children<dim, spacedim>(
                  dof_cell_iterator);

            const auto &fe = dof_handler.get_fe(fe_index);
            Assert(fe.n_dofs_per_cell() > 0,
                   ExcMessage(
                     "Cannot coarsen onto a FiniteElement with no DoFs."));
            AssertDimension(dof_cell_iterator->n_children(),
                            children_values.size());

            const auto child_iterators = dof_cell_iterator->child_iterators();
            const unsigned int n_children_with_fe_nothing =
              std::count_if(child_iterators.begin(),
                            child_iterators.end(),
                            [](const auto &child_cell) {
                              return child_cell->get_fe().n_dofs_per_cell() ==
                                     0;
                            });

            Assert(
              n_children_with_fe_nothing == 0 ||
                n_children_with_fe_nothing == dof_cell_iterator->n_children(),
              ExcMessage(
                "Coarsening is only supported for parent cells where either all"
                " or none of the child cells are FE_Nothing."));

            // in case all children are FE_Nothing there is nothing to
            // interpolate and we just return the first entry from the children
            // values (containing invalid entries)
            if (n_children_with_fe_nothing == dof_cell_iterator->n_children())
              {
                return children_values[0];
              }

            const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
            Vector<Number>     tmp(dofs_per_cell);
            Vector<Number>     interpolated_values(dofs_per_cell);

            // Otherwise, perform the actual interpolation here. Due to the
            // assert above, we know that all child cells have data to
            // interpolate.
            for (unsigned int child = 0;
                 child < dof_cell_iterator->n_children();
                 ++child)
              {
                // interpolate the previously stored values on a child to the
                // mother cell
                fe.get_restriction_matrix(child,
                                          dof_cell_iterator->refinement_case())
                  .vmult(tmp, children_values[child]);

                // and add up or set them in the output vector
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  if (fe.restriction_is_additive(i))
                    interpolated_values(i) += tmp(i);
                  else if (tmp(i) != Number())
                    interpolated_values(i) = tmp(i);
              }

            return interpolated_values;
          };

        cell_data_transfer = std::make_unique<
          CellDataTransfer<dim, spacedim, std::vector<Vector<Number>>>>(
          dynamic_cast<
            dealii::parallel::distributed::Triangulation<dim, spacedim> &>(
            const_cast<dealii::Triangulation<dim, spacedim> &>(
              dof_handler.get_triangulation())),
          false,
          &dealii::AdaptationStrategies::Refinement::
            preserve<dim, spacedim, Vector<Number>>,
          coarsening_strategy);
      }



      template <int dim, typename VectorType, int spacedim>
      void
      FieldTransfer<dim, VectorType, spacedim>::
        prepare_for_coarsening_and_refinement(
          const VectorType  &in,
          const unsigned int fe_nothing_index)
      {
        const unsigned int dofs_per_cell =
          dof_handler.get_fe_collection().max_dofs_per_cell();

        Vector<Number> cell_solution(dofs_per_cell);
        Vector<Number> dummy_cell_solution(dofs_per_cell);

        for (auto &val : dummy_cell_solution)
          {
            val = std::numeric_limits<Number>::infinity();
          }

        in.update_ghost_values();

        std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
        for (const auto &cell : dof_handler.active_cell_iterators())
          {
            if ((cell->is_locally_owned()) &&
                (cell->active_fe_index() != fe_nothing_index))
              {
                cell->get_dof_values(in, cell_solution);
                data_to_transfer.push_back(cell_solution);
              }
            else
              {
                data_to_transfer.push_back(dummy_cell_solution);
              }
          }

        cell_data_transfer->prepare_for_coarsening_and_refinement(
          data_to_transfer);
      }



      template <int dim, typename VectorType, int spacedim>
      void
      FieldTransfer<dim, VectorType, spacedim>::interpolate(
        const Number                    &new_value,
        const AffineConstraints<Number> &affine_constraints,
        VectorType                      &out)
      {
        const unsigned int dofs_per_cell =
          dof_handler.get_fe_collection().max_dofs_per_cell();
        std::vector<Vector<Number>> transferred_data(
          dof_handler.get_triangulation().n_active_cells(),
          Vector<Number>(dofs_per_cell));

        cell_data_transfer->unpack(transferred_data);

        // Free the memory allocated by data_to_transfer
        data_to_transfer.clear();

        for (unsigned int i = 0; i < out.locally_owned_size(); ++i)
          out.local_element(i) = std::numeric_limits<Number>::infinity();

        unsigned int cell_i = 0;
        for (const auto &cell : dof_handler.active_cell_iterators())
          {
            if ((cell->is_locally_owned()) &&
                (transferred_data[cell_i][0] !=
                 std::numeric_limits<Number>::infinity()))
              {
                cell->set_dof_values(transferred_data[cell_i], out);
              }
            ++cell_i;
          }


        // Communicate the results.
        out.compress(VectorOperation::insert);

        // Treat hanging nodes
        std::vector<types::global_dof_index> dof_indices;
        std::vector<types::global_dof_index> dofs_map;
        std::vector<std::vector<std::pair<types::global_dof_index, Number>>>
                            constraint_lines;
        std::vector<Number> constraint_values;
        IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
        for (auto constrained_dof : locally_owned_dofs)
          if (affine_constraints.is_constrained(constrained_dof))
            {
              auto *constraint =
                affine_constraints.get_constraint_entries(constrained_dof);
              const unsigned int line_size = constraint->size();
              bool               add_line  = false;
              for (unsigned int i = 0; i < line_size; ++i)
                {
                  types::global_dof_index constraining_dof =
                    (*constraint)[i].first;
                  // If one of the constraining value is infinity, we need to
                  // reverse the relationship
                  if (out[constraining_dof] ==
                      std::numeric_limits<Number>::infinity())
                    {
                      add_line = true;
                      break;
                    }
                }

              if (add_line)
                {
                  std::vector<std::pair<types::global_dof_index, Number>> line;
                  Number val = out[constrained_dof];

                  for (unsigned int i = 0; i < line_size; ++i)
                    {
                      types::global_dof_index constraining_dof =
                        (*constraint)[i].first;

                      if (out[constraining_dof] ==
                          std::numeric_limits<Number>::infinity())
                        {
                          auto constraining_dof_map_it =
                            std::find(dofs_map.begin(),
                                      dofs_map.end(),
                                      constraining_dof);
                          if (constraining_dof_map_it == dofs_map.end())
                            {
                              dofs_map.push_back(constraining_dof);
                            }
                          line.push_back((*constraint)[i]);
                        }
                      else
                        {
                          val -=
                            out[constraining_dof] * (*constraint)[i].second;
                        }
                    }
                  constraint_lines.push_back(line);
                  constraint_values.push_back(val);
                }
            }

        // Build a constraint matrix that we invert
        const unsigned int n_rows = constraint_lines.size();
        if (n_rows > 0)
          {
            const unsigned int               n_cols = dofs_map.size();
            dealii::LAPACKFullMatrix<Number> constraints_matrix(n_rows, n_cols);
            dealii::Vector<Number>           constraint_values_vec(n_rows);

            for (unsigned int i = 0; i < n_rows; ++i)
              {
                for (unsigned int j = 0; j < n_cols; ++j)
                  {
                    if (j < constraint_lines[i].size())
                      {
                        auto constraint_it =
                          std::find(dofs_map.begin(),
                                    dofs_map.end(),
                                    constraint_lines[i][j].first);
                        constraints_matrix(i,
                                           constraint_it - dofs_map.begin()) =
                          constraint_lines[i][j].second;
                      }
                  }

                constraint_values_vec[i] = constraint_values[i];
              }

            constraints_matrix.compute_inverse_svd();

            dealii::Vector<Number> new_constrained_values(n_cols);
            constraints_matrix.vmult(new_constrained_values,
                                     constraint_values_vec);

            for (unsigned int i = 0; i < n_cols; ++i)
              {
                out[dofs_map[i]] = new_constrained_values[i];
              }
          }

        // Set the value to the newly create DoFs.
        std::for_each(out.begin(), out.end(), [&](Number &val) {
          if (val == std::numeric_limits<Number>::infinity())
            {
              val = new_value;
            }
        });
      }
    } // namespace experimental
  }   // namespace distributed
} // namespace parallel

// explicit instantiations
#  include "distributed/field_transfer.inst"

DEAL_II_NAMESPACE_CLOSE

#endif

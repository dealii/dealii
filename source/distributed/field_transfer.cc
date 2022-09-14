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
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

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
        cell_data_transfer = std::make_unique<
          CellDataTransfer<dim, spacedim, std::vector<Vector<Number>>>>(
          dynamic_cast<
            dealii::parallel::distributed::Triangulation<dim, spacedim> &>(
            const_cast<dealii::Triangulation<dim, spacedim> &>(
              dof_handler.get_triangulation())));
      }



      template <int dim, typename VectorType, int spacedim>
      void
      FieldTransfer<dim, VectorType, spacedim>::
        prepare_for_coarsening_and_refinement(
          const VectorType & in,
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
        const Number &                   new_value,
        const AffineConstraints<Number> &affine_constraints,
        VectorType &                     out)
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
        for (auto const &cell : dof_handler.active_cell_iterators())
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
        out.compress(dealii::VectorOperation::min);

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
#  include "field_transfer.inst"

DEAL_II_NAMESPACE_CLOSE

#endif

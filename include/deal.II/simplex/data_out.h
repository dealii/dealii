// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

#ifndef dealii_simplex_data_out_h
#define dealii_simplex_data_out_h

#include <deal.II/base/config.h>

#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/lac/la_parallel_vector.h>

DEAL_II_NAMESPACE_OPEN

namespace Simplex
{
  /**
   * This class is the main class to provide output of data described by
   * finite element fields defined on a collection of simplex cells.
   */
  class DataOut
  {
  public:
    /**
     * Write data to @p stream in Vtk.
     *
     * @note Currently, only working for the vector types Vector and
     *   LinearAlgebra::distributed::Vector.
     */
    template <int dim, int spacedim, typename VectorType, typename StreamType>
    void static write_vtk(const Mapping<dim, spacedim> &   mapping,
                          const DoFHandler<dim, spacedim> &dof_handler,
                          const VectorType &               vector,
                          const std::string &              label,
                          StreamType &                     stream)
    {
      Assert(dof_handler.hp_capability_enabled == false, ExcNotImplemented());

      IndexSet is_local, is_ghost;

      if (auto vector_ =
            dynamic_cast<const typename LinearAlgebra::distributed::Vector<
              typename VectorType::value_type> *>(&vector))
        {
          is_local = vector_->get_partitioner()->locally_owned_range();
          is_ghost = vector_->get_partitioner()->ghost_indices();
        }
      else
        {
          is_local = IndexSet(vector.size());
          is_local.add_range(0, vector.size());
          is_ghost = IndexSet(vector.size());
        }

      const unsigned int n_dofs = is_local.n_elements() + is_ghost.n_elements();

      const auto &       fe            = dof_handler.get_fe();
      const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

      const Quadrature<dim> quad(fe.get_unit_support_points());

      const UpdateFlags flag = update_quadrature_points;

      FEValues<dim, spacedim> fe_values(mapping, fe, quad, flag);

      std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

      std::vector<Point<spacedim>> all_points(n_dofs);
      unsigned int                 n_locally_owned_cells = 0;

      const auto global_to_local = [&](const types::global_dof_index index) {
        if (is_local.is_element(index))
          return is_local.index_within_set(index);

        if (is_ghost.is_element(index))
          return is_ghost.index_within_set(index) + is_local.n_elements();

        Assert(false, ExcNotImplemented());
      };

      // prepare points
      for (const auto &cell : dof_handler.cell_iterators())
        {
          if (!cell->is_locally_owned())
            continue;

          n_locally_owned_cells++;

          fe_values.reinit(cell);

          cell->get_dof_indices(dof_indices);

          const auto &points = fe_values.get_quadrature_points();

          for (unsigned int i = 0; i < dofs_per_cell; i++)
            all_points[global_to_local(dof_indices[i])] = points[i];
        }

      stream << "# vtk DataFile Version 2.0" << std::endl;
      stream << "Cube example" << std::endl;
      stream << "ASCII" << std::endl;
      stream << "DATASET UNSTRUCTURED_GRID" << std::endl;

      stream << "POINTS " << all_points.size() << " float" << std::endl;
      for (const auto &point : all_points)
        {
          for (int d = 0; d < spacedim; ++d)
            stream << point[d] << " ";
          for (int d = spacedim; d < 3; ++d)
            stream << 0.0 << " ";
          stream << std::endl;
        }

      stream << "CELLS " << n_locally_owned_cells << " "
             << n_locally_owned_cells * (dofs_per_cell + 1) << std::endl;

      for (const auto &cell : dof_handler.cell_iterators())
        {
          if (!cell->is_locally_owned())
            continue;

          FEValues<dim, spacedim> fe_values(mapping, fe, quad, flag);
          fe_values.reinit(cell);

          cell->get_dof_indices(dof_indices);

          stream << dofs_per_cell << " ";

          for (unsigned int i = 0; i < dofs_per_cell; i++)
            stream << global_to_local(dof_indices[i]) << " ";

          stream << std::endl;
        }

      stream << "CELL_TYPES " << n_locally_owned_cells << std::endl;

      const auto cell_type = [](const auto dofs_per_cell) {
        if (dim == 2 && dofs_per_cell == 3)
          return 5;
        if (dim == 2 && dofs_per_cell == 6)
          return 22;
        if (dim == 3 && dofs_per_cell == 4)
          return 10;
        if (dim == 3 && dofs_per_cell == 10)
          return 24;

        Assert(false, ExcNotImplemented());
      }(dofs_per_cell);

      for (unsigned int cell = 0; cell < n_locally_owned_cells; cell++)
        stream << cell_type << std::endl;

      stream << "POINT_DATA " << n_dofs << std::endl;
      stream << "SCALARS " << label << " double 1" << std::endl;
      stream << "LOOKUP_TABLE default" << std::endl;

      for (unsigned int i = 0; i < n_dofs; ++i)
        if (auto vector_ =
              dynamic_cast<const typename LinearAlgebra::distributed::Vector<
                typename VectorType::value_type> *>(&vector))
          stream << vector_->local_element(i) << std::endl;
        else
          stream << vector[i] << std::endl;
    }
  };
} // namespace Simplex

DEAL_II_NAMESPACE_CLOSE

#endif

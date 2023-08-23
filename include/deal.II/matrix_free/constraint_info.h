// ---------------------------------------------------------------------
//
// Copyright (C) 2022 - 2023 by the deal.II authors
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


#ifndef dealii_matrix_free_constraint_info_h
#define dealii_matrix_free_constraint_info_h


#include <deal.II/base/config.h>

#include <deal.II/base/floating_point_comparator.h>

#include <deal.II/matrix_free/dof_info.h>
#include <deal.II/matrix_free/evaluation_template_factory.h>
#include <deal.II/matrix_free/hanging_nodes_internal.h>
#include <deal.II/matrix_free/mapping_info.h>

#include <limits>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MatrixFreeFunctions
  {
    /**
     * A struct that takes entries describing a constraint and puts them into
     * a sorted list where duplicates are filtered out
     */
    template <typename Number>
    struct ConstraintValues
    {
      ConstraintValues();

      /**
       * This function inserts some constrained entries to the collection of
       * all values. It stores the (reordered) numbering of the dofs
       * (according to the ordering that matches with the function) in
       * new_indices, and returns the storage position the double array for
       * access later on.
       */
      template <typename number2>
      unsigned short
      insert_entries(
        const std::vector<std::pair<types::global_dof_index, number2>>
          &entries);

      /**
       * Return the memory consumption of the allocated memory in this class.
       */
      std::size_t
      memory_consumption() const;

      std::vector<std::pair<types::global_dof_index, double>>
                                           constraint_entries;
      std::vector<types::global_dof_index> constraint_indices;

      std::pair<std::vector<Number>, types::global_dof_index> next_constraint;
      std::map<std::vector<Number>,
               types::global_dof_index,
               FloatingPointComparator<Number>>
        constraints;
    };



    /**
     * A helper class to apply constraints in matrix-free loops in
     * user code. It combines constraint related functionalties from
     * MatrixFree and FEEvaluation.
     */
    template <int dim, typename Number>
    class ConstraintInfo
    {
    public:
      /**
       * Version 1: indices are extracted from DoFCellAccessor and
       * constraints are resolved with the help of AffineConstraints.
       */
      void
      reinit(const DoFHandler<dim> &dof_handler,
             const unsigned int     n_cells,
             const bool             use_fast_hanging_node_algorithm = true);

      void
      read_dof_indices(
        const unsigned int                                    cell_no,
        const unsigned int                                    mg_level,
        const TriaIterator<DoFCellAccessor<dim, dim, false>> &cell,
        const dealii::AffineConstraints<typename Number::value_type>
                                                                 &constraints,
        const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner);

      /**
       * Version 2: no constraints, indices are user-provided.
       */
      void
      reinit(const unsigned int n_cells);

      void
      read_dof_indices(
        const unsigned int                                        cell_no,
        const std::vector<types::global_dof_index>               &dof_indices,
        const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner);

      void
      finalize();

      template <typename T, typename VectorType>
      void
      read_write_operation(const T           &operation,
                           VectorType        &global_vector,
                           Number            *local_vector,
                           const unsigned int first_cell,
                           const unsigned int n_cells,
                           const unsigned int n_dofs_per_cell,
                           const bool         apply_constraints) const;

      void
      apply_hanging_node_constraints(
        const unsigned int     first_cell,
        const unsigned int     n_lanes_filled,
        const bool             transpose,
        AlignedVector<Number> &evaluation_data_coarse) const;

      /**
       * Return the memory consumption of the allocated memory in this class.
       */
      std::size_t
      memory_consumption() const;

    private:
      // for setup
      ConstraintValues<double>               constraint_values;
      std::vector<std::vector<unsigned int>> dof_indices_per_cell;
      std::vector<std::vector<unsigned int>> plain_dof_indices_per_cell;
      std::vector<std::vector<std::pair<unsigned short, unsigned short>>>
        constraint_indicator_per_cell;

      std::unique_ptr<HangingNodes<dim>>     hanging_nodes;
      std::vector<std::vector<unsigned int>> lexicographic_numbering;

    public:
      // for read_write_operation()
      std::vector<unsigned int> dof_indices;
      std::vector<std::pair<unsigned short, unsigned short>>
                                                         constraint_indicator;
      std::vector<std::pair<unsigned int, unsigned int>> row_starts;

      std::vector<unsigned int> plain_dof_indices;
      std::vector<unsigned int> row_starts_plain_indices;

      // for constraint_pool_begin/end()
      std::vector<typename Number::value_type> constraint_pool_data;
      std::vector<unsigned int>                constraint_pool_row_index;

      std::vector<ShapeInfo<Number>>          shape_infos;
      std::vector<compressed_constraint_kind> hanging_node_constraint_masks;
      std::vector<unsigned int>               active_fe_indices;

    private:
      inline const typename Number::value_type *
      constraint_pool_begin(const unsigned int row) const;

      inline const typename Number::value_type *
      constraint_pool_end(const unsigned int row) const;
    };



    // ------------------------- inline functions --------------------------

    // NOLINTNEXTLINE(modernize-use-equals-default)
    template <typename Number>
    ConstraintValues<Number>::ConstraintValues()
      : constraints(FloatingPointComparator<Number>(
          1. * std::numeric_limits<double>::epsilon() * 1024.))
    {}



    template <typename Number>
    template <typename number2>
    unsigned short
    ConstraintValues<Number>::insert_entries(
      const std::vector<std::pair<types::global_dof_index, number2>> &entries)
    {
      next_constraint.first.resize(entries.size());
      if (entries.size() > 0)
        {
          constraint_indices.resize(entries.size());
          // Use assign so that values for nonmatching Number / number2 are
          // converted:
          constraint_entries.assign(entries.begin(), entries.end());
          std::sort(constraint_entries.begin(),
                    constraint_entries.end(),
                    [](const std::pair<types::global_dof_index, double> &p1,
                       const std::pair<types::global_dof_index, double> &p2) {
                      return p1.second < p2.second;
                    });
          for (types::global_dof_index j = 0; j < constraint_entries.size();
               j++)
            {
              // copy the indices of the constraint entries after sorting.
              constraint_indices[j] = constraint_entries[j].first;

              // one_constraint takes the weights of the constraint
              next_constraint.first[j] = constraint_entries[j].second;
            }
        }

      // check whether or not constraint is already in pool. the initial
      // implementation computed a hash value based on the truncated array (to
      // given accuracy around 1e-13) in order to easily detect different
      // arrays and then made a fine-grained check when the hash values were
      // equal. this was quite lengthy and now we use a std::map with a
      // user-defined comparator to compare floating point arrays to a
      // tolerance 1e-13.
      types::global_dof_index insert_position = numbers::invalid_dof_index;
      const auto position = constraints.find(next_constraint.first);
      if (position != constraints.end())
        insert_position = position->second;
      else
        {
          next_constraint.second = constraints.size();
          constraints.insert(next_constraint);
          insert_position = next_constraint.second;
        }

      // we want to store the result as a short variable, so we have to make
      // sure that the result does not exceed the limits when casting.
      Assert(insert_position < (1 << (8 * sizeof(unsigned short))),
             ExcInternalError());
      return static_cast<unsigned short>(insert_position);
    }



    template <int dim, typename Number>
    inline void
    ConstraintInfo<dim, Number>::reinit(
      const DoFHandler<dim> &dof_handler,
      const unsigned int     n_cells,
      const bool             use_fast_hanging_node_algorithm)
    {
      this->dof_indices_per_cell.resize(n_cells);
      this->plain_dof_indices_per_cell.resize(n_cells);
      this->constraint_indicator_per_cell.resize(n_cells);

      // note: has_hanging_nodes() is a global operatrion
      const bool has_hanging_nodes =
        dof_handler.get_triangulation().has_hanging_nodes();

      if (use_fast_hanging_node_algorithm && has_hanging_nodes)
        {
          hanging_nodes = std::make_unique<HangingNodes<dim>>(
            dof_handler.get_triangulation());

          hanging_node_constraint_masks.resize(n_cells);
        }

      const auto &fes = dof_handler.get_fe_collection();
      lexicographic_numbering.resize(fes.size());
      shape_infos.resize(fes.size());

      for (unsigned int i = 0; i < fes.size(); ++i)
        {
          if (fes[i].reference_cell().is_hyper_cube())
            {
              const Quadrature<1> dummy_quadrature(
                std::vector<Point<1>>(1, Point<1>()));
              shape_infos[i].reinit(dummy_quadrature, fes[i], 0);
            }
          else
            {
              const auto dummy_quadrature =
                fes[i].reference_cell().template get_gauss_type_quadrature<dim>(
                  1);
              shape_infos[i].reinit(dummy_quadrature, fes[i], 0);
            }

          lexicographic_numbering[i] = shape_infos[i].lexicographic_numbering;
        }
      active_fe_indices.resize(n_cells);
    }



    template <int dim, typename Number>
    inline void
    ConstraintInfo<dim, Number>::reinit(const unsigned int n_cells)
    {
      this->dof_indices_per_cell.resize(n_cells);
      this->plain_dof_indices_per_cell.resize(0);
      this->constraint_indicator_per_cell.resize(n_cells);
    }



    template <int dim, typename Number>
    inline void
    ConstraintInfo<dim, Number>::read_dof_indices(
      const unsigned int                                            cell_no,
      const unsigned int                                            mg_level,
      const TriaIterator<DoFCellAccessor<dim, dim, false>>         &cell,
      const dealii::AffineConstraints<typename Number::value_type> &constraints,
      const std::shared_ptr<const Utilities::MPI::Partitioner>     &partitioner)
    {
      std::vector<types::global_dof_index> local_dof_indices(
        cell->get_fe().n_dofs_per_cell());
      std::vector<types::global_dof_index> local_dof_indices_lex(
        cell->get_fe().n_dofs_per_cell());

      if (mg_level == numbers::invalid_unsigned_int)
        cell->get_dof_indices(local_dof_indices);
      else
        cell->get_mg_dof_indices(local_dof_indices);

      {
        AssertIndexRange(cell->active_fe_index(), shape_infos.size());

        const auto &lexicographic_numbering =
          shape_infos[cell->active_fe_index()].lexicographic_numbering;

        AssertDimension(lexicographic_numbering.size(),
                        local_dof_indices.size());

        for (unsigned int i = 0; i < cell->get_fe().n_dofs_per_cell(); ++i)
          local_dof_indices_lex[i] =
            local_dof_indices[lexicographic_numbering[i]];
      }

      std::pair<unsigned short, unsigned short> constraint_iterator(0, 0);

      AssertIndexRange(cell_no, this->constraint_indicator_per_cell.size());
      AssertIndexRange(cell_no, this->dof_indices_per_cell.size());
      AssertIndexRange(cell_no, this->plain_dof_indices_per_cell.size());
      AssertIndexRange(cell_no, this->plain_dof_indices_per_cell.size());

      auto &constraint_indicator = this->constraint_indicator_per_cell[cell_no];
      auto &dof_indices          = this->dof_indices_per_cell[cell_no];
      auto &plain_dof_indices    = this->plain_dof_indices_per_cell[cell_no];

      AssertDimension(constraint_indicator_per_cell[cell_no].size(), 0);
      AssertDimension(dof_indices_per_cell[cell_no].size(), 0);
      AssertDimension(plain_dof_indices_per_cell[cell_no].size(), 0);

      const auto global_to_local =
        [&](const types::global_dof_index global_index) -> unsigned int {
        if (partitioner)
          return partitioner->global_to_local(global_index);
        else
          return global_index;
      };

      // plain indices
      plain_dof_indices.resize(local_dof_indices_lex.size());
      for (unsigned int i = 0; i < local_dof_indices_lex.size(); ++i)
        plain_dof_indices[i] = global_to_local(local_dof_indices_lex[i]);

      if (hanging_nodes)
        {
          AssertIndexRange(cell_no, this->hanging_node_constraint_masks.size());
          AssertIndexRange(cell_no, this->active_fe_indices.size());

          std::vector<ConstraintKinds> mask(cell->get_fe().n_components());
          hanging_nodes->setup_constraints(
            cell, {}, lexicographic_numbering, local_dof_indices_lex, mask);

          hanging_node_constraint_masks[cell_no] = compress(mask[0], dim);
          active_fe_indices[cell_no]             = cell->active_fe_index();
        }

      for (auto current_dof : local_dof_indices_lex)
        {
          const auto *entries_ptr =
            constraints.get_constraint_entries(current_dof);

          // dof is constrained
          if (entries_ptr != nullptr)
            {
              const auto                   &entries   = *entries_ptr;
              const types::global_dof_index n_entries = entries.size();
              if (n_entries == 1 &&
                  std::abs(entries[0].second - 1.) <
                    100 * std::numeric_limits<double>::epsilon())
                {
                  current_dof = entries[0].first;
                  goto no_constraint;
                }

              constraint_indicator.push_back(constraint_iterator);
              constraint_indicator.back().second =
                constraint_values.insert_entries(entries);

              // reset constraint iterator for next round
              constraint_iterator.first = 0;

              if (n_entries > 0)
                {
                  const std::vector<types::global_dof_index>
                    &constraint_indices = constraint_values.constraint_indices;
                  for (unsigned int j = 0; j < n_entries; ++j)
                    {
                      dof_indices.push_back(
                        global_to_local(constraint_indices[j]));
                    }
                }
            }
          else
            {
            no_constraint:
              dof_indices.push_back(global_to_local(current_dof));

              // make sure constraint_iterator.first is always within the
              // bounds of unsigned short
              Assert(constraint_iterator.first <
                       (1 << (8 * sizeof(unsigned short))) - 1,
                     ExcInternalError());
              constraint_iterator.first++;
            }
        }
    }



    template <int dim, typename Number>
    inline void
    ConstraintInfo<dim, Number>::read_dof_indices(
      const unsigned int                          cell_no,
      const std::vector<types::global_dof_index> &local_dof_indices_lex,
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner)
    {
      const auto global_to_local =
        [&](const types::global_dof_index global_index) -> unsigned int {
        if (partitioner)
          return partitioner->global_to_local(global_index);
        else
          return global_index;
      };

      std::pair<unsigned short, unsigned short> constraint_iterator(0, 0);

      auto &constraint_indicator = this->constraint_indicator_per_cell[cell_no];
      auto &dof_indices          = this->dof_indices_per_cell[cell_no];

      for (const auto current_dof : local_dof_indices_lex)
        {
          // dof is constrained
          if (current_dof == numbers::invalid_dof_index)
            {
              const std::vector<
                std::pair<types::global_dof_index, typename Number::value_type>>
                entries = {};

              constraint_indicator.push_back(constraint_iterator);
              constraint_indicator.back().second =
                constraint_values.insert_entries(entries);

              constraint_iterator.first = 0;
            }
          else
            {
              dof_indices.push_back(global_to_local(current_dof));

              // make sure constraint_iterator.first is always within the
              // bounds of unsigned short
              Assert(constraint_iterator.first <
                       (1 << (8 * sizeof(unsigned short))) - 1,
                     ExcInternalError());
              constraint_iterator.first++;
            }
        }
    }



    template <int dim, typename Number>
    inline void
    ConstraintInfo<dim, Number>::finalize()
    {
      this->dof_indices          = {};
      this->plain_dof_indices    = {};
      this->constraint_indicator = {};

      this->row_starts = {};
      this->row_starts.emplace_back(0, 0);

      if (plain_dof_indices_per_cell.empty() == false)
        {
          this->row_starts_plain_indices = {};
          this->row_starts_plain_indices.emplace_back(0);
        }

      for (unsigned int i = 0; i < dof_indices_per_cell.size(); ++i)
        {
          this->dof_indices.insert(this->dof_indices.end(),
                                   dof_indices_per_cell[i].begin(),
                                   dof_indices_per_cell[i].end());
          this->constraint_indicator.insert(
            this->constraint_indicator.end(),
            constraint_indicator_per_cell[i].begin(),
            constraint_indicator_per_cell[i].end());

          this->row_starts.emplace_back(this->dof_indices.size(),
                                        this->constraint_indicator.size());

          if (plain_dof_indices_per_cell.empty() == false)
            {
              this->plain_dof_indices.insert(
                this->plain_dof_indices.end(),
                plain_dof_indices_per_cell[i].begin(),
                plain_dof_indices_per_cell[i].end());

              this->row_starts_plain_indices.emplace_back(
                this->plain_dof_indices.size());
            }
        }

      std::vector<const std::vector<double> *> constraints(
        constraint_values.constraints.size());
      unsigned int length = 0;
      for (const auto &it : constraint_values.constraints)
        {
          AssertIndexRange(it.second, constraints.size());
          constraints[it.second] = &it.first;
          length += it.first.size();
        }

      constraint_pool_data.clear();
      constraint_pool_data.reserve(length);
      constraint_pool_row_index.reserve(constraint_values.constraints.size() +
                                        1);
      constraint_pool_row_index.resize(1, 0);

      for (const auto &constraint : constraints)
        {
          Assert(constraint != nullptr, ExcInternalError());
          constraint_pool_data.insert(constraint_pool_data.end(),
                                      constraint->begin(),
                                      constraint->end());
          constraint_pool_row_index.push_back(constraint_pool_data.size());
        }

      AssertDimension(constraint_pool_data.size(), length);

      dof_indices_per_cell.clear();
      constraint_indicator_per_cell.clear();

      if (hanging_nodes &&
          std::all_of(hanging_node_constraint_masks.begin(),
                      hanging_node_constraint_masks.end(),
                      [](const auto i) {
                        return i == unconstrained_compressed_constraint_kind;
                      }))
        hanging_node_constraint_masks.clear();
    }



    template <int dim, typename Number>
    template <typename T, typename VectorType>
    inline void
    ConstraintInfo<dim, Number>::read_write_operation(
      const T           &operation,
      VectorType        &global_vector,
      Number            *local_vector,
      const unsigned int first_cell,
      const unsigned int n_cells,
      const unsigned int n_dofs_per_cell,
      const bool         apply_constraints) const
    {
      if ((row_starts_plain_indices.empty() == false) &&
          (apply_constraints == false))
        {
          for (unsigned int v = 0; v < n_cells; ++v)
            {
              const unsigned int  cell_index = first_cell + v;
              const unsigned int *dof_indices =
                this->plain_dof_indices.data() +
                this->row_starts_plain_indices[cell_index];

              for (unsigned int i = 0; i < n_dofs_per_cell; ++dof_indices, ++i)
                operation.process_dof(*dof_indices,
                                      global_vector,
                                      local_vector[i][v]);
            }

          return;
        }

      for (unsigned int v = 0; v < n_cells; ++v)
        {
          const unsigned int  cell_index = first_cell + v;
          const unsigned int *dof_indices =
            this->dof_indices.data() + this->row_starts[cell_index].first;
          unsigned int index_indicators = this->row_starts[cell_index].second;
          unsigned int next_index_indicators =
            this->row_starts[cell_index + 1].second;

          unsigned int ind_local = 0;
          for (; index_indicators != next_index_indicators; ++index_indicators)
            {
              const std::pair<unsigned short, unsigned short> indicator =
                this->constraint_indicator[index_indicators];

              // run through values up to next constraint
              for (unsigned int j = 0; j < indicator.first; ++j)
                operation.process_dof(dof_indices[j],
                                      global_vector,
                                      local_vector[ind_local + j][v]);

              ind_local += indicator.first;
              dof_indices += indicator.first;

              // constrained case: build the local value as a linear
              // combination of the global value according to constraints
              typename Number::value_type value;
              operation.pre_constraints(local_vector[ind_local][v], value);

              const typename Number::value_type *data_val =
                this->constraint_pool_begin(indicator.second);
              const typename Number::value_type *end_pool =
                this->constraint_pool_end(indicator.second);
              for (; data_val != end_pool; ++data_val, ++dof_indices)
                operation.process_constraint(*dof_indices,
                                             *data_val,
                                             global_vector,
                                             value);

              operation.post_constraints(value, local_vector[ind_local][v]);
              ind_local++;
            }

          AssertIndexRange(ind_local, n_dofs_per_cell + 1);

          for (; ind_local < n_dofs_per_cell; ++dof_indices, ++ind_local)
            operation.process_dof(*dof_indices,
                                  global_vector,
                                  local_vector[ind_local][v]);
        }
    }



    template <int dim, typename Number>
    inline void
    ConstraintInfo<dim, Number>::apply_hanging_node_constraints(
      const unsigned int     first_cell,
      const unsigned int     n_lanes_filled,
      const bool             transpose,
      AlignedVector<Number> &evaluation_data_coarse) const
    {
      if (hanging_node_constraint_masks.empty())
        return;

      std::array<MatrixFreeFunctions::compressed_constraint_kind,
                 Number::size()>
        constraint_mask;

      bool hn_available = false;

      for (unsigned int v = 0; v < n_lanes_filled; ++v)
        {
          const auto mask = hanging_node_constraint_masks[first_cell + v];

          constraint_mask[v] = mask;

          hn_available |= (mask != unconstrained_compressed_constraint_kind);
        }

      if (hn_available == true)
        {
          for (unsigned int v = n_lanes_filled; v < Number::size(); ++v)
            constraint_mask[v] = unconstrained_compressed_constraint_kind;

          for (unsigned int i = 1; i < n_lanes_filled; ++i)
            AssertDimension(active_fe_indices[first_cell],
                            active_fe_indices[first_cell + i]);

          const auto &shape_info = shape_infos[active_fe_indices[first_cell]];

          dealii::internal::FEEvaluationHangingNodesFactory<
            dim,
            typename Number::value_type,
            Number>::apply(shape_info.n_components,
                           shape_info.data.front().fe_degree,
                           shape_info,
                           transpose,
                           constraint_mask,
                           evaluation_data_coarse.begin());
        }
    }



    template <int dim, typename Number>
    inline const typename Number::value_type *
    ConstraintInfo<dim, Number>::constraint_pool_begin(
      const unsigned int row) const
    {
      AssertIndexRange(row, constraint_pool_row_index.size() - 1);
      return constraint_pool_data.empty() ?
               nullptr :
               constraint_pool_data.data() + constraint_pool_row_index[row];
    }



    template <int dim, typename Number>
    inline const typename Number::value_type *
    ConstraintInfo<dim, Number>::constraint_pool_end(
      const unsigned int row) const
    {
      AssertIndexRange(row, constraint_pool_row_index.size() - 1);
      return constraint_pool_data.empty() ?
               nullptr :
               constraint_pool_data.data() + constraint_pool_row_index[row + 1];
    }



    template <int dim, typename Number>
    inline std::size_t
    ConstraintInfo<dim, Number>::memory_consumption() const
    {
      std::size_t size = 0;

      size += MemoryConsumption::memory_consumption(constraint_values);
      size += MemoryConsumption::memory_consumption(dof_indices_per_cell);
      size += MemoryConsumption::memory_consumption(plain_dof_indices_per_cell);
      size +=
        MemoryConsumption::memory_consumption(constraint_indicator_per_cell);

      if (hanging_nodes)
        size += MemoryConsumption::memory_consumption(*hanging_nodes);

      size += MemoryConsumption::memory_consumption(lexicographic_numbering);
      size += MemoryConsumption::memory_consumption(dof_indices);
      size += MemoryConsumption::memory_consumption(constraint_indicator);
      size += MemoryConsumption::memory_consumption(row_starts);
      size += MemoryConsumption::memory_consumption(plain_dof_indices);
      size += MemoryConsumption::memory_consumption(row_starts_plain_indices);
      size += MemoryConsumption::memory_consumption(constraint_pool_data);
      size += MemoryConsumption::memory_consumption(constraint_pool_row_index);
      size += MemoryConsumption::memory_consumption(shape_infos);
      size +=
        MemoryConsumption::memory_consumption(hanging_node_constraint_masks);
      size += MemoryConsumption::memory_consumption(active_fe_indices);

      return size;
    }



    template <typename Number>
    inline std::size_t
    ConstraintValues<Number>::memory_consumption() const
    {
      std::size_t size = 0;

      size += MemoryConsumption::memory_consumption(constraint_entries);
      size += MemoryConsumption::memory_consumption(constraint_indices);

      // TODO: map does not have memory_consumption()
      // size += MemoryConsumption::memory_consumption(constraints);

      return size;
    }

  } // namespace MatrixFreeFunctions
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif

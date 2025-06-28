// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/evaluation_kernels.h>
#include <deal.II/matrix_free/tensor_product_kernels.h>

#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_internal.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.templates.h>

#include <algorithm>

DEAL_II_NAMESPACE_OPEN


template <int dim, typename Number>
MGTransferMatrixFree<dim, Number>::MGTransferMatrixFree()
  : fe_degree(0)
  , element_is_continuous(false)
  , n_components(0)
  , n_child_cell_dofs(0)
{}



template <int dim, typename Number>
MGTransferMatrixFree<dim, Number>::MGTransferMatrixFree(
  const MGConstrainedDoFs &mg_c)
  : fe_degree(0)
  , element_is_continuous(false)
  , n_components(0)
  , n_child_cell_dofs(0)
{
  this->mg_constrained_dofs = &mg_c;
}



template <int dim, typename Number>
void
MGTransferMatrixFree<dim, Number>::initialize_constraints(
  const MGConstrainedDoFs &mg_c)
{
  this->mg_constrained_dofs = &mg_c;
}



template <int dim, typename Number>
void
MGTransferMatrixFree<dim, Number>::clear()
{
  this->MGLevelGlobalTransfer<
    LinearAlgebra::distributed::Vector<Number>>::clear();
  fe_degree             = 0;
  element_is_continuous = false;
  n_components          = 0;
  n_child_cell_dofs     = 0;
  level_dof_indices.clear();
  parent_child_connect.clear();
  dirichlet_indices.clear();
  n_owned_level_cells.clear();
  prolongation_matrix_1d.clear();
  evaluation_data.clear();
  weights_on_refined.clear();
}



template <int dim, typename Number>
void
MGTransferMatrixFree<dim, Number>::build(
  const DoFHandler<dim> &dof_handler,
  const std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
    &external_partitioners)
{
  Assert(dof_handler.has_level_dofs(),
         ExcMessage(
           "The underlying DoFHandler object has not had its "
           "distribute_mg_dofs() function called, but this is a prerequisite "
           "for multigrid transfers. You will need to call this function, "
           "probably close to where you already call distribute_dofs()."));
  if (external_partitioners.size() > 0)
    {
      Assert(this->initialize_dof_vector == nullptr,
             ExcMessage("An initialize_dof_vector function has already "
                        "been registered in the constructor!"));

      this->initialize_dof_vector =
        [external_partitioners](
          const unsigned int                          level,
          LinearAlgebra::distributed::Vector<Number> &vec) {
          vec.reinit(external_partitioners[level]);
        };
    }

  this->fill_and_communicate_copy_indices(dof_handler);

  vector_partitioners.resize(0,
                             dof_handler.get_triangulation().n_global_levels() -
                               1);
  for (unsigned int level = 0; level <= this->ghosted_level_vector.max_level();
       ++level)
    vector_partitioners[level] =
      this->ghosted_level_vector[level].get_partitioner();

  std::vector<std::vector<Number>> weights_unvectorized;

  internal::MGTransfer::ElementInfo<Number> elem_info;

  internal::MGTransfer::setup_transfer<dim, Number>(
    dof_handler,
    this->mg_constrained_dofs,
    external_partitioners,
    elem_info,
    level_dof_indices,
    parent_child_connect,
    n_owned_level_cells,
    dirichlet_indices,
    weights_unvectorized,
    this->copy_indices_global_mine,
    vector_partitioners);

  // reinit the ghosted level vector in case its partitioner differs from the
  // one the user intends to use externally
  this->ghosted_level_vector.resize(0, vector_partitioners.max_level());
  for (unsigned int level = 0; level <= vector_partitioners.max_level();
       ++level)
    if (external_partitioners.size() == vector_partitioners.max_level() + 1 &&
        external_partitioners[level].get() == vector_partitioners[level].get())
      this->ghosted_level_vector[level].reinit(0);
    else
      this->ghosted_level_vector[level].reinit(vector_partitioners[level]);

  // unpack element info data
  fe_degree             = elem_info.fe_degree;
  element_is_continuous = elem_info.element_is_continuous;
  n_components          = elem_info.n_components;
  n_child_cell_dofs     = elem_info.n_child_cell_dofs;

  // duplicate and put into vectorized array
  prolongation_matrix_1d.resize(elem_info.prolongation_matrix_1d.size());
  for (unsigned int i = 0; i < elem_info.prolongation_matrix_1d.size(); ++i)
    prolongation_matrix_1d[i] = elem_info.prolongation_matrix_1d[i];

  // reshuffle into aligned vector of vectorized arrays
  const unsigned int vec_size = VectorizedArray<Number>::size();
  const unsigned int n_levels =
    dof_handler.get_triangulation().n_global_levels();

  const unsigned int n_weights_per_cell = Utilities::fixed_power<dim>(3);
  weights_on_refined.resize(n_levels - 1);
  for (unsigned int level = 1; level < n_levels; ++level)
    {
      weights_on_refined[level - 1].resize(
        ((n_owned_level_cells[level - 1] + vec_size - 1) / vec_size) *
        n_weights_per_cell);

      for (unsigned int c = 0; c < n_owned_level_cells[level - 1]; ++c)
        {
          const unsigned int comp = c / vec_size;
          const unsigned int v    = c % vec_size;
          for (unsigned int i = 0; i < n_weights_per_cell; ++i)
            {
              weights_on_refined[level - 1][comp * n_weights_per_cell + i][v] =
                weights_unvectorized[level - 1][c * n_weights_per_cell + i];
            }
        }
    }

  evaluation_data.resize(n_child_cell_dofs);
}



template <int dim, typename Number>
void
MGTransferMatrixFree<dim, Number>::build(
  const DoFHandler<dim> &dof_handler,
  const std::function<void(const unsigned int,
                           LinearAlgebra::distributed::Vector<Number> &)>
    &initialize_dof_vector)
{
  if (initialize_dof_vector)
    {
      const unsigned int n_levels =
        dof_handler.get_triangulation().n_global_levels();

      std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
        external_partitioners(n_levels);

      for (unsigned int level = 0; level < n_levels; ++level)
        {
          LinearAlgebra::distributed::Vector<Number> vector;
          initialize_dof_vector(level, vector);
          external_partitioners[level] = vector.get_partitioner();
        }

      build(dof_handler, external_partitioners);
    }
  else
    {
      build(dof_handler);
    }
}



template <int dim, typename Number>
void
MGTransferMatrixFree<dim, Number>::prolongate(
  const unsigned int                                to_level,
  LinearAlgebra::distributed::Vector<Number>       &dst,
  const LinearAlgebra::distributed::Vector<Number> &src) const
{
  dst = 0;
  prolongate_and_add(to_level, dst, src);
}



template <int dim, typename Number>
void
MGTransferMatrixFree<dim, Number>::prolongate_and_add(
  const unsigned int                                to_level,
  LinearAlgebra::distributed::Vector<Number>       &dst,
  const LinearAlgebra::distributed::Vector<Number> &src) const
{
  Assert((to_level >= 1) && (to_level <= level_dof_indices.size()),
         ExcIndexRange(to_level, 1, level_dof_indices.size() + 1));

  const bool src_inplace = src.get_partitioner().get() ==
                           this->vector_partitioners[to_level - 1].get();
  if (src_inplace == false)
    {
      if (this->ghosted_level_vector[to_level - 1].get_partitioner().get() !=
          this->vector_partitioners[to_level - 1].get())
        this->ghosted_level_vector[to_level - 1].reinit(
          this->vector_partitioners[to_level - 1]);
      this->ghosted_level_vector[to_level - 1].copy_locally_owned_data_from(
        src);
    }

  const bool dst_inplace =
    dst.get_partitioner().get() == this->vector_partitioners[to_level].get();
  if (dst_inplace == false)
    {
      if (this->ghosted_level_vector[to_level].get_partitioner().get() !=
          this->vector_partitioners[to_level].get())
        this->ghosted_level_vector[to_level].reinit(
          this->vector_partitioners[to_level]);
      AssertDimension(this->ghosted_level_vector[to_level].locally_owned_size(),
                      dst.locally_owned_size());
      this->ghosted_level_vector[to_level] = 0.;
    }

  const LinearAlgebra::distributed::Vector<Number> &src_vec =
    src_inplace ? src : this->ghosted_level_vector[to_level - 1];
  LinearAlgebra::distributed::Vector<Number> &dst_vec =
    dst_inplace ? dst : this->ghosted_level_vector[to_level];

  src_vec.update_ghost_values();
  // the implementation in do_prolongate_add is templated in the degree of the
  // element (for efficiency reasons), so we need to find the appropriate
  // kernel here...
  if (fe_degree == 0)
    do_prolongate_add<0>(to_level, dst_vec, src_vec);
  else if (fe_degree == 1)
    do_prolongate_add<1>(to_level, dst_vec, src_vec);
  else if (fe_degree == 2)
    do_prolongate_add<2>(to_level, dst_vec, src_vec);
  else if (fe_degree == 3)
    do_prolongate_add<3>(to_level, dst_vec, src_vec);
  else if (fe_degree == 4)
    do_prolongate_add<4>(to_level, dst_vec, src_vec);
  else if (fe_degree == 5)
    do_prolongate_add<5>(to_level, dst_vec, src_vec);
  else if (fe_degree == 6)
    do_prolongate_add<6>(to_level, dst_vec, src_vec);
  else if (fe_degree == 7)
    do_prolongate_add<7>(to_level, dst_vec, src_vec);
  else if (fe_degree == 8)
    do_prolongate_add<8>(to_level, dst_vec, src_vec);
  else if (fe_degree == 9)
    do_prolongate_add<9>(to_level, dst_vec, src_vec);
  else if (fe_degree == 10)
    do_prolongate_add<10>(to_level, dst_vec, src_vec);
  else
    do_prolongate_add<-1>(to_level, dst_vec, src_vec);

  dst_vec.compress(VectorOperation::add);
  if (dst_inplace == false)
    dst += dst_vec;

  if (src_inplace == true)
    src.zero_out_ghost_values();
}



template <int dim, typename Number>
void
MGTransferMatrixFree<dim, Number>::restrict_and_add(
  const unsigned int                                from_level,
  LinearAlgebra::distributed::Vector<Number>       &dst,
  const LinearAlgebra::distributed::Vector<Number> &src) const
{
  Assert((from_level >= 1) && (from_level <= level_dof_indices.size()),
         ExcIndexRange(from_level, 1, level_dof_indices.size() + 1));

  const bool src_inplace =
    src.get_partitioner().get() == this->vector_partitioners[from_level].get();
  if (src_inplace == false)
    {
      if (this->ghosted_level_vector[from_level].get_partitioner().get() !=
          this->vector_partitioners[from_level].get())
        this->ghosted_level_vector[from_level].reinit(
          this->vector_partitioners[from_level]);
      this->ghosted_level_vector[from_level].copy_locally_owned_data_from(src);
    }

  const bool dst_inplace = dst.get_partitioner().get() ==
                           this->vector_partitioners[from_level - 1].get();
  if (dst_inplace == false)
    {
      if (this->ghosted_level_vector[from_level - 1].get_partitioner().get() !=
          this->vector_partitioners[from_level - 1].get())
        this->ghosted_level_vector[from_level - 1].reinit(
          this->vector_partitioners[from_level - 1]);
      AssertDimension(
        this->ghosted_level_vector[from_level - 1].locally_owned_size(),
        dst.locally_owned_size());
      this->ghosted_level_vector[from_level - 1] = 0.;
    }

  const LinearAlgebra::distributed::Vector<Number> &src_vec =
    src_inplace ? src : this->ghosted_level_vector[from_level];
  LinearAlgebra::distributed::Vector<Number> &dst_vec =
    dst_inplace ? dst : this->ghosted_level_vector[from_level - 1];

  src_vec.update_ghost_values();

  if (fe_degree == 0)
    do_restrict_add<0>(from_level, dst_vec, src_vec);
  else if (fe_degree == 1)
    do_restrict_add<1>(from_level, dst_vec, src_vec);
  else if (fe_degree == 2)
    do_restrict_add<2>(from_level, dst_vec, src_vec);
  else if (fe_degree == 3)
    do_restrict_add<3>(from_level, dst_vec, src_vec);
  else if (fe_degree == 4)
    do_restrict_add<4>(from_level, dst_vec, src_vec);
  else if (fe_degree == 5)
    do_restrict_add<5>(from_level, dst_vec, src_vec);
  else if (fe_degree == 6)
    do_restrict_add<6>(from_level, dst_vec, src_vec);
  else if (fe_degree == 7)
    do_restrict_add<7>(from_level, dst_vec, src_vec);
  else if (fe_degree == 8)
    do_restrict_add<8>(from_level, dst_vec, src_vec);
  else if (fe_degree == 9)
    do_restrict_add<9>(from_level, dst_vec, src_vec);
  else if (fe_degree == 10)
    do_restrict_add<10>(from_level, dst_vec, src_vec);
  else
    // go to the non-templated version of the evaluator
    do_restrict_add<-1>(from_level, dst_vec, src_vec);

  dst_vec.compress(VectorOperation::add);
  if (dst_inplace == false)
    dst += dst_vec;

  if (src_inplace == true)
    src.zero_out_ghost_values();
}



template <int dim, typename Number>
template <int degree>
void
MGTransferMatrixFree<dim, Number>::do_prolongate_add(
  const unsigned int                                to_level,
  LinearAlgebra::distributed::Vector<Number>       &dst,
  const LinearAlgebra::distributed::Vector<Number> &src) const
{
  const unsigned int vec_size        = VectorizedArray<Number>::size();
  const unsigned int degree_size     = (degree > -1 ? degree : fe_degree) + 1;
  const unsigned int n_child_dofs_1d = 2 * degree_size - element_is_continuous;
  const unsigned int n_scalar_cell_dofs =
    Utilities::fixed_power<dim>(n_child_dofs_1d);
  constexpr unsigned int three_to_dim = Utilities::pow(3, dim);

  // If we have user defined MG constraints, we must create
  // a non-const, ghosted version of the source vector to distribute
  // constraints.
  const LinearAlgebra::distributed::Vector<Number> *to_use = &src;
  LinearAlgebra::distributed::Vector<Number>        new_src;
  if (this->mg_constrained_dofs != nullptr &&
      this->mg_constrained_dofs->get_user_constraint_matrix(to_level - 1)
          .get_local_lines()
          .size() > 0)
    {
      LinearAlgebra::distributed::Vector<Number> copy_src(src);

      // Distribute any user defined constraints
      this->mg_constrained_dofs->get_user_constraint_matrix(to_level - 1)
        .distribute(copy_src);

      // Re-initialize new ghosted vector with correct constraints
      new_src.reinit(copy_src);
      new_src = copy_src;
      new_src.update_ghost_values();

      to_use = &new_src;
    }

  for (unsigned int cell = 0; cell < n_owned_level_cells[to_level - 1];
       cell += vec_size)
    {
      const unsigned int n_lanes =
        (cell + vec_size > n_owned_level_cells[to_level - 1]) ?
          (n_owned_level_cells[to_level - 1] - cell) :
          vec_size;

      // read from source vector
      for (unsigned int v = 0; v < n_lanes; ++v)
        {
          const unsigned int shift =
            internal::MGTransfer::compute_shift_within_children<dim>(
              parent_child_connect[to_level - 1][cell + v].second,
              fe_degree + 1 - element_is_continuous,
              fe_degree);
          const unsigned int *indices =
            &level_dof_indices[to_level - 1]
                              [parent_child_connect[to_level - 1][cell + v]
                                   .first *
                                 n_child_cell_dofs +
                               shift];
          for (unsigned int c = 0, m = 0; c < n_components; ++c)
            {
              for (unsigned int k = 0; k < (dim > 2 ? degree_size : 1); ++k)
                for (unsigned int j = 0; j < (dim > 1 ? degree_size : 1); ++j)
                  for (unsigned int i = 0; i < degree_size; ++i, ++m)
                    evaluation_data[m][v] = to_use->local_element(
                      indices[c * n_scalar_cell_dofs +
                              k * n_child_dofs_1d * n_child_dofs_1d +
                              j * n_child_dofs_1d + i]);

              // apply Dirichlet boundary conditions on parent cell
              for (std::vector<unsigned short>::const_iterator i =
                     dirichlet_indices[to_level - 1][cell + v].begin();
                   i != dirichlet_indices[to_level - 1][cell + v].end();
                   ++i)
                evaluation_data[*i][v] = 0.;
            }
        }

      AssertDimension(prolongation_matrix_1d.size(),
                      degree_size * n_child_dofs_1d);
      // perform tensorized operation
      if (element_is_continuous)
        {
          // must go through the components backwards because we want to write
          // the output to the same array as the input
          for (int c = n_components - 1; c >= 0; --c)
            internal::FEEvaluationImplBasisChange<
              internal::evaluate_general,
              internal::EvaluatorQuantity::value,
              dim,
              degree + 1,
              2 * degree + 1>::do_forward(1,
                                          prolongation_matrix_1d,
                                          evaluation_data.begin() +
                                            c * Utilities::fixed_power<dim>(
                                                  degree_size),
                                          evaluation_data.begin() +
                                            c * n_scalar_cell_dofs,
                                          fe_degree + 1,
                                          2 * fe_degree + 1);
          internal::weight_fe_q_dofs_by_entity<dim,
                                               degree != -1 ? 2 * degree + 1 :
                                                              -1,
                                               VectorizedArray<Number>>(
            &weights_on_refined[to_level - 1][(cell / vec_size) * three_to_dim],
            n_components,
            2 * fe_degree + 1,
            evaluation_data.begin());
        }
      else
        {
          for (int c = n_components - 1; c >= 0; --c)
            internal::FEEvaluationImplBasisChange<
              internal::evaluate_general,
              internal::EvaluatorQuantity::value,
              dim,
              degree + 1,
              2 * degree + 2>::do_forward(1,
                                          prolongation_matrix_1d,
                                          evaluation_data.begin() +
                                            c * Utilities::fixed_power<dim>(
                                                  degree_size),
                                          evaluation_data.begin() +
                                            c * n_scalar_cell_dofs,
                                          fe_degree + 1,
                                          2 * fe_degree + 2);
        }

      // write into dst vector
      const unsigned int *indices =
        &level_dof_indices[to_level][cell * n_child_cell_dofs];
      for (unsigned int v = 0; v < n_lanes; ++v)
        {
          for (unsigned int i = 0; i < n_child_cell_dofs; ++i)
            dst.local_element(indices[i]) += evaluation_data[i][v];
          indices += n_child_cell_dofs;
        }
    }
}



template <int dim, typename Number>
template <int degree>
void
MGTransferMatrixFree<dim, Number>::do_restrict_add(
  const unsigned int                                from_level,
  LinearAlgebra::distributed::Vector<Number>       &dst,
  const LinearAlgebra::distributed::Vector<Number> &src) const
{
  const unsigned int vec_size        = VectorizedArray<Number>::size();
  const unsigned int degree_size     = (degree > -1 ? degree : fe_degree) + 1;
  const unsigned int n_child_dofs_1d = 2 * degree_size - element_is_continuous;
  const unsigned int n_scalar_cell_dofs =
    Utilities::fixed_power<dim>(n_child_dofs_1d);
  constexpr unsigned int three_to_dim = Utilities::pow(3, dim);

  for (unsigned int cell = 0; cell < n_owned_level_cells[from_level - 1];
       cell += vec_size)
    {
      const unsigned int n_lanes =
        (cell + vec_size > n_owned_level_cells[from_level - 1]) ?
          (n_owned_level_cells[from_level - 1] - cell) :
          vec_size;

      // read from source vector
      {
        const unsigned int *indices =
          &level_dof_indices[from_level][cell * n_child_cell_dofs];
        for (unsigned int v = 0; v < n_lanes; ++v)
          {
            for (unsigned int i = 0; i < n_child_cell_dofs; ++i)
              evaluation_data[i][v] = src.local_element(indices[i]);
            indices += n_child_cell_dofs;
          }
      }

      AssertDimension(prolongation_matrix_1d.size(),
                      degree_size * n_child_dofs_1d);
      // perform tensorized operation
      if (element_is_continuous)
        {
          internal::weight_fe_q_dofs_by_entity<dim,
                                               degree != -1 ? 2 * degree + 1 :
                                                              -1,
                                               VectorizedArray<Number>>(
            &weights_on_refined[from_level - 1]
                               [(cell / vec_size) * three_to_dim],
            n_components,
            2 * fe_degree + 1,
            evaluation_data.data());
          for (unsigned int c = 0; c < n_components; ++c)
            internal::FEEvaluationImplBasisChange<
              internal::evaluate_general,
              internal::EvaluatorQuantity::value,
              dim,
              degree + 1,
              2 * degree + 1>::do_backward(1,
                                           prolongation_matrix_1d,
                                           false,
                                           evaluation_data.begin() +
                                             c * n_scalar_cell_dofs,
                                           evaluation_data.begin() +
                                             c * Utilities::fixed_power<dim>(
                                                   degree_size),
                                           fe_degree + 1,
                                           2 * fe_degree + 1);
        }
      else
        {
          for (unsigned int c = 0; c < n_components; ++c)
            internal::FEEvaluationImplBasisChange<
              internal::evaluate_general,
              internal::EvaluatorQuantity::value,
              dim,
              degree + 1,
              2 * degree + 2>::do_backward(1,
                                           prolongation_matrix_1d,
                                           false,
                                           evaluation_data.begin() +
                                             c * n_scalar_cell_dofs,
                                           evaluation_data.begin() +
                                             c * Utilities::fixed_power<dim>(
                                                   degree_size),
                                           fe_degree + 1,
                                           2 * fe_degree + 2);
        }

      // write into dst vector
      for (unsigned int v = 0; v < n_lanes; ++v)
        {
          const unsigned int shift =
            internal::MGTransfer::compute_shift_within_children<dim>(
              parent_child_connect[from_level - 1][cell + v].second,
              fe_degree + 1 - element_is_continuous,
              fe_degree);
          AssertIndexRange(
            parent_child_connect[from_level - 1][cell + v].first *
                n_child_cell_dofs +
              n_child_cell_dofs - 1,
            level_dof_indices[from_level - 1].size());
          const unsigned int *indices =
            &level_dof_indices[from_level - 1]
                              [parent_child_connect[from_level - 1][cell + v]
                                   .first *
                                 n_child_cell_dofs +
                               shift];
          for (unsigned int c = 0, m = 0; c < n_components; ++c)
            {
              // apply Dirichlet boundary conditions on parent cell
              for (std::vector<unsigned short>::const_iterator i =
                     dirichlet_indices[from_level - 1][cell + v].begin();
                   i != dirichlet_indices[from_level - 1][cell + v].end();
                   ++i)
                evaluation_data[*i][v] = 0.;

              for (unsigned int k = 0; k < (dim > 2 ? degree_size : 1); ++k)
                for (unsigned int j = 0; j < (dim > 1 ? degree_size : 1); ++j)
                  for (unsigned int i = 0; i < degree_size; ++i, ++m)
                    dst.local_element(
                      indices[c * n_scalar_cell_dofs +
                              k * n_child_dofs_1d * n_child_dofs_1d +
                              j * n_child_dofs_1d + i]) +=
                      evaluation_data[m][v];
            }
        }
    }
}



template <int dim, typename Number>
std::size_t
MGTransferMatrixFree<dim, Number>::memory_consumption() const
{
  std::size_t memory = MGLevelGlobalTransfer<
    LinearAlgebra::distributed::Vector<Number>>::memory_consumption();
  memory += MemoryConsumption::memory_consumption(level_dof_indices);
  memory += MemoryConsumption::memory_consumption(parent_child_connect);
  memory += MemoryConsumption::memory_consumption(n_owned_level_cells);
  memory += MemoryConsumption::memory_consumption(prolongation_matrix_1d);
  memory += MemoryConsumption::memory_consumption(evaluation_data);
  memory += MemoryConsumption::memory_consumption(weights_on_refined);
  memory += MemoryConsumption::memory_consumption(dirichlet_indices);
  return memory;
}



template <int dim, typename Number>
MGTransferBlockMatrixFree<dim, Number>::MGTransferBlockMatrixFree(
  const MGConstrainedDoFs &mg_c)
  : MGTransferBlockMatrixFreeBase<dim,
                                  Number,
                                  MGTransferMatrixFree<dim, Number>>(true)
{
  this->matrix_free_transfer_vector.emplace_back(mg_c);
}



template <int dim, typename Number>
MGTransferBlockMatrixFree<dim, Number>::MGTransferBlockMatrixFree(
  const std::vector<MGConstrainedDoFs> &mg_c)
  : MGTransferBlockMatrixFreeBase<dim,
                                  Number,
                                  MGTransferMatrixFree<dim, Number>>(false)
{
  for (const auto &constrained_block_dofs : mg_c)
    this->matrix_free_transfer_vector.emplace_back(constrained_block_dofs);
}



template <int dim, typename Number>
void
MGTransferBlockMatrixFree<dim, Number>::initialize_constraints(
  const MGConstrainedDoFs &mg_c)
{
  Assert(this->same_for_all,
         ExcMessage("This object was initialized with support for usage with "
                    "one DoFHandler for each block, but this method assumes "
                    "that the same DoFHandler is used for all the blocks!"));
  AssertDimension(this->matrix_free_transfer_vector.size(), 1);

  this->matrix_free_transfer_vector[0].initialize_constraints(mg_c);
}



template <int dim, typename Number>
void
MGTransferBlockMatrixFree<dim, Number>::initialize_constraints(
  const std::vector<MGConstrainedDoFs> &mg_c)
{
  Assert(!this->same_for_all,
         ExcMessage("This object was initialized with support for using "
                    "the same DoFHandler for all the blocks, but this "
                    "method assumes that there is a separate DoFHandler "
                    "for each block!"));
  AssertDimension(this->matrix_free_transfer_vector.size(), mg_c.size());

  for (unsigned int i = 0; i < mg_c.size(); ++i)
    this->matrix_free_transfer_vector[i].initialize_constraints(mg_c[i]);
}



template <int dim, typename Number>
void
MGTransferBlockMatrixFree<dim, Number>::build(
  const DoFHandler<dim> &dof_handler)
{
  AssertDimension(this->matrix_free_transfer_vector.size(), 1);
  this->matrix_free_transfer_vector[0].build(dof_handler);
}



template <int dim, typename Number>
void
MGTransferBlockMatrixFree<dim, Number>::build(
  const std::vector<const DoFHandler<dim> *> &dof_handler)
{
  AssertDimension(this->matrix_free_transfer_vector.size(), dof_handler.size());
  for (unsigned int i = 0; i < dof_handler.size(); ++i)
    this->matrix_free_transfer_vector[i].build(*dof_handler[i]);
}



template <int dim, typename Number>
std::size_t
MGTransferBlockMatrixFree<dim, Number>::memory_consumption() const
{
  std::size_t total_memory_consumption = 0;
  for (const auto &el : this->matrix_free_transfer_vector)
    total_memory_consumption += el.memory_consumption();
  return total_memory_consumption;
}



template <int dim, typename Number>
const MGTransferMatrixFree<dim, Number> &
MGTransferBlockMatrixFree<dim, Number>::get_matrix_free_transfer(
  const unsigned int b) const
{
  return matrix_free_transfer_vector[b];
}



template <int dim, typename Number>
void
MGTransferBlockMatrixFree<dim, Number>::clear()
{
  this->matrix_free_transfer_vector.clear();
}



// explicit instantiation
#include "multigrid/mg_transfer_matrix_free.inst"


DEAL_II_NAMESPACE_CLOSE

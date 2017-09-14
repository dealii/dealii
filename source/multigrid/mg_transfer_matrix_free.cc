// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
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


#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/mg_transfer_internal.h>

#include <deal.II/matrix_free/tensor_product_kernels.h>

#include <algorithm>

DEAL_II_NAMESPACE_OPEN


template <int dim, typename Number>
MGTransferMatrixFree<dim,Number>::MGTransferMatrixFree ()
  :
  fe_degree(0),
  element_is_continuous(false),
  n_components(0),
  n_child_cell_dofs(0)
{}



template <int dim, typename Number>
MGTransferMatrixFree<dim,Number>::MGTransferMatrixFree (const MGConstrainedDoFs &mg_c)
  :
  fe_degree(0),
  element_is_continuous(false),
  n_components(0),
  n_child_cell_dofs(0)
{
  this->mg_constrained_dofs = &mg_c;
}



template <int dim, typename Number>
void MGTransferMatrixFree<dim,Number>::initialize_constraints
(const MGConstrainedDoFs &mg_c)
{
  this->mg_constrained_dofs = &mg_c;
}



template <int dim, typename Number>
void MGTransferMatrixFree<dim,Number>::clear ()
{
  this->MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number> >::clear();
  fe_degree = 0;
  element_is_continuous = false;
  n_components = 0;
  n_child_cell_dofs = 0;
  level_dof_indices.clear();
  parent_child_connect.clear();
  dirichlet_indices.clear();
  n_owned_level_cells.clear();
  prolongation_matrix_1d.clear();
  evaluation_data.clear();
  weights_on_refined.clear();
}


template <int dim, typename Number>
void MGTransferMatrixFree<dim,Number>::build
(const DoFHandler<dim,dim>  &mg_dof)
{
  this->fill_and_communicate_copy_indices(mg_dof);

  std::vector<std::vector<Number> > weights_unvectorized;

  internal::MGTransfer::ElementInfo<Number> elem_info;

  internal::MGTransfer::setup_transfer<dim,Number>(mg_dof,
                                                   this->mg_constrained_dofs,
                                                   elem_info,
                                                   level_dof_indices,
                                                   parent_child_connect,
                                                   n_owned_level_cells,
                                                   dirichlet_indices,
                                                   weights_unvectorized,
                                                   this->copy_indices_global_mine,
                                                   this->ghosted_level_vector);
  // unpack element info data
  fe_degree                = elem_info.fe_degree;
  element_is_continuous    = elem_info.element_is_continuous;
  n_components             = elem_info.n_components;
  n_child_cell_dofs        = elem_info.n_child_cell_dofs;

  // duplicate and put into vectorized array
  prolongation_matrix_1d.resize(elem_info.prolongation_matrix_1d.size());
  for (unsigned int i=0; i<elem_info.prolongation_matrix_1d.size(); i++)
    prolongation_matrix_1d[i] = elem_info.prolongation_matrix_1d[i];

  // reshuffle into aligned vector of vectorized arrays
  const unsigned int vec_size = VectorizedArray<Number>::n_array_elements;
  const unsigned int n_levels = mg_dof.get_triangulation().n_global_levels();

  const unsigned int n_weights_per_cell = Utilities::fixed_power<dim>(3);
  weights_on_refined.resize(n_levels-1);
  for (unsigned int level = 1; level<n_levels; ++level)
    {
      weights_on_refined[level-1].resize(((n_owned_level_cells[level-1]+vec_size-1)/vec_size)*n_weights_per_cell);

      for (unsigned int c=0; c<n_owned_level_cells[level-1]; ++c)
        {
          const unsigned int comp = c/vec_size;
          const unsigned int v = c%vec_size;
          for (unsigned int i = 0; i<n_weights_per_cell; ++i)
            {

              weights_on_refined[level-1][comp*n_weights_per_cell+i][v] = weights_unvectorized[level-1][c*n_weights_per_cell+i];
            }
        }
    }

  evaluation_data.resize(3*n_child_cell_dofs);
}



template <int dim, typename Number>
void MGTransferMatrixFree<dim,Number>
::prolongate (const unsigned int                           to_level,
              LinearAlgebra::distributed::Vector<Number>       &dst,
              const LinearAlgebra::distributed::Vector<Number> &src) const
{
  Assert ((to_level >= 1) && (to_level<=level_dof_indices.size()),
          ExcIndexRange (to_level, 1, level_dof_indices.size()+1));

  AssertDimension(this->ghosted_level_vector[to_level].local_size(),
                  dst.local_size());
  AssertDimension(this->ghosted_level_vector[to_level-1].local_size(),
                  src.local_size());

  this->ghosted_level_vector[to_level-1] = src;
  this->ghosted_level_vector[to_level-1].update_ghost_values();
  this->ghosted_level_vector[to_level] = 0.;

  // the implementation in do_prolongate_add is templated in the degree of the
  // element (for efficiency reasons), so we need to find the appropriate
  // kernel here...
  if (fe_degree == 0)
    do_prolongate_add<0>(to_level, this->ghosted_level_vector[to_level],
                         this->ghosted_level_vector[to_level-1]);
  else if (fe_degree == 1)
    do_prolongate_add<1>(to_level, this->ghosted_level_vector[to_level],
                         this->ghosted_level_vector[to_level-1]);
  else if (fe_degree == 2)
    do_prolongate_add<2>(to_level, this->ghosted_level_vector[to_level],
                         this->ghosted_level_vector[to_level-1]);
  else if (fe_degree == 3)
    do_prolongate_add<3>(to_level, this->ghosted_level_vector[to_level],
                         this->ghosted_level_vector[to_level-1]);
  else if (fe_degree == 4)
    do_prolongate_add<4>(to_level, this->ghosted_level_vector[to_level],
                         this->ghosted_level_vector[to_level-1]);
  else if (fe_degree == 5)
    do_prolongate_add<5>(to_level, this->ghosted_level_vector[to_level],
                         this->ghosted_level_vector[to_level-1]);
  else if (fe_degree == 6)
    do_prolongate_add<6>(to_level, this->ghosted_level_vector[to_level],
                         this->ghosted_level_vector[to_level-1]);
  else if (fe_degree == 7)
    do_prolongate_add<7>(to_level, this->ghosted_level_vector[to_level],
                         this->ghosted_level_vector[to_level-1]);
  else if (fe_degree == 8)
    do_prolongate_add<8>(to_level, this->ghosted_level_vector[to_level],
                         this->ghosted_level_vector[to_level-1]);
  else if (fe_degree == 9)
    do_prolongate_add<9>(to_level, this->ghosted_level_vector[to_level],
                         this->ghosted_level_vector[to_level-1]);
  else if (fe_degree == 10)
    do_prolongate_add<10>(to_level, this->ghosted_level_vector[to_level],
                          this->ghosted_level_vector[to_level-1]);
  else
    do_prolongate_add<-1>(to_level, this->ghosted_level_vector[to_level],
                          this->ghosted_level_vector[to_level-1]);

  this->ghosted_level_vector[to_level].compress(VectorOperation::add);
  dst = this->ghosted_level_vector[to_level];
}



template <int dim, typename Number>
void MGTransferMatrixFree<dim,Number>
::restrict_and_add (const unsigned int                           from_level,
                    LinearAlgebra::distributed::Vector<Number>       &dst,
                    const LinearAlgebra::distributed::Vector<Number> &src) const
{
  Assert ((from_level >= 1) && (from_level<=level_dof_indices.size()),
          ExcIndexRange (from_level, 1, level_dof_indices.size()+1));

  AssertDimension(this->ghosted_level_vector[from_level].local_size(),
                  src.local_size());
  AssertDimension(this->ghosted_level_vector[from_level-1].local_size(),
                  dst.local_size());

  this->ghosted_level_vector[from_level] = src;
  this->ghosted_level_vector[from_level].update_ghost_values();
  this->ghosted_level_vector[from_level-1] = 0.;

  if (fe_degree == 0)
    do_restrict_add<0>(from_level, this->ghosted_level_vector[from_level-1],
                       this->ghosted_level_vector[from_level]);
  else if (fe_degree == 1)
    do_restrict_add<1>(from_level, this->ghosted_level_vector[from_level-1],
                       this->ghosted_level_vector[from_level]);
  else if (fe_degree == 2)
    do_restrict_add<2>(from_level, this->ghosted_level_vector[from_level-1],
                       this->ghosted_level_vector[from_level]);
  else if (fe_degree == 3)
    do_restrict_add<3>(from_level, this->ghosted_level_vector[from_level-1],
                       this->ghosted_level_vector[from_level]);
  else if (fe_degree == 4)
    do_restrict_add<4>(from_level, this->ghosted_level_vector[from_level-1],
                       this->ghosted_level_vector[from_level]);
  else if (fe_degree == 5)
    do_restrict_add<5>(from_level, this->ghosted_level_vector[from_level-1],
                       this->ghosted_level_vector[from_level]);
  else if (fe_degree == 6)
    do_restrict_add<6>(from_level, this->ghosted_level_vector[from_level-1],
                       this->ghosted_level_vector[from_level]);
  else if (fe_degree == 7)
    do_restrict_add<7>(from_level, this->ghosted_level_vector[from_level-1],
                       this->ghosted_level_vector[from_level]);
  else if (fe_degree == 8)
    do_restrict_add<8>(from_level, this->ghosted_level_vector[from_level-1],
                       this->ghosted_level_vector[from_level]);
  else if (fe_degree == 9)
    do_restrict_add<9>(from_level, this->ghosted_level_vector[from_level-1],
                       this->ghosted_level_vector[from_level]);
  else if (fe_degree == 10)
    do_restrict_add<10>(from_level, this->ghosted_level_vector[from_level-1],
                        this->ghosted_level_vector[from_level]);
  else
    // go to the non-templated version of the evaluator
    do_restrict_add<-1>(from_level, this->ghosted_level_vector[from_level-1],
                        this->ghosted_level_vector[from_level]);

  this->ghosted_level_vector[from_level-1].compress(VectorOperation::add);
  dst += this->ghosted_level_vector[from_level-1];
}



namespace
{
  template <int dim, typename Eval, typename Number, bool prolongate>
  void
  perform_tensorized_op(const Eval &evaluator,
                        const unsigned int n_child_cell_dofs,
                        const unsigned int n_components,
                        AlignedVector<VectorizedArray<Number> > &evaluation_data)
  {
    if (Eval::n_q_points != numbers::invalid_unsigned_int)
      AssertDimension(n_components * Eval::n_q_points, n_child_cell_dofs);
    VectorizedArray<Number> *t0 = &evaluation_data[0];
    VectorizedArray<Number> *t1 = &evaluation_data[n_child_cell_dofs];
    VectorizedArray<Number> *t2 = &evaluation_data[2*n_child_cell_dofs];

    for (unsigned int c=0; c<n_components; ++c)
      {
        // for the prolongate case, we go from dofs (living on the parent cell) to
        // quads (living on all children) in the FEEvaluation terminology
        if (dim == 1)
          evaluator.template values<0,prolongate,false>(t0, t2);
        else if (dim == 2)
          {
            evaluator.template values<0,prolongate,false>(t0, t1);
            evaluator.template values<1,prolongate,false>(t1, t2);
          }
        else if (dim == 3)
          {
            evaluator.template values<0,prolongate,false>(t0, t2);
            evaluator.template values<1,prolongate,false>(t2, t1);
            evaluator.template values<2,prolongate,false>(t1, t2);
          }
        else
          Assert(false, ExcNotImplemented());
        if (prolongate)
          {
            t0 += Eval::dofs_per_cell;
            t2 += Eval::n_q_points;
          }
        else
          {
            t0 += Eval::n_q_points;
            t2 += Eval::dofs_per_cell;
          }
      }
  }

  template <int dim, int degree, typename Number>
  void weight_dofs_on_child (const VectorizedArray<Number> *weights,
                             const unsigned int n_components,
                             const unsigned int fe_degree,
                             VectorizedArray<Number> *data)
  {
    Assert(fe_degree > 0, ExcNotImplemented());
    Assert(fe_degree < 100, ExcNotImplemented());
    const int loop_length = degree != -1 ? 2*degree+1 : 2*fe_degree+1;
    unsigned int degree_to_3 [100];
    degree_to_3[0] = 0;
    for (int i=1; i<loop_length-1; ++i)
      degree_to_3[i] = 1;
    degree_to_3[loop_length-1] = 2;
    for (unsigned int c=0; c<n_components; ++c)
      for (int k=0; k<(dim>2 ? loop_length : 1); ++k)
        for (int j=0; j<(dim>1 ? loop_length : 1); ++j)
          {
            const unsigned int shift = 9*degree_to_3[k] + 3*degree_to_3[j];
            data[0] *= weights[shift];
            // loop bound as int avoids compiler warnings in case loop_length
            // == 1 (polynomial degree 0)
            for (int i=1; i<loop_length-1; ++i)
              data[i] *= weights[shift+1];
            data[loop_length-1] *= weights[shift+2];
            data += loop_length;
          }
  }
}



template <int dim, typename Number>
template <int degree>
void MGTransferMatrixFree<dim,Number>
::do_prolongate_add (const unsigned int                           to_level,
                     LinearAlgebra::distributed::Vector<Number>       &dst,
                     const LinearAlgebra::distributed::Vector<Number> &src) const
{
  const unsigned int vec_size = VectorizedArray<Number>::n_array_elements;
  const unsigned int degree_size = (degree > -1 ? degree : fe_degree) + 1;
  const unsigned int n_child_dofs_1d = 2*degree_size - element_is_continuous;
  const unsigned int n_scalar_cell_dofs = Utilities::fixed_power<dim>(n_child_dofs_1d);
  const unsigned int three_to_dim = Utilities::fixed_int_power<3,dim>::value;

  for (unsigned int cell=0; cell < n_owned_level_cells[to_level-1];
       cell += vec_size)
    {
      const unsigned int n_chunks = cell+vec_size > n_owned_level_cells[to_level-1] ?
                                    n_owned_level_cells[to_level-1] - cell : vec_size;

      // read from source vector
      for (unsigned int v=0; v<n_chunks; ++v)
        {
          const unsigned int shift = internal::MGTransfer::compute_shift_within_children<dim>
                                     (parent_child_connect[to_level-1][cell+v].second,
                                      fe_degree+1-element_is_continuous, fe_degree);
          const unsigned int *indices = &level_dof_indices[to_level-1][parent_child_connect[to_level-1][cell+v].first*n_child_cell_dofs+shift];
          for (unsigned int c=0, m=0; c<n_components; ++c)
            {
              for (unsigned int k=0; k<(dim>2 ? degree_size : 1); ++k)
                for (unsigned int j=0; j<(dim>1 ? degree_size : 1); ++j)
                  for (unsigned int i=0; i<degree_size; ++i, ++m)
                    evaluation_data[m][v] =
                      src.local_element(indices[c*n_scalar_cell_dofs +
                                                k*n_child_dofs_1d*n_child_dofs_1d+
                                                j*n_child_dofs_1d+i]);

              // apply Dirichlet boundary conditions on parent cell
              for (std::vector<unsigned short>::const_iterator i=dirichlet_indices[to_level-1][cell+v].begin(); i!=dirichlet_indices[to_level-1][cell+v].end(); ++i)
                evaluation_data[*i][v] = 0.;
            }
        }

      AssertDimension(prolongation_matrix_1d.size(),
                      degree_size * n_child_dofs_1d);
      // perform tensorized operation
      if (element_is_continuous)
        {
          typedef internal::EvaluatorTensorProduct<internal::evaluate_general,dim,degree,degree!=-1 ? 2*degree+1 : 0,VectorizedArray<Number> > Evaluator;
          Evaluator evaluator(prolongation_matrix_1d,
                              prolongation_matrix_1d,
                              prolongation_matrix_1d,
                              fe_degree,
                              2*fe_degree+1);
          perform_tensorized_op<dim,Evaluator,Number,true>(evaluator,
                                                           n_child_cell_dofs,
                                                           n_components,
                                                           evaluation_data);
          weight_dofs_on_child<dim,degree,Number>(&weights_on_refined[to_level-1][(cell/vec_size)*three_to_dim],
                                                  n_components, fe_degree,
                                                  &evaluation_data[2*n_child_cell_dofs]);
        }
      else
        {
          typedef internal::EvaluatorTensorProduct<internal::evaluate_general,dim,degree,2*(degree+1),VectorizedArray<Number> > Evaluator;
          Evaluator evaluator(prolongation_matrix_1d,
                              prolongation_matrix_1d,
                              prolongation_matrix_1d,
                              fe_degree,
                              2*(fe_degree+1));
          perform_tensorized_op<dim,Evaluator,Number,true>(evaluator,
                                                           n_child_cell_dofs,
                                                           n_components,
                                                           evaluation_data);
        }

      // write into dst vector
      const unsigned int *indices = &level_dof_indices[to_level][cell*
                                                                 n_child_cell_dofs];
      for (unsigned int v=0; v<n_chunks; ++v)
        {
          for (unsigned int i=0; i<n_child_cell_dofs; ++i)
            dst.local_element(indices[i]) += evaluation_data[2*n_child_cell_dofs+i][v];
          indices += n_child_cell_dofs;
        }
    }
}



template <int dim, typename Number>
template <int degree>
void MGTransferMatrixFree<dim,Number>
::do_restrict_add (const unsigned int                           from_level,
                   LinearAlgebra::distributed::Vector<Number>       &dst,
                   const LinearAlgebra::distributed::Vector<Number> &src) const
{
  const unsigned int vec_size = VectorizedArray<Number>::n_array_elements;
  const unsigned int degree_size = (degree > -1 ? degree : fe_degree) + 1;
  const unsigned int n_child_dofs_1d = 2*degree_size - element_is_continuous;
  const unsigned int n_scalar_cell_dofs = Utilities::fixed_power<dim>(n_child_dofs_1d);
  const unsigned int three_to_dim = Utilities::fixed_int_power<3,dim>::value;

  for (unsigned int cell=0; cell < n_owned_level_cells[from_level-1];
       cell += vec_size)
    {
      const unsigned int n_chunks = cell+vec_size > n_owned_level_cells[from_level-1] ?
                                    n_owned_level_cells[from_level-1] - cell : vec_size;

      // read from source vector
      {
        const unsigned int *indices = &level_dof_indices[from_level][cell*
                                      n_child_cell_dofs];
        for (unsigned int v=0; v<n_chunks; ++v)
          {
            for (unsigned int i=0; i<n_child_cell_dofs; ++i)
              evaluation_data[i][v] = src.local_element(indices[i]);
            indices += n_child_cell_dofs;
          }
      }

      AssertDimension(prolongation_matrix_1d.size(),
                      degree_size * n_child_dofs_1d);
      // perform tensorized operation
      if (element_is_continuous)
        {
          typedef internal::EvaluatorTensorProduct<internal::evaluate_general,dim,degree,degree!=-1 ? 2*degree+1 : 0,VectorizedArray<Number> > Evaluator;
          Evaluator evaluator(prolongation_matrix_1d,
                              prolongation_matrix_1d,
                              prolongation_matrix_1d,
                              fe_degree,
                              2*fe_degree+1);
          weight_dofs_on_child<dim,degree,Number>(&weights_on_refined[from_level-1][(cell/vec_size)*three_to_dim],
                                                  n_components, fe_degree,
                                                  &evaluation_data[0]);
          perform_tensorized_op<dim,Evaluator,Number,false>(evaluator,
                                                            n_child_cell_dofs,
                                                            n_components,
                                                            evaluation_data);
        }
      else
        {
          typedef internal::EvaluatorTensorProduct<internal::evaluate_general,dim,degree,2*(degree+1),VectorizedArray<Number> > Evaluator;
          Evaluator evaluator(prolongation_matrix_1d,
                              prolongation_matrix_1d,
                              prolongation_matrix_1d,
                              fe_degree,
                              2*(fe_degree+1));
          perform_tensorized_op<dim,Evaluator,Number,false>(evaluator,
                                                            n_child_cell_dofs,
                                                            n_components,
                                                            evaluation_data);
        }

      // write into dst vector
      for (unsigned int v=0; v<n_chunks; ++v)
        {
          const unsigned int shift = internal::MGTransfer::compute_shift_within_children<dim>
                                     (parent_child_connect[from_level-1][cell+v].second,
                                      fe_degree+1-element_is_continuous, fe_degree);
          AssertIndexRange(parent_child_connect[from_level-1][cell+v].first*
                           n_child_cell_dofs+n_child_cell_dofs-1,
                           level_dof_indices[from_level-1].size());
          const unsigned int *indices = &level_dof_indices[from_level-1][parent_child_connect[from_level-1][cell+v].first*n_child_cell_dofs+shift];
          for (unsigned int c=0, m=0; c<n_components; ++c)
            {
              // apply Dirichlet boundary conditions on parent cell
              for (std::vector<unsigned short>::const_iterator i=dirichlet_indices[from_level-1][cell+v].begin(); i!=dirichlet_indices[from_level-1][cell+v].end(); ++i)
                evaluation_data[2*n_child_cell_dofs+(*i)][v] = 0.;

              for (unsigned int k=0; k<(dim>2 ? degree_size : 1); ++k)
                for (unsigned int j=0; j<(dim>1 ? degree_size : 1); ++j)
                  for (unsigned int i=0; i<degree_size; ++i, ++m)
                    dst.local_element(indices[c*n_scalar_cell_dofs +
                                              k*n_child_dofs_1d*n_child_dofs_1d+
                                              j*n_child_dofs_1d+i])
                    += evaluation_data[2*n_child_cell_dofs+m][v];
            }
        }
    }
}



template <int dim, typename Number>
std::size_t
MGTransferMatrixFree<dim,Number>::memory_consumption() const
{
  std::size_t memory = MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number> >::memory_consumption();
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
MGTransferBlockMatrixFree<dim,Number>::MGTransferBlockMatrixFree (const MGConstrainedDoFs &mg_c)
  :
  matrix_free_transfer(mg_c)
{
}



template <int dim, typename Number>
void MGTransferBlockMatrixFree<dim,Number>::initialize_constraints
(const MGConstrainedDoFs &mg_c)
{
  matrix_free_transfer.initialize_constraints(mg_c);
}



template <int dim, typename Number>
void MGTransferBlockMatrixFree<dim,Number>::clear ()
{
  matrix_free_transfer.clear();
}



template <int dim, typename Number>
void MGTransferBlockMatrixFree<dim,Number>::build
(const DoFHandler<dim,dim>  &mg_dof)
{
  matrix_free_transfer.build(mg_dof);
}



template <int dim, typename Number>
void MGTransferBlockMatrixFree<dim,Number>
::prolongate (const unsigned int                           to_level,
              LinearAlgebra::distributed::BlockVector<Number>       &dst,
              const LinearAlgebra::distributed::BlockVector<Number> &src) const
{
  AssertDimension(dst.n_blocks(), src.n_blocks());
  const unsigned int n_blocks = dst.n_blocks();
  for (unsigned int b = 0; b < n_blocks; ++b)
    matrix_free_transfer.prolongate(to_level, dst.block(b), src.block(b));
}



template <int dim, typename Number>
void MGTransferBlockMatrixFree<dim,Number>
::restrict_and_add (const unsigned int                           from_level,
                    LinearAlgebra::distributed::BlockVector<Number>       &dst,
                    const LinearAlgebra::distributed::BlockVector<Number> &src) const
{
  AssertDimension(dst.n_blocks(), src.n_blocks());
  const unsigned int n_blocks = dst.n_blocks();
  for (unsigned int b = 0; b < n_blocks; ++b)
    matrix_free_transfer.restrict_and_add(from_level, dst.block(b), src.block(b));
}



template <int dim, typename Number>
std::size_t
MGTransferBlockMatrixFree<dim,Number>::memory_consumption() const
{
  return matrix_free_transfer.memory_consumption();
}



// explicit instantiation
#include "mg_transfer_matrix_free.inst"


DEAL_II_NAMESPACE_CLOSE

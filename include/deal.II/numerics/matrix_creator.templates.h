// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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

#ifndef dealii_matrix_creator_templates_h
#define dealii_matrix_creator_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/matrix_tools.h>

#ifdef DEAL_II_WITH_PETSC
#  include <deal.II/lac/petsc_block_sparse_matrix.h>
#  include <deal.II/lac/petsc_sparse_matrix.h>
#  include <deal.II/lac/petsc_vector.h>
#endif

#ifdef DEAL_II_WITH_TRILINOS
#  include <deal.II/lac/trilinos_block_sparse_matrix.h>
#  include <deal.II/lac/trilinos_parallel_block_vector.h>
#  include <deal.II/lac/trilinos_sparse_matrix.h>
#  include <deal.II/lac/trilinos_vector.h>
#endif


#include <algorithm>
#include <cmath>
#include <set>


DEAL_II_NAMESPACE_OPEN
namespace MatrixCreator
{
  namespace internal
  {
    namespace AssemblerData
    {
      template <int dim, int spacedim, typename number>
      struct Scratch
      {
        Scratch(const ::dealii::hp::FECollection<dim, spacedim> &fe,
                const UpdateFlags                                update_flags,
                const Function<spacedim, number> *               coefficient,
                const Function<spacedim, number> *               rhs_function,
                const ::dealii::hp::QCollection<dim> &           quadrature,
                const ::dealii::hp::MappingCollection<dim, spacedim> &mapping)
          : fe_collection(fe)
          , quadrature_collection(quadrature)
          , mapping_collection(mapping)
          , x_fe_values(mapping_collection,
                        fe_collection,
                        quadrature_collection,
                        update_flags)
          , coefficient_values(quadrature_collection.max_n_quadrature_points())
          , coefficient_vector_values(
              quadrature_collection.max_n_quadrature_points(),
              dealii::Vector<number>(fe_collection.n_components()))
          , rhs_values(quadrature_collection.max_n_quadrature_points())
          , rhs_vector_values(quadrature_collection.max_n_quadrature_points(),
                              dealii::Vector<number>(
                                fe_collection.n_components()))
          , coefficient(coefficient)
          , rhs_function(rhs_function)
          , update_flags(update_flags)
        {}

        Scratch(const Scratch &data)
          : fe_collection(data.fe_collection)
          , quadrature_collection(data.quadrature_collection)
          , mapping_collection(data.mapping_collection)
          , x_fe_values(mapping_collection,
                        fe_collection,
                        quadrature_collection,
                        data.update_flags)
          , coefficient_values(data.coefficient_values)
          , coefficient_vector_values(data.coefficient_vector_values)
          , rhs_values(data.rhs_values)
          , rhs_vector_values(data.rhs_vector_values)
          , coefficient(data.coefficient)
          , rhs_function(data.rhs_function)
          , update_flags(data.update_flags)
        {}

        Scratch &
        operator=(const Scratch &)
        {
          Assert(false, ExcNotImplemented());
          return *this;
        }


        const ::dealii::hp::FECollection<dim, spacedim> &fe_collection;
        const ::dealii::hp::QCollection<dim> &           quadrature_collection;
        const ::dealii::hp::MappingCollection<dim, spacedim>
          &mapping_collection;

        ::dealii::hp::FEValues<dim, spacedim> x_fe_values;

        std::vector<number>                 coefficient_values;
        std::vector<dealii::Vector<number>> coefficient_vector_values;
        std::vector<number>                 rhs_values;
        std::vector<dealii::Vector<number>> rhs_vector_values;

        const Function<spacedim, number> *coefficient;
        const Function<spacedim, number> *rhs_function;

        const UpdateFlags update_flags;
      };


      template <typename number>
      struct CopyData
      {
        std::vector<types::global_dof_index> dof_indices;
        FullMatrix<number>                   cell_matrix;
        dealii::Vector<number>               cell_rhs;
        const AffineConstraints<number> *    constraints;
      };
    } // namespace AssemblerData


    template <int dim, int spacedim, typename CellIterator, typename number>
    void
    mass_assembler(
      const CellIterator &cell,
      MatrixCreator::internal::AssemblerData::Scratch<dim, spacedim, number>
        &                                                       data,
      MatrixCreator::internal::AssemblerData::CopyData<number> &copy_data)
    {
      data.x_fe_values.reinit(cell);
      const FEValues<dim, spacedim> &fe_values =
        data.x_fe_values.get_present_fe_values();

      const unsigned int dofs_per_cell       = fe_values.dofs_per_cell,
                         n_q_points          = fe_values.n_quadrature_points;
      const FiniteElement<dim, spacedim> &fe = fe_values.get_fe();
      const unsigned int                  n_components = fe.n_components();

      Assert(data.rhs_function == nullptr ||
               data.rhs_function->n_components == 1 ||
               data.rhs_function->n_components == n_components,
             ::dealii::MatrixCreator::ExcComponentMismatch());
      Assert(data.coefficient == nullptr ||
               data.coefficient->n_components == 1 ||
               data.coefficient->n_components == n_components,
             ::dealii::MatrixCreator::ExcComponentMismatch());

      copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
      copy_data.cell_rhs.reinit(dofs_per_cell);

      copy_data.dof_indices.resize(dofs_per_cell);
      cell->get_dof_indices(copy_data.dof_indices);

      const bool use_rhs_function = data.rhs_function != nullptr;
      if (use_rhs_function)
        {
          if (data.rhs_function->n_components == 1)
            {
              data.rhs_values.resize(n_q_points);
              data.rhs_function->value_list(fe_values.get_quadrature_points(),
                                            data.rhs_values);
            }
          else
            {
              data.rhs_vector_values.resize(
                n_q_points, dealii::Vector<number>(n_components));
              data.rhs_function->vector_value_list(
                fe_values.get_quadrature_points(), data.rhs_vector_values);
            }
        }

      const bool use_coefficient = data.coefficient != nullptr;
      if (use_coefficient)
        {
          if (data.coefficient->n_components == 1)
            {
              data.coefficient_values.resize(n_q_points);
              data.coefficient->value_list(fe_values.get_quadrature_points(),
                                           data.coefficient_values);
            }
          else
            {
              data.coefficient_vector_values.resize(
                n_q_points, dealii::Vector<number>(n_components));
              data.coefficient->vector_value_list(
                fe_values.get_quadrature_points(),
                data.coefficient_vector_values);
            }
        }


      const std::vector<double> &JxW = fe_values.get_JxW_values();
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        if (fe.is_primitive())
          {
            const unsigned int component_i =
              fe.system_to_component_index(i).first;
            const double *phi_i = &fe_values.shape_value(i, 0);

            // use symmetry in the mass matrix here:
            // just need to calculate the diagonal
            // and half of the elements above the
            // diagonal
            for (unsigned int j = i; j < dofs_per_cell; ++j)
              if ((n_components == 1) ||
                  (fe.system_to_component_index(j).first == component_i))
                {
                  const double *phi_j    = &fe_values.shape_value(j, 0);
                  number        add_data = 0;
                  if (use_coefficient)
                    {
                      if (data.coefficient->n_components == 1)
                        for (unsigned int point = 0; point < n_q_points;
                             ++point)
                          add_data +=
                            (data.coefficient_values[point] * phi_i[point] *
                             phi_j[point] * JxW[point]);
                      else
                        for (unsigned int point = 0; point < n_q_points;
                             ++point)
                          add_data +=
                            (data.coefficient_vector_values[point](
                               component_i) *
                             phi_i[point] * phi_j[point] * JxW[point]);
                    }
                  else
                    for (unsigned int point = 0; point < n_q_points; ++point)
                      add_data += phi_i[point] * phi_j[point] * JxW[point];

                  // this is even ok for i==j, since then
                  // we just write the same value twice.
                  copy_data.cell_matrix(i, j) = add_data;
                  copy_data.cell_matrix(j, i) = add_data;
                }

            if (use_rhs_function)
              {
                number add_data = 0;
                if (data.rhs_function->n_components == 1)
                  for (unsigned int point = 0; point < n_q_points; ++point)
                    add_data +=
                      data.rhs_values[point] * phi_i[point] * JxW[point];
                else
                  for (unsigned int point = 0; point < n_q_points; ++point)
                    add_data += data.rhs_vector_values[point](component_i) *
                                phi_i[point] * JxW[point];
                copy_data.cell_rhs(i) = add_data;
              }
          }
        else
          {
            // non-primitive vector-valued FE, using
            // symmetry again
            for (unsigned int j = i; j < dofs_per_cell; ++j)
              {
                number add_data = 0;
                for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
                  if (fe.get_nonzero_components(i)[comp_i] &&
                      fe.get_nonzero_components(j)[comp_i])
                    {
                      if (use_coefficient)
                        {
                          if (data.coefficient->n_components == 1)
                            for (unsigned int point = 0; point < n_q_points;
                                 ++point)
                              add_data +=
                                (data.coefficient_values[point] *
                                 fe_values.shape_value_component(i,
                                                                 point,
                                                                 comp_i) *
                                 fe_values.shape_value_component(j,
                                                                 point,
                                                                 comp_i) *
                                 JxW[point]);
                          else
                            for (unsigned int point = 0; point < n_q_points;
                                 ++point)
                              add_data +=
                                (data.coefficient_vector_values[point](comp_i) *
                                 fe_values.shape_value_component(i,
                                                                 point,
                                                                 comp_i) *
                                 fe_values.shape_value_component(j,
                                                                 point,
                                                                 comp_i) *
                                 JxW[point]);
                        }
                      else
                        for (unsigned int point = 0; point < n_q_points;
                             ++point)
                          add_data +=
                            fe_values.shape_value_component(i, point, comp_i) *
                            fe_values.shape_value_component(j, point, comp_i) *
                            JxW[point];
                    }

                copy_data.cell_matrix(i, j) = add_data;
                copy_data.cell_matrix(j, i) = add_data;
              }

            if (use_rhs_function)
              {
                number add_data = 0;
                for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
                  if (fe.get_nonzero_components(i)[comp_i])
                    {
                      if (data.rhs_function->n_components == 1)
                        for (unsigned int point = 0; point < n_q_points;
                             ++point)
                          add_data +=
                            data.rhs_values[point] *
                            fe_values.shape_value_component(i, point, comp_i) *
                            JxW[point];
                      else
                        for (unsigned int point = 0; point < n_q_points;
                             ++point)
                          add_data +=
                            data.rhs_vector_values[point](comp_i) *
                            fe_values.shape_value_component(i, point, comp_i) *
                            JxW[point];
                    }
                copy_data.cell_rhs(i) = add_data;
              }
          }
    }



    template <int dim, int spacedim, typename CellIterator>
    void
    laplace_assembler(
      const CellIterator &cell,
      MatrixCreator::internal::AssemblerData::Scratch<dim, spacedim, double>
        &                                                       data,
      MatrixCreator::internal::AssemblerData::CopyData<double> &copy_data)
    {
      data.x_fe_values.reinit(cell);
      const FEValues<dim, spacedim> &fe_values =
        data.x_fe_values.get_present_fe_values();

      const unsigned int dofs_per_cell       = fe_values.dofs_per_cell,
                         n_q_points          = fe_values.n_quadrature_points;
      const FiniteElement<dim, spacedim> &fe = fe_values.get_fe();
      const unsigned int                  n_components = fe.n_components();

      Assert(data.rhs_function == nullptr ||
               data.rhs_function->n_components == 1 ||
               data.rhs_function->n_components == n_components,
             ::dealii::MatrixCreator::ExcComponentMismatch());
      Assert(data.coefficient == nullptr ||
               data.coefficient->n_components == 1 ||
               data.coefficient->n_components == n_components,
             ::dealii::MatrixCreator::ExcComponentMismatch());

      copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
      copy_data.cell_rhs.reinit(dofs_per_cell);
      copy_data.dof_indices.resize(dofs_per_cell);
      cell->get_dof_indices(copy_data.dof_indices);


      const bool use_rhs_function = data.rhs_function != nullptr;
      if (use_rhs_function)
        {
          if (data.rhs_function->n_components == 1)
            {
              data.rhs_values.resize(n_q_points);
              data.rhs_function->value_list(fe_values.get_quadrature_points(),
                                            data.rhs_values);
            }
          else
            {
              data.rhs_vector_values.resize(
                n_q_points, dealii::Vector<double>(n_components));
              data.rhs_function->vector_value_list(
                fe_values.get_quadrature_points(), data.rhs_vector_values);
            }
        }

      const bool use_coefficient = data.coefficient != nullptr;
      if (use_coefficient)
        {
          if (data.coefficient->n_components == 1)
            {
              data.coefficient_values.resize(n_q_points);
              data.coefficient->value_list(fe_values.get_quadrature_points(),
                                           data.coefficient_values);
            }
          else
            {
              data.coefficient_vector_values.resize(
                n_q_points, dealii::Vector<double>(n_components));
              data.coefficient->vector_value_list(
                fe_values.get_quadrature_points(),
                data.coefficient_vector_values);
            }
        }


      const std::vector<double> &JxW = fe_values.get_JxW_values();
      double                     add_data;
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        if (fe.is_primitive())
          {
            const unsigned int component_i =
              fe.system_to_component_index(i).first;
            const Tensor<1, spacedim> *grad_phi_i = &fe_values.shape_grad(i, 0);

            // can use symmetry
            for (unsigned int j = i; j < dofs_per_cell; ++j)
              if ((n_components == 1) ||
                  (fe.system_to_component_index(j).first == component_i))
                {
                  const Tensor<1, spacedim> *grad_phi_j =
                    &fe_values.shape_grad(j, 0);
                  add_data = 0;
                  if (use_coefficient)
                    {
                      if (data.coefficient->n_components == 1)
                        for (unsigned int point = 0; point < n_q_points;
                             ++point)
                          add_data +=
                            ((grad_phi_i[point] * grad_phi_j[point]) *
                             JxW[point] * data.coefficient_values[point]);
                      else
                        for (unsigned int point = 0; point < n_q_points;
                             ++point)
                          add_data += ((grad_phi_i[point] * grad_phi_j[point]) *
                                       JxW[point] *
                                       data.coefficient_vector_values[point](
                                         component_i));
                    }
                  else
                    for (unsigned int point = 0; point < n_q_points; ++point)
                      add_data +=
                        (grad_phi_i[point] * grad_phi_j[point]) * JxW[point];

                  copy_data.cell_matrix(i, j) = add_data;
                  copy_data.cell_matrix(j, i) = add_data;
                }

            if (use_rhs_function)
              {
                const double *phi_i = &fe_values.shape_value(i, 0);
                add_data            = 0;
                if (data.rhs_function->n_components == 1)
                  for (unsigned int point = 0; point < n_q_points; ++point)
                    add_data +=
                      phi_i[point] * JxW[point] * data.rhs_values[point];
                else
                  for (unsigned int point = 0; point < n_q_points; ++point)
                    add_data += phi_i[point] * JxW[point] *
                                data.rhs_vector_values[point](component_i);
                copy_data.cell_rhs(i) = add_data;
              }
          }
        else
          {
            // non-primitive vector-valued FE
            for (unsigned int j = i; j < dofs_per_cell; ++j)
              {
                add_data = 0;
                for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
                  if (fe.get_nonzero_components(i)[comp_i] &&
                      fe.get_nonzero_components(j)[comp_i])
                    {
                      if (use_coefficient)
                        {
                          if (data.coefficient->n_components == 1)
                            for (unsigned int point = 0; point < n_q_points;
                                 ++point)
                              add_data +=
                                ((fe_values.shape_grad_component(i,
                                                                 point,
                                                                 comp_i) *
                                  fe_values.shape_grad_component(j,
                                                                 point,
                                                                 comp_i)) *
                                 JxW[point] * data.coefficient_values[point]);
                          else
                            for (unsigned int point = 0; point < n_q_points;
                                 ++point)
                              add_data +=
                                ((fe_values.shape_grad_component(i,
                                                                 point,
                                                                 comp_i) *
                                  fe_values.shape_grad_component(j,
                                                                 point,
                                                                 comp_i)) *
                                 JxW[point] *
                                 data.coefficient_vector_values[point](comp_i));
                        }
                      else
                        for (unsigned int point = 0; point < n_q_points;
                             ++point)
                          add_data +=
                            (fe_values.shape_grad_component(i, point, comp_i) *
                             fe_values.shape_grad_component(j, point, comp_i)) *
                            JxW[point];
                    }

                copy_data.cell_matrix(i, j) = add_data;
                copy_data.cell_matrix(j, i) = add_data;
              }

            if (use_rhs_function)
              {
                add_data = 0;
                for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
                  if (fe.get_nonzero_components(i)[comp_i])
                    {
                      if (data.rhs_function->n_components == 1)
                        for (unsigned int point = 0; point < n_q_points;
                             ++point)
                          add_data +=
                            fe_values.shape_value_component(i, point, comp_i) *
                            JxW[point] * data.rhs_values[point];
                      else
                        for (unsigned int point = 0; point < n_q_points;
                             ++point)
                          add_data +=
                            fe_values.shape_value_component(i, point, comp_i) *
                            JxW[point] * data.rhs_vector_values[point](comp_i);
                    }
                copy_data.cell_rhs(i) = add_data;
              }
          }
    }



    template <typename number, typename MatrixType, typename VectorType>
    void
    copy_local_to_global(const AssemblerData::CopyData<number> &data,
                         MatrixType *                           matrix,
                         VectorType *                           right_hand_side)
    {
      const unsigned int dofs_per_cell = data.dof_indices.size();
      (void)dofs_per_cell;

      Assert(data.cell_matrix.m() == dofs_per_cell, ExcInternalError());
      Assert(data.cell_matrix.n() == dofs_per_cell, ExcInternalError());
      Assert((right_hand_side == nullptr) ||
               (data.cell_rhs.size() == dofs_per_cell),
             ExcInternalError());

      if (right_hand_side != nullptr)
        data.constraints->distribute_local_to_global(data.cell_matrix,
                                                     data.cell_rhs,
                                                     data.dof_indices,
                                                     *matrix,
                                                     *right_hand_side);
      else
        data.constraints->distribute_local_to_global(data.cell_matrix,
                                                     data.dof_indices,
                                                     *matrix);
    }



    namespace AssemblerBoundary
    {
      struct Scratch
      {
        Scratch() = default;
      };

      template <int dim, int spacedim, typename number>
      struct CopyData
      {
        CopyData();

        unsigned int                                             dofs_per_cell;
        std::vector<types::global_dof_index>                     dofs;
        std::vector<std::vector<bool>>                           dof_is_on_face;
        typename DoFHandler<dim, spacedim>::active_cell_iterator cell;
        std::vector<FullMatrix<number>>                          cell_matrix;
        std::vector<Vector<number>>                              cell_vector;
      };


      template <int dim, int spacedim, typename number>
      CopyData<dim, spacedim, number>::CopyData()
        : dofs_per_cell(numbers::invalid_unsigned_int)
      {}


    } // namespace AssemblerBoundary
  }   // namespace internal
} // namespace MatrixCreator


namespace MatrixCreator
{
  template <int dim, int spacedim, typename number>
  void
  create_mass_matrix(const Mapping<dim, spacedim> &          mapping,
                     const DoFHandler<dim, spacedim> &       dof,
                     const Quadrature<dim> &                 q,
                     SparseMatrix<number> &                  matrix,
                     const Function<spacedim, number> *const coefficient,
                     const AffineConstraints<number> &       constraints)
  {
    Assert(matrix.m() == dof.n_dofs(),
           ExcDimensionMismatch(matrix.m(), dof.n_dofs()));
    Assert(matrix.n() == dof.n_dofs(),
           ExcDimensionMismatch(matrix.n(), dof.n_dofs()));

    hp::FECollection<dim, spacedim>      fe_collection(dof.get_fe());
    hp::QCollection<dim>                 q_collection(q);
    hp::MappingCollection<dim, spacedim> mapping_collection(mapping);
    MatrixCreator::internal::AssemblerData::Scratch<dim, spacedim, number>
      assembler_data(fe_collection,
                     update_values | update_JxW_values |
                       (coefficient != nullptr ? update_quadrature_points :
                                                 UpdateFlags(0)),
                     coefficient,
                     /*rhs_function=*/nullptr,
                     q_collection,
                     mapping_collection);

    MatrixCreator::internal::AssemblerData::CopyData<number> copy_data;
    copy_data.cell_matrix.reinit(
      assembler_data.fe_collection.max_dofs_per_cell(),
      assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.cell_rhs.reinit(assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.dof_indices.resize(
      assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.constraints = &constraints;

    WorkStream::run(
      dof.begin_active(),
      static_cast<typename DoFHandler<dim, spacedim>::active_cell_iterator>(
        dof.end()),
      &MatrixCreator::internal::mass_assembler<
        dim,
        spacedim,
        typename DoFHandler<dim, spacedim>::active_cell_iterator,
        number>,
      [&matrix](const internal::AssemblerData::CopyData<number> &data) {
        MatrixCreator::internal::copy_local_to_global(
          data, &matrix, static_cast<Vector<number> *>(nullptr));
      },
      assembler_data,
      copy_data);
  }



  template <int dim, int spacedim, typename number>
  void
  create_mass_matrix(const DoFHandler<dim, spacedim> &       dof,
                     const Quadrature<dim> &                 q,
                     SparseMatrix<number> &                  matrix,
                     const Function<spacedim, number> *const coefficient,
                     const AffineConstraints<number> &       constraints)
  {
    create_mass_matrix(get_default_linear_mapping(dof.get_triangulation()),
                       dof,
                       q,
                       matrix,
                       coefficient,
                       constraints);
  }



  template <int dim, int spacedim, typename number>
  void
  create_mass_matrix(const Mapping<dim, spacedim> &          mapping,
                     const DoFHandler<dim, spacedim> &       dof,
                     const Quadrature<dim> &                 q,
                     SparseMatrix<number> &                  matrix,
                     const Function<spacedim, number> &      rhs,
                     Vector<number> &                        rhs_vector,
                     const Function<spacedim, number> *const coefficient,
                     const AffineConstraints<number> &       constraints)
  {
    Assert(matrix.m() == dof.n_dofs(),
           ExcDimensionMismatch(matrix.m(), dof.n_dofs()));
    Assert(matrix.n() == dof.n_dofs(),
           ExcDimensionMismatch(matrix.n(), dof.n_dofs()));

    hp::FECollection<dim, spacedim>      fe_collection(dof.get_fe());
    hp::QCollection<dim>                 q_collection(q);
    hp::MappingCollection<dim, spacedim> mapping_collection(mapping);
    MatrixCreator::internal::AssemblerData::Scratch<dim, spacedim, number>
                                                             assembler_data(fe_collection,
                     update_values | update_JxW_values |
                       update_quadrature_points,
                     coefficient,
                     &rhs,
                     q_collection,
                     mapping_collection);
    MatrixCreator::internal::AssemblerData::CopyData<number> copy_data;
    copy_data.cell_matrix.reinit(
      assembler_data.fe_collection.max_dofs_per_cell(),
      assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.cell_rhs.reinit(assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.dof_indices.resize(
      assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.constraints = &constraints;

    WorkStream::run(
      dof.begin_active(),
      static_cast<typename DoFHandler<dim, spacedim>::active_cell_iterator>(
        dof.end()),
      &MatrixCreator::internal::mass_assembler<
        dim,
        spacedim,
        typename DoFHandler<dim, spacedim>::active_cell_iterator,
        number>,
      [&matrix,
       &rhs_vector](const internal::AssemblerData::CopyData<number> &data) {
        MatrixCreator::internal::copy_local_to_global(data,
                                                      &matrix,
                                                      &rhs_vector);
      },
      assembler_data,
      copy_data);
  }



  template <int dim, int spacedim, typename number>
  void
  create_mass_matrix(const DoFHandler<dim, spacedim> &       dof,
                     const Quadrature<dim> &                 q,
                     SparseMatrix<number> &                  matrix,
                     const Function<spacedim, number> &      rhs,
                     Vector<number> &                        rhs_vector,
                     const Function<spacedim, number> *const coefficient,
                     const AffineConstraints<number> &       constraints)
  {
    create_mass_matrix(get_default_linear_mapping(dof.get_triangulation()),
                       dof,
                       q,
                       matrix,
                       rhs,
                       rhs_vector,
                       coefficient,
                       constraints);
  }



  template <int dim, int spacedim, typename number>
  void
  create_mass_matrix(const hp::MappingCollection<dim, spacedim> &mapping,
                     const DoFHandler<dim, spacedim> &           dof,
                     const hp::QCollection<dim> &                q,
                     SparseMatrix<number> &                      matrix,
                     const Function<spacedim, number> *const     coefficient,
                     const AffineConstraints<number> &           constraints)
  {
    Assert(matrix.m() == dof.n_dofs(),
           ExcDimensionMismatch(matrix.m(), dof.n_dofs()));
    Assert(matrix.n() == dof.n_dofs(),
           ExcDimensionMismatch(matrix.n(), dof.n_dofs()));

    MatrixCreator::internal::AssemblerData::Scratch<dim, spacedim, number>
                                                             assembler_data(dof.get_fe_collection(),
                     update_values | update_JxW_values |
                       (coefficient != nullptr ? update_quadrature_points :
                                                 UpdateFlags(0)),
                     coefficient,
                     /*rhs_function=*/nullptr,
                     q,
                     mapping);
    MatrixCreator::internal::AssemblerData::CopyData<number> copy_data;
    copy_data.cell_matrix.reinit(
      assembler_data.fe_collection.max_dofs_per_cell(),
      assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.cell_rhs.reinit(assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.dof_indices.resize(
      assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.constraints = &constraints;

    WorkStream::run(
      dof.begin_active(),
      static_cast<typename DoFHandler<dim, spacedim>::active_cell_iterator>(
        dof.end()),
      &MatrixCreator::internal::mass_assembler<
        dim,
        spacedim,
        typename DoFHandler<dim, spacedim>::active_cell_iterator,
        number>,
      [&matrix](const internal::AssemblerData::CopyData<number> &data) {
        MatrixCreator::internal::copy_local_to_global(
          data, &matrix, static_cast<Vector<number> *>(nullptr));
      },
      assembler_data,
      copy_data);
  }



  template <int dim, int spacedim, typename number>
  void
  create_mass_matrix(const DoFHandler<dim, spacedim> &       dof,
                     const hp::QCollection<dim> &            q,
                     SparseMatrix<number> &                  matrix,
                     const Function<spacedim, number> *const coefficient,
                     const AffineConstraints<number> &       constraints)
  {
    create_mass_matrix(hp::StaticMappingQ1<dim, spacedim>::mapping_collection,
                       dof,
                       q,
                       matrix,
                       coefficient,
                       constraints);
  }



  template <int dim, int spacedim, typename number>
  void
  create_mass_matrix(const hp::MappingCollection<dim, spacedim> &mapping,
                     const DoFHandler<dim, spacedim> &           dof,
                     const hp::QCollection<dim> &                q,
                     SparseMatrix<number> &                      matrix,
                     const Function<spacedim, number> &          rhs,
                     Vector<number> &                            rhs_vector,
                     const Function<spacedim, number> *const     coefficient,
                     const AffineConstraints<number> &           constraints)
  {
    Assert(matrix.m() == dof.n_dofs(),
           ExcDimensionMismatch(matrix.m(), dof.n_dofs()));
    Assert(matrix.n() == dof.n_dofs(),
           ExcDimensionMismatch(matrix.n(), dof.n_dofs()));

    MatrixCreator::internal::AssemblerData::Scratch<dim, spacedim, number>
                                                             assembler_data(dof.get_fe_collection(),
                     update_values | update_JxW_values |
                       update_quadrature_points,
                     coefficient,
                     &rhs,
                     q,
                     mapping);
    MatrixCreator::internal::AssemblerData::CopyData<number> copy_data;
    copy_data.cell_matrix.reinit(
      assembler_data.fe_collection.max_dofs_per_cell(),
      assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.cell_rhs.reinit(assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.dof_indices.resize(
      assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.constraints = &constraints;

    WorkStream::run(
      dof.begin_active(),
      static_cast<typename DoFHandler<dim, spacedim>::active_cell_iterator>(
        dof.end()),
      &MatrixCreator::internal::mass_assembler<
        dim,
        spacedim,
        typename DoFHandler<dim, spacedim>::active_cell_iterator,
        number>,
      [&matrix,
       &rhs_vector](const internal::AssemblerData::CopyData<number> &data) {
        MatrixCreator::internal::copy_local_to_global(data,
                                                      &matrix,
                                                      &rhs_vector);
      },
      assembler_data,
      copy_data);
  }



  template <int dim, int spacedim, typename number>
  void
  create_mass_matrix(const DoFHandler<dim, spacedim> &       dof,
                     const hp::QCollection<dim> &            q,
                     SparseMatrix<number> &                  matrix,
                     const Function<spacedim, number> &      rhs,
                     Vector<number> &                        rhs_vector,
                     const Function<spacedim, number> *const coefficient,
                     const AffineConstraints<number> &       constraints)
  {
    create_mass_matrix(hp::StaticMappingQ1<dim, spacedim>::mapping_collection,
                       dof,
                       q,
                       matrix,
                       rhs,
                       rhs_vector,
                       coefficient,
                       constraints);
  }



  namespace internal
  {
    template <int dim, int spacedim, typename number>
    void static inline create_boundary_mass_matrix_1(
      typename DoFHandler<dim, spacedim>::active_cell_iterator const &cell,
      MatrixCreator::internal::AssemblerBoundary::Scratch const &,
      MatrixCreator::internal::AssemblerBoundary::
        CopyData<dim, spacedim, number> & copy_data,
      Mapping<dim, spacedim> const &      mapping,
      FiniteElement<dim, spacedim> const &fe,
      Quadrature<dim - 1> const &         q,
      std::map<types::boundary_id, const Function<spacedim, number> *> const
        &                                     boundary_functions,
      Function<spacedim, number> const *const coefficient,
      std::vector<unsigned int> const &       component_mapping)

    {
      // Most assertions for this function are in the calling function
      // before creating threads.
      const unsigned int n_components = fe.n_components();
      const unsigned int n_function_components =
        boundary_functions.begin()->second->n_components;
      const bool fe_is_system    = (n_components != 1);
      const bool fe_is_primitive = fe.is_primitive();

      copy_data.cell          = cell;
      copy_data.dofs_per_cell = fe.n_dofs_per_cell();

      UpdateFlags update_flags =
        UpdateFlags(update_values | update_JxW_values | update_normal_vectors |
                    update_quadrature_points);
      FEFaceValues<dim, spacedim> fe_values(mapping, fe, q, update_flags);

      // two variables for the coefficient, one for the two cases
      // indicated in the name
      std::vector<number> coefficient_values(fe_values.n_quadrature_points, 1.);
      std::vector<Vector<number>> coefficient_vector_values(
        fe_values.n_quadrature_points, Vector<number>(n_components));
      const bool coefficient_is_vector =
        (coefficient != nullptr && coefficient->n_components != 1);

      std::vector<number> rhs_values_scalar(fe_values.n_quadrature_points);
      std::vector<Vector<number>> rhs_values_system(
        fe_values.n_quadrature_points, Vector<number>(n_function_components));

      copy_data.dofs.resize(copy_data.dofs_per_cell);
      cell->get_dof_indices(copy_data.dofs);

      std::vector<types::global_dof_index> dofs_on_face_vector;

      // Because CopyData objects are reused and emplace_back is
      // used, dof_is_on_face, cell_matrix, and cell_vector must be
      // cleared before they are reused
      copy_data.dof_is_on_face.clear();
      copy_data.cell_matrix.clear();
      copy_data.cell_vector.clear();

      for (const unsigned int face : cell->face_indices())
        // check if this face is on that part of the boundary we are
        // interested in
        if (boundary_functions.find(cell->face(face)->boundary_id()) !=
            boundary_functions.end())
          {
            copy_data.cell_matrix.emplace_back(copy_data.dofs_per_cell,
                                               copy_data.dofs_per_cell);
            copy_data.cell_vector.emplace_back(copy_data.dofs_per_cell);
            fe_values.reinit(cell, face);

            if (fe_is_system)
              // FE has several components
              {
                boundary_functions.find(cell->face(face)->boundary_id())
                  ->second->vector_value_list(fe_values.get_quadrature_points(),
                                              rhs_values_system);

                if (coefficient_is_vector)
                  // If coefficient is vector valued, fill all
                  // components
                  coefficient->vector_value_list(
                    fe_values.get_quadrature_points(),
                    coefficient_vector_values);
                else
                  {
                    // If a scalar function is given, update the
                    // values, if not, use the default one set in the
                    // constructor above
                    if (coefficient != nullptr)
                      coefficient->value_list(fe_values.get_quadrature_points(),
                                              coefficient_values);
                    // Copy scalar values into vector
                    for (unsigned int point = 0;
                         point < fe_values.n_quadrature_points;
                         ++point)
                      coefficient_vector_values[point] =
                        coefficient_values[point];
                  }

                // Special treatment for Hdiv elements,
                // where only normal components should be projected.
                // For Hdiv we need to compute (u dot n, v dot n) which
                // can be done as sum over dim components of
                // u[c] * n[c] * v[c] * n[c] = u[c] * v[c] *
                // normal_adjustment[c] Same approach does not work for Hcurl,
                // so we throw an exception. Default value 1.0 allows for use
                // with non Hdiv elements
                std::vector<std::vector<double>> normal_adjustment(
                  fe_values.n_quadrature_points,
                  std::vector<double>(n_components, 1.));

                for (unsigned int comp = 0; comp < n_components; ++comp)
                  {
                    const FiniteElement<dim, spacedim> &base =
                      fe.base_element(fe.component_to_base_index(comp).first);
                    const unsigned int bcomp =
                      fe.component_to_base_index(comp).second;

                    if (!base.conforms(FiniteElementData<dim>::H1) &&
                        base.conforms(FiniteElementData<dim>::Hdiv) &&
                        fe_is_primitive)
                      Assert(false, ExcNotImplemented());

                    if (!base.conforms(FiniteElementData<dim>::H1) &&
                        base.conforms(FiniteElementData<dim>::Hcurl))
                      Assert(false, ExcNotImplemented());

                    if (!base.conforms(FiniteElementData<dim>::H1) &&
                        base.conforms(FiniteElementData<dim>::Hdiv))
                      for (unsigned int point = 0;
                           point < fe_values.n_quadrature_points;
                           ++point)
                        normal_adjustment[point][comp] =
                          fe_values.normal_vector(point)[bcomp] *
                          fe_values.normal_vector(point)[bcomp];
                  }

                for (unsigned int point = 0;
                     point < fe_values.n_quadrature_points;
                     ++point)
                  {
                    const double weight = fe_values.JxW(point);
                    for (unsigned int i = 0; i < fe_values.dofs_per_cell; ++i)
                      if (fe_is_primitive)
                        {
                          for (unsigned int j = 0; j < fe_values.dofs_per_cell;
                               ++j)
                            {
                              if (fe.system_to_component_index(j).first ==
                                  fe.system_to_component_index(i).first)
                                {
                                  copy_data.cell_matrix.back()(i, j) +=
                                    coefficient_vector_values[point](
                                      fe.system_to_component_index(i).first) *
                                    weight * fe_values.shape_value(j, point) *
                                    fe_values.shape_value(i, point);
                                }
                            }
                          copy_data.cell_vector.back()(i) +=
                            rhs_values_system[point](
                              component_mapping[fe.system_to_component_index(i)
                                                  .first]) *
                            fe_values.shape_value(i, point) * weight;
                        }
                      else
                        {
                          for (unsigned int comp = 0; comp < n_components;
                               ++comp)
                            {
                              for (unsigned int j = 0;
                                   j < fe_values.dofs_per_cell;
                                   ++j)
                                copy_data.cell_matrix.back()(i, j) +=
                                  coefficient_vector_values[point](comp) *
                                  fe_values.shape_value_component(j,
                                                                  point,
                                                                  comp) *
                                  fe_values.shape_value_component(i,
                                                                  point,
                                                                  comp) *
                                  normal_adjustment[point][comp] * weight;
                              copy_data.cell_vector.back()(i) +=
                                rhs_values_system[point](
                                  component_mapping[comp]) *
                                fe_values.shape_value_component(i,
                                                                point,
                                                                comp) *
                                normal_adjustment[point][comp] * weight;
                            }
                        }
                  }
              }
            else
              // FE is a scalar one
              {
                boundary_functions.find(cell->face(face)->boundary_id())
                  ->second->value_list(fe_values.get_quadrature_points(),
                                       rhs_values_scalar);

                if (coefficient != nullptr)
                  coefficient->value_list(fe_values.get_quadrature_points(),
                                          coefficient_values);
                for (unsigned int point = 0;
                     point < fe_values.n_quadrature_points;
                     ++point)
                  {
                    const double weight = fe_values.JxW(point);
                    for (unsigned int i = 0; i < fe_values.dofs_per_cell; ++i)
                      {
                        const double v = fe_values.shape_value(i, point);
                        for (unsigned int j = 0; j < fe_values.dofs_per_cell;
                             ++j)
                          {
                            const double u = fe_values.shape_value(j, point);
                            copy_data.cell_matrix.back()(i, j) +=
                              (coefficient_values[point] * u * v * weight);
                          }
                        copy_data.cell_vector.back()(i) +=
                          rhs_values_scalar[point] * v * weight;
                      }
                  }
              }


            dofs_on_face_vector.resize(fe.n_dofs_per_face(face));
            cell->face(face)->get_dof_indices(dofs_on_face_vector);
            // for each dof on the cell, have a flag whether it is on
            // the face
            copy_data.dof_is_on_face.emplace_back(copy_data.dofs_per_cell);
            // check for each of the dofs on this cell whether it is
            // on the face
            for (unsigned int i = 0; i < copy_data.dofs_per_cell; ++i)
              copy_data.dof_is_on_face.back()[i] =
                (std::find(dofs_on_face_vector.begin(),
                           dofs_on_face_vector.end(),
                           copy_data.dofs[i]) != dofs_on_face_vector.end());
          }
    }

    template <int dim, int spacedim, typename number>
    void
    copy_boundary_mass_matrix_1(
      MatrixCreator::internal::AssemblerBoundary::
        CopyData<dim, spacedim, number> const &copy_data,
      std::map<types::boundary_id, const Function<spacedim, number> *> const
        &                                         boundary_functions,
      std::vector<types::global_dof_index> const &dof_to_boundary_mapping,
      SparseMatrix<number> &                      matrix,
      Vector<number> &                            rhs_vector)
    {
      // now transfer cell matrix and vector to the whole boundary matrix
      //
      // in the following: dof[i] holds the global index of the i-th degree of
      // freedom on the present cell. If it is also a dof on the boundary, it
      // must be a nonzero entry in the dof_to_boundary_mapping and then
      // the boundary index of this dof is dof_to_boundary_mapping[dof[i]].
      //
      // if dof[i] is not on the boundary, it should be zero on the boundary
      // therefore on all quadrature points and finally all of its
      // entries in the cell matrix and vector should be zero. If not, we
      // throw an error (note: because of the evaluation of the shape
      // functions only up to machine precision, the term "must be zero"
      // really should mean: "should be very small". since this is only an
      // assertion and not part of the code, we may choose "very small"
      // quite arbitrarily)
      //
      // the main problem here is that the matrix or vector entry should also
      // be zero if the degree of freedom dof[i] is on the boundary, but not
      // on the present face, i.e. on another face of the same cell also
      // on the boundary. We can therefore not rely on the
      // dof_to_boundary_mapping[dof[i]] being !=-1, we really have to
      // determine whether dof[i] is a dof on the present face. We do so
      // by getting the dofs on the face into @p{dofs_on_face_vector},
      // a vector as always. Usually, searching in a vector is
      // inefficient, so we copy the dofs into a set, which enables binary
      // searches.
      unsigned int pos(0);
      for (const unsigned int face : copy_data.cell->face_indices())
        {
          // check if this face is on that part of
          // the boundary we are interested in
          if (boundary_functions.find(
                copy_data.cell->face(face)->boundary_id()) !=
              boundary_functions.end())
            {
              for (unsigned int i = 0; i < copy_data.dofs_per_cell; ++i)
                {
                  if (copy_data.dof_is_on_face[pos][i] &&
                      dof_to_boundary_mapping[copy_data.dofs[i]] !=
                        numbers::invalid_dof_index)
                    {
                      for (unsigned int j = 0; j < copy_data.dofs_per_cell; ++j)
                        if (copy_data.dof_is_on_face[pos][j] &&
                            dof_to_boundary_mapping[copy_data.dofs[j]] !=
                              numbers::invalid_dof_index)
                          {
                            AssertIsFinite(copy_data.cell_matrix[pos](i, j));
                            matrix.add(
                              dof_to_boundary_mapping[copy_data.dofs[i]],
                              dof_to_boundary_mapping[copy_data.dofs[j]],
                              copy_data.cell_matrix[pos](i, j));
                          }
                      AssertIsFinite(copy_data.cell_vector[pos](i));
                      rhs_vector(dof_to_boundary_mapping[copy_data.dofs[i]]) +=
                        copy_data.cell_vector[pos](i);
                    }
                }
              ++pos;
            }
        }
    }


    template <>
    void inline create_boundary_mass_matrix_1<1, 3, float>(
      DoFHandler<1, 3>::active_cell_iterator const & /*cell*/,
      MatrixCreator::internal::AssemblerBoundary::Scratch const &,
      MatrixCreator::internal::AssemblerBoundary::CopyData<1, 3, float>
        & /*copy_data*/,
      Mapping<1, 3> const &,
      FiniteElement<1, 3> const &,
      Quadrature<0> const &,
      std::map<types::boundary_id, const Function<3, float> *> const
        & /*boundary_functions*/,
      Function<3, float> const *const /*coefficient*/,
      std::vector<unsigned int> const & /*component_mapping*/)
    {
      Assert(false, ExcNotImplemented());
    }

    template <>
    void inline create_boundary_mass_matrix_1<1, 3, double>(
      DoFHandler<1, 3>::active_cell_iterator const & /*cell*/,
      MatrixCreator::internal::AssemblerBoundary::Scratch const &,
      MatrixCreator::internal::AssemblerBoundary::CopyData<1, 3, double>
        & /*copy_data*/,
      Mapping<1, 3> const &,
      FiniteElement<1, 3> const &,
      Quadrature<0> const &,
      std::map<types::boundary_id, const Function<3, double> *> const
        & /*boundary_functions*/,
      Function<3, double> const *const /*coefficient*/,
      std::vector<unsigned int> const & /*component_mapping*/)
    {
      Assert(false, ExcNotImplemented());
    }

  } // namespace internal



  template <int dim, int spacedim, typename number>
  void
  create_boundary_mass_matrix(
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof,
    const Quadrature<dim - 1> &      q,
    SparseMatrix<number> &           matrix,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                     boundary_functions,
    Vector<number> &                        rhs_vector,
    std::vector<types::global_dof_index> &  dof_to_boundary_mapping,
    const Function<spacedim, number> *const coefficient,
    std::vector<unsigned int>               component_mapping)
  {
    // what would that be in 1d? the
    // identity matrix on the boundary
    // dofs?
    if (dim == 1)
      {
        Assert(false, ExcNotImplemented());
        return;
      }

    const FiniteElement<dim, spacedim> &fe           = dof.get_fe();
    const unsigned int                  n_components = fe.n_components();

    Assert(matrix.n() == dof.n_boundary_dofs(boundary_functions),
           ExcInternalError());
    Assert(matrix.n() == matrix.m(), ExcInternalError());
    Assert(matrix.n() == rhs_vector.size(), ExcInternalError());
    Assert(boundary_functions.size() != 0, ExcInternalError());
    Assert(dof_to_boundary_mapping.size() == dof.n_dofs(), ExcInternalError());
    Assert(coefficient == nullptr || coefficient->n_components == 1 ||
             coefficient->n_components == n_components,
           ExcComponentMismatch());

    if (component_mapping.size() == 0)
      {
        AssertDimension(n_components,
                        boundary_functions.begin()->second->n_components);
        for (unsigned int i = 0; i < n_components; ++i)
          component_mapping.push_back(i);
      }
    else
      AssertDimension(n_components, component_mapping.size());

    MatrixCreator::internal::AssemblerBoundary::Scratch scratch;
    MatrixCreator::internal::AssemblerBoundary::CopyData<dim, spacedim, number>
      copy_data;

    WorkStream::run(
      dof.begin_active(),
      dof.end(),
      [&mapping, &fe, &q, &boundary_functions, coefficient, &component_mapping](
        typename DoFHandler<dim, spacedim>::active_cell_iterator const &cell,
        MatrixCreator::internal::AssemblerBoundary::Scratch const &scratch_data,
        MatrixCreator::internal::AssemblerBoundary::
          CopyData<dim, spacedim, number> &copy_data) {
        internal::create_boundary_mass_matrix_1(cell,
                                                scratch_data,
                                                copy_data,
                                                mapping,
                                                fe,
                                                q,
                                                boundary_functions,
                                                coefficient,
                                                component_mapping);
      },
      [&boundary_functions, &dof_to_boundary_mapping, &matrix, &rhs_vector](
        MatrixCreator::internal::AssemblerBoundary::
          CopyData<dim, spacedim, number> const &copy_data) {
        internal::copy_boundary_mass_matrix_1(copy_data,
                                              boundary_functions,
                                              dof_to_boundary_mapping,
                                              matrix,
                                              rhs_vector);
      },
      scratch,
      copy_data);
  }



  namespace internal
  {
    template <int dim, int spacedim, typename number>
    void
    create_hp_boundary_mass_matrix_1(
      typename DoFHandler<dim, spacedim>::active_cell_iterator const &cell,
      MatrixCreator::internal::AssemblerBoundary::Scratch const &,
      MatrixCreator::internal::AssemblerBoundary ::
        CopyData<dim, spacedim, number> &         copy_data,
      hp::MappingCollection<dim, spacedim> const &mapping,
      hp::FECollection<dim, spacedim> const &     fe_collection,
      hp::QCollection<dim - 1> const &            q,
      const std::map<types::boundary_id, const Function<spacedim, number> *>
        &                                     boundary_functions,
      Function<spacedim, number> const *const coefficient,
      std::vector<unsigned int> const &       component_mapping)
    {
      const unsigned int n_components = fe_collection.n_components();
      const unsigned int n_function_components =
        boundary_functions.begin()->second->n_components;
      const FiniteElement<dim, spacedim> &fe              = cell->get_fe();
      const bool                          fe_is_system    = (n_components != 1);
      const bool                          fe_is_primitive = fe.is_primitive();

      copy_data.cell          = cell;
      copy_data.dofs_per_cell = fe.n_dofs_per_cell();
      copy_data.dofs.resize(copy_data.dofs_per_cell);
      cell->get_dof_indices(copy_data.dofs);


      UpdateFlags update_flags =
        UpdateFlags(update_values | update_JxW_values | update_normal_vectors |
                    update_quadrature_points);
      hp::FEFaceValues<dim, spacedim> x_fe_values(mapping,
                                                  fe_collection,
                                                  q,
                                                  update_flags);

      // two variables for the coefficient,
      // one for the two cases indicated in
      // the name
      std::vector<number>         coefficient_values;
      std::vector<Vector<number>> coefficient_vector_values;

      const bool coefficient_is_vector =
        (coefficient != nullptr && coefficient->n_components != 1);

      std::vector<number>         rhs_values_scalar;
      std::vector<Vector<number>> rhs_values_system;

      std::vector<types::global_dof_index> dofs_on_face_vector;

      copy_data.dofs.resize(copy_data.dofs_per_cell);
      cell->get_dof_indices(copy_data.dofs);

      // Because CopyData objects are reused and that push_back is used,
      // dof_is_on_face, cell_matrix, and cell_vector must be cleared before
      // they are reused
      copy_data.dof_is_on_face.clear();
      copy_data.cell_matrix.clear();
      copy_data.cell_vector.clear();


      for (const unsigned int face : cell->face_indices())
        // check if this face is on that part of
        // the boundary we are interested in
        if (boundary_functions.find(cell->face(face)->boundary_id()) !=
            boundary_functions.end())
          {
            x_fe_values.reinit(cell, face);

            const FEFaceValues<dim, spacedim> &fe_values =
              x_fe_values.get_present_fe_values();

            copy_data.cell_matrix.emplace_back(copy_data.dofs_per_cell,
                                               copy_data.dofs_per_cell);
            copy_data.cell_vector.emplace_back(copy_data.dofs_per_cell);

            if (fe_is_system)
              // FE has several components
              {
                coefficient_vector_values.resize(fe_values.n_quadrature_points,
                                                 Vector<number>(n_components));
                coefficient_values.resize(fe_values.n_quadrature_points, 1.);

                rhs_values_system.resize(fe_values.n_quadrature_points,
                                         Vector<number>(n_function_components));
                boundary_functions.find(cell->face(face)->boundary_id())
                  ->second->vector_value_list(fe_values.get_quadrature_points(),
                                              rhs_values_system);
                if (coefficient_is_vector)
                  // In case coefficient is vector-valued, fill
                  // all components
                  coefficient->vector_value_list(
                    fe_values.get_quadrature_points(),
                    coefficient_vector_values);
                else
                  // In case the scalar function is given, update the
                  // values, if not - use the default (1.0)
                  {
                    if (coefficient != nullptr)
                      coefficient->value_list(fe_values.get_quadrature_points(),
                                              coefficient_values);

                    for (unsigned int point = 0;
                         point < fe_values.n_quadrature_points;
                         ++point)
                      coefficient_vector_values[point] =
                        coefficient_values[point];
                  }

                // Special treatment for Hdiv elements,
                // where only normal components should be projected.
                // For Hdiv we need to compute (u dot n, v dot n) which
                // can be done as sum over dim components of
                // u[c] * n[c] * v[c] * n[c] = u[c] * v[c] *
                // normal_adjustment[c] Same approach does not work for Hcurl,
                // so we throw an exception. Default value 1.0 allows for use
                // with non Hdiv elements
                std::vector<std::vector<double>> normal_adjustment(
                  fe_values.n_quadrature_points,
                  std::vector<double>(n_components, 1.));

                for (unsigned int comp = 0; comp < n_components; ++comp)
                  {
                    const FiniteElement<dim, spacedim> &base =
                      fe.base_element(fe.component_to_base_index(comp).first);
                    const unsigned int bcomp =
                      fe.component_to_base_index(comp).second;

                    if (!base.conforms(FiniteElementData<dim>::H1) &&
                        base.conforms(FiniteElementData<dim>::Hdiv) &&
                        fe_is_primitive)
                      Assert(false, ExcNotImplemented());

                    if (!base.conforms(FiniteElementData<dim>::H1) &&
                        base.conforms(FiniteElementData<dim>::Hcurl))
                      Assert(false, ExcNotImplemented());

                    if (!base.conforms(FiniteElementData<dim>::H1) &&
                        base.conforms(FiniteElementData<dim>::Hdiv))
                      for (unsigned int point = 0;
                           point < fe_values.n_quadrature_points;
                           ++point)
                        normal_adjustment[point][comp] =
                          fe_values.normal_vector(point)[bcomp] *
                          fe_values.normal_vector(point)[bcomp];
                  }

                for (unsigned int point = 0;
                     point < fe_values.n_quadrature_points;
                     ++point)
                  {
                    const double weight = fe_values.JxW(point);
                    for (unsigned int i = 0; i < fe_values.dofs_per_cell; ++i)
                      if (fe_is_primitive)
                        {
                          for (unsigned int j = 0; j < fe_values.dofs_per_cell;
                               ++j)
                            {
                              if (fe.system_to_component_index(i).first ==
                                  fe.system_to_component_index(j).first)
                                {
                                  copy_data.cell_matrix.back()(i, j) +=
                                    coefficient_vector_values[point](
                                      fe.system_to_component_index(i).first) *
                                    fe_values.shape_value(i, point) *
                                    fe_values.shape_value(j, point) * weight;
                                }
                            }

                          copy_data.cell_vector.back()(i) +=
                            rhs_values_system[point](
                              component_mapping[fe.system_to_component_index(i)
                                                  .first]) *
                            fe_values.shape_value(i, point) * weight;
                        }
                      else
                        {
                          for (unsigned int comp = 0; comp < n_components;
                               ++comp)
                            {
                              for (unsigned int j = 0;
                                   j < fe_values.dofs_per_cell;
                                   ++j)
                                copy_data.cell_matrix.back()(i, j) +=
                                  coefficient_vector_values[point](comp) *
                                  fe_values.shape_value_component(i,
                                                                  point,
                                                                  comp) *
                                  fe_values.shape_value_component(j,
                                                                  point,
                                                                  comp) *
                                  normal_adjustment[point][comp] * weight;
                              copy_data.cell_vector.back()(i) +=
                                rhs_values_system[point](
                                  component_mapping[comp]) *
                                fe_values.shape_value_component(i,
                                                                point,
                                                                comp) *
                                normal_adjustment[point][comp] * weight;
                            }
                        }
                  }
              }
            else
              // FE is scalar
              {
                coefficient_values.resize(fe_values.n_quadrature_points, 1.);
                rhs_values_scalar.resize(fe_values.n_quadrature_points);
                boundary_functions.find(cell->face(face)->boundary_id())
                  ->second->value_list(fe_values.get_quadrature_points(),
                                       rhs_values_scalar);

                if (coefficient != nullptr)
                  {
                    coefficient_values.resize(fe_values.n_quadrature_points);
                    coefficient->value_list(fe_values.get_quadrature_points(),
                                            coefficient_values);
                  }

                for (unsigned int point = 0;
                     point < fe_values.n_quadrature_points;
                     ++point)
                  {
                    const double weight = fe_values.JxW(point);
                    for (unsigned int i = 0; i < fe_values.dofs_per_cell; ++i)
                      {
                        const double v = fe_values.shape_value(i, point);
                        for (unsigned int j = 0; j < fe_values.dofs_per_cell;
                             ++j)
                          {
                            const double u = fe_values.shape_value(j, point);
                            copy_data.cell_matrix.back()(i, j) +=
                              (coefficient_values[point] * v * u * weight);
                          }
                        copy_data.cell_vector.back()(i) +=
                          rhs_values_scalar[point] * v * weight;
                      }
                  }
              }

            dofs_on_face_vector.resize(fe.n_dofs_per_face(face));
            cell->face(face)->get_dof_indices(dofs_on_face_vector,
                                              cell->active_fe_index());
            // for each dof on the cell, have a
            // flag whether it is on the face
            copy_data.dof_is_on_face.emplace_back(copy_data.dofs_per_cell);
            // check for each of the dofs on this cell
            // whether it is on the face
            for (unsigned int i = 0; i < copy_data.dofs_per_cell; ++i)
              copy_data.dof_is_on_face.back()[i] =
                (std::find(dofs_on_face_vector.begin(),
                           dofs_on_face_vector.end(),
                           copy_data.dofs[i]) != dofs_on_face_vector.end());
          }
    }


    template <int dim, int spacedim, typename number>
    void
    copy_hp_boundary_mass_matrix_1(
      MatrixCreator::internal::AssemblerBoundary ::
        CopyData<dim, spacedim, number> const &copy_data,
      std::map<types::boundary_id, const Function<spacedim, number> *> const
        &                                         boundary_functions,
      std::vector<types::global_dof_index> const &dof_to_boundary_mapping,
      SparseMatrix<number> &                      matrix,
      Vector<number> &                            rhs_vector)
    {
      // now transfer cell matrix and vector to the whole boundary matrix
      //
      // in the following: dof[i] holds the  global index of the i-th degree of
      // freedom on the present cell. If it is also a dof on the boundary, it
      // must be a nonzero entry in the dof_to_boundary_mapping and then
      // the boundary index of this dof is dof_to_boundary_mapping[dof[i]].
      //
      // if dof[i] is not on the boundary, it should be zero on the boundary
      // therefore on all quadrature points and finally all of its
      // entries in the cell matrix and vector should be zero. If not, we
      // throw an error (note: because of the evaluation of the shape
      // functions only up to machine precision, the term "must be zero"
      // really should mean: "should be very small". since this is only an
      // assertion and not part of the code, we may choose "very small"
      // quite arbitrarily)
      //
      // the main problem here is that the matrix or vector entry should also
      // be zero if the degree of freedom dof[i] is on the boundary, but not
      // on the present face, i.e. on another face of the same cell also
      // on the boundary. We can therefore not rely on the
      // dof_to_boundary_mapping[dof[i]] being !=-1, we really have to
      // determine whether dof[i] is a dof on the present face. We do so
      // by getting the dofs on the face into @p{dofs_on_face_vector},
      // a vector as always. Usually, searching in a vector is
      // inefficient, so we copy the dofs into a set, which enables binary
      // searches.
      unsigned int pos(0);
      for (const unsigned int face : copy_data.cell->face_indices())
        {
          // check if this face is on that part of
          // the boundary we are interested in
          if (boundary_functions.find(
                copy_data.cell->face(face)->boundary_id()) !=
              boundary_functions.end())
            {
#ifdef DEBUG
              // in debug mode: compute an element in the matrix which is
              // guaranteed to belong to a boundary dof. We do this to check
              // that the entries in the cell matrix are guaranteed to be zero
              // if the respective dof is not on the boundary. Since because of
              // round-off, the actual value of the matrix entry may be
              // only close to zero, we assert that it is small relative to an
              // element which is guaranteed to be nonzero. (absolute smallness
              // does not suffice since the size of the domain scales in here)
              //
              // for this purpose we seek the diagonal of the matrix, where
              // there must be an element belonging to the boundary. we take the
              // maximum diagonal entry.
              types::global_dof_index max_element = 0;
              for (const auto index : dof_to_boundary_mapping)
                if ((index != numbers::invalid_dof_index) &&
                    (index > max_element))
                  max_element = index;
              Assert(max_element == matrix.n() - 1, ExcInternalError());

              double max_diag_entry = 0;
              for (unsigned int i = 0; i < copy_data.dofs_per_cell; ++i)
                if (std::abs(copy_data.cell_matrix[pos](i, i)) > max_diag_entry)
                  max_diag_entry = std::abs(copy_data.cell_matrix[pos](i, i));
#endif

              for (unsigned int i = 0; i < copy_data.dofs_per_cell; ++i)
                {
                  if (copy_data.dof_is_on_face[pos][i] &&
                      dof_to_boundary_mapping[copy_data.dofs[i]] !=
                        numbers::invalid_dof_index)
                    {
                      for (unsigned int j = 0; j < copy_data.dofs_per_cell; ++j)
                        if (copy_data.dof_is_on_face[pos][j] &&
                            dof_to_boundary_mapping[copy_data.dofs[j]] !=
                              numbers::invalid_dof_index)
                          {
                            AssertIsFinite(copy_data.cell_matrix[pos](i, j));
                            matrix.add(
                              dof_to_boundary_mapping[copy_data.dofs[i]],
                              dof_to_boundary_mapping[copy_data.dofs[j]],
                              copy_data.cell_matrix[pos](i, j));
                          }
                      AssertIsFinite(copy_data.cell_vector[pos](i));
                      rhs_vector(dof_to_boundary_mapping[copy_data.dofs[i]]) +=
                        copy_data.cell_vector[pos](i);
                    }
                }
              ++pos;
            }
        }
    }
  } // namespace internal



  template <int dim, int spacedim, typename number>
  void
  create_boundary_mass_matrix(
    const DoFHandler<dim, spacedim> &dof,
    const Quadrature<dim - 1> &      q,
    SparseMatrix<number> &           matrix,
    const std::map<types::boundary_id, const Function<spacedim, number> *> &rhs,
    Vector<number> &                        rhs_vector,
    std::vector<types::global_dof_index> &  dof_to_boundary_mapping,
    const Function<spacedim, number> *const a,
    std::vector<unsigned int>               component_mapping)
  {
    create_boundary_mass_matrix(get_default_linear_mapping(
                                  dof.get_triangulation()),
                                dof,
                                q,
                                matrix,
                                rhs,
                                rhs_vector,
                                dof_to_boundary_mapping,
                                a,
                                component_mapping);
  }



  template <int dim, int spacedim, typename number>
  void
  create_boundary_mass_matrix(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim> &           dof,
    const hp::QCollection<dim - 1> &            q,
    SparseMatrix<number> &                      matrix,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                     boundary_functions,
    Vector<number> &                        rhs_vector,
    std::vector<types::global_dof_index> &  dof_to_boundary_mapping,
    const Function<spacedim, number> *const coefficient,
    std::vector<unsigned int>               component_mapping)
  {
    // what would that be in 1d? the
    // identity matrix on the boundary
    // dofs?
    if (dim == 1)
      {
        Assert(false, ExcNotImplemented());
        return;
      }

    const hp::FECollection<dim, spacedim> &fe_collection =
      dof.get_fe_collection();
    const unsigned int n_components = fe_collection.n_components();

    Assert(matrix.n() == dof.n_boundary_dofs(boundary_functions),
           ExcInternalError());
    Assert(matrix.n() == matrix.m(), ExcInternalError());
    Assert(matrix.n() == rhs_vector.size(), ExcInternalError());
    Assert(boundary_functions.size() != 0, ExcInternalError());
    Assert(dof_to_boundary_mapping.size() == dof.n_dofs(), ExcInternalError());
    Assert(coefficient == nullptr || coefficient->n_components == 1 ||
             coefficient->n_components == n_components,
           ExcComponentMismatch());

    if (component_mapping.size() == 0)
      {
        AssertDimension(n_components,
                        boundary_functions.begin()->second->n_components);
        for (unsigned int i = 0; i < n_components; ++i)
          component_mapping.push_back(i);
      }
    else
      AssertDimension(n_components, component_mapping.size());

    MatrixCreator::internal::AssemblerBoundary::Scratch scratch;
    MatrixCreator::internal::AssemblerBoundary::CopyData<dim, spacedim, number>
      copy_data;

    WorkStream::run(
      dof.begin_active(),
      dof.end(),
      [&mapping,
       &fe_collection,
       &q,
       &boundary_functions,
       coefficient,
       &component_mapping](
        typename DoFHandler<dim, spacedim>::active_cell_iterator const &cell,
        MatrixCreator::internal::AssemblerBoundary::Scratch const &scratch_data,
        MatrixCreator::internal::AssemblerBoundary ::
          CopyData<dim, spacedim, number> &copy_data) {
        internal::create_hp_boundary_mass_matrix_1(cell,
                                                   scratch_data,
                                                   copy_data,
                                                   mapping,
                                                   fe_collection,
                                                   q,
                                                   boundary_functions,
                                                   coefficient,
                                                   component_mapping);
      },
      [&boundary_functions, &dof_to_boundary_mapping, &matrix, &rhs_vector](
        MatrixCreator::internal::AssemblerBoundary ::
          CopyData<dim, spacedim, number> const &copy_data) {
        internal::copy_hp_boundary_mass_matrix_1(copy_data,
                                                 boundary_functions,
                                                 dof_to_boundary_mapping,
                                                 matrix,
                                                 rhs_vector);
      },
      scratch,
      copy_data);
  }



  template <int dim, int spacedim, typename number>
  void
  create_boundary_mass_matrix(
    const DoFHandler<dim, spacedim> &dof,
    const hp::QCollection<dim - 1> & q,
    SparseMatrix<number> &           matrix,
    const std::map<types::boundary_id, const Function<spacedim, number> *> &rhs,
    Vector<number> &                        rhs_vector,
    std::vector<types::global_dof_index> &  dof_to_boundary_mapping,
    const Function<spacedim, number> *const a,
    std::vector<unsigned int>               component_mapping)
  {
    create_boundary_mass_matrix(
      hp::StaticMappingQ1<dim, spacedim>::mapping_collection,
      dof,
      q,
      matrix,
      rhs,
      rhs_vector,
      dof_to_boundary_mapping,
      a,
      component_mapping);
  }



  template <int dim, int spacedim>
  void
  create_laplace_matrix(const Mapping<dim, spacedim> &   mapping,
                        const DoFHandler<dim, spacedim> &dof,
                        const Quadrature<dim> &          q,
                        SparseMatrix<double> &           matrix,
                        const Function<spacedim> *const  coefficient,
                        const AffineConstraints<double> &constraints)
  {
    Assert(matrix.m() == dof.n_dofs(),
           ExcDimensionMismatch(matrix.m(), dof.n_dofs()));
    Assert(matrix.n() == dof.n_dofs(),
           ExcDimensionMismatch(matrix.n(), dof.n_dofs()));

    hp::FECollection<dim, spacedim>      fe_collection(dof.get_fe());
    hp::QCollection<dim>                 q_collection(q);
    hp::MappingCollection<dim, spacedim> mapping_collection(mapping);
    MatrixCreator::internal::AssemblerData::Scratch<dim, spacedim, double>
                                                             assembler_data(fe_collection,
                     update_gradients | update_JxW_values |
                       (coefficient != nullptr ? update_quadrature_points :
                                                 UpdateFlags(0)),
                     coefficient,
                     /*rhs_function=*/nullptr,
                     q_collection,
                     mapping_collection);
    MatrixCreator::internal::AssemblerData::CopyData<double> copy_data;
    copy_data.cell_matrix.reinit(
      assembler_data.fe_collection.max_dofs_per_cell(),
      assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.cell_rhs.reinit(assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.dof_indices.resize(
      assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.constraints = &constraints;

    WorkStream::run(
      dof.begin_active(),
      static_cast<typename DoFHandler<dim, spacedim>::active_cell_iterator>(
        dof.end()),
      &MatrixCreator::internal::laplace_assembler<
        dim,
        spacedim,
        typename DoFHandler<dim, spacedim>::active_cell_iterator>,
      [&matrix](const internal::AssemblerData::CopyData<double> &data) {
        MatrixCreator::internal::copy_local_to_global(
          data, &matrix, static_cast<Vector<double> *>(nullptr));
      },
      assembler_data,
      copy_data);
  }



  template <int dim, int spacedim>
  void
  create_laplace_matrix(const DoFHandler<dim, spacedim> &dof,
                        const Quadrature<dim> &          q,
                        SparseMatrix<double> &           matrix,
                        const Function<spacedim> *const  coefficient,
                        const AffineConstraints<double> &constraints)
  {
    create_laplace_matrix(get_default_linear_mapping(dof.get_triangulation()),
                          dof,
                          q,
                          matrix,
                          coefficient,
                          constraints);
  }



  template <int dim, int spacedim>
  void
  create_laplace_matrix(const Mapping<dim, spacedim> &   mapping,
                        const DoFHandler<dim, spacedim> &dof,
                        const Quadrature<dim> &          q,
                        SparseMatrix<double> &           matrix,
                        const Function<spacedim> &       rhs,
                        Vector<double> &                 rhs_vector,
                        const Function<spacedim> *const  coefficient,
                        const AffineConstraints<double> &constraints)
  {
    Assert(matrix.m() == dof.n_dofs(),
           ExcDimensionMismatch(matrix.m(), dof.n_dofs()));
    Assert(matrix.n() == dof.n_dofs(),
           ExcDimensionMismatch(matrix.n(), dof.n_dofs()));

    hp::FECollection<dim, spacedim>      fe_collection(dof.get_fe());
    hp::QCollection<dim>                 q_collection(q);
    hp::MappingCollection<dim, spacedim> mapping_collection(mapping);
    MatrixCreator::internal::AssemblerData::Scratch<dim, spacedim, double>
                                                             assembler_data(fe_collection,
                     update_gradients | update_values | update_JxW_values |
                       update_quadrature_points,
                     coefficient,
                     &rhs,
                     q_collection,
                     mapping_collection);
    MatrixCreator::internal::AssemblerData::CopyData<double> copy_data;
    copy_data.cell_matrix.reinit(
      assembler_data.fe_collection.max_dofs_per_cell(),
      assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.cell_rhs.reinit(assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.dof_indices.resize(
      assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.constraints = &constraints;

    WorkStream::run(
      dof.begin_active(),
      static_cast<typename DoFHandler<dim, spacedim>::active_cell_iterator>(
        dof.end()),
      &MatrixCreator::internal::laplace_assembler<
        dim,
        spacedim,
        typename DoFHandler<dim, spacedim>::active_cell_iterator>,
      [&matrix,
       &rhs_vector](const internal::AssemblerData::CopyData<double> &data) {
        MatrixCreator::internal::copy_local_to_global(data,
                                                      &matrix,
                                                      &rhs_vector);
      },
      assembler_data,
      copy_data);
  }



  template <int dim, int spacedim>
  void
  create_laplace_matrix(const DoFHandler<dim, spacedim> &dof,
                        const Quadrature<dim> &          q,
                        SparseMatrix<double> &           matrix,
                        const Function<spacedim> &       rhs,
                        Vector<double> &                 rhs_vector,
                        const Function<spacedim> *const  coefficient,
                        const AffineConstraints<double> &constraints)
  {
    create_laplace_matrix(get_default_linear_mapping(dof.get_triangulation()),
                          dof,
                          q,
                          matrix,
                          rhs,
                          rhs_vector,
                          coefficient,
                          constraints);
  }



  template <int dim, int spacedim>
  void
  create_laplace_matrix(const hp::MappingCollection<dim, spacedim> &mapping,
                        const DoFHandler<dim, spacedim> &           dof,
                        const hp::QCollection<dim> &                q,
                        SparseMatrix<double> &                      matrix,
                        const Function<spacedim> *const             coefficient,
                        const AffineConstraints<double> &           constraints)
  {
    Assert(matrix.m() == dof.n_dofs(),
           ExcDimensionMismatch(matrix.m(), dof.n_dofs()));
    Assert(matrix.n() == dof.n_dofs(),
           ExcDimensionMismatch(matrix.n(), dof.n_dofs()));

    MatrixCreator::internal::AssemblerData::Scratch<dim, spacedim, double>
                                                             assembler_data(dof.get_fe_collection(),
                     update_gradients | update_JxW_values |
                       (coefficient != nullptr ? update_quadrature_points :
                                                 UpdateFlags(0)),
                     coefficient,
                     /*rhs_function=*/nullptr,
                     q,
                     mapping);
    MatrixCreator::internal::AssemblerData::CopyData<double> copy_data;
    copy_data.cell_matrix.reinit(
      assembler_data.fe_collection.max_dofs_per_cell(),
      assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.cell_rhs.reinit(assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.dof_indices.resize(
      assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.constraints = &constraints;

    WorkStream::run(
      dof.begin_active(),
      static_cast<typename DoFHandler<dim, spacedim>::active_cell_iterator>(
        dof.end()),
      &MatrixCreator::internal::laplace_assembler<
        dim,
        spacedim,
        typename DoFHandler<dim, spacedim>::active_cell_iterator>,
      [&matrix](const internal::AssemblerData::CopyData<double> &data) {
        MatrixCreator::internal::copy_local_to_global(
          data, &matrix, static_cast<Vector<double> *>(nullptr));
      },
      assembler_data,
      copy_data);
  }



  template <int dim, int spacedim>
  void
  create_laplace_matrix(const DoFHandler<dim, spacedim> &dof,
                        const hp::QCollection<dim> &     q,
                        SparseMatrix<double> &           matrix,
                        const Function<spacedim> *const  coefficient,
                        const AffineConstraints<double> &constraints)
  {
    create_laplace_matrix(
      hp::StaticMappingQ1<dim, spacedim>::mapping_collection,
      dof,
      q,
      matrix,
      coefficient,
      constraints);
  }



  template <int dim, int spacedim>
  void
  create_laplace_matrix(const hp::MappingCollection<dim, spacedim> &mapping,
                        const DoFHandler<dim, spacedim> &           dof,
                        const hp::QCollection<dim> &                q,
                        SparseMatrix<double> &                      matrix,
                        const Function<spacedim> &                  rhs,
                        Vector<double> &                            rhs_vector,
                        const Function<spacedim> *const             coefficient,
                        const AffineConstraints<double> &           constraints)
  {
    Assert(matrix.m() == dof.n_dofs(),
           ExcDimensionMismatch(matrix.m(), dof.n_dofs()));
    Assert(matrix.n() == dof.n_dofs(),
           ExcDimensionMismatch(matrix.n(), dof.n_dofs()));

    MatrixCreator::internal::AssemblerData::Scratch<dim, spacedim, double>
                                                             assembler_data(dof.get_fe_collection(),
                     update_gradients | update_values | update_JxW_values |
                       update_quadrature_points,
                     coefficient,
                     &rhs,
                     q,
                     mapping);
    MatrixCreator::internal::AssemblerData::CopyData<double> copy_data;
    copy_data.cell_matrix.reinit(
      assembler_data.fe_collection.max_dofs_per_cell(),
      assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.cell_rhs.reinit(assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.dof_indices.resize(
      assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.constraints = &constraints;

    WorkStream::run(
      dof.begin_active(),
      static_cast<typename DoFHandler<dim, spacedim>::active_cell_iterator>(
        dof.end()),
      &MatrixCreator::internal::laplace_assembler<
        dim,
        spacedim,
        typename DoFHandler<dim, spacedim>::active_cell_iterator>,
      [&matrix,
       &rhs_vector](const internal::AssemblerData::CopyData<double> &data) {
        MatrixCreator::internal::copy_local_to_global(data,
                                                      &matrix,
                                                      &rhs_vector);
      },
      assembler_data,
      copy_data);
  }



  template <int dim, int spacedim>
  void
  create_laplace_matrix(const DoFHandler<dim, spacedim> &dof,
                        const hp::QCollection<dim> &     q,
                        SparseMatrix<double> &           matrix,
                        const Function<spacedim> &       rhs,
                        Vector<double> &                 rhs_vector,
                        const Function<spacedim> *const  coefficient,
                        const AffineConstraints<double> &constraints)
  {
    create_laplace_matrix(
      hp::StaticMappingQ1<dim, spacedim>::mapping_collection,
      dof,
      q,
      matrix,
      rhs,
      rhs_vector,
      coefficient,
      constraints);
  }

} // namespace MatrixCreator

DEAL_II_NAMESPACE_CLOSE

#endif

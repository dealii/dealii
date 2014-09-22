// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2014 by the deal.II authors
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

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>

#ifdef DEAL_II_WITH_PETSC
#  include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#  include <deal.II/lac/petsc_sparse_matrix.h>
#  include <deal.II/lac/petsc_parallel_vector.h>
#  include <deal.II/lac/petsc_vector.h>
#  include <deal.II/lac/petsc_parallel_block_sparse_matrix.h>
#endif

#ifdef DEAL_II_WITH_TRILINOS
#  include <deal.II/lac/trilinos_sparse_matrix.h>
#  include <deal.II/lac/trilinos_vector.h>
#  include <deal.II/lac/trilinos_block_sparse_matrix.h>
#  include <deal.II/lac/trilinos_block_vector.h>
#endif

#include <algorithm>


#include <algorithm>
#include <set>
#include <cmath>


DEAL_II_NAMESPACE_OPEN


namespace MatrixCreator
{
  namespace internal
  {
    /**
     * Convenience abbreviation for
     * pairs of DoF handler cell
     * iterators. This type works
     * just like a
     * <tt>std::pair<iterator,iterator></tt>
     * but is templatized on the
     * dof handler that shouls be used.
     */
    template <typename DH>
    struct IteratorRange
    {
      /**
       * Typedef for the iterator type.
       */
      typedef typename DH::active_cell_iterator active_cell_iterator;

      /**
       * Abbreviation for a pair of
       * iterators.
       */
      typedef std::pair<active_cell_iterator,active_cell_iterator> iterator_pair;

      /**
       * Constructor. Initialize
       * the two values by the
       * given values.
       */
      IteratorRange (const active_cell_iterator &first,
                     const active_cell_iterator &second);

      /**
       * Constructor taking a pair
       * of values for
       * initialization.
       */
      IteratorRange (const iterator_pair &ip);

      /**
       * Pair of iterators denoting
       * a half-open range.
       */
      active_cell_iterator first, second;
    };




    template <typename DH>
    inline
    IteratorRange<DH>::
    IteratorRange (const active_cell_iterator &first,
                   const active_cell_iterator &second)
      :
      first (first),
      second (second)
    {}



    template <typename DH>
    inline
    IteratorRange<DH>::IteratorRange (const iterator_pair &ip)
      :
      first (ip.first),
      second (ip.second)
    {}



    namespace AssemblerData
    {
      template <int dim,
                int spacedim>
      struct Scratch
      {
        Scratch (const ::dealii::hp::FECollection<dim,spacedim> &fe,
                 const UpdateFlags         update_flags,
                 const Function<spacedim> *coefficient,
                 const Function<spacedim> *rhs_function,
                 const ::dealii::hp::QCollection<dim> &quadrature,
                 const ::dealii::hp::MappingCollection<dim,spacedim> &mapping)
          :
          fe_collection (fe),
          quadrature_collection (quadrature),
          mapping_collection (mapping),
          x_fe_values (mapping_collection,
                       fe_collection,
                       quadrature_collection,
                       update_flags),
          coefficient_values(quadrature_collection.max_n_quadrature_points()),
          coefficient_vector_values (quadrature_collection.max_n_quadrature_points(),
                                     dealii::Vector<double> (fe_collection.n_components())),
          rhs_values(quadrature_collection.max_n_quadrature_points()),
          rhs_vector_values(quadrature_collection.max_n_quadrature_points(),
                            dealii::Vector<double> (fe_collection.n_components())),
          coefficient (coefficient),
          rhs_function (rhs_function),
          update_flags (update_flags)
        {}

        Scratch (const Scratch &data)
          :
          fe_collection (data.fe_collection),
          quadrature_collection (data.quadrature_collection),
          mapping_collection (data.mapping_collection),
          x_fe_values (mapping_collection,
                       fe_collection,
                       quadrature_collection,
                       data.update_flags),
          coefficient_values (data.coefficient_values),
          coefficient_vector_values (data.coefficient_vector_values),
          rhs_values (data.rhs_values),
          rhs_vector_values (data.rhs_vector_values),
          coefficient (data.coefficient),
          rhs_function (data.rhs_function),
          update_flags (data.update_flags)
        {}

        const ::dealii::hp::FECollection<dim,spacedim>      &fe_collection;
        const ::dealii::hp::QCollection<dim>                &quadrature_collection;
        const ::dealii::hp::MappingCollection<dim,spacedim> &mapping_collection;

        ::dealii::hp::FEValues<dim,spacedim> x_fe_values;

        std::vector<double>                  coefficient_values;
        std::vector<dealii::Vector<double> > coefficient_vector_values;
        std::vector<double>                  rhs_values;
        std::vector<dealii::Vector<double> > rhs_vector_values;

        const Function<spacedim>   *coefficient;
        const Function<spacedim>   *rhs_function;

        const UpdateFlags update_flags;
      };


      struct CopyData
      {
        std::vector<types::global_dof_index> dof_indices;
        FullMatrix<double>        cell_matrix;
        dealii::Vector<double>    cell_rhs;
        const ConstraintMatrix   *constraints;
      };
    }


    template <int dim,
              int spacedim,
              typename CellIterator>
    void mass_assembler (const CellIterator &cell,
                         MatrixCreator::internal::AssemblerData::Scratch<dim,spacedim> &data,
                         MatrixCreator::internal::AssemblerData::CopyData              &copy_data)
    {
      data.x_fe_values.reinit (cell);
      const FEValues<dim,spacedim> &fe_values = data.x_fe_values.get_present_fe_values ();

      const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                         n_q_points    = fe_values.n_quadrature_points;
      const FiniteElement<dim,spacedim> &fe  = fe_values.get_fe();
      const unsigned int n_components  = fe.n_components();

      Assert(data.rhs_function == 0 ||
             data.rhs_function->n_components==1 ||
             data.rhs_function->n_components==n_components,
             ::dealii::MatrixCreator::ExcComponentMismatch());
      Assert(data.coefficient == 0 ||
             data.coefficient->n_components==1 ||
             data.coefficient->n_components==n_components,
             ::dealii::MatrixCreator::ExcComponentMismatch());

      copy_data.cell_matrix.reinit (dofs_per_cell, dofs_per_cell);
      copy_data.cell_rhs.reinit (dofs_per_cell);

      copy_data.dof_indices.resize (dofs_per_cell);
      cell->get_dof_indices (copy_data.dof_indices);

      const bool use_rhs_function = data.rhs_function != 0;
      if (use_rhs_function)
        {
          if (data.rhs_function->n_components==1)
            {
              data.rhs_values.resize (n_q_points);
              data.rhs_function->value_list (fe_values.get_quadrature_points(),
                                             data.rhs_values);
            }
          else
            {
              data.rhs_vector_values.resize (n_q_points,
                                             dealii::Vector<double>(n_components));
              data.rhs_function->vector_value_list (fe_values.get_quadrature_points(),
                                                    data.rhs_vector_values);
            }
        }

      const bool use_coefficient = data.coefficient != 0;
      if (use_coefficient)
        {
          if (data.coefficient->n_components==1)
            {
              data.coefficient_values.resize (n_q_points);
              data.coefficient->value_list (fe_values.get_quadrature_points(),
                                            data.coefficient_values);
            }
          else
            {
              data.coefficient_vector_values.resize (n_q_points,
                                                     dealii::Vector<double>(n_components));
              data.coefficient->vector_value_list (fe_values.get_quadrature_points(),
                                                   data.coefficient_vector_values);
            }
        }


      double add_data;
      const std::vector<double> &JxW = fe_values.get_JxW_values();
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        if (fe.is_primitive ())
          {
            const unsigned int component_i =
              fe.system_to_component_index(i).first;
            const double *phi_i = &fe_values.shape_value(i,0);
            add_data = 0;

            // use symmetry in the mass matrix here:
            // just need to calculate the diagonal
            // and half of the elements above the
            // diagonal
            for (unsigned int j=i; j<dofs_per_cell; ++j)
              if ((n_components==1) ||
                  (fe.system_to_component_index(j).first ==
                   component_i))
                {
                  const double *phi_j = &fe_values.shape_value(j,0);
                  add_data = 0;
                  if (use_coefficient)
                    {
                      if (data.coefficient->n_components==1)
                        for (unsigned int point=0; point<n_q_points; ++point)
                          add_data += (phi_i[point] * phi_j[point] * JxW[point] *
                                       data.coefficient_values[point]);
                      else
                        for (unsigned int point=0; point<n_q_points; ++point)
                          add_data += (phi_i[point] * phi_j[point] * JxW[point] *
                                       data.coefficient_vector_values[point](component_i));
                    }
                  else
                    for (unsigned int point=0; point<n_q_points; ++point)
                      add_data += phi_i[point] * phi_j[point] * JxW[point];

                  // this is even ok for i==j, since then
                  // we just write the same value twice.
                  copy_data.cell_matrix(i,j) = add_data;
                  copy_data.cell_matrix(j,i) = add_data;
                }

            if (use_rhs_function)
              {
                add_data = 0;
                if (data.rhs_function->n_components==1)
                  for (unsigned int point=0; point<n_q_points; ++point)
                    add_data += phi_i[point] * JxW[point] *
                                data.rhs_values[point];
                else
                  for (unsigned int point=0; point<n_q_points; ++point)
                    add_data += phi_i[point] * JxW[point] *
                                data.rhs_vector_values[point](component_i);
                copy_data.cell_rhs(i) = add_data;
              }
          }
        else
          {
            // non-primitive vector-valued FE, using
            // symmetry again
            for (unsigned int j=i; j<dofs_per_cell; ++j)
              {
                add_data = 0;
                for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
                  if (fe.get_nonzero_components(i)[comp_i] &&
                      fe.get_nonzero_components(j)[comp_i])
                    {
                      if (use_coefficient)
                        {
                          if (data.coefficient->n_components==1)
                            for (unsigned int point=0; point<n_q_points; ++point)
                              add_data += (fe_values.shape_value_component(i,point,comp_i) *
                                           fe_values.shape_value_component(j,point,comp_i) *
                                           JxW[point] *
                                           data.coefficient_values[point]);
                          else
                            for (unsigned int point=0; point<n_q_points; ++point)
                              add_data += (fe_values.shape_value_component(i,point,comp_i) *
                                           fe_values.shape_value_component(j,point,comp_i) *
                                           JxW[point] *
                                           data.coefficient_vector_values[point](comp_i));
                        }
                      else
                        for (unsigned int point=0; point<n_q_points; ++point)
                          add_data += fe_values.shape_value_component(i,point,comp_i) *
                                      fe_values.shape_value_component(j,point,comp_i) * JxW[point];
                    }

                copy_data.cell_matrix(i,j) = add_data;
                copy_data.cell_matrix(j,i) = add_data;
              }

            if (use_rhs_function)
              {
                add_data = 0;
                for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
                  if (fe.get_nonzero_components(i)[comp_i])
                    {
                      if (data.rhs_function->n_components==1)
                        for (unsigned int point=0; point<n_q_points; ++point)
                          add_data += fe_values.shape_value_component(i,point,comp_i) *
                                      JxW[point] * data.rhs_values[point];
                      else
                        for (unsigned int point=0; point<n_q_points; ++point)
                          add_data += fe_values.shape_value_component(i,point,comp_i) *
                                      JxW[point] * data.rhs_vector_values[point](comp_i);
                    }
                copy_data.cell_rhs(i) = add_data;
              }
          }
    }



    template <int dim,
              int spacedim,
              typename CellIterator>
    void laplace_assembler (const CellIterator &cell,
                            MatrixCreator::internal::AssemblerData::Scratch<dim,spacedim> &data,
                            MatrixCreator::internal::AssemblerData::CopyData              &copy_data)
    {
      data.x_fe_values.reinit (cell);
      const FEValues<dim,spacedim> &fe_values = data.x_fe_values.get_present_fe_values ();

      const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                         n_q_points    = fe_values.n_quadrature_points;
      const FiniteElement<dim,spacedim>    &fe  = fe_values.get_fe();
      const unsigned int n_components  = fe.n_components();

      Assert(data.rhs_function == 0 ||
             data.rhs_function->n_components==1 ||
             data.rhs_function->n_components==n_components,
             ::dealii::MatrixCreator::ExcComponentMismatch());
      Assert(data.coefficient == 0 ||
             data.coefficient->n_components==1 ||
             data.coefficient->n_components==n_components,
             ::dealii::MatrixCreator::ExcComponentMismatch());

      copy_data.cell_matrix.reinit (dofs_per_cell, dofs_per_cell);
      copy_data.cell_rhs.reinit (dofs_per_cell);
      copy_data.dof_indices.resize (dofs_per_cell);
      cell->get_dof_indices (copy_data.dof_indices);


      const bool use_rhs_function = data.rhs_function != 0;
      if (use_rhs_function)
        {
          if (data.rhs_function->n_components==1)
            {
              data.rhs_values.resize (n_q_points);
              data.rhs_function->value_list (fe_values.get_quadrature_points(),
                                             data.rhs_values);
            }
          else
            {
              data.rhs_vector_values.resize (n_q_points,
                                             dealii::Vector<double>(n_components));
              data.rhs_function->vector_value_list (fe_values.get_quadrature_points(),
                                                    data.rhs_vector_values);
            }
        }

      const bool use_coefficient = data.coefficient != 0;
      if (use_coefficient)
        {
          if (data.coefficient->n_components==1)
            {
              data.coefficient_values.resize (n_q_points);
              data.coefficient->value_list (fe_values.get_quadrature_points(),
                                            data.coefficient_values);
            }
          else
            {
              data.coefficient_vector_values.resize (n_q_points,
                                                     dealii::Vector<double>(n_components));
              data.coefficient->vector_value_list (fe_values.get_quadrature_points(),
                                                   data.coefficient_vector_values);
            }
        }


      const std::vector<double> &JxW = fe_values.get_JxW_values();
      double add_data;
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        if (fe.is_primitive ())
          {
            const unsigned int component_i =
              fe.system_to_component_index(i).first;
            const Tensor<1,spacedim> *grad_phi_i =
              &fe_values.shape_grad(i,0);

            // can use symmetry
            for (unsigned int j=i; j<dofs_per_cell; ++j)
              if ((n_components==1) ||
                  (fe.system_to_component_index(j).first ==
                   component_i))
                {
                  const Tensor<1,spacedim> *grad_phi_j =
                    & fe_values.shape_grad(j,0);
                  add_data = 0;
                  if (use_coefficient)
                    {
                      if (data.coefficient->n_components==1)
                        for (unsigned int point=0; point<n_q_points; ++point)
                          add_data += ((grad_phi_i[point]*grad_phi_j[point]) *
                                       JxW[point] *
                                       data.coefficient_values[point]);
                      else
                        for (unsigned int point=0; point<n_q_points; ++point)
                          add_data += ((grad_phi_i[point]*grad_phi_j[point]) *
                                       JxW[point] *
                                       data.coefficient_vector_values[point](component_i));
                    }
                  else
                    for (unsigned int point=0; point<n_q_points; ++point)
                      add_data += (grad_phi_i[point]*grad_phi_j[point]) *
                                  JxW[point];

                  copy_data.cell_matrix(i,j) = add_data;
                  copy_data.cell_matrix(j,i) = add_data;
                }

            if (use_rhs_function)
              {
                const double *phi_i = &fe_values.shape_value(i,0);
                add_data = 0;
                if (data.rhs_function->n_components==1)
                  for (unsigned int point=0; point<n_q_points; ++point)
                    add_data += phi_i[point] * JxW[point] *
                                data.rhs_values[point];
                else
                  for (unsigned int point=0; point<n_q_points; ++point)
                    add_data += phi_i[point] * JxW[point] *
                                data.rhs_vector_values[point](component_i);
                copy_data.cell_rhs(i) = add_data;
              }
          }
        else
          {
            // non-primitive vector-valued FE
            for (unsigned int j=i; j<dofs_per_cell; ++j)
              {
                add_data = 0;
                for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
                  if (fe.get_nonzero_components(i)[comp_i] &&
                      fe.get_nonzero_components(j)[comp_i])
                    {
                      if (use_coefficient)
                        {
                          if (data.coefficient->n_components==1)
                            for (unsigned int point=0; point<n_q_points; ++point)
                              add_data += ((fe_values.shape_grad_component(i,point,comp_i) *
                                            fe_values.shape_grad_component(j,point,comp_i)) *
                                           JxW[point] *
                                           data.coefficient_values[point]);
                          else
                            for (unsigned int point=0; point<n_q_points; ++point)
                              add_data += ((fe_values.shape_grad_component(i,point,comp_i) *
                                            fe_values.shape_grad_component(j,point,comp_i)) *
                                           JxW[point] *
                                           data.coefficient_vector_values[point](comp_i));
                        }
                      else
                        for (unsigned int point=0; point<n_q_points; ++point)
                          add_data += (fe_values.shape_grad_component(i,point,comp_i) *
                                       fe_values.shape_grad_component(j,point,comp_i)) *
                                      JxW[point];
                    }

                copy_data.cell_matrix(i,j) = add_data;
                copy_data.cell_matrix(j,i) = add_data;
              }

            if (use_rhs_function)
              {
                add_data = 0;
                for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
                  if (fe.get_nonzero_components(i)[comp_i])
                    {
                      if (data.rhs_function->n_components==1)
                        for (unsigned int point=0; point<n_q_points; ++point)
                          add_data += fe_values.shape_value_component(i,point,comp_i) *
                                      JxW[point] * data.rhs_values[point];
                      else
                        for (unsigned int point=0; point<n_q_points; ++point)
                          add_data += fe_values.shape_value_component(i,point,comp_i) *
                                      JxW[point] * data.rhs_vector_values[point](comp_i);
                    }
                copy_data.cell_rhs(i) = add_data;
              }
          }
    }



    template <typename MatrixType,
              typename VectorType>
    void copy_local_to_global (const AssemblerData::CopyData &data,
                               MatrixType *matrix,
                               VectorType *right_hand_side)
    {
      const unsigned int dofs_per_cell = data.dof_indices.size();

      Assert (data.cell_matrix.m() == dofs_per_cell,
              ExcInternalError());
      Assert (data.cell_matrix.n() == dofs_per_cell,
              ExcInternalError());
      Assert ((right_hand_side == 0)
              ||
              (data.cell_rhs.size() == dofs_per_cell),
              ExcInternalError());

      if (right_hand_side != 0)
        data.constraints->distribute_local_to_global(data.cell_matrix,
                                                     data.cell_rhs,
                                                     data.dof_indices,
                                                     *matrix, *right_hand_side);
      else
        data.constraints->distribute_local_to_global(data.cell_matrix,
                                                     data.dof_indices,
                                                     *matrix);
    }



    namespace AssemblerBoundary
    {
      struct Scratch
      {
        Scratch() {}
      };

      template <typename DH>
      struct CopyData
      {
        CopyData() {};

        CopyData(CopyData const &data);

        unsigned int dofs_per_cell;
        std::vector<types::global_dof_index> dofs;
        std::vector<std::vector<bool> > dof_is_on_face;
        typename DH::active_cell_iterator cell;
        std::vector<FullMatrix<double> > cell_matrix;
        std::vector<Vector<double> > cell_vector;
      };

      template <typename DH>
      CopyData<DH>::CopyData(CopyData const &data) :
        dofs_per_cell(data.dofs_per_cell),
        dofs(data.dofs),
        dof_is_on_face(data.dof_is_on_face),
        cell(data.cell),
        cell_matrix(data.cell_matrix),
        cell_vector(data.cell_vector)
      {}
    }
  }
}


namespace MatrixCreator
{

  template <int dim, typename number, int spacedim>
  void create_mass_matrix (const Mapping<dim,spacedim>       &mapping,
                           const DoFHandler<dim,spacedim>    &dof,
                           const Quadrature<dim>    &q,
                           SparseMatrix<number>     &matrix,
                           const Function<spacedim> *const coefficient,
                           const ConstraintMatrix   &constraints)
  {
    Assert (matrix.m() == dof.n_dofs(),
            ExcDimensionMismatch (matrix.m(), dof.n_dofs()));
    Assert (matrix.n() == dof.n_dofs(),
            ExcDimensionMismatch (matrix.n(), dof.n_dofs()));

    hp::FECollection<dim,spacedim>      fe_collection (dof.get_fe());
    hp::QCollection<dim>                q_collection (q);
    hp::MappingCollection<dim,spacedim> mapping_collection (mapping);
    MatrixCreator::internal::AssemblerData::Scratch<dim, spacedim>
    assembler_data (fe_collection,
                    update_values | update_JxW_values |
                    (coefficient != 0 ? update_quadrature_points : UpdateFlags(0)),
                    coefficient, /*rhs_function=*/0,
                    q_collection, mapping_collection);

    MatrixCreator::internal::AssemblerData::CopyData copy_data;
    copy_data.cell_matrix.reinit (assembler_data.fe_collection.max_dofs_per_cell(),
                                  assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.cell_rhs.reinit (assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.dof_indices.resize (assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.constraints = &constraints;

    WorkStream::run (dof.begin_active(),
                     static_cast<typename DoFHandler<dim,spacedim>::active_cell_iterator>(dof.end()),
                     &MatrixCreator::internal::mass_assembler<dim, spacedim, typename DoFHandler<dim,spacedim>::active_cell_iterator>,
                     std_cxx11::bind (&MatrixCreator::internal::
                                      copy_local_to_global<SparseMatrix<number>, Vector<double> >,
                                      std_cxx11::_1, &matrix, (Vector<double> *)0),
                     assembler_data,
                     copy_data);
  }



  template <int dim, typename number, int spacedim>
  void create_mass_matrix (const DoFHandler<dim,spacedim>    &dof,
                           const Quadrature<dim>    &q,
                           SparseMatrix<number>     &matrix,
                           const Function<spacedim> *const coefficient,
                           const ConstraintMatrix   &constraints)
  {
    create_mass_matrix(StaticMappingQ1<dim,spacedim>::mapping, dof,
                       q, matrix, coefficient, constraints);
  }



  template <int dim, typename number, int spacedim>
  void create_mass_matrix (const Mapping<dim,spacedim>       &mapping,
                           const DoFHandler<dim,spacedim>    &dof,
                           const Quadrature<dim>    &q,
                           SparseMatrix<number>     &matrix,
                           const Function<spacedim>      &rhs,
                           Vector<double>           &rhs_vector,
                           const Function<spacedim> *const coefficient,
                           const ConstraintMatrix   &constraints)
  {
    Assert (matrix.m() == dof.n_dofs(),
            ExcDimensionMismatch (matrix.m(), dof.n_dofs()));
    Assert (matrix.n() == dof.n_dofs(),
            ExcDimensionMismatch (matrix.n(), dof.n_dofs()));

    hp::FECollection<dim,spacedim>      fe_collection (dof.get_fe());
    hp::QCollection<dim>                q_collection (q);
    hp::MappingCollection<dim,spacedim> mapping_collection (mapping);
    MatrixCreator::internal::AssemblerData::Scratch<dim, spacedim>
    assembler_data (fe_collection,
                    update_values |
                    update_JxW_values | update_quadrature_points,
                    coefficient, &rhs,
                    q_collection, mapping_collection);
    MatrixCreator::internal::AssemblerData::CopyData copy_data;
    copy_data.cell_matrix.reinit (assembler_data.fe_collection.max_dofs_per_cell(),
                                  assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.cell_rhs.reinit (assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.dof_indices.resize (assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.constraints = &constraints;

    WorkStream::run (dof.begin_active(),
                     static_cast<typename DoFHandler<dim,spacedim>::active_cell_iterator>(dof.end()),
                     &MatrixCreator::internal::mass_assembler<dim, spacedim, typename DoFHandler<dim,spacedim>::active_cell_iterator>,
                     std_cxx11::bind(&MatrixCreator::internal::
                                     copy_local_to_global<SparseMatrix<number>, Vector<double> >,
                                     std_cxx11::_1, &matrix, &rhs_vector),
                     assembler_data,
                     copy_data);
  }



  template <int dim, typename number, int spacedim>
  void create_mass_matrix (const DoFHandler<dim,spacedim>    &dof,
                           const Quadrature<dim>    &q,
                           SparseMatrix<number>     &matrix,
                           const Function<spacedim>      &rhs,
                           Vector<double>           &rhs_vector,
                           const Function<spacedim> *const coefficient,
                           const ConstraintMatrix   &constraints)
  {
    create_mass_matrix(StaticMappingQ1<dim,spacedim>::mapping,
                       dof, q, matrix, rhs, rhs_vector, coefficient,
                       constraints);
  }



  template <int dim, typename number, int spacedim>
  void create_mass_matrix (const hp::MappingCollection<dim,spacedim> &mapping,
                           const hp::DoFHandler<dim,spacedim>    &dof,
                           const hp::QCollection<dim>    &q,
                           SparseMatrix<number>     &matrix,
                           const Function<spacedim> *const coefficient,
                           const ConstraintMatrix   &constraints)
  {
    Assert (matrix.m() == dof.n_dofs(),
            ExcDimensionMismatch (matrix.m(), dof.n_dofs()));
    Assert (matrix.n() == dof.n_dofs(),
            ExcDimensionMismatch (matrix.n(), dof.n_dofs()));

    MatrixCreator::internal::AssemblerData::Scratch<dim, spacedim>
    assembler_data (dof.get_fe(),
                    update_values | update_JxW_values |
                    (coefficient != 0 ? update_quadrature_points : UpdateFlags(0)),
                    coefficient, /*rhs_function=*/0,
                    q, mapping);
    MatrixCreator::internal::AssemblerData::CopyData copy_data;
    copy_data.cell_matrix.reinit (assembler_data.fe_collection.max_dofs_per_cell(),
                                  assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.cell_rhs.reinit (assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.dof_indices.resize (assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.constraints = &constraints;

    WorkStream::run (dof.begin_active(),
                     static_cast<typename hp::DoFHandler<dim,spacedim>::active_cell_iterator>(dof.end()),
                     &MatrixCreator::internal::mass_assembler<dim, spacedim, typename hp::DoFHandler<dim,spacedim>::active_cell_iterator>,
                     std_cxx11::bind (&MatrixCreator::internal::
                                      copy_local_to_global<SparseMatrix<number>, Vector<double> >,
                                      std_cxx11::_1, &matrix, (Vector<double> *)0),
                     assembler_data,
                     copy_data);
  }



  template <int dim, typename number, int spacedim>
  void create_mass_matrix (const hp::DoFHandler<dim,spacedim> &dof,
                           const hp::QCollection<dim> &q,
                           SparseMatrix<number>     &matrix,
                           const Function<spacedim> *const coefficient,
                           const ConstraintMatrix   &constraints)
  {
    create_mass_matrix(hp::StaticMappingQ1<dim,spacedim>::mapping_collection,
                       dof, q, matrix, coefficient, constraints);
  }



  template <int dim, typename number, int spacedim>
  void create_mass_matrix (const hp::MappingCollection<dim,spacedim> &mapping,
                           const hp::DoFHandler<dim,spacedim> &dof,
                           const hp::QCollection<dim>    &q,
                           SparseMatrix<number>     &matrix,
                           const Function<spacedim>      &rhs,
                           Vector<double>           &rhs_vector,
                           const Function<spacedim> *const coefficient,
                           const ConstraintMatrix   &constraints)
  {
    Assert (matrix.m() == dof.n_dofs(),
            ExcDimensionMismatch (matrix.m(), dof.n_dofs()));
    Assert (matrix.n() == dof.n_dofs(),
            ExcDimensionMismatch (matrix.n(), dof.n_dofs()));

    MatrixCreator::internal::AssemblerData::Scratch<dim, spacedim>
    assembler_data (dof.get_fe(),
                    update_values |
                    update_JxW_values | update_quadrature_points,
                    coefficient, &rhs,
                    q, mapping);
    MatrixCreator::internal::AssemblerData::CopyData copy_data;
    copy_data.cell_matrix.reinit (assembler_data.fe_collection.max_dofs_per_cell(),
                                  assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.cell_rhs.reinit (assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.dof_indices.resize (assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.constraints = &constraints;

    WorkStream::run (dof.begin_active(),
                     static_cast<typename hp::DoFHandler<dim,spacedim>::active_cell_iterator>(dof.end()),
                     &MatrixCreator::internal::mass_assembler<dim, spacedim, typename hp::DoFHandler<dim,spacedim>::active_cell_iterator>,
                     std_cxx11::bind (&MatrixCreator::internal::
                                      copy_local_to_global<SparseMatrix<number>, Vector<double> >,
                                      std_cxx11::_1, &matrix, &rhs_vector),
                     assembler_data,
                     copy_data);
  }



  template <int dim, typename number, int spacedim>
  void create_mass_matrix (const hp::DoFHandler<dim,spacedim> &dof,
                           const hp::QCollection<dim> &q,
                           SparseMatrix<number>     &matrix,
                           const Function<spacedim> &rhs,
                           Vector<double>           &rhs_vector,
                           const Function<spacedim> *const coefficient,
                           const ConstraintMatrix   &constraints)
  {
    create_mass_matrix(hp::StaticMappingQ1<dim,spacedim>::mapping_collection, dof, q,
                       matrix, rhs, rhs_vector, coefficient, constraints);
  }



  namespace
  {
    template <int dim, int spacedim>
    void
    create_boundary_mass_matrix_1 (typename DoFHandler<dim,spacedim>::active_cell_iterator const &cell,
                                   MatrixCreator::internal::AssemblerBoundary::Scratch const &scratch,
                                   MatrixCreator::internal::AssemblerBoundary::CopyData<DoFHandler<dim,
                                   spacedim> > &copy_data,
                                   Mapping<dim, spacedim> const &mapping,
                                   FiniteElement<dim,spacedim> const &fe,
                                   Quadrature<dim-1> const &q,
                                   typename FunctionMap<spacedim>::type const &boundary_functions,
                                   Function<spacedim> const *const coefficient,
                                   std::vector<unsigned int> const &component_mapping)

    {
      // All assertions for this function
      // are in the calling function
      // before creating threads.
      const unsigned int n_components = fe.n_components();
      const unsigned int n_function_components = boundary_functions.begin()->second->n_components;
      const bool         fe_is_system = (n_components != 1);
      const bool         fe_is_primitive = fe.is_primitive();

      const unsigned int dofs_per_face = fe.dofs_per_face;

      copy_data.cell = cell;
      copy_data.dofs_per_cell = fe.dofs_per_cell;

      UpdateFlags update_flags = UpdateFlags (update_values     |
                                              update_JxW_values |
                                              update_normal_vectors |
                                              update_quadrature_points);
      FEFaceValues<dim,spacedim> fe_values (mapping, fe, q, update_flags);

      // two variables for the coefficient,
      // one for the two cases indicated in
      // the name
      std::vector<double>          coefficient_values (fe_values.n_quadrature_points, 1.);
      std::vector<Vector<double> > coefficient_vector_values (fe_values.n_quadrature_points,
                                                              Vector<double>(n_components));
      const bool coefficient_is_vector = (coefficient != 0 && coefficient->n_components != 1)
                                         ? true : false;

      std::vector<double>          rhs_values_scalar (fe_values.n_quadrature_points);
      std::vector<Vector<double> > rhs_values_system (fe_values.n_quadrature_points,
                                                      Vector<double>(n_function_components));

      copy_data.dofs.resize(copy_data.dofs_per_cell);
      cell->get_dof_indices (copy_data.dofs);

      std::vector<types::global_dof_index> dofs_on_face_vector (dofs_per_face);

      // Because CopyData objects are reused and that push_back is used,
      // dof_is_on_face, cell_matrix, and cell_vector must be cleared before
      // they are reused
      copy_data.dof_is_on_face.clear();
      copy_data.cell_matrix.clear();
      copy_data.cell_vector.clear();

      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
        // check if this face is on that part of
        // the boundary we are interested in
        if (boundary_functions.find(cell->face(face)->boundary_indicator()) !=
            boundary_functions.end())
          {
            copy_data.cell_matrix.push_back(FullMatrix<double> (copy_data.dofs_per_cell,
                                                                copy_data.dofs_per_cell));
            copy_data.cell_vector.push_back(Vector<double> (copy_data.dofs_per_cell));
            fe_values.reinit (cell, face);

            if (fe_is_system)
              // FE has several components
              {
                boundary_functions.find(cell->face(face)->boundary_indicator())
                ->second->vector_value_list (fe_values.get_quadrature_points(),
                                             rhs_values_system);

                if (coefficient_is_vector)
                  // If coefficient is
                  // vector valued, fill
                  // all components
                  coefficient->vector_value_list (fe_values.get_quadrature_points(),
                                                  coefficient_vector_values);
                else
                  {
                    // If a scalar
                    // function is
                    // given, update
                    // the values, if
                    // not, use the
                    // default one set
                    // in the
                    // constructor above
                    if (coefficient != 0)
                      coefficient->value_list (fe_values.get_quadrature_points(),
                                               coefficient_values);
                    // Copy scalar
                    // values into vector
                    for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
                      coefficient_vector_values[point] = coefficient_values[point];
                  }

                // Special treatment
                // for Hdiv and Hcurl
                // elements, where only
                // the normal or
                // tangential component
                // should be projected.
                std::vector<std::vector<double> > normal_adjustment(fe_values.n_quadrature_points,
                                                                    std::vector<double>(n_components, 1.));

                for (unsigned int comp = 0; comp<n_components; ++comp)
                  {
                    const FiniteElement<dim,spacedim> &base = fe.base_element(fe.component_to_base_index(comp).first);
                    const unsigned int bcomp = fe.component_to_base_index(comp).second;

                    if (!base.conforms(FiniteElementData<dim>::H1) &&
                        base.conforms(FiniteElementData<dim>::Hdiv))
                      for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
                        normal_adjustment[point][comp] = fe_values.normal_vector(point)(bcomp)
                                                         * fe_values.normal_vector(point)(bcomp);
                  }

                for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
                  {
                    const double weight = fe_values.JxW(point);
                    for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i)
                      if (fe_is_primitive)
                        {
                          for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
                            {
                              if (fe.system_to_component_index(j).first
                                  == fe.system_to_component_index(i).first)
                                {
                                  copy_data.cell_matrix.back()(i,j)
                                  += weight
                                     * fe_values.shape_value(j,point)
                                     * fe_values.shape_value(i,point)
                                     * coefficient_vector_values[point](fe.system_to_component_index(i).first);
                                }
                            }
                          copy_data.cell_vector.back()(i) += fe_values.shape_value(i,point)
                                                             * rhs_values_system[point](component_mapping[fe.system_to_component_index(i).first])
                                                             * weight;
                        }
                      else
                        {
                          for (unsigned int comp=0; comp<n_components; ++comp)
                            {
                              for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
                                copy_data.cell_matrix.back()(i,j)
                                += fe_values.shape_value_component(j,point,comp)
                                   * fe_values.shape_value_component(i,point,comp)
                                   * normal_adjustment[point][comp]
                                   * weight * coefficient_vector_values[point](comp);
                              copy_data.cell_vector.back()(i) += fe_values.shape_value_component(i,point,comp) *
                                                                 rhs_values_system[point](component_mapping[comp])
                                                                 * normal_adjustment[point][comp]
                                                                 * weight;
                            }
                        }
                  }
              }
            else
              // FE is a scalar one
              {
                boundary_functions.find(cell->face(face)->boundary_indicator())
                ->second->value_list (fe_values.get_quadrature_points(), rhs_values_scalar);

                if (coefficient != 0)
                  coefficient->value_list (fe_values.get_quadrature_points(),
                                           coefficient_values);
                for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
                  {
                    const double weight = fe_values.JxW(point);
                    for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i)
                      {
                        const double v = fe_values.shape_value(i,point);
                        for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
                          {
                            const double u = fe_values.shape_value(j,point);
                            copy_data.cell_matrix.back()(i,j) += (u*v*weight*coefficient_values[point]);
                          }
                        copy_data.cell_vector.back()(i) += v * rhs_values_scalar[point] *weight;
                      }
                  }
              }


            cell->face(face)->get_dof_indices (dofs_on_face_vector);
            // for each dof on the cell, have a
            // flag whether it is on the face
            copy_data.dof_is_on_face.push_back(std::vector<bool> (copy_data.dofs_per_cell));
            // check for each of the dofs on this cell
            // whether it is on the face
            for (unsigned int i=0; i<copy_data.dofs_per_cell; ++i)
              copy_data.dof_is_on_face.back()[i] = (std::find(dofs_on_face_vector.begin(),
                                                              dofs_on_face_vector.end(),
                                                              copy_data.dofs[i])
                                                    !=
                                                    dofs_on_face_vector.end());
          }
    }

    template <int dim,int spacedim>
    void copy_boundary_mass_matrix_1(MatrixCreator::internal::AssemblerBoundary::CopyData<DoFHandler<dim,
                                     spacedim> > const &copy_data,
                                     typename FunctionMap<spacedim>::type const &boundary_functions,
                                     std::vector<types::global_dof_index> const &dof_to_boundary_mapping,
                                     SparseMatrix<double> &matrix,
                                     Vector<double> &rhs_vector)
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
      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
        {
          // check if this face is on that part of
          // the boundary we are interested in
          if (boundary_functions.find(copy_data.cell->face(face)->boundary_indicator()) !=
              boundary_functions.end())
            {
              for (unsigned int i=0; i<copy_data.dofs_per_cell; ++i)
                {
                  if (copy_data.dof_is_on_face[pos][i] &&
                      dof_to_boundary_mapping[copy_data.dofs[i]] != numbers::invalid_dof_index)
                    {
                      for (unsigned int j=0; j<copy_data.dofs_per_cell; ++j)
                        if (copy_data.dof_is_on_face[pos][j] &&
                            dof_to_boundary_mapping[copy_data.dofs[j]] != numbers::invalid_dof_index)
                          {
                            Assert(numbers::is_finite(copy_data.cell_matrix[pos](i,j)),
                                   ExcNumberNotFinite());
                            matrix.add(dof_to_boundary_mapping[copy_data.dofs[i]],
                                       dof_to_boundary_mapping[copy_data.dofs[j]],
                                       copy_data.cell_matrix[pos](i,j));
                          }
                      Assert(numbers::is_finite(copy_data.cell_vector[pos](i)), ExcNumberNotFinite());
                      rhs_vector(dof_to_boundary_mapping[copy_data.dofs[i]]) += copy_data.cell_vector[pos](i);
                    }
                }
              ++pos;
            }
        }
    }


    template <>
    void
    create_boundary_mass_matrix_1<1,3> (DoFHandler<1,3>::active_cell_iterator const &cell,
                                        MatrixCreator::internal::AssemblerBoundary::Scratch const
                                        &scratch,
                                        MatrixCreator::internal::AssemblerBoundary::CopyData<DoFHandler<1,
                                        3> > &copy_data,
                                        Mapping<1,3> const &mapping,
                                        FiniteElement<1,3> const &fe,
                                        Quadrature<0> const &q,
                                        FunctionMap<3>::type const &boundary_functions,
                                        Function<3> const *const coefficient,
                                        std::vector<unsigned int> const &component_mapping)
    {
      Assert(false,ExcNotImplemented());
    }
  }




  template <int dim, int spacedim>
  void
  create_boundary_mass_matrix (const Mapping<dim, spacedim>  &mapping,
                               const DoFHandler<dim,spacedim> &dof,
                               const Quadrature<dim-1>  &q,
                               SparseMatrix<double>  &matrix,
                               const typename FunctionMap<spacedim>::type  &boundary_functions,
                               Vector<double>            &rhs_vector,
                               std::vector<types::global_dof_index> &dof_to_boundary_mapping,
                               const Function<spacedim> *const coefficient,
                               std::vector<unsigned int> component_mapping)
  {
    // what would that be in 1d? the
    // identity matrix on the boundary
    // dofs?
    if (dim == 1)
      {
        Assert (false, ExcNotImplemented());
        return;
      }

    const FiniteElement<dim,spacedim> &fe = dof.get_fe();
    const unsigned int n_components  = fe.n_components();

    Assert (matrix.n() == dof.n_boundary_dofs(boundary_functions),
            ExcInternalError());
    Assert (matrix.n() == matrix.m(), ExcInternalError());
    Assert (matrix.n() == rhs_vector.size(), ExcInternalError());
    Assert (boundary_functions.size() != 0, ExcInternalError());
    Assert (dof_to_boundary_mapping.size() == dof.n_dofs(),
            ExcInternalError());
    Assert (coefficient ==0 ||
            coefficient->n_components==1 ||
            coefficient->n_components==n_components, ExcComponentMismatch());

    if (component_mapping.size() == 0)
      {
        AssertDimension (n_components, boundary_functions.begin()->second->n_components);
        for (unsigned int i=0; i<n_components; ++i)
          component_mapping.push_back(i);
      }
    else
      AssertDimension (n_components, component_mapping.size());

    MatrixCreator::internal::AssemblerBoundary::Scratch scratch;
    MatrixCreator::internal::AssemblerBoundary::CopyData<DoFHandler<dim,spacedim> > copy_data;

    WorkStream::run(dof.begin_active(),dof.end(),
                    static_cast<std_cxx11::function<void (typename DoFHandler<dim,spacedim>::active_cell_iterator
                                                          const &,MatrixCreator::internal::AssemblerBoundary::Scratch const &,
                                                          MatrixCreator::internal::AssemblerBoundary::CopyData<DoFHandler<dim,spacedim> > &)> >
                    (std_cxx11::bind(&create_boundary_mass_matrix_1<dim,spacedim>,std_cxx11::_1,std_cxx11::_2,
                                     std_cxx11::_3,
                                     std_cxx11::cref(mapping),std_cxx11::cref(fe),std_cxx11::cref(q),
                                     std_cxx11::cref(boundary_functions),coefficient,
                                     std_cxx11::cref(component_mapping))),
                    static_cast<std_cxx11::function<void (MatrixCreator::internal::AssemblerBoundary
                                                          ::CopyData<DoFHandler<dim,spacedim> > const &)> > (std_cxx11::bind(
                                                                &copy_boundary_mass_matrix_1<dim,spacedim>,
                                                                std_cxx11::_1,
                                                                std_cxx11::cref(boundary_functions),
                                                                std_cxx11::cref(dof_to_boundary_mapping),
                                                                std_cxx11::ref(matrix),
                                                                std_cxx11::ref(rhs_vector))),
                    scratch,
                    copy_data);
  }



  namespace
  {

    template <int dim, int spacedim>
    void
    create_hp_boundary_mass_matrix_1 (typename hp::DoFHandler<dim,spacedim>::active_cell_iterator const
                                      &cell,
                                      MatrixCreator::internal::AssemblerBoundary::Scratch const &scratch,
                                      MatrixCreator::internal::AssemblerBoundary
                                      ::CopyData<hp::DoFHandler<dim,spacedim> > &copy_data,
                                      hp::MappingCollection<dim,spacedim> const &mapping,
                                      hp::FECollection<dim,spacedim> const &fe_collection,
                                      hp::QCollection<dim-1> const &q,
                                      const typename FunctionMap<spacedim>::type &boundary_functions,
                                      Function<spacedim> const *const coefficient,
                                      std::vector<unsigned int> const &component_mapping)
    {
      const unsigned int n_components  = fe_collection.n_components();
      const unsigned int n_function_components = boundary_functions.begin()->second->n_components;
      const bool         fe_is_system  = (n_components != 1);
      const FiniteElement<dim,spacedim> &fe = cell->get_fe();
      const unsigned int dofs_per_face = fe.dofs_per_face;

      copy_data.cell = cell;
      copy_data.dofs_per_cell = fe.dofs_per_cell;
      copy_data.dofs.resize(copy_data.dofs_per_cell);
      cell->get_dof_indices (copy_data.dofs);


      UpdateFlags update_flags = UpdateFlags (update_values     |
                                              update_JxW_values |
                                              update_quadrature_points);
      hp::FEFaceValues<dim,spacedim> x_fe_values (mapping, fe_collection, q, update_flags);

      // two variables for the coefficient,
      // one for the two cases indicated in
      // the name
      std::vector<double>          coefficient_values;
      std::vector<Vector<double> > coefficient_vector_values;

      std::vector<double>          rhs_values_scalar;
      std::vector<Vector<double> > rhs_values_system;

      std::vector<types::global_dof_index> dofs_on_face_vector (dofs_per_face);

      copy_data.dofs.resize(copy_data.dofs_per_cell);
      cell->get_dof_indices (copy_data.dofs);

      // Because CopyData objects are reused and that push_back is used,
      // dof_is_on_face, cell_matrix, and cell_vector must be cleared before
      // they are reused
      copy_data.dof_is_on_face.clear();
      copy_data.cell_matrix.clear();
      copy_data.cell_vector.clear();


      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
        // check if this face is on that part of
        // the boundary we are interested in
        if (boundary_functions.find(cell->face(face)->boundary_indicator()) !=
            boundary_functions.end())
          {
            x_fe_values.reinit (cell, face);

            const FEFaceValues<dim,spacedim> &fe_values = x_fe_values.get_present_fe_values ();

            copy_data.cell_matrix.push_back(FullMatrix<double> (copy_data.dofs_per_cell,
                                                                copy_data.dofs_per_cell));
            copy_data.cell_vector.push_back(Vector<double> (copy_data.dofs_per_cell));

            if (fe_is_system)
              // FE has several components
              {
                rhs_values_system.resize (fe_values.n_quadrature_points,
                                          Vector<double>(n_function_components));
                boundary_functions.find(cell->face(face)->boundary_indicator())
                ->second->vector_value_list (fe_values.get_quadrature_points(),
                                             rhs_values_system);

                if (coefficient != 0)
                  {
                    if (coefficient->n_components==1)
                      {
                        coefficient_values.resize (fe_values.n_quadrature_points);
                        coefficient->value_list (fe_values.get_quadrature_points(),
                                                 coefficient_values);
                        for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
                          {
                            const double weight = fe_values.JxW(point);
                            for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i)
                              {
                                const double v = fe_values.shape_value(i,point);
                                for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
                                  if (fe.system_to_component_index(i).first ==
                                      fe.system_to_component_index(j).first)
                                    {
                                      const double u = fe_values.shape_value(j,point);
                                      copy_data.cell_matrix.back()(i,j)
                                      += (u * v * weight * coefficient_values[point]);
                                    }

                                copy_data.cell_vector.back()(i) += v *
                                                                   rhs_values_system[point](
                                                                     component_mapping[fe.system_to_component_index(i).first]) * weight;
                              }
                          }
                      }
                    else
                      {
                        coefficient_vector_values.resize (fe_values.n_quadrature_points,
                                                          Vector<double>(n_components));
                        coefficient->vector_value_list (fe_values.get_quadrature_points(),
                                                        coefficient_vector_values);
                        for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
                          {
                            const double weight = fe_values.JxW(point);
                            for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i)
                              {
                                const double v = fe_values.shape_value(i,point);
                                const unsigned int component_i=
                                  fe.system_to_component_index(i).first;
                                for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
                                  if (fe.system_to_component_index(j).first ==
                                      component_i)
                                    {
                                      const double u = fe_values.shape_value(j,point);
                                      copy_data.cell_matrix.back()(i,j) +=
                                        (u * v * weight * coefficient_vector_values[point](component_i));
                                    }
                                copy_data.cell_vector.back()(i) += v *
                                                                   rhs_values_system[point](component_mapping[component_i]) * weight;
                              }
                          }
                      }
                  }
                else  //      if (coefficient == 0)
                  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
                    {
                      const double weight = fe_values.JxW(point);
                      for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i)
                        {
                          const double v = fe_values.shape_value(i,point);
                          for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
                            if (fe.system_to_component_index(i).first ==
                                fe.system_to_component_index(j).first)
                              {
                                const double u = fe_values.shape_value(j,point);
                                copy_data.cell_matrix.back()(i,j) += (u * v * weight);
                              }
                          copy_data.cell_vector.back()(i) += v *
                                                             rhs_values_system[point](
                                                               fe.system_to_component_index(i).first) *
                                                             weight;
                        }
                    }
              }
            else
              // FE is a scalar one
              {
                rhs_values_scalar.resize (fe_values.n_quadrature_points);
                boundary_functions.find(cell->face(face)->boundary_indicator())
                ->second->value_list (fe_values.get_quadrature_points(), rhs_values_scalar);

                if (coefficient != 0)
                  {
                    coefficient_values.resize (fe_values.n_quadrature_points);
                    coefficient->value_list (fe_values.get_quadrature_points(),
                                             coefficient_values);
                    for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
                      {
                        const double weight = fe_values.JxW(point);
                        for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i)
                          {
                            const double v = fe_values.shape_value(i,point);
                            for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
                              {
                                const double u = fe_values.shape_value(j,point);
                                copy_data.cell_matrix.back()(i,j) += (u * v * weight *
                                                                      coefficient_values[point]);
                              }
                            copy_data.cell_vector.back()(i) += v * rhs_values_scalar[point] *weight;
                          }
                      }
                  }
                else
                  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
                    {
                      const double weight = fe_values.JxW(point);
                      for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i)
                        {
                          const double v = fe_values.shape_value(i,point);
                          for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
                            {
                              const double u = fe_values.shape_value(j,point);
                              copy_data.cell_matrix.back()(i,j) += (u * v * weight);
                            }
                          copy_data.cell_vector.back()(i) += v * rhs_values_scalar[point] * weight;
                        }
                    }
              }

            cell->face(face)->get_dof_indices (dofs_on_face_vector,
                                               cell->active_fe_index());
            // for each dof on the cell, have a
            // flag whether it is on the face
            copy_data.dof_is_on_face.push_back(std::vector<bool> (copy_data.dofs_per_cell));
            // check for each of the dofs on this cell
            // whether it is on the face
            for (unsigned int i=0; i<copy_data.dofs_per_cell; ++i)
              copy_data.dof_is_on_face.back()[i] = (std::find(dofs_on_face_vector.begin(),
                                                              dofs_on_face_vector.end(),
                                                              copy_data.dofs[i])
                                                    !=
                                                    dofs_on_face_vector.end());
          }
    }



    template <int dim,int spacedim>
    void copy_hp_boundary_mass_matrix_1(MatrixCreator::internal::AssemblerBoundary
                                        ::CopyData<hp::DoFHandler<dim,spacedim> > const &copy_data,
                                        typename FunctionMap<spacedim>::type const &boundary_functions,
                                        std::vector<types::global_dof_index> const &dof_to_boundary_mapping,
                                        SparseMatrix<double> &matrix,
                                        Vector<double> &rhs_vector)
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
      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
        {
          // check if this face is on that part of
          // the boundary we are interested in
          if (boundary_functions.find(copy_data.cell->face(face)->boundary_indicator()) !=
              boundary_functions.end())
            {
#ifdef DEBUG
              // in debug mode: compute an element in the matrix which is
              // guaranteed to belong to a boundary dof. We do this to check that the
              // entries in the cell matrix are guaranteed to be zero if the
              // respective dof is not on the boundary. Since because of
              // round-off, the actual value of the matrix entry may be
              // only close to zero, we assert that it is small relative to an element
              // which is guaranteed to be nonzero. (absolute smallness does not
              // suffice since the size of the domain scales in here)
              //
              // for this purpose we seek the diagonal of the matrix, where there
              // must be an element belonging to the boundary. we take the maximum
              // diagonal entry.
              types::global_dof_index max_element = static_cast<types::global_dof_index>(0);
              for (std::vector<types::global_dof_index>::const_iterator i=dof_to_boundary_mapping.begin();
                   i!=dof_to_boundary_mapping.end(); ++i)
                if ((*i != hp::DoFHandler<dim,spacedim>::invalid_dof_index) &&
                    (*i > max_element))
                  max_element = *i;
              Assert (max_element  == matrix.n()-1, ExcInternalError());

              double max_diag_entry = 0;
              for (unsigned int i=0; i<copy_data.dofs_per_cell; ++i)
                if (std::fabs(copy_data.cell_matrix[pos](i,i)) > max_diag_entry)
                  max_diag_entry = std::fabs(copy_data.cell_matrix[pos](i,i));
#endif

              for (unsigned int i=0; i<copy_data.dofs_per_cell; ++i)
                for (unsigned int j=0; j<copy_data.dofs_per_cell; ++j)
                  {
                    if (copy_data.dof_is_on_face[pos][i] && copy_data.dof_is_on_face[pos][j])
                      matrix.add(dof_to_boundary_mapping[copy_data.dofs[i]],
                                 dof_to_boundary_mapping[copy_data.dofs[j]],
                                 copy_data.cell_matrix[pos](i,j));
                    else
                      {
                        // assume that all shape functions that are nonzero on the boundary
                        // are also listed in the @p{dof_to_boundary} mapping. if that
                        // is not the case, then the boundary mass matrix does not
                        // make that much sense anyway, as it only contains entries for
                        // parts of the functions living on the boundary
                        //
                        // these, we may compare here for relative smallness of all
                        // entries in the local matrix which are not taken over to
                        // the global one
                        Assert (std::fabs(copy_data.cell_matrix[pos](i,j)) <= 1e-10 * max_diag_entry,
                                ExcInternalError ());
                      }
                  }

              for (unsigned int j=0; j<copy_data.dofs_per_cell; ++j)
                if (copy_data.dof_is_on_face[pos][j])
                  rhs_vector(dof_to_boundary_mapping[copy_data.dofs[j]]) += copy_data.cell_vector[pos](j);
                else
                  {
                    // compare here for relative
                    // smallness
                    Assert (std::fabs(copy_data.cell_vector[pos](j)) <= 1e-10 * max_diag_entry,
                            ExcInternalError());
                  }
              ++pos;
            }
        }
    }
  }



  template <int dim, int spacedim>
  void create_boundary_mass_matrix (const DoFHandler<dim,spacedim>     &dof,
                                    const Quadrature<dim-1>   &q,
                                    SparseMatrix<double>      &matrix,
                                    const typename FunctionMap<spacedim>::type &rhs,
                                    Vector<double>            &rhs_vector,
                                    std::vector<types::global_dof_index> &dof_to_boundary_mapping,
                                    const Function<spacedim> *const a,
                                    std::vector<unsigned int> component_mapping)
  {
    create_boundary_mass_matrix(StaticMappingQ1<dim,spacedim>::mapping, dof, q,
                                matrix,rhs, rhs_vector, dof_to_boundary_mapping, a, component_mapping);
  }



  template <int dim, int spacedim>
  void
  create_boundary_mass_matrix (const hp::MappingCollection<dim,spacedim>        &mapping,
                               const hp::DoFHandler<dim,spacedim>     &dof,
                               const hp::QCollection<dim-1>   &q,
                               SparseMatrix<double>      &matrix,
                               const typename FunctionMap<spacedim>::type         &boundary_functions,
                               Vector<double>            &rhs_vector,
                               std::vector<types::global_dof_index> &dof_to_boundary_mapping,
                               const Function<spacedim> *const coefficient,
                               std::vector<unsigned int> component_mapping)
  {
    // what would that be in 1d? the
    // identity matrix on the boundary
    // dofs?
    if (dim == 1)
      {
        Assert (false, ExcNotImplemented());
        return;
      }

    const hp::FECollection<dim,spacedim> &fe_collection = dof.get_fe();
    const unsigned int n_components  = fe_collection.n_components();

    Assert (matrix.n() == dof.n_boundary_dofs(boundary_functions),
            ExcInternalError());
    Assert (matrix.n() == matrix.m(), ExcInternalError());
    Assert (matrix.n() == rhs_vector.size(), ExcInternalError());
    Assert (boundary_functions.size() != 0, ExcInternalError());
    Assert (dof_to_boundary_mapping.size() == dof.n_dofs(),
            ExcInternalError());
    Assert (coefficient ==0 ||
            coefficient->n_components==1 ||
            coefficient->n_components==n_components, ExcComponentMismatch());

    if (component_mapping.size() == 0)
      {
        AssertDimension (n_components, boundary_functions.begin()->second->n_components);
        for (unsigned int i=0; i<n_components; ++i)
          component_mapping.push_back(i);
      }
    else
      AssertDimension (n_components, component_mapping.size());

    MatrixCreator::internal::AssemblerBoundary::Scratch scratch;
    MatrixCreator::internal::AssemblerBoundary::CopyData<hp::DoFHandler<dim,spacedim> > copy_data;

    WorkStream::run(dof.begin_active(),dof.end(),
                    static_cast<std_cxx11::function<void (typename hp::DoFHandler<dim,spacedim>::active_cell_iterator
                                                          const &,MatrixCreator::internal::AssemblerBoundary::Scratch const &,
                                                          MatrixCreator::internal::AssemblerBoundary::CopyData<hp::DoFHandler<dim,spacedim> > &)> >
                    (std_cxx11::bind( &create_hp_boundary_mass_matrix_1<dim,spacedim>,std_cxx11::_1,std_cxx11::_2,
                                      std_cxx11::_3,
                                      std_cxx11::cref(mapping),std_cxx11::cref(fe_collection),std_cxx11::cref(q),
                                      std_cxx11::cref(boundary_functions),coefficient,
                                      std_cxx11::cref(component_mapping))),
                    static_cast<std_cxx11::function<void (MatrixCreator::internal::AssemblerBoundary
                                                          ::CopyData<hp::DoFHandler<dim,spacedim> > const &)> > (
                      std_cxx11::bind( &copy_hp_boundary_mass_matrix_1<dim,spacedim>,
                                       std_cxx11::_1,
                                       std_cxx11::cref(boundary_functions),
                                       std_cxx11::cref(dof_to_boundary_mapping),
                                       std_cxx11::ref(matrix),
                                       std_cxx11::ref(rhs_vector))),
                    scratch,
                    copy_data);
  }




  template <int dim, int spacedim>
  void create_boundary_mass_matrix (const hp::DoFHandler<dim,spacedim>     &dof,
                                    const hp::QCollection<dim-1>   &q,
                                    SparseMatrix<double>      &matrix,
                                    const typename FunctionMap<spacedim>::type &rhs,
                                    Vector<double>            &rhs_vector,
                                    std::vector<types::global_dof_index> &dof_to_boundary_mapping,
                                    const Function<spacedim> *const a,
                                    std::vector<unsigned int> component_mapping)
  {
    create_boundary_mass_matrix(hp::StaticMappingQ1<dim,spacedim>::mapping_collection, dof, q,
                                matrix,rhs, rhs_vector, dof_to_boundary_mapping, a, component_mapping);
  }



  template <int dim, int spacedim>
  void create_laplace_matrix (const Mapping<dim, spacedim>       &mapping,
                              const DoFHandler<dim,spacedim>    &dof,
                              const Quadrature<dim>    &q,
                              SparseMatrix<double>     &matrix,
                              const Function<spacedim> *const coefficient,
                              const ConstraintMatrix   &constraints)
  {
    Assert (matrix.m() == dof.n_dofs(),
            ExcDimensionMismatch (matrix.m(), dof.n_dofs()));
    Assert (matrix.n() == dof.n_dofs(),
            ExcDimensionMismatch (matrix.n(), dof.n_dofs()));

    hp::FECollection<dim,spacedim>      fe_collection (dof.get_fe());
    hp::QCollection<dim>                q_collection (q);
    hp::MappingCollection<dim,spacedim> mapping_collection (mapping);
    MatrixCreator::internal::AssemblerData::Scratch<dim, spacedim>
    assembler_data (fe_collection,
                    update_gradients  | update_JxW_values |
                    (coefficient != 0 ? update_quadrature_points : UpdateFlags(0)),
                    coefficient, /*rhs_function=*/0,
                    q_collection, mapping_collection);
    MatrixCreator::internal::AssemblerData::CopyData copy_data;
    copy_data.cell_matrix.reinit (assembler_data.fe_collection.max_dofs_per_cell(),
                                  assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.cell_rhs.reinit (assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.dof_indices.resize (assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.constraints = &constraints;

    WorkStream::run (dof.begin_active(),
                     static_cast<typename DoFHandler<dim,spacedim>::active_cell_iterator>(dof.end()),
                     &MatrixCreator::internal::laplace_assembler<dim, spacedim, typename DoFHandler<dim,spacedim>::active_cell_iterator>,
                     std_cxx11::bind (&MatrixCreator::internal::
                                      copy_local_to_global<SparseMatrix<double>, Vector<double> >,
                                      std_cxx11::_1,
                                      &matrix,
                                      (Vector<double> *)NULL),
                     assembler_data,
                     copy_data);
  }



  template <int dim, int spacedim>
  void create_laplace_matrix (const DoFHandler<dim,spacedim>    &dof,
                              const Quadrature<dim>    &q,
                              SparseMatrix<double>     &matrix,
                              const Function<spacedim> *const coefficient,
                              const ConstraintMatrix   &constraints)
  {
    create_laplace_matrix(StaticMappingQ1<dim,spacedim>::mapping, dof, q, matrix, coefficient);
  }



  template <int dim, int spacedim>
  void create_laplace_matrix (const Mapping<dim, spacedim>       &mapping,
                              const DoFHandler<dim,spacedim>    &dof,
                              const Quadrature<dim>    &q,
                              SparseMatrix<double>     &matrix,
                              const Function<spacedim>      &rhs,
                              Vector<double>           &rhs_vector,
                              const Function<spacedim> *const coefficient,
                              const ConstraintMatrix   &constraints)
  {
    Assert (matrix.m() == dof.n_dofs(),
            ExcDimensionMismatch (matrix.m(), dof.n_dofs()));
    Assert (matrix.n() == dof.n_dofs(),
            ExcDimensionMismatch (matrix.n(), dof.n_dofs()));

    hp::FECollection<dim,spacedim>      fe_collection (dof.get_fe());
    hp::QCollection<dim>                q_collection (q);
    hp::MappingCollection<dim,spacedim> mapping_collection (mapping);
    MatrixCreator::internal::AssemblerData::Scratch<dim, spacedim>
    assembler_data (fe_collection,
                    update_gradients  | update_values |
                    update_JxW_values | update_quadrature_points,
                    coefficient, &rhs,
                    q_collection, mapping_collection);
    MatrixCreator::internal::AssemblerData::CopyData copy_data;
    copy_data.cell_matrix.reinit (assembler_data.fe_collection.max_dofs_per_cell(),
                                  assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.cell_rhs.reinit (assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.dof_indices.resize (assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.constraints = &constraints;

    WorkStream::run (dof.begin_active(),
                     static_cast<typename DoFHandler<dim,spacedim>::active_cell_iterator>(dof.end()),
                     &MatrixCreator::internal::laplace_assembler<dim, spacedim, typename DoFHandler<dim,spacedim>::active_cell_iterator>,
                     std_cxx11::bind (&MatrixCreator::internal::
                                      copy_local_to_global<SparseMatrix<double>, Vector<double> >,
                                      std_cxx11::_1,
                                      &matrix,
                                      &rhs_vector),
                     assembler_data,
                     copy_data);
  }



  template <int dim, int spacedim>
  void create_laplace_matrix (const DoFHandler<dim,spacedim>    &dof,
                              const Quadrature<dim>    &q,
                              SparseMatrix<double>     &matrix,
                              const Function<spacedim>      &rhs,
                              Vector<double>           &rhs_vector,
                              const Function<spacedim> *const coefficient,
                              const ConstraintMatrix   &constraints)
  {
    create_laplace_matrix(StaticMappingQ1<dim,spacedim>::mapping, dof, q,
                          matrix, rhs, rhs_vector, coefficient);
  }



  template <int dim, int spacedim>
  void create_laplace_matrix (const hp::MappingCollection<dim,spacedim>       &mapping,
                              const hp::DoFHandler<dim,spacedim>    &dof,
                              const hp::QCollection<dim>    &q,
                              SparseMatrix<double>     &matrix,
                              const Function<spacedim> *const coefficient,
                              const ConstraintMatrix   &constraints)
  {
    Assert (matrix.m() == dof.n_dofs(),
            ExcDimensionMismatch (matrix.m(), dof.n_dofs()));
    Assert (matrix.n() == dof.n_dofs(),
            ExcDimensionMismatch (matrix.n(), dof.n_dofs()));

    MatrixCreator::internal::AssemblerData::Scratch<dim, spacedim>
    assembler_data (dof.get_fe(),
                    update_gradients  | update_JxW_values |
                    (coefficient != 0 ? update_quadrature_points : UpdateFlags(0)),
                    coefficient, /*rhs_function=*/0,
                    q, mapping);
    MatrixCreator::internal::AssemblerData::CopyData copy_data;
    copy_data.cell_matrix.reinit (assembler_data.fe_collection.max_dofs_per_cell(),
                                  assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.cell_rhs.reinit (assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.dof_indices.resize (assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.constraints = &constraints;

    WorkStream::run (dof.begin_active(),
                     static_cast<typename hp::DoFHandler<dim,spacedim>::active_cell_iterator>(dof.end()),
                     &MatrixCreator::internal::laplace_assembler<dim, spacedim, typename hp::DoFHandler<dim,spacedim>::active_cell_iterator>,
                     std_cxx11::bind (&MatrixCreator::internal::
                                      copy_local_to_global<SparseMatrix<double>, Vector<double> >,
                                      std_cxx11::_1,
                                      &matrix,
                                      (Vector<double> *)0),
                     assembler_data,
                     copy_data);
  }



  template <int dim, int spacedim>
  void create_laplace_matrix (const hp::DoFHandler<dim,spacedim>    &dof,
                              const hp::QCollection<dim>    &q,
                              SparseMatrix<double>     &matrix,
                              const Function<spacedim> *const coefficient,
                              const ConstraintMatrix   &constraints)
  {
    create_laplace_matrix(hp::StaticMappingQ1<dim,spacedim>::mapping_collection, dof, q, matrix, coefficient);
  }



  template <int dim, int spacedim>
  void create_laplace_matrix (const hp::MappingCollection<dim,spacedim>       &mapping,
                              const hp::DoFHandler<dim,spacedim>    &dof,
                              const hp::QCollection<dim>    &q,
                              SparseMatrix<double>     &matrix,
                              const Function<spacedim>      &rhs,
                              Vector<double>           &rhs_vector,
                              const Function<spacedim> *const coefficient,
                              const ConstraintMatrix   &constraints)
  {
    Assert (matrix.m() == dof.n_dofs(),
            ExcDimensionMismatch (matrix.m(), dof.n_dofs()));
    Assert (matrix.n() == dof.n_dofs(),
            ExcDimensionMismatch (matrix.n(), dof.n_dofs()));

    MatrixCreator::internal::AssemblerData::Scratch<dim, spacedim>
    assembler_data (dof.get_fe(),
                    update_gradients  | update_values |
                    update_JxW_values | update_quadrature_points,
                    coefficient, &rhs,
                    q, mapping);
    MatrixCreator::internal::AssemblerData::CopyData copy_data;
    copy_data.cell_matrix.reinit (assembler_data.fe_collection.max_dofs_per_cell(),
                                  assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.cell_rhs.reinit (assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.dof_indices.resize (assembler_data.fe_collection.max_dofs_per_cell());
    copy_data.constraints = &constraints;

    WorkStream::run (dof.begin_active(),
                     static_cast<typename hp::DoFHandler<dim,spacedim>::active_cell_iterator>(dof.end()),
                     &MatrixCreator::internal::laplace_assembler<dim, spacedim, typename hp::DoFHandler<dim,spacedim>::active_cell_iterator>,
                     std_cxx11::bind (&MatrixCreator::internal::
                                      copy_local_to_global<SparseMatrix<double>, Vector<double> >,
                                      std_cxx11::_1,
                                      &matrix,
                                      &rhs_vector),
                     assembler_data,
                     copy_data);
  }



  template <int dim, int spacedim>
  void create_laplace_matrix (const hp::DoFHandler<dim,spacedim>    &dof,
                              const hp::QCollection<dim>    &q,
                              SparseMatrix<double>     &matrix,
                              const Function<spacedim>      &rhs,
                              Vector<double>           &rhs_vector,
                              const Function<spacedim> *const coefficient,
                              const ConstraintMatrix   &constraints)
  {
    create_laplace_matrix(hp::StaticMappingQ1<dim,spacedim>::mapping_collection, dof, q,
                          matrix, rhs, rhs_vector, coefficient);
  }

}  // namespace MatrixCreator


namespace MatrixTools
{
  namespace
  {
    template <typename Iterator>
    bool column_less_than(const typename Iterator::value_type p,
                          const unsigned int column)
    {
      return (p.column() < column);
    }
  }

//TODO:[WB] I don't think that the optimized storage of diagonals is needed (GK)
  template <typename number>
  void
  apply_boundary_values (const std::map<types::global_dof_index,double> &boundary_values,
                         SparseMatrix<number>  &matrix,
                         Vector<number>   &solution,
                         Vector<number>   &right_hand_side,
                         const bool        eliminate_columns)
  {
    Assert (matrix.n() == right_hand_side.size(),
            ExcDimensionMismatch(matrix.n(), right_hand_side.size()));
    Assert (matrix.n() == solution.size(),
            ExcDimensionMismatch(matrix.n(), solution.size()));
    Assert (matrix.n() == matrix.m(),
            ExcDimensionMismatch(matrix.n(), matrix.m()));

    // if no boundary values are to be applied
    // simply return
    if (boundary_values.size() == 0)
      return;


    const types::global_dof_index n_dofs = matrix.m();

    // if a diagonal entry is zero
    // later, then we use another
    // number instead. take it to be
    // the first nonzero diagonal
    // element of the matrix, or 1 if
    // there is no such thing
    number first_nonzero_diagonal_entry = 1;
    for (unsigned int i=0; i<n_dofs; ++i)
      if (matrix.diag_element(i) != 0)
        {
          first_nonzero_diagonal_entry = matrix.diag_element(i);
          break;
        }


    std::map<types::global_dof_index,double>::const_iterator dof  = boundary_values.begin(),
                                                             endd = boundary_values.end();
    for (; dof != endd; ++dof)
      {
        Assert (dof->first < n_dofs, ExcInternalError());

        const types::global_dof_index dof_number = dof->first;
        // for each boundary dof:

        // set entries of this line to zero except for the diagonal
        // entry
        for (typename SparseMatrix<number>::iterator
             p = matrix.begin(dof_number);
             p != matrix.end(dof_number); ++p)
          if (p->column() != dof_number)
            p->value() = 0.;

        // set right hand side to
        // wanted value: if main diagonal
        // entry nonzero, don't touch it
        // and scale rhs accordingly. If
        // zero, take the first main
        // diagonal entry we can find, or
        // one if no nonzero main diagonal
        // element exists. Normally, however,
        // the main diagonal entry should
        // not be zero.
        //
        // store the new rhs entry to make
        // the gauss step more efficient
        number new_rhs;
        if (matrix.diag_element(dof_number) != 0.0)
          {
            new_rhs = dof->second * matrix.diag_element(dof_number);
            right_hand_side(dof_number) = new_rhs;
          }
        else
          {
            matrix.set (dof_number, dof_number,
                        first_nonzero_diagonal_entry);
            new_rhs = dof->second * first_nonzero_diagonal_entry;
            right_hand_side(dof_number) = new_rhs;
          }


        // if the user wants to have
        // the symmetry of the matrix
        // preserved, and if the
        // sparsity pattern is
        // symmetric, then do a Gauss
        // elimination step with the
        // present row
        if (eliminate_columns)
          {
            // store the only nonzero entry
            // of this line for the Gauss
            // elimination step
            const number diagonal_entry = matrix.diag_element(dof_number);

            // we have to loop over all rows of the matrix which have
            // a nonzero entry in the column which we work in
            // presently. if the sparsity pattern is symmetric, then
            // we can get the positions of these rows cheaply by
            // looking at the nonzero column numbers of the present
            // row. we need not look at the first entry of each row,
            // since that is the diagonal element and thus the present
            // row
            for (typename SparseMatrix<number>::iterator
                 q = matrix.begin(dof_number)+1;
                 q != matrix.end(dof_number); ++q)
              {
                const types::global_dof_index row = q->column();

                // find the position of
                // element
                // (row,dof_number)
                bool (*comp)(const typename SparseMatrix<number>::iterator::value_type p,
                             const unsigned int column)
                  = &column_less_than<typename SparseMatrix<number>::iterator>;
                const typename SparseMatrix<number>::iterator
                p = Utilities::lower_bound(matrix.begin(row)+1,
                                           matrix.end(row),
                                           dof_number,
                                           comp);

                // check whether this line has an entry in the
                // regarding column (check for ==dof_number and !=
                // next_row, since if row==dof_number-1, *p is a
                // past-the-end pointer but points to dof_number
                // anyway...)
                //
                // there should be such an entry! we know this because
                // we have assumed that the sparsity pattern is
                // symmetric and we only walk over those rows for
                // which the current row has a column entry
                Assert ((p != matrix.end(row))
                        &&
                        (p->column() == dof_number),
                        ExcMessage("This function is trying to access an element of the "
                                   "matrix that doesn't seem to exist. Are you using a "
                                   "nonsymmetric sparsity pattern? If so, you are not "
                                   "allowed to set the eliminate_column argument of this "
                                   "function, see the documentation."));

                // correct right hand side
                right_hand_side(row) -= p->value() /
                                        diagonal_entry * new_rhs;

                // set matrix entry to zero
                p->value() = 0.;
              }
          }

        // preset solution vector
        solution(dof_number) = dof->second;
      }
  }



  template <typename number>
  void
  apply_boundary_values (const std::map<types::global_dof_index,double> &boundary_values,
                         BlockSparseMatrix<number>  &matrix,
                         BlockVector<number>   &solution,
                         BlockVector<number>   &right_hand_side,
                         const bool             eliminate_columns)
  {
    const unsigned int blocks = matrix.n_block_rows();

    Assert (matrix.n() == right_hand_side.size(),
            ExcDimensionMismatch(matrix.n(), right_hand_side.size()));
    Assert (matrix.n() == solution.size(),
            ExcDimensionMismatch(matrix.n(), solution.size()));
    Assert (matrix.n_block_rows() == matrix.n_block_cols(),
            ExcNotQuadratic());
    Assert (matrix.get_sparsity_pattern().get_row_indices() ==
            matrix.get_sparsity_pattern().get_column_indices(),
            ExcNotQuadratic());
    Assert (matrix.get_sparsity_pattern().get_column_indices() ==
            solution.get_block_indices (),
            ExcBlocksDontMatch ());
    Assert (matrix.get_sparsity_pattern().get_row_indices() ==
            right_hand_side.get_block_indices (),
            ExcBlocksDontMatch ());

    // if no boundary values are to be applied
    // simply return
    if (boundary_values.size() == 0)
      return;


    const types::global_dof_index n_dofs = matrix.m();

    // if a diagonal entry is zero
    // later, then we use another
    // number instead. take it to be
    // the first nonzero diagonal
    // element of the matrix, or 1 if
    // there is no such thing
    number first_nonzero_diagonal_entry = 0;
    for (unsigned int diag_block=0; diag_block<blocks; ++diag_block)
      {
        for (unsigned int i=0; i<matrix.block(diag_block,diag_block).n(); ++i)
          if (matrix.block(diag_block,diag_block).diag_element(i) != 0)
            {
              first_nonzero_diagonal_entry
                = matrix.block(diag_block,diag_block).diag_element(i);
              break;
            }
        // check whether we have found
        // something in the present
        // block
        if (first_nonzero_diagonal_entry != 0)
          break;
      }
    // nothing found on all diagonal
    // blocks? if so, use 1.0 instead
    if (first_nonzero_diagonal_entry == 0)
      first_nonzero_diagonal_entry = 1;


    std::map<types::global_dof_index,double>::const_iterator dof  = boundary_values.begin(),
                                                             endd = boundary_values.end();
    const BlockSparsityPattern &
    sparsity_pattern = matrix.get_sparsity_pattern();

    // pointer to the mapping between
    // global and block indices. since
    // the row and column mappings are
    // equal, store a pointer on only
    // one of them
    const BlockIndices &
    index_mapping = sparsity_pattern.get_column_indices();

    // now loop over all boundary dofs
    for (; dof != endd; ++dof)
      {
        Assert (dof->first < n_dofs, ExcInternalError());

        // get global index and index
        // in the block in which this
        // dof is located
        const types::global_dof_index dof_number = dof->first;
        const std::pair<unsigned int,types::global_dof_index>
        block_index = index_mapping.global_to_local (dof_number);

        // for each boundary dof:

        // set entries of this line
        // to zero except for the diagonal
        // entry. Note that the diagonal
        // entry is always the first one
        // in a row for square matrices
        for (unsigned int block_col=0; block_col<blocks; ++block_col)
          for (typename SparseMatrix<number>::iterator
               p = (block_col == block_index.first ?
                    matrix.block(block_index.first,block_col).begin(block_index.second) + 1 :
                    matrix.block(block_index.first,block_col).begin(block_index.second));
               p != matrix.block(block_index.first,block_col).end(block_index.second);
               ++p)
            p->value() = 0;

        // set right hand side to
        // wanted value: if main diagonal
        // entry nonzero, don't touch it
        // and scale rhs accordingly. If
        // zero, take the first main
        // diagonal entry we can find, or
        // one if no nonzero main diagonal
        // element exists. Normally, however,
        // the main diagonal entry should
        // not be zero.
        //
        // store the new rhs entry to make
        // the gauss step more efficient
        number new_rhs;
        if (matrix.block(block_index.first, block_index.first)
            .diag_element(block_index.second) != 0.0)
          new_rhs = dof->second *
                    matrix.block(block_index.first, block_index.first)
                    .diag_element(block_index.second);
        else
          {
            matrix.block(block_index.first, block_index.first)
            .diag_element(block_index.second)
              = first_nonzero_diagonal_entry;
            new_rhs = dof->second * first_nonzero_diagonal_entry;
          }
        right_hand_side.block(block_index.first)(block_index.second)
          = new_rhs;


        // if the user wants to have
        // the symmetry of the matrix
        // preserved, and if the
        // sparsity pattern is
        // symmetric, then do a Gauss
        // elimination step with the
        // present row. this is a
        // little more complicated for
        // block matrices.
        if (eliminate_columns)
          {
            // store the only nonzero entry
            // of this line for the Gauss
            // elimination step
            const number diagonal_entry
              = matrix.block(block_index.first,block_index.first)
                .diag_element(block_index.second);

            // we have to loop over all
            // rows of the matrix which
            // have a nonzero entry in
            // the column which we work
            // in presently. if the
            // sparsity pattern is
            // symmetric, then we can
            // get the positions of
            // these rows cheaply by
            // looking at the nonzero
            // column numbers of the
            // present row.
            //
            // note that if we check
            // whether row @p{row} in
            // block (r,c) is non-zero,
            // then we have to check
            // for the existence of
            // column @p{row} in block
            // (c,r), i.e. of the
            // transpose block
            for (unsigned int block_row=0; block_row<blocks; ++block_row)
              {
                // get pointers to the sparsity patterns of this block and of
                // the transpose one
                const SparsityPattern &this_sparsity
                  = sparsity_pattern.block (block_row, block_index.first);

                SparseMatrix<number> &this_matrix
                  = matrix.block(block_row, block_index.first);
                SparseMatrix<number> &transpose_matrix
                  = matrix.block(block_index.first, block_row);

                // traverse the row of the transpose block to find the
                // interesting rows in the present block.  don't use the
                // diagonal element of the diagonal block
                for (typename SparseMatrix<number>::iterator
                     q = (block_index.first == block_row ?
                          transpose_matrix.begin(block_index.second)+1 :
                          transpose_matrix.begin(block_index.second));
                     q != transpose_matrix.end(block_index.second);
                     ++q)
                  {
                    // get the number of the column in this row in which a
                    // nonzero entry is. this is also the row of the transpose
                    // block which has an entry in the interesting row
                    const types::global_dof_index row = q->column();

                    // find the position of element (row,dof_number) in this
                    // block (not in the transpose one). note that we have to
                    // take care of special cases with square sub-matrices
                    bool (*comp)(typename SparseMatrix<number>::iterator::value_type p,
                                 const unsigned int column)
                      = &column_less_than<typename SparseMatrix<number>::iterator>;

                    typename SparseMatrix<number>::iterator p = this_matrix.end();

                    if (this_sparsity.n_rows() == this_sparsity.n_cols())
                      {
                        if (this_matrix.begin(row)->column()
                            ==
                            block_index.second)
                          p = this_matrix.begin(row);
                        else
                          p = Utilities::lower_bound(this_matrix.begin(row)+1,
                                                     this_matrix.end(row),
                                                     block_index.second,
                                                     comp);
                      }
                    else
                      p = Utilities::lower_bound(this_matrix.begin(row),
                                                 this_matrix.end(row),
                                                 block_index.second,
                                                 comp);

                    // check whether this line has an entry in the
                    // regarding column (check for ==dof_number and !=
                    // next_row, since if row==dof_number-1, *p is a
                    // past-the-end pointer but points to dof_number
                    // anyway...)
                    //
                    // there should be such an entry! we know this because
                    // we have assumed that the sparsity pattern is
                    // symmetric and we only walk over those rows for
                    // which the current row has a column entry
                    Assert ((p->column() == block_index.second) &&
                            (p != this_matrix.end(row)),
                            ExcInternalError());

                    // correct right hand side
                    right_hand_side.block(block_row)(row)
                    -= p->value() /
                       diagonal_entry * new_rhs;

                    // set matrix entry to zero
                    p->value() = 0.;
                  }
              }
          }

        // preset solution vector
        solution.block(block_index.first)(block_index.second) = dof->second;
      }
  }



#ifdef DEAL_II_WITH_PETSC

  namespace internal
  {
    namespace PETScWrappers
    {
      template <typename PETScMatrix, typename PETScVector>
      void
      apply_boundary_values (const std::map<types::global_dof_index,double> &boundary_values,
                             PETScMatrix      &matrix,
                             PETScVector      &solution,
                             PETScVector      &right_hand_side,
                             const bool        eliminate_columns)
      {
        Assert (eliminate_columns == false, ExcNotImplemented());

        Assert (matrix.n() == right_hand_side.size(),
                ExcDimensionMismatch(matrix.n(), right_hand_side.size()));
        Assert (matrix.n() == solution.size(),
                ExcDimensionMismatch(matrix.n(), solution.size()));

        // if no boundary values are to be applied, then
        // jump straight to the compress() calls that we still have
        // to perform because they are collective operations
        if (boundary_values.size() > 0)
          {
            const std::pair<types::global_dof_index, types::global_dof_index> local_range
              = matrix.local_range();
            Assert (local_range == right_hand_side.local_range(),
                    ExcInternalError());
            Assert (local_range == solution.local_range(),
                    ExcInternalError());

            // determine the first nonzero diagonal
            // entry from within the part of the
            // matrix that we can see. if we can't
            // find such an entry, take one
            PetscScalar average_nonzero_diagonal_entry = 1;
            for (types::global_dof_index i=local_range.first; i<local_range.second; ++i)
              if (matrix.diag_element(i) != 0)
                {
                  average_nonzero_diagonal_entry = std::fabs(matrix.diag_element(i));
                  break;
                }

            // figure out which rows of the matrix we
            // have to eliminate on this processor
            std::vector<types::global_dof_index> constrained_rows;
            for (std::map<types::global_dof_index,double>::const_iterator
                 dof  = boundary_values.begin();
                 dof != boundary_values.end();
                 ++dof)
              if ((dof->first >= local_range.first) &&
                  (dof->first < local_range.second))
                constrained_rows.push_back (dof->first);

            // then eliminate these rows and set
            // their diagonal entry to what we have
            // determined above. note that for petsc
            // matrices interleaving read with write
            // operations is very expensive. thus, we
            // here always replace the diagonal
            // element, rather than first checking
            // whether it is nonzero and in that case
            // preserving it. this is different from
            // the case of deal.II sparse matrices
            // treated in the other functions.
            matrix.clear_rows (constrained_rows, average_nonzero_diagonal_entry);

            std::vector<types::global_dof_index> indices;
            std::vector<PetscScalar>  solution_values;
            for (std::map<types::global_dof_index,double>::const_iterator
                 dof  = boundary_values.begin();
                 dof != boundary_values.end();
                 ++dof)
              if ((dof->first >= local_range.first) &&
                  (dof->first < local_range.second))
                {
                  indices.push_back (dof->first);
                  solution_values.push_back (dof->second);
                }
            solution.set (indices, solution_values);

            // now also set appropriate values for
            // the rhs
            for (unsigned int i=0; i<solution_values.size(); ++i)
              solution_values[i] *= average_nonzero_diagonal_entry;

            right_hand_side.set (indices, solution_values);
          }
        else
          {
            // clear_rows() is a collective operation so we still have to call
            // it:
            std::vector<types::global_dof_index> constrained_rows;
            matrix.clear_rows (constrained_rows, 1.);
          }

        // clean up
        matrix.compress ();
        solution.compress (VectorOperation::insert);
        right_hand_side.compress (VectorOperation::insert);
      }
    }
  }



  void
  apply_boundary_values (const std::map<types::global_dof_index,double> &boundary_values,
                         PETScWrappers::SparseMatrix   &matrix,
                         PETScWrappers::Vector   &solution,
                         PETScWrappers::Vector   &right_hand_side,
                         const bool        eliminate_columns)
  {
    // simply redirect to the generic function
    // used for both petsc matrix types
    internal::PETScWrappers::apply_boundary_values (boundary_values, matrix, solution,
                                                    right_hand_side, eliminate_columns);
  }



  void
  apply_boundary_values (const std::map<types::global_dof_index,double> &boundary_values,
                         PETScWrappers::MPI::SparseMatrix   &matrix,
                         PETScWrappers::MPI::Vector   &solution,
                         PETScWrappers::MPI::Vector   &right_hand_side,
                         const bool        eliminate_columns)
  {
    // simply redirect to the generic function
    // used for both petsc matrix types
    internal::PETScWrappers::apply_boundary_values (boundary_values, matrix, solution,
                                                    right_hand_side, eliminate_columns);

    // compress the matrix once we're done
    matrix.compress ();
  }


  void
  apply_boundary_values (const std::map<types::global_dof_index,double>  &boundary_values,
                         PETScWrappers::MPI::BlockSparseMatrix &matrix,
                         PETScWrappers::MPI::BlockVector        &solution,
                         PETScWrappers::MPI::BlockVector        &right_hand_side,
                         const bool                            eliminate_columns)
  {
    Assert (matrix.n() == right_hand_side.size(),
            ExcDimensionMismatch(matrix.n(), right_hand_side.size()));
    Assert (matrix.n() == solution.size(),
            ExcDimensionMismatch(matrix.n(), solution.size()));
    Assert (matrix.n_block_rows() == matrix.n_block_cols(),
            ExcNotQuadratic());

    const unsigned int n_blocks = matrix.n_block_rows();

    matrix.compress();

    // We need to find the subdivision
    // into blocks for the boundary values.
    // To this end, generate a vector of
    // maps with the respective indices.
    std::vector<std::map<dealii::types::global_dof_index,double> > block_boundary_values(n_blocks);
    {
      int block = 0;
      dealii::types::global_dof_index offset = 0;
      for (std::map<types::global_dof_index,double>::const_iterator
           dof  = boundary_values.begin();
           dof != boundary_values.end();
           ++dof)
        {
          if (dof->first >= matrix.block(block,0).m() + offset)
            {
              offset += matrix.block(block,0).m();
              block++;
            }
          const types::global_dof_index index = dof->first - offset;
          block_boundary_values[block].insert(std::pair<types::global_dof_index, double> (index,dof->second));
        }
    }

    // Now call the non-block variants on
    // the diagonal subblocks and the
    // solution/rhs.
    for (unsigned int block=0; block<n_blocks; ++block)
      internal::PETScWrappers::apply_boundary_values(block_boundary_values[block],
                                                     matrix.block(block,block),
                                                     solution.block(block),
                                                     right_hand_side.block(block),
                                                     eliminate_columns);

    // Finally, we need to do something
    // about the off-diagonal matrices. This
    // is luckily not difficult. Just clear
    // the whole row.
    for (unsigned int block_m=0; block_m<n_blocks; ++block_m)
      {
        const std::pair<types::global_dof_index, types::global_dof_index> local_range
          = matrix.block(block_m,0).local_range();

        std::vector<types::global_dof_index> constrained_rows;
        for (std::map<types::global_dof_index,double>::const_iterator
             dof  = block_boundary_values[block_m].begin();
             dof != block_boundary_values[block_m].end();
             ++dof)
          if ((dof->first >= local_range.first) &&
              (dof->first < local_range.second))
            constrained_rows.push_back (dof->first);

        for (unsigned int block_n=0; block_n<n_blocks; ++block_n)
          if (block_m != block_n)
            matrix.block(block_m,block_n).clear_rows(constrained_rows);
      }
  }

#endif



#ifdef DEAL_II_WITH_TRILINOS

  namespace internal
  {
    namespace TrilinosWrappers
    {
      template <typename TrilinosMatrix, typename TrilinosVector>
      void
      apply_boundary_values (const std::map<types::global_dof_index,double> &boundary_values,
                             TrilinosMatrix      &matrix,
                             TrilinosVector      &solution,
                             TrilinosVector      &right_hand_side,
                             const bool           eliminate_columns)
      {
        Assert (eliminate_columns == false, ExcNotImplemented());

        Assert (matrix.n() == right_hand_side.size(),
                ExcDimensionMismatch(matrix.n(), right_hand_side.size()));
        Assert (matrix.n() == solution.size(),
                ExcDimensionMismatch(matrix.m(), solution.size()));

        // if no boundary values are to be applied, then
        // jump straight to the compress() calls that we still have
        // to perform because they are collective operations
        if (boundary_values.size() > 0)
          {
            const std::pair<types::global_dof_index, types::global_dof_index> local_range
              = matrix.local_range();
            Assert (local_range == right_hand_side.local_range(),
                    ExcInternalError());
            Assert (local_range == solution.local_range(),
                    ExcInternalError());

            // we have to read and write from this
            // matrix (in this order). this will only
            // work if we compress the matrix first,
            // done here
            matrix.compress ();

            // determine the first nonzero diagonal
            // entry from within the part of the
            // matrix that we can see. if we can't
            // find such an entry, take one
            TrilinosScalar average_nonzero_diagonal_entry = 1;
            for (types::global_dof_index i=local_range.first; i<local_range.second; ++i)
              if (matrix.diag_element(i) != 0)
                {
                  average_nonzero_diagonal_entry = std::fabs(matrix.diag_element(i));
                  break;
                }

            // figure out which rows of the matrix we
            // have to eliminate on this processor
            std::vector<types::global_dof_index> constrained_rows;
            for (std::map<types::global_dof_index,double>::const_iterator
                 dof  = boundary_values.begin();
                 dof != boundary_values.end();
                 ++dof)
              if ((dof->first >= local_range.first) &&
                  (dof->first < local_range.second))
                constrained_rows.push_back (dof->first);

            // then eliminate these rows and
            // set their diagonal entry to
            // what we have determined
            // above. if the value already is
            // nonzero, it will be preserved,
            // in accordance with the basic
            // matrix classes in deal.II.
            matrix.clear_rows (constrained_rows, average_nonzero_diagonal_entry);

            std::vector<types::global_dof_index> indices;
            std::vector<TrilinosScalar>  solution_values;
            for (std::map<types::global_dof_index,double>::const_iterator
                 dof  = boundary_values.begin();
                 dof != boundary_values.end();
                 ++dof)
              if ((dof->first >= local_range.first) &&
                  (dof->first < local_range.second))
                {
                  indices.push_back (dof->first);
                  solution_values.push_back (dof->second);
                }
            solution.set (indices, solution_values);

            // now also set appropriate
            // values for the rhs
            for (unsigned int i=0; i<solution_values.size(); ++i)
              solution_values[i] *= matrix.diag_element(indices[i]);

            right_hand_side.set (indices, solution_values);
          }
        else
          {
            // clear_rows() is a collective operation so we still have to call
            // it:
            std::vector<types::global_dof_index> constrained_rows;
            matrix.clear_rows (constrained_rows, 1.);
          }

        // clean up
        matrix.compress ();
        solution.compress (VectorOperation::insert);
        right_hand_side.compress (VectorOperation::insert);
      }



      template <typename TrilinosMatrix, typename TrilinosBlockVector>
      void
      apply_block_boundary_values (const std::map<types::global_dof_index,double> &boundary_values,
                                   TrilinosMatrix      &matrix,
                                   TrilinosBlockVector &solution,
                                   TrilinosBlockVector &right_hand_side,
                                   const bool          eliminate_columns)
      {
        Assert (eliminate_columns == false, ExcNotImplemented());

        Assert (matrix.n() == right_hand_side.size(),
                ExcDimensionMismatch(matrix.n(), right_hand_side.size()));
        Assert (matrix.n() == solution.size(),
                ExcDimensionMismatch(matrix.n(), solution.size()));
        Assert (matrix.n_block_rows() == matrix.n_block_cols(),
                ExcNotQuadratic());

        const unsigned int n_blocks = matrix.n_block_rows();

        matrix.compress();

        // We need to find the subdivision
        // into blocks for the boundary values.
        // To this end, generate a vector of
        // maps with the respective indices.
        std::vector<std::map<types::global_dof_index,double> > block_boundary_values(n_blocks);
        {
          int block=0;
          types::global_dof_index offset = 0;
          for (std::map<types::global_dof_index,double>::const_iterator
               dof  = boundary_values.begin();
               dof != boundary_values.end();
               ++dof)
            {
              if (dof->first >= matrix.block(block,0).m() + offset)
                {
                  offset += matrix.block(block,0).m();
                  block++;
                }
              const types::global_dof_index index = dof->first - offset;
              block_boundary_values[block].insert(
                std::pair<types::global_dof_index, double> (index,dof->second));
            }
        }

        // Now call the non-block variants on
        // the diagonal subblocks and the
        // solution/rhs.
        for (unsigned int block=0; block<n_blocks; ++block)
          TrilinosWrappers::apply_boundary_values(block_boundary_values[block],
                                                  matrix.block(block,block),
                                                  solution.block(block),
                                                  right_hand_side.block(block),
                                                  eliminate_columns);

        // Finally, we need to do something
        // about the off-diagonal matrices. This
        // is luckily not difficult. Just clear
        // the whole row.
        for (unsigned int block_m=0; block_m<n_blocks; ++block_m)
          {
            const std::pair<types::global_dof_index, types::global_dof_index> local_range
              = matrix.block(block_m,0).local_range();

            std::vector<types::global_dof_index> constrained_rows;
            for (std::map<types::global_dof_index,double>::const_iterator
                 dof  = block_boundary_values[block_m].begin();
                 dof != block_boundary_values[block_m].end();
                 ++dof)
              if ((dof->first >= local_range.first) &&
                  (dof->first < local_range.second))
                constrained_rows.push_back (dof->first);

            for (unsigned int block_n=0; block_n<n_blocks; ++block_n)
              if (block_m != block_n)
                matrix.block(block_m,block_n).clear_rows(constrained_rows);
          }
      }
    }
  }




  void
  apply_boundary_values (const std::map<types::global_dof_index,double> &boundary_values,
                         TrilinosWrappers::SparseMatrix   &matrix,
                         TrilinosWrappers::Vector         &solution,
                         TrilinosWrappers::Vector         &right_hand_side,
                         const bool        eliminate_columns)
  {
    // simply redirect to the generic function
    // used for both trilinos matrix types
    internal::TrilinosWrappers::apply_boundary_values (boundary_values, matrix, solution,
                                                       right_hand_side, eliminate_columns);
  }



  void
  apply_boundary_values (const std::map<types::global_dof_index,double> &boundary_values,
                         TrilinosWrappers::SparseMatrix   &matrix,
                         TrilinosWrappers::MPI::Vector    &solution,
                         TrilinosWrappers::MPI::Vector    &right_hand_side,
                         const bool        eliminate_columns)
  {
    // simply redirect to the generic function
    // used for both trilinos matrix types
    internal::TrilinosWrappers::apply_boundary_values (boundary_values, matrix, solution,
                                                       right_hand_side, eliminate_columns);
  }



  void
  apply_boundary_values (const std::map<types::global_dof_index,double>  &boundary_values,
                         TrilinosWrappers::BlockSparseMatrix &matrix,
                         TrilinosWrappers::BlockVector        &solution,
                         TrilinosWrappers::BlockVector        &right_hand_side,
                         const bool                            eliminate_columns)
  {
    internal::TrilinosWrappers::apply_block_boundary_values (boundary_values, matrix,
                                                             solution, right_hand_side,
                                                             eliminate_columns);
  }



  void
  apply_boundary_values (const std::map<types::global_dof_index,double>  &boundary_values,
                         TrilinosWrappers::BlockSparseMatrix &matrix,
                         TrilinosWrappers::MPI::BlockVector   &solution,
                         TrilinosWrappers::MPI::BlockVector   &right_hand_side,
                         const bool                            eliminate_columns)
  {
    internal::TrilinosWrappers::apply_block_boundary_values (boundary_values, matrix,
                                                             solution, right_hand_side,
                                                             eliminate_columns);
  }

#endif



  void
  local_apply_boundary_values (const std::map<types::global_dof_index,double> &boundary_values,
                               const std::vector<types::global_dof_index> &local_dof_indices,
                               FullMatrix<double> &local_matrix,
                               Vector<double>     &local_rhs,
                               const bool          eliminate_columns)
  {
    Assert (local_dof_indices.size() == local_matrix.m(),
            ExcDimensionMismatch(local_dof_indices.size(),
                                 local_matrix.m()));
    Assert (local_dof_indices.size() == local_matrix.n(),
            ExcDimensionMismatch(local_dof_indices.size(),
                                 local_matrix.n()));
    Assert (local_dof_indices.size() == local_rhs.size(),
            ExcDimensionMismatch(local_dof_indices.size(),
                                 local_rhs.size()));

    // if there is nothing to do, then exit
    // right away
    if (boundary_values.size() == 0)
      return;

    // otherwise traverse all the dofs used in
    // the local matrices and vectors and see
    // what's there to do

    // if we need to treat an entry, then we
    // set the diagonal entry to its absolute
    // value. if it is zero, we used to set it
    // to one, which is a really terrible
    // choice that can lead to hours of
    // searching for bugs in programs (I
    // experienced this :-( ) if the matrix
    // entries are otherwise very large. this
    // is so since iterative solvers would
    // simply not correct boundary nodes for
    // their correct values since the residual
    // contributions of their rows of the
    // linear system is almost zero if the
    // diagonal entry is one. thus, set it to
    // the average absolute value of the
    // nonzero diagonal elements.
    //
    // we only compute this value lazily the
    // first time we need it.
    double average_diagonal = 0;
    const unsigned int n_local_dofs = local_dof_indices.size();
    for (unsigned int i=0; i<n_local_dofs; ++i)
      {
        const std::map<types::global_dof_index, double>::const_iterator
        boundary_value = boundary_values.find (local_dof_indices[i]);
        if (boundary_value != boundary_values.end())
          {
            // remove this row, except for the
            // diagonal element
            for (unsigned int j=0; j<n_local_dofs; ++j)
              if (i != j)
                local_matrix(i,j) = 0;

            // replace diagonal entry by its
            // absolute value to make sure that
            // everything remains positive, or
            // by the average diagonal value if
            // zero
            if (local_matrix(i,i) == 0.)
              {
                // if average diagonal hasn't
                // yet been computed, do so now
                if (average_diagonal == 0.)
                  {
                    unsigned int nonzero_diagonals = 0;
                    for (unsigned int k=0; k<n_local_dofs; ++k)
                      if (local_matrix(k,k) != 0.)
                        {
                          average_diagonal += std::fabs(local_matrix(k,k));
                          ++nonzero_diagonals;
                        }
                    if (nonzero_diagonals != 0)
                      average_diagonal /= nonzero_diagonals;
                    else
                      average_diagonal = 0;
                  }

                // only if all diagonal entries
                // are zero, then resort to the
                // last measure: choose one
                if (average_diagonal == 0.)
                  average_diagonal = 1.;

                local_matrix(i,i) = average_diagonal;
              }
            else
              local_matrix(i,i) = std::fabs(local_matrix(i,i));

            // and replace rhs entry by correct
            // value
            local_rhs(i) = local_matrix(i,i) * boundary_value->second;

            // finally do the elimination step
            // if requested
            if (eliminate_columns == true)
              {
                for (unsigned int row=0; row<n_local_dofs; ++row)
                  if (row != i)
                    {
                      local_rhs(row) -= local_matrix(row,i) * boundary_value->second;
                      local_matrix(row,i) = 0;
                    }
              }
          }
      }
  }
}



// explicit instantiations
#include "matrix_tools.inst"

namespace MatrixTools
{
  template
  void
  apply_boundary_values<double> (const std::map<types::global_dof_index,double> &boundary_values,
                                 SparseMatrix<double>  &matrix,
                                 Vector<double>   &solution,
                                 Vector<double>   &right_hand_side,
                                 const bool        eliminate_columns);
  template
  void
  apply_boundary_values<float> (const std::map<types::global_dof_index,double> &boundary_values,
                                SparseMatrix<float>  &matrix,
                                Vector<float>   &solution,
                                Vector<float>   &right_hand_side,
                                const bool        eliminate_columns);

  template
  void
  apply_boundary_values<double> (const std::map<types::global_dof_index,double> &boundary_values,
                                 BlockSparseMatrix<double>  &matrix,
                                 BlockVector<double>   &solution,
                                 BlockVector<double>   &right_hand_side,
                                 const bool        eliminate_columns);
  template
  void
  apply_boundary_values<float> (const std::map<types::global_dof_index,double> &boundary_values,
                                BlockSparseMatrix<float>  &matrix,
                                BlockVector<float>   &solution,
                                BlockVector<float>   &right_hand_side,
                                const bool        eliminate_columns);
}

DEAL_II_NAMESPACE_CLOSE

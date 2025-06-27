// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_vector_tools_mean_value_templates_h
#define dealii_vector_tools_mean_value_templates_h

#include <deal.II/distributed/tria_base.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/filtered_iterator.h>

#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/vector_tools_mean_value.h>


DEAL_II_NAMESPACE_OPEN

namespace VectorTools
{
  namespace internal
  {
    template <typename VectorType>
    std::enable_if_t<dealii::is_serial_vector<VectorType>::value == true>
    subtract_mean_value(VectorType &v, const std::vector<bool> &p_select)
    {
      if (p_select.empty())
        {
          // In case of an empty boolean mask operate on the whole vector:
          v.add(-v.mean_value());
        }
      else
        {
          const unsigned int n = v.size();

          Assert(p_select.size() == n,
                 ExcDimensionMismatch(p_select.size(), n));

          typename VectorType::value_type s       = 0.;
          unsigned int                    counter = 0;
          for (unsigned int i = 0; i < n; ++i)
            if (p_select[i])
              {
                typename VectorType::value_type vi = v(i);
                s += vi;
                ++counter;
              }
          // Error out if we have not constrained anything. Note that in this
          // case the vector v is always nonempty.
          Assert(n == 0 || counter > 0,
                 ComponentMask::ExcNoComponentSelected());

          s /= counter;

          for (unsigned int i = 0; i < n; ++i)
            if (p_select[i])
              v(i) -= s;
        }
    }



    template <typename VectorType>
    std::enable_if_t<dealii::is_serial_vector<VectorType>::value == false>
    subtract_mean_value(VectorType &v, const std::vector<bool> &p_select)
    {
      Assert(p_select.empty(), ExcNotImplemented());
      // In case of an empty boolean mask operate on the whole vector:
      v.add(-v.mean_value());
    }
  } // namespace internal


  template <typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void subtract_mean_value(VectorType &v, const std::vector<bool> &p_select)
  {
    internal::subtract_mean_value(v, p_select);
  }

  namespace internal
  {
    template <typename Number>
    void
    set_possibly_complex_number(const double r, const double, Number &n)
    {
      n = r;
    }



    template <typename Type>
    void
    set_possibly_complex_number(const double        r,
                                const double        i,
                                std::complex<Type> &n)
    {
      n = std::complex<Type>(r, i);
    }
  } // namespace internal



  template <typename VectorType, int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void add_constant(VectorType                           &solution,
                    const DoFHandler<dim, spacedim>      &dof_handler,
                    const unsigned int                    component,
                    const typename VectorType::value_type constant_adjustment)
  {
    Assert(dof_handler.has_hp_capabilities() == false, ExcNotImplemented());

    AssertDimension(solution.size(), dof_handler.n_dofs());
    AssertIndexRange(component, dof_handler.get_fe().n_components());

    const FiniteElement<dim, spacedim> &fe_system = dof_handler.get_fe();
    const FiniteElement<dim, spacedim> &fe = fe_system.get_sub_fe(component, 1);

    if ((dynamic_cast<const FE_DGP<dim, spacedim> *>(&fe) != nullptr))
      {
        // The FE to modify is an FE_DGP, which is not a nodal
        // element. The first shape function of a DGP element happens
        // to be the constant function, so we just have to adjust the
        // corresponding DoF on each cell:
        std::vector<types::global_dof_index> local_dof_indices(
          dof_handler.get_fe().dofs_per_cell);

        for (const auto &cell : dof_handler.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              cell->get_dof_indices(local_dof_indices);
              const unsigned int first_pressure_dof =
                fe_system.component_to_system_index(component, 0);

              // Make sure that this DoF is really owned by the
              // current processor:
              Assert(dof_handler.locally_owned_dofs().is_element(
                       local_dof_indices[first_pressure_dof]),
                     ExcInternalError());

              // Then adjust its value:
              solution(local_dof_indices[first_pressure_dof]) +=
                constant_adjustment;
            }

        solution.compress(VectorOperation::add);
      }
    else if ((dynamic_cast<const FE_Q<dim, spacedim> *>(&fe) != nullptr) ||
             (dynamic_cast<const FE_SimplexP<dim, spacedim> *>(&fe) != nullptr))
      {
        // We need to make sure to not touch DoFs shared between cells more
        // than once. Instead of counting or limiting the number of times
        // we touch an individual vector entry, we make a copy of the vector
        // and use that as the input and add the constant to it. This way
        // it does not matter how often an individual entry is touched.

        VectorType copy(solution);
        copy = solution;

        std::vector<types::global_dof_index> local_dof_indices(
          dof_handler.get_fe().dofs_per_cell);

        for (const auto &cell : dof_handler.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              cell->get_dof_indices(local_dof_indices);
              for (unsigned i = 0; i < dof_handler.get_fe().dofs_per_cell; ++i)
                {
                  if (!fe_system.is_primitive(i))
                    continue;

                  const auto component_and_index =
                    fe_system.system_to_component_index(i);
                  if (component_and_index.first == component)
                    {
                      const types::global_dof_index idx = local_dof_indices[i];
                      // Make sure that this DoF is really owned by the
                      // current processor:
                      if (!dof_handler.locally_owned_dofs().is_element(idx))
                        continue;

                      // Then adjust its value:
                      solution(idx) =
                        typename VectorType::value_type(copy(idx)) +
                        constant_adjustment;
                    }
                }
            }

        solution.compress(VectorOperation::insert);
      }
    else if ((dynamic_cast<const FE_DGQ<dim, spacedim> *>(&fe) != nullptr) ||
             (dynamic_cast<const FE_SimplexDGP<dim, spacedim> *>(&fe) !=
              nullptr))
      {
        // Add the constant to every single shape function per cell

        std::vector<types::global_dof_index> local_dof_indices(
          dof_handler.get_fe().dofs_per_cell);

        for (const auto &cell : dof_handler.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              cell->get_dof_indices(local_dof_indices);
              for (unsigned i = 0; i < dof_handler.get_fe().dofs_per_cell; ++i)
                {
                  if (!fe_system.is_primitive(i))
                    continue;

                  const auto component_and_index =
                    fe_system.system_to_component_index(i);
                  if (component_and_index.first == component)
                    {
                      // Make sure that this DoF is really owned by the
                      // current processor:
                      Assert(dof_handler.locally_owned_dofs().is_element(
                               local_dof_indices[i]),
                             ExcInternalError());

                      // Then adjust its value:
                      solution(local_dof_indices[i]) += constant_adjustment;
                    }
                }
            }

        solution.compress(VectorOperation::add);
      }
    else
      AssertThrow(false, ExcNotImplemented());
  }



#ifdef DEAL_II_WITH_TRILINOS
  template <int dim, int spacedim>
  void
  add_constant(LinearAlgebra::EpetraWrappers::Vector &,
               const DoFHandler<dim, spacedim> &,
               const unsigned int,
               const double)
  {
    // TODO: no vector access using operator()
    AssertThrow(false, ExcNotImplemented());
  }



#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
  template <int dim, int spacedim, typename ValueType, typename MemorySpace>
  void
  add_constant(LinearAlgebra::TpetraWrappers::Vector<ValueType, MemorySpace> &,
               const DoFHandler<dim, spacedim> &,
               const unsigned int,
               const ValueType)
  {
    // TODO: no vector access using operator()
    AssertThrow(false, ExcNotImplemented());
  }
#  endif
#endif



  template <int dim, typename Number, int spacedim>
  Number
  compute_mean_value(
    const hp::MappingCollection<dim, spacedim> &mapping_collection,
    const DoFHandler<dim, spacedim>            &dof,
    const hp::QCollection<dim>                 &q_collection,
    const ReadVector<Number>                   &v,
    const unsigned int                          component)
  {
    const hp::FECollection<dim, spacedim> &fe_collection =
      dof.get_fe_collection();
    const unsigned int n_components = fe_collection.n_components();

    AssertDimension(v.size(), dof.n_dofs());
    AssertIndexRange(component, n_components);

    hp::FEValues<dim, spacedim> fe_values_collection(
      mapping_collection,
      fe_collection,
      q_collection,
      UpdateFlags(update_JxW_values | update_values));

    std::vector<Vector<Number>> values;

    Number                                            mean = Number();
    typename numbers::NumberTraits<Number>::real_type area = 0.;
    // Compute mean value
    for (const auto &cell :
         dof.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
      {
        fe_values_collection.reinit(cell);
        const FEValues<dim, spacedim> &fe_values =
          fe_values_collection.get_present_fe_values();

        values.resize(fe_values.n_quadrature_points,
                      Vector<Number>(n_components));
        fe_values.get_function_values(v, values);
        for (unsigned int k = 0; k < fe_values.n_quadrature_points; ++k)
          {
            mean += fe_values.JxW(k) * values[k](component);
            area += fe_values.JxW(k);
          }
      }

#ifdef DEAL_II_WITH_MPI
    // if this was a distributed DoFHandler, we need to do the reduction
    // over the entire domain
    if (const parallel::TriangulationBase<dim, spacedim> *p_triangulation =
          dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
            &dof.get_triangulation()))
      {
        // The type used to store the elements of the global vector may be a
        // real or a complex number. Do the global reduction always with real
        // and imaginary types so that we don't have to distinguish, and to
        // this end just copy everything into a complex number and, later,
        // back into the original data type.
        std::complex<double> mean_double = mean;
        double my_values[3] = {mean_double.real(), mean_double.imag(), area};
        double global_values[3];

        const int ierr = MPI_Allreduce(my_values,
                                       global_values,
                                       3,
                                       MPI_DOUBLE,
                                       MPI_SUM,
                                       p_triangulation->get_mpi_communicator());
        AssertThrowMPI(ierr);

        internal::set_possibly_complex_number(global_values[0],
                                              global_values[1],
                                              mean);
        area = global_values[2];
      }
#endif

    return (mean / area);
  }


  template <int dim, typename Number, int spacedim>
  Number
  compute_mean_value(const Mapping<dim, spacedim>    &mapping,
                     const DoFHandler<dim, spacedim> &dof,
                     const Quadrature<dim>           &quadrature,
                     const ReadVector<Number>        &v,
                     const unsigned int               component)
  {
    return compute_mean_value(hp::MappingCollection<dim, spacedim>(mapping),
                              dof,
                              hp::QCollection<dim>(quadrature),
                              v,
                              component);
  }


  template <int dim, typename Number, int spacedim>
  Number
  compute_mean_value(const DoFHandler<dim, spacedim> &dof,
                     const Quadrature<dim>           &quadrature,
                     const ReadVector<Number>        &v,
                     const unsigned int               component)
  {
    return compute_mean_value(get_default_linear_mapping(
                                dof.get_triangulation()),
                              dof,
                              quadrature,
                              v,
                              component);
  }
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_mean_value_templates_h

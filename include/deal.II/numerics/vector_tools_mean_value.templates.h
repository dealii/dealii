// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_vector_tools_mean_value_templates_h
#define dealii_vector_tools_mean_value_templates_h


#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_vector.h>
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
    typename std::enable_if<dealii::is_serial_vector<VectorType>::value ==
                            true>::type
    subtract_mean_value(VectorType &v, const std::vector<bool> &p_select)
    {
      if (p_select.size() == 0)
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
    typename std::enable_if<dealii::is_serial_vector<VectorType>::value ==
                            false>::type
    subtract_mean_value(VectorType &v, const std::vector<bool> &p_select)
    {
      (void)p_select;
      Assert(p_select.size() == 0, ExcNotImplemented());
      // In case of an empty boolean mask operate on the whole vector:
      v.add(-v.mean_value());
    }
  } // namespace internal


  template <typename VectorType>
  void
  subtract_mean_value(VectorType &v, const std::vector<bool> &p_select)
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

  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  compute_mean_value(const Mapping<dim, spacedim> &   mapping,
                     const DoFHandler<dim, spacedim> &dof,
                     const Quadrature<dim> &          quadrature,
                     const VectorType &               v,
                     const unsigned int               component)
  {
    using Number = typename VectorType::value_type;
    Assert(v.size() == dof.n_dofs(),
           ExcDimensionMismatch(v.size(), dof.n_dofs()));
    AssertIndexRange(component, dof.get_fe(0).n_components());

    FEValues<dim, spacedim> fe(mapping,
                               dof.get_fe(),
                               quadrature,
                               UpdateFlags(update_JxW_values | update_values));

    std::vector<Vector<Number>> values(
      quadrature.size(), Vector<Number>(dof.get_fe(0).n_components()));

    Number                                            mean = Number();
    typename numbers::NumberTraits<Number>::real_type area = 0.;
    // Compute mean value
    for (const auto &cell : dof.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe.reinit(cell);
          fe.get_function_values(v, values);
          for (unsigned int k = 0; k < quadrature.size(); ++k)
            {
              mean += fe.JxW(k) * values[k](component);
              area += fe.JxW(k);
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
                                       p_triangulation->get_communicator());
        AssertThrowMPI(ierr);

        internal::set_possibly_complex_number(global_values[0],
                                              global_values[1],
                                              mean);
        area = global_values[2];
      }
#endif

    return (mean / area);
  }


  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  compute_mean_value(const DoFHandler<dim, spacedim> &dof,
                     const Quadrature<dim> &          quadrature,
                     const VectorType &               v,
                     const unsigned int               component)
  {
    return compute_mean_value(
      StaticMappingQ1<dim, spacedim>::mapping, dof, quadrature, v, component);
  }
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_mean_value_templates_h

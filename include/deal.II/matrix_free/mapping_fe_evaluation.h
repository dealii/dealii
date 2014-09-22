// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2014 by the deal.II authors
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


#ifndef __deal2__matrix_free_mapping_fe_evaluation_h
#define __deal2__matrix_free_mapping_fe_evaluation_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/base/aligned_vector.h>
#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/mapping_info.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_nothing.h>


DEAL_II_NAMESPACE_OPEN


/**
 * This class provides evaluated mapping information using standard deal.II
 * information in a form that FEEvaluation and friends can use for vectorized
 * access. Since no vectorization over cells is available with the
 * DoFHandler/Triangulation cell iterators, the interface to FEEvaluation's
 * vectorization model is to use @p VectorizedArray::n_array_element copies of
 * the same element. This interface is thus primarily useful for evaluating
 * several operators on the same cell, e.g., when assembling cell matrices.
 *
 * As opposed to the standard Mapping classes in deal.II, this class does not
 * actually provide a boundary description that can be used to evaluate the
 * geometry, but it rather provides the evaluated geometry from a given
 * deal.II mapping (as passed to the constructor of this class) in a form
 * accessible to FEEvaluation.
 *
 * @author Martin Kronbichler, 2014
 */
template <int dim, typename Number=double>
class MappingFEEvaluation : public Subscriptor
{
public:
  /**
   * Constructor, similar to FEValues. Since this class only evaluates the
   * geometry, no finite element has to be specified and the simplest element,
   * FE_Nothing, is used internally for the underlying FEValues object.
   */
  MappingFEEvaluation (const Mapping<dim> &mapping,
                       const Quadrature<1> &quadrature,
                       const UpdateFlags update_flags);

  /**
   * Constructor. Instead of providing a mapping, use MappingQ1.
   */
  MappingFEEvaluation (const Quadrature<1> &quadrature,
                       const UpdateFlags update_flags);

  /**
   * Initialize with the given cell iterator. This works for all type of
   * iterators that can be converted into a Triangulation::cell_iterator
   * (DoFHandler::cell_iterator or Triangulation::cell_iterator).
   */
  template <typename ITERATOR>
  void reinit(ITERATOR cell);

  /**
   * Return a triangulation iterator to the current cell.
   */
  typename Triangulation<dim>::cell_iterator get_cell () const;

  /**
   * Return a reference to the underlying FEValues object that evaluates
   * certain quantities (only mapping-related ones like Jacobians or mapped
   * quadrature points are accessible, as no finite element data is actually
   * used).
   */
  const FEValues<dim> &get_fe_values () const;

  /**
   * Return a vector of inverse transpose Jacobians. For compatibility with
   * FEEvaluation, it returns tensors of vectorized arrays, even though all
   * components are equal.
   */
  const AlignedVector<Tensor<2,dim,VectorizedArray<Number> > > &
  get_inverse_jacobians() const;

  /**
   * Return a vector of quadrature weights times the Jacobian determinant
   * (JxW). For compatibility with FEEvaluation, it returns tensors of
   * vectorized arrays, even though all components are equal.
   */
  const AlignedVector<VectorizedArray<Number> > &
  get_JxW_values() const;

  /**
   * Return a vector of quadrature points in real space on the given
   * cell. For compatibility with FEEvaluation, it returns tensors of
   * vectorized arrays, even though all components are equal.
   */
  const AlignedVector<Point<dim,VectorizedArray<Number> > > &
  get_quadrature_points() const;

  /**
   * Return a vector of quadrature points in real space on the given
   * cell. For compatibility with FEEvaluation, it returns tensors of
   * vectorized arrays, even though all components are equal.
   */
  const AlignedVector<Tensor<1,dim,VectorizedArray<Number> > > &
  get_normal_vectors() const;

  /**
   * Return a reference to 1D quadrature underlying this object.
   */
  const Quadrature<1> &
  get_quadrature () const;

private:
  /**
   * Dummy finite element object necessary for initializing the FEValues
   * object.
   */
  FE_Nothing<dim> fe_dummy;

  /**
   * An underlying FEValues object that performs the (scalar) evaluation.
   */
  FEValues<dim> fe_values;

  /**
   * Get 1D quadrature formula to be used for reinitializing shape info.
   */
  const Quadrature<1> quadrature_1d;

  /**
   * Inverse Jacobians, stored in vectorized array form.
   */
  AlignedVector<Tensor<2,dim,VectorizedArray<Number> > > inverse_jacobians;

  /**
   * Stored Jacobian determinants and quadrature weights
   */
  AlignedVector<VectorizedArray<Number> > jxw_values;

  /**
   * Stored quadrature points
   */
  AlignedVector<Point<dim,VectorizedArray<Number> > > quadrature_points;

  /**
   * Stored normal vectors (for face integration)
   */
  AlignedVector<Tensor<1,dim,VectorizedArray<Number> > > normal_vectors;
};


/*----------------------- Inline functions ----------------------------------*/

template <int dim, typename Number>
inline
MappingFEEvaluation<dim,Number>::MappingFEEvaluation (const Mapping<dim> &mapping,
                                                      const Quadrature<1> &quadrature,
                                                      const UpdateFlags update_flags)
  :
  fe_values(mapping, fe_dummy, Quadrature<dim>(quadrature),
            internal::MatrixFreeFunctions::MappingInfo<dim,Number>::compute_update_flags(update_flags)),
  quadrature_1d(quadrature),
  inverse_jacobians(fe_values.get_quadrature().size()),
  jxw_values(fe_values.get_quadrature().size()),
  quadrature_points(fe_values.get_quadrature().size()),
  normal_vectors(fe_values.get_quadrature().size())
{
  Assert(!(fe_values.get_update_flags() & update_jacobian_grads),
         ExcNotImplemented());
}



template <int dim, typename Number>
inline
MappingFEEvaluation<dim,Number>::MappingFEEvaluation (const Quadrature<1> &quadrature,
                                                      const UpdateFlags update_flags)
  :
  fe_values(fe_dummy, Quadrature<dim>(quadrature),
            internal::MatrixFreeFunctions::MappingInfo<dim,Number>::compute_update_flags(update_flags)),
  quadrature_1d(quadrature),
  inverse_jacobians(fe_values.get_quadrature().size()),
  jxw_values(fe_values.get_quadrature().size()),
  quadrature_points(fe_values.get_quadrature().size()),
  normal_vectors(fe_values.get_quadrature().size())
{
  Assert(!(fe_values.get_update_flags() & update_jacobian_grads),
         ExcNotImplemented());
}



template <int dim, typename Number>
template <typename ITERATOR>
inline
void
MappingFEEvaluation<dim,Number>::reinit(ITERATOR cell)
{
  fe_values.reinit(static_cast<typename Triangulation<dim>::cell_iterator>(cell));
  for (unsigned int q=0; q<fe_values.get_quadrature().size(); ++q)
    {
      if (fe_values.get_update_flags() & update_inverse_jacobians)
        for (unsigned int d=0; d<dim; ++d)
          for (unsigned int e=0; e<dim; ++e)
            inverse_jacobians[q][d][e] = fe_values.inverse_jacobian(q)[e][d];
      if (fe_values.get_update_flags() & update_quadrature_points)
        for (unsigned int d=0; d<dim; ++d)
          quadrature_points[q][d] = fe_values.quadrature_point(q)[d];
      if (fe_values.get_update_flags() & update_normal_vectors)
        for (unsigned int d=0; d<dim; ++d)
          normal_vectors[q][d] = fe_values.normal_vector(q)[d];
      if (fe_values.get_update_flags() & update_JxW_values)
        jxw_values[q] = fe_values.JxW(q);
    }
}



template <int dim, typename Number>
inline
typename Triangulation<dim>::cell_iterator
MappingFEEvaluation<dim,Number>::get_cell() const
{
  return fe_values.get_cell();
}



template <int dim, typename Number>
inline
const FEValues<dim> &
MappingFEEvaluation<dim,Number>::get_fe_values() const
{
  return fe_values;
}



template <int dim, typename Number>
inline
const AlignedVector<Tensor<2,dim,VectorizedArray<Number> > > &
MappingFEEvaluation<dim,Number>::get_inverse_jacobians() const
{
  return inverse_jacobians;
}



template <int dim, typename Number>
inline
const AlignedVector<Tensor<1,dim,VectorizedArray<Number> > > &
MappingFEEvaluation<dim,Number>::get_normal_vectors() const
{
  return normal_vectors;
}



template <int dim, typename Number>
inline
const AlignedVector<Point<dim,VectorizedArray<Number> > > &
MappingFEEvaluation<dim,Number>::get_quadrature_points() const
{
  return quadrature_points;
}



template <int dim, typename Number>
inline
const AlignedVector<VectorizedArray<Number> > &
MappingFEEvaluation<dim,Number>::get_JxW_values() const
{
  return jxw_values;
}



template <int dim, typename Number>
inline
const Quadrature<1> &
MappingFEEvaluation<dim,Number>::get_quadrature() const
{
  return quadrature_1d;
}


#ifndef DOXYGEN


#endif  // ifndef DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif

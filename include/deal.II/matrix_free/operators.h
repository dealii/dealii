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


#ifndef __deal2__matrix_free_operators_h
#define __deal2__matrix_free_operators_h


#include <deal.II/base/exceptions.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/matrix_free/fe_evaluation.h>


DEAL_II_NAMESPACE_OPEN


namespace MatrixFreeOperators
{
  /**
   * This class implements the operation of the action of the inverse of a mass
   * matrix on an element for the special case of an evaluation object with as
   * many quadrature points as there are cell degrees of freedom. It uses
   * algorithms from FEEvaluation and produces the exact mass matrix for DGQ
   * elements. This algorithm uses tensor products of inverse 1D shape matrices
   * over quadrature points, so the inverse operation is exactly as expensive as
   * applying the forward operator on each cell. Of course, for continuous
   * finite elements this operation does not produce the inverse of a mass
   * operation as the coupling between the elements cannot be considered by this
   * operation.
   *
   * The equation may contain variable coefficients, so the user is required to
   * provide an array for the inverse of the local coefficient (this class
   * provide a helper method 'fill_inverse_JxW_values' to get the inverse of a
   * constant-coefficient operator).
   *
   * @author Martin Kronbichler, 2014
   */
  template <int dim, int fe_degree, int n_components = 1, typename Number = double>
  class CellwiseInverseMassMatrix
  {
  public:
    /**
     * Constructor. Initializes the shape information from the ShapeInfo field
     * in the class FEEval.
     */
    CellwiseInverseMassMatrix (const FEEvaluationBase<dim,n_components,Number> &fe_eval);

    /**
     * Applies the inverse mass matrix operation on an input array. It is
     * assumed that the passed input and output arrays are of correct size,
     * namely FEEval::dofs_per_cell * n_components long. The inverse of the
     * local coefficient (also containing the inverse JxW values) must be passed
     * as first argument. Passing more than one component in the coefficient is
     * allowed.
     */
    void apply(const AlignedVector<VectorizedArray<Number> > &inverse_coefficient,
               const unsigned int             n_actual_components,
               const VectorizedArray<Number> *in_array,
               VectorizedArray<Number>       *out_array) const;

    /**
     * Fills the given array with the inverse of the JxW values, i.e., a mass
     * matrix with coefficient 1. Non-unit coefficients must be multiplied (in
     * inverse form) to this array.
     */
    void fill_inverse_JxW_values(AlignedVector<VectorizedArray<Number> > &inverse_jxw) const;

  private:
    /**
     * A reference to the FEEvaluation object for getting the JxW_values.
     */
    const FEEvaluationBase<dim,n_components,Number> &fe_eval;

    /**
     * A structure to hold inverse shape functions
     */
    AlignedVector<VectorizedArray<Number> > inverse_shape;
  };



  // ------------------------------------ inline functions ---------------------

  template <int dim, int fe_degree, int n_components, typename Number>
  inline
  CellwiseInverseMassMatrix<dim,fe_degree,n_components,Number>
  ::CellwiseInverseMassMatrix (const FEEvaluationBase<dim,n_components,Number> &fe_eval)
    :
    fe_eval (fe_eval)
  {
    FullMatrix<double> shapes_1d(fe_degree+1, fe_degree+1);
    for (unsigned int i=0, c=0; i<shapes_1d.m(); ++i)
      for (unsigned int j=0; j<shapes_1d.n(); ++j, ++c)
        shapes_1d(i,j) = fe_eval.get_shape_info().shape_values_number[c];
    shapes_1d.gauss_jordan();
    const unsigned int stride = (fe_degree+2)/2;
    inverse_shape.resize(stride*(fe_degree+1));
    for (unsigned int i=0; i<stride; ++i)
      for (unsigned int q=0; q<(fe_degree+2)/2; ++q)
        {
          inverse_shape[i*stride+q] =
            0.5 * (shapes_1d(i,q) + shapes_1d(i,fe_degree-q));
          inverse_shape[(fe_degree-i)*stride+q] =
            0.5 * (shapes_1d(i,q) - shapes_1d(i,fe_degree-q));
        }
    if (fe_degree % 2 == 0)
      for (unsigned int q=0; q<(fe_degree+2)/2; ++q)
        inverse_shape[fe_degree/2*stride+q] = shapes_1d(fe_degree/2,q);
  }



  template <int dim, int fe_degree, int n_components, typename Number>
  inline
  void
  CellwiseInverseMassMatrix<dim,fe_degree,n_components,Number>
  ::fill_inverse_JxW_values(AlignedVector<VectorizedArray<Number> > &inverse_jxw) const
  {
    const unsigned int dofs_per_cell = Utilities::fixed_int_power<fe_degree+1,dim>::value;
    Assert(inverse_jxw.size() > 0 &&
           inverse_jxw.size() % dofs_per_cell == 0,
           ExcMessage("Expected diagonal to be a multiple of scalar dof per cells"));

    // temporarily reduce size of inverse_jxw to dofs_per_cell to get JxW values
    // from fe_eval (will not reallocate any memory)
    const unsigned int previous_size = inverse_jxw.size();
    inverse_jxw.resize(dofs_per_cell);
    fe_eval.fill_JxW_values(inverse_jxw);

    // invert
    inverse_jxw.resize_fast(previous_size);
    for (unsigned int q=0; q<dofs_per_cell; ++q)
      inverse_jxw[q] = 1./inverse_jxw[q];
    // copy values to rest of vector
    for (unsigned int q=dofs_per_cell; q<previous_size; )
      for (unsigned int i=0; i<dofs_per_cell; ++i, ++q)
        inverse_jxw[q] = inverse_jxw[i];
  }



  template <int dim, int fe_degree, int n_components, typename Number>
  inline
  void
  CellwiseInverseMassMatrix<dim,fe_degree,n_components,Number>
  ::apply(const AlignedVector<VectorizedArray<Number> > &inverse_coefficients,
          const unsigned int             n_actual_components,
          const VectorizedArray<Number> *in_array,
          VectorizedArray<Number>       *out_array) const
  {
    const unsigned int dofs_per_cell = Utilities::fixed_int_power<fe_degree+1,dim>::value;
    Assert(inverse_coefficients.size() > 0 &&
           inverse_coefficients.size() % dofs_per_cell == 0,
           ExcMessage("Expected diagonal to be a multiple of scalar dof per cells"));
    if (inverse_coefficients.size() != dofs_per_cell)
      AssertDimension(n_actual_components * dofs_per_cell, inverse_coefficients.size());

    Assert(dim == 2 || dim == 3, ExcNotImplemented());

    internal::EvaluatorTensorProduct<internal::evaluate_evenodd,dim,fe_degree,
             fe_degree+1, VectorizedArray<Number> >
             evaluator(inverse_shape, inverse_shape, inverse_shape);

    const unsigned int shift_coefficient =
      inverse_coefficients.size() > dofs_per_cell ? dofs_per_cell : 0;
    const VectorizedArray<Number> *inv_coefficient = &inverse_coefficients[0];
    VectorizedArray<Number> temp_data_field[dofs_per_cell];
    for (unsigned int d=0; d<n_actual_components; ++d)
      {
        const VectorizedArray<Number> *in = in_array+d*dofs_per_cell;
        VectorizedArray<Number> *out = out_array+d*dofs_per_cell;
        // Need to select 'apply' method with hessian slot because values
        // assume symmetries that do not exist in the inverse shapes
        evaluator.template hessians<0,false,false> (in, temp_data_field);
        evaluator.template hessians<1,false,false> (temp_data_field, out);

        if (dim == 3)
          {
            evaluator.template hessians<2,false,false> (out, temp_data_field);
            for (unsigned int q=0; q<dofs_per_cell; ++q)
              temp_data_field[q] *= inv_coefficient[q];
            evaluator.template hessians<2,true,false> (temp_data_field, out);
          }
        else if (dim == 2)
          for (unsigned int q=0; q<dofs_per_cell; ++q)
            out[q] *= inv_coefficient[q];

        evaluator.template hessians<1,true,false>(out, temp_data_field);
        evaluator.template hessians<0,true,false>(temp_data_field, out);

        inv_coefficient += shift_coefficient;
      }
  }

} // end of namespace MatrixFreeOperators


DEAL_II_NAMESPACE_CLOSE

#endif

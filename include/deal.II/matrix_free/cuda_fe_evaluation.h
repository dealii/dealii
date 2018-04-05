// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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

#ifndef dealii_cuda_fe_evaluation_h
#define dealii_cuda_fe_evaluation_h

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/cuda_vector.h>
#include <deal.II/matrix_free/cuda_matrix_free.h>
#include <deal.II/matrix_free/cuda_matrix_free.templates.h>
#include <deal.II/matrix_free/cuda_tensor_product_kernels.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Namespace for the CUDA wrappers
 */
namespace CUDAWrappers
{
  namespace internal
  {
    template <int dim, int fe_degree, bool transpose, typename Number>
    __device__ void resolve_hanging_nodes_shmem(Number *values, const unsigned
                                                int constr)
    {
      //TODO
    }
  }



  /**
   * This class provides all the functions necessary to evaluate functions at
   * quadrature points and cell integrations. In functionality, this class is
   * similar to FEValues<dim>.
   *
   * This class has five template arguments:
   *
   * @ptaram dim Dimension in which this class is to be used
   *
   * @tparam fe_degree Degree of the tensor prodict finite element with fe_degree+1
   * degrees of freedom per coordinate direction
   *
   * @tparam n_q_points_1d Number of points in the quadrature formular in 1D,
   * defaults to fe_degree+1
   *
   * @tparam n_components Number of vector components when solving a system of
   * PDEs. If the same operation is applied to several components of a PDE (e.g.
   * a vector Laplace equation), they can be applied simultaneously with one call
   * (and often more efficiently). Defaults to 1
   *
   * @tparam Number Number format, @p double or @p float. Defaults to @p
   * double.
   *
   * @ingroup CUDAWrappers
   *
   * @author Karl Ljungkvist, Bruno Turcksin, 2016
   */
  template <int dim, int fe_degree, int n_q_points_1d = fe_degree+1,
            int n_components_ = 1, typename Number = double>
  class FEEvaluation
  {
  public:
    typedef Number                                  value_type;
    typedef Tensor<1,dim,Number>                    gradient_type;
    typedef typename MatrixFree<dim, Number>::Data  data_type;
    static constexpr unsigned int dimension       = dim;
    static constexpr unsigned int n_components    = n_components_;
    static constexpr unsigned int n_q_points      = Utilities::pow(n_q_points_1d, dim);
    static constexpr unsigned int tensor_dofs_per_cell = Utilities::pow(fe_degree + 1, dim);

    /**
     * Constructor.
     */
    __device__ FEEvaluation(int cell_id,
                            const data_type *data,
                            SharedData<dim,Number> *shdata);

    /**
     * For the vector @p src, read out the values on the degrees of freedom of
     * the current cell, and store them internally. Similar functionality as
     * the function DoFAccessor::get_interpolated_dof_values when no
     * constraints are present, but it also includes constraints from hanging
     * nodes, so once can see it as a similar function to
     * ConstraintMatrix::read_dof_valuess as well.
     */
    __device__ void read_dof_values(const Number *src);

    /**
     * Take the value stored internally on dof values of the current cell and
     * sum them into the vector @p dst. The function also applies constraints
     * during the write operation. The functionality is hence similar to the
     * function ConstraintMatrix::distribute_local_to_global.
     */
    __device__  void distribute_local_to_global(Number *dst) const;

    /**
     * Evaluate the function values and the gradients of the FE function given
     * at the DoF values in the input vector at the quadrature points on the
     * unit cell. The function arguments specify which parts shall actually be
     * computed. This function needs to be called before the functions
     * @p get_value() or @p get_gradient() give useful information.
     */
    __device__ void evaluate(const bool evaluate_val,
                             const bool evaluate_grad);

    /**
     * This function takes the values and/or gradients that are stored on
     * quadrature points, tests them by all the basis functions/gradients on
     * the cell and performs the cell integration. The two function arguments
     * @p integrate_val and @p integrate_grad are used to enable/disable some
     * of the values or the gradients.
     */
    __device__ void integrate(const bool integrate_val,
                              const bool integrate_grad);

    /**
     * Return the value of a finite element function at quadrature point
     * number @p q_point after a call to @p evalue(true,...).
     */
    __device__ value_type get_value(const unsigned int q_point) const;

    /**
     * Write a value to the field containing the values on quadrature points
     * with component @p q_point. Access to the same fields as through @p
     * get_value(), This specifies the value which is tested by all basis
     * function on the current cell and integrated over.
     */
    __device__ void submit_value(const value_type &val_in,
                                 const unsigned int q_point);

    /**
     * Return the gradient of a finite element function at quadrature point
     * number @p q_point after a call to @p evaluate(...,true).
     */
    __device__ gradient_type get_gradient(const unsigned int q_point) const;

    /**
     * Write a contribution that is tested by the gradient to the field
     * containing the values on quadrature points with component @p q_point
     */
    __device__ void submit_gradient(const gradient_type &grad_in,
                                    const unsigned int q_point);

    /**
     * Apply the function @p func on every quadrature point.
     */
    template <typename functor>
    __device__ void apply_quad_point_operations(const functor &func);

  private:
    unsigned int *local_to_global;
    unsigned int n_cells;
    unsigned int padding_length;

    const unsigned int constraint_mask;

    Number *inv_jac;
    Number *JxW;

    // Internal buffer
    Number *values;
    Number *gradients[dim];
  };



  template <int dim, int fe_degree, int n_q_points_1d, int n_components_,
            typename Number>
  __device__
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
  FEEvaluation(int cell_id,
               const data_type *data,
               SharedData<dim,Number> *shdata)
    :
    n_cells(data->n_cells),
    padding_length(data->padding_length),
    constraint_mask(data->constraint_mask[cell_id]),
    values(shdata->values)
  {
    local_to_global = data->local_to_global + padding_length*cell_id;
    inv_jac = data->inv_jacobian + padding_length*cell_id;
    JxW = data->JxW + padding_length*cell_id;

    for (unsigned int i=0; i < dim; ++i)
      gradients[i] = shdata->gradients[i];
  }



  template <int dim, int fe_degree, int n_q_points_1d, int n_components_,
            typename Number>
  __device__ void
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
  read_dof_values(const Number *src)
  {
    static_assert(n_components_ == 1, "This function only supports FE with one \
                  components");
    const unsigned int idx = (threadIdx.x%n_q_points_1d)
                             +(dim>1 ? threadIdx.y : 0)*n_q_points_1d
                             +(dim>2 ? threadIdx.z : 0)*n_q_points_1d*n_q_points_1d;

    const unsigned int src_idx = local_to_global[idx];
    // Use the read-only data cache.
    values[idx] = __ldg(&src[src_idx]);

    if (constraint_mask)
      internal::resolve_hanging_nodes_shmem<dim,fe_degree,false>(values,
                                                                 constraint_mask);

    __syncthreads();
  }



  template <int dim, int fe_degree, int n_q_points_1d, int n_components_,
            typename Number>
  __device__ void
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
  distribute_local_to_global(Number *dst) const
  {
    static_assert(n_components_ == 1, "This function only supports FE with one \
                  components");
    if (constraint_mask)
      internal::resolve_hanging_nodes_shmem<dim,fe_degree,true>(values,
                                                                constraint_mask);


    const unsigned int idx = (threadIdx.x%n_q_points_1d)
                             + (dim>1 ? threadIdx.y : 0) * n_q_points_1d
                             + (dim>2 ? threadIdx.z : 0) * n_q_points_1d * n_q_points_1d;
    const unsigned int destination_idx = local_to_global[idx];

    dst[destination_idx] += values[idx];
  }



  template <int dim, int fe_degree, int n_q_points_1d, int n_components_,
            typename Number>
  __device__ void
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
  evaluate(const bool evaluate_val, const bool evaluate_grad)
  {
    // First evaluate the gradients because it requires values that will be
    // changed if evaluate_val is true
    internal::EvaluatorTensorProduct<internal::EvaluatorVariant::evaluate_general,
             dim, fe_degree,n_q_points_1d, Number> evaluator_tensor_product;
    if (evaluate_grad == true)
      {
        evaluator_tensor_product.gradient_at_quad_pts(values, gradients);
        __syncthreads();
      }

    if (evaluate_val == true)
      {
        evaluator_tensor_product.value_at_quad_pts(values);
        __syncthreads();
      }
  }



  template <int dim, int fe_degree, int n_q_points_1d, int n_components_,
            typename Number>
  __device__ void
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
  integrate(const bool integrate_val, const bool integrate_grad)
  {
    internal::EvaluatorTensorProduct<internal::EvaluatorVariant::evaluate_general,
             dim, fe_degree,n_q_points_1d, Number> evaluator_tensor_product;
    if (integrate_val == true)
      {
        evaluator_tensor_product.integrate_value(values);
        __syncthreads();
        if (integrate_grad == true)
          {
            evaluator_tensor_product.integrate_gradient<true>(values, gradients);
            __syncthreads();
          }
      }
    else if (integrate_grad == true)
      {
        evaluator_tensor_product.integrate_gradient<false>(values, gradients);
        __syncthreads();
      }
  }



  template <int dim, int fe_degree, int n_q_points_1d, int n_components_,
            typename Number>
  __device__
  typename FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::value_type
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
  get_value(const unsigned int q_point) const
  {
    return values[q_point];
  }



  template <int dim, int fe_degree, int n_q_points_1d, int n_components_,
            typename Number>
  __device__ void
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
  submit_value(const value_type &val_in, const unsigned int q_point)
  {
    values[q_point] = val_in * JxW[q_point];
  }



  template <int dim, int fe_degree, int n_q_points_1d, int n_components_,
            typename Number>
  __device__
  typename FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::gradient_type
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
  get_gradient(const unsigned int q_point) const
  {
    static_assert(n_components_ == 1, "This function only supports FE with one \
                  components");
    // TODO optimize if the mesh is uniform
    const Number *inv_jacobian = &inv_jac[q_point];
    gradient_type grad;
    for (int d_1=0; d_1<dim; ++d_1)
      {
        Number tmp = 0.;
        for (int d_2=0; d_2<dim; ++d_2)
          tmp += inv_jacobian[padding_length*n_cells*(dim*d_2+d_1)] *
                 gradients[d_2][q_point];
        grad[d_1] = tmp;
      }

    return grad;
  }



  template <int dim, int fe_degree, int n_q_points_1d, int n_components_,
            typename Number>
  __device__ void
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
  submit_gradient(const gradient_type &grad_in, const unsigned int q_point)
  {
    // TODO optimize if the mesh is uniform
    const Number *inv_jacobian = &inv_jac[q_point];
    for (int d_1=0; d_1<dim; ++d_1)
      {
        Number tmp = 0.;
        for (int d_2=0; d_2<dim; ++d_2)
          tmp += inv_jacobian[n_cells*padding_length*(dim*d_1+d_2)] *
                 grad_in[d_2];
        gradients[d_1][q_point] = tmp * JxW[q_point];
      }
  }



  template <int dim, int fe_degree, int n_q_points_1d, int n_components_,
            typename Number>
  template <typename functor>
  __device__ void
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
  apply_quad_point_operations(const functor &func)
  {
    const unsigned int q_point = (dim == 1 ? threadIdx.x%n_q_points_1d :
                                  dim == 2 ? threadIdx.x%n_q_points_1d + n_q_points_1d *threadIdx.y :
                                  threadIdx.x%n_q_points_1d + n_q_points_1d * (threadIdx.y +
                                      n_q_points_1d*threadIdx.z));
    func(this, q_point);

    __syncthreads();
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif

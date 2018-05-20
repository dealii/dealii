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

#ifndef dealii_cuda_matrix_free_h
#define dealii_cuda_matrix_free_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_CUDA

#  include <cuda_runtime_api.h>
#  include <deal.II/base/quadrature.h>
#  include <deal.II/base/tensor.h>
#  include <deal.II/dofs/dof_handler.h>
#  include <deal.II/fe/fe_update_flags.h>
#  include <deal.II/fe/mapping.h>
#  include <deal.II/fe/mapping_q1.h>
#  include <deal.II/lac/constraint_matrix.h>
#  include <deal.II/lac/cuda_vector.h>

DEAL_II_NAMESPACE_OPEN

namespace CUDAWrappers
{
  // forward declaration
  namespace internal
  {
    template <int dim, typename Number>
    class ReinitHelper;
  }

  /**
   * This class collects all the data that is stored for the matrix free
   * implementation. The storage scheme is tailored towards several loops
   * performed with the same data, i.e., typically doing many matrix-vector
   * products or residual computations on the same mesh.
   *
   * This class does not implement any operations involving finite element basis
   * functions, i.e., regarding the operation performed on the cells. For these
   * operations, the class FEEvaluation is designed to use the data collected in
   * this class.
   *
   * This class implements a loop over all cells (cell_loop()). This loop is
   * scheduled in such a way that cells that cells that share degrees of freedom
   * are not worked on simultaneously, which implies that it is possible to
   * write to vectors in parallel without having to explicitly synchronize
   * access to these vectors and matrices. This class does not implement any
   * shape values, all it does is to cache the respective data. To implement
   * finite element operations, use the class CUDAWrappers::FEEvalutation.
   *
   * This class traverse the cells in a different order than the usual
   * Triangulation class in deal.II.
   *
   * @note Only float and double are supported.
   *
   * @ingroup CUDAWrappers
   */
  template <int dim, typename Number = double>
  class MatrixFree : public Subscriptor
  {
  public:
    typedef Tensor<2, dim, Tensor<1, dim, Number>> jacobian_type;
    // TODO this should really be a CUDAWrappers::Point
    typedef Tensor<1, dim, Number> point_type;

    // Use Number2 so we don't hide the template parameter Number
    template <typename Number2>
    using CUDAVector = ::dealii::LinearAlgebra::CUDAWrappers::Vector<Number2>;

    /**
     * Parallelization scheme used: parallel_in_elem (parallelism at the level
     * of degrees of freedom) or parallel_over_elem (parallelism at the level of
     * cells)
     */
    enum ParallelizationScheme
    {
      parallel_in_elem,
      parallel_over_elem
    };

    struct AdditionalData
    {
      AdditionalData(const ParallelizationScheme parallelization_scheme
                     = parallel_in_elem,
                     const UpdateFlags mapping_update_flags
                     = update_gradients | update_JxW_values)
        : parallelization_scheme(parallelization_scheme),
          mapping_update_flags(mapping_update_flags)
      {}

      /**
       * Number of colors created by the graph coloring algorithm.
       */
      unsigned int n_colors;
      /**
       * Parallelization scheme used, parallelization over degrees of freedom or
       * over cells.
       */
      ParallelizationScheme parallelization_scheme;
      /**
       * This flag is used to determine which quantities should be cached. This
       * class can cache data needed for gradient computations (inverse
       * Jacobians), Jacobian determinants (JxW), quadrature points as well as
       * data for Hessians (derivative of Jacobians). By default, only data for
       * gradients and Jacobian determinants times quadrature weights, JxW, are
       * cached. If quadrature points of second derivatives are needed, they
       * must be specified by this field.
       */
      UpdateFlags mapping_update_flags;
    };

    /**
     * Structure which is passed to the kernel. It is used to pass all the
     * necessary information from the CPU to the GPU.
     */
    struct Data
    {
      point_type*   q_points;
      unsigned int* local_to_global;
      Number*       inv_jacobian;
      Number*       JxW;
      unsigned int  n_cells;
      unsigned int  padding_length;
      unsigned int  row_start;
      unsigned int* constraint_mask;
    };

    /**
     * Default constructor.
     */
    MatrixFree();

    unsigned int
    get_padding_length() const;

    /**
     * Extracts the information needed to perform loops over cells. The
     * DoFHandler and ConstraintMatrix describe the layout of degrees of
     * freedom, the DoFHandler and the mapping describe the transformation from
     * unit to real cell, and the finite element underlying the DoFHandler
     * together with the quadrature formula describe the local operations.
     */
    void
    reinit(const Mapping<dim>&     mapping,
           const DoFHandler<dim>&  dof_handler,
           const ConstraintMatrix& constraints,
           const Quadrature<1>&    quad,
           const AdditionalData    additional_data = AdditionalData());

    /**
     * Initializes the data structures. Same as above but using a Q1 mapping.
     */
    void
    reinit(const DoFHandler<dim>&  dof_handler,
           const ConstraintMatrix& constraints,
           const Quadrature<1>&    quad,
           const AdditionalData    AdditionalData = AdditionalData());

    /**
     * Return the Data structure associated with @p color.
     */
    Data
    get_data(unsigned int color) const;

    /**
     * This method runs the loop over all cells and apply the local operation on
     * each element in parallel. @p func is a functor which is appplied on each color.
     */
    template <typename functor>
    void
    cell_loop(const functor&            func,
              const CUDAVector<Number>& src,
              CUDAVector<Number>&       dst) const;

    void
    copy_constrained_values(const CUDAVector<Number>& src,
                            CUDAVector<Number>&       dst) const;

    void
    set_constrained_values(const Number val, CUDAVector<Number>& dst) const;

    /**
     * Free all the memory allocated.
     */
    void
    free();

    /**
     * Return an approximation of the memory consumption of this class in bytes.
     */
    std::size_t
    memory_consumption() const;

  private:
    /**
     * Parallelization scheme used, parallelization over degrees of freedom or
     * over cells.
     */
    ParallelizationScheme parallelization_scheme;
    /**
     * Degree of the finite element used.
     */
    unsigned int fe_degree;
    /**
     * Number of degrees of freedom per cell.
     */
    unsigned int dofs_per_cell;
    /**
     * Number of constrained degrees of freedom.
     */
    unsigned int n_constrained_dofs;
    /**
     * Number of quadrature points per cells.
     */
    unsigned int q_points_per_cell;
    /**
     * Number of colors produced by the graph coloring algorithm.
     */
    unsigned int n_colors;
    /**
     * Number of cells in each color.
     */
    std::vector<unsigned int> n_cells;
    /**
     * Vector of pointers to the quadrature points associated to the cells of
     * each color.
     */
    std::vector<point_type*> q_points;
    /**
     * Map the position in the local vector to the position in the global
     * vector.
     */
    std::vector<unsigned int*> local_to_global;
    /**
     * Vector of pointer to the inverse Jacobian associated to the cells of each
     * color.
     */
    std::vector<Number*> inv_jacobian;
    /**
     * Vector of pointer to the Jacobian time the weights associated to the
     * cells of each color.
     */
    std::vector<Number*> JxW;

    // Constraints
    unsigned int*              constrained_dofs;
    std::vector<unsigned int*> constraint_mask;
    /**
     * Grid dimensions associated to the different colors. The grid dimensions
     * are used to launch the CUDA kernels.
     */
    std::vector<dim3> grid_dim;
    /**
     * Block dimensions associated to the different colors. The block dimensions
     * are used to launch the CUDA kernels.
     */
    std::vector<dim3> block_dim;

    // Parallelization parameter
    unsigned int cells_per_block;
    dim3         constraint_grid_dim;
    dim3         constraint_block_dim;

    unsigned int              padding_length;
    std::vector<unsigned int> row_start;

    friend class internal::ReinitHelper<dim, Number>;
  };

  // TODO find a better place to put these things
  // Structure to pass the shared memory into a general user function.
  template <int dim, typename Number>
  struct SharedData
  {
    __device__
    SharedData(Number* vd, Number* gq[dim])
      : values(vd)
    {
      for(int d = 0; d < dim; ++d)
        gradients[d] = gq[d];
    }

    Number* values;
    Number* gradients[dim];
  };

  // This function determines the number of cells per block, possibly at compile
  // time
  // TODO this function should be rewritten using meta-programming
  __host__ __device__ constexpr unsigned int
           cells_per_block_shmem(int dim, int fe_degree)
  {
    return dim == 2 ?
             (fe_degree == 1 ?
                32 :
                fe_degree == 2 ? 8 :
                                 fe_degree == 3 ? 4 : fe_degree == 4 ? 4 : 1) :
             dim == 3 ? (fe_degree == 1 ? 8 : fe_degree == 2 ? 2 : 1) : 1;
  }
} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE

#endif

#endif

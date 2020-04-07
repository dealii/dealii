// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2019 by the deal.II authors
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


#ifndef dealii_matrix_free_operators_h
#define dealii_matrix_free_operators_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>


DEAL_II_NAMESPACE_OPEN


namespace MatrixFreeOperators
{
  namespace BlockHelper
  {
    // workaroud for unifying non-block vector and block vector implementations
    // a non-block vector has one block and the only subblock is the vector
    // itself
    template <typename VectorType>
    typename std::enable_if<IsBlockVector<VectorType>::value,
                            unsigned int>::type
    n_blocks(const VectorType &vector)
    {
      return vector.n_blocks();
    }

    template <typename VectorType>
    typename std::enable_if<!IsBlockVector<VectorType>::value,
                            unsigned int>::type
    n_blocks(const VectorType &)
    {
      return 1;
    }

    template <typename VectorType>
    typename std::enable_if<IsBlockVector<VectorType>::value,
                            typename VectorType::BlockType &>::type
    subblock(VectorType &vector, unsigned int block_no)
    {
      return vector.block(block_no);
    }

    template <typename VectorType>
    typename std::enable_if<IsBlockVector<VectorType>::value,
                            const typename VectorType::BlockType &>::type
    subblock(const VectorType &vector, unsigned int block_no)
    {
      AssertIndexRange(block_no, vector.n_blocks());
      return vector.block(block_no);
    }

    template <typename VectorType>
    typename std::enable_if<!IsBlockVector<VectorType>::value,
                            VectorType &>::type
    subblock(VectorType &vector, unsigned int)
    {
      return vector;
    }

    template <typename VectorType>
    typename std::enable_if<!IsBlockVector<VectorType>::value,
                            const VectorType &>::type
    subblock(const VectorType &vector, unsigned int)
    {
      return vector;
    }

    template <typename VectorType>
    typename std::enable_if<IsBlockVector<VectorType>::value, void>::type
    collect_sizes(VectorType &vector)
    {
      vector.collect_sizes();
    }

    template <typename VectorType>
    typename std::enable_if<!IsBlockVector<VectorType>::value, void>::type
    collect_sizes(const VectorType &)
    {}
  } // namespace BlockHelper

  /**
   * Abstract base class for matrix-free operators which can be used both at
   * the finest mesh or at a certain level in geometric multigrid.
   *
   * A derived class has to implement apply_add() method as well as
   * compute_diagonal() to initialize the protected member
   * inverse_diagonal_entries and/or diagonal_entries. In case of a
   * non-symmetric operator, Tapply_add() should be additionally implemented.
   *
   * Currently, the only supported vectors are
   * LinearAlgebra::distributed::Vector and
   * LinearAlgebra::distributed::BlockVector.
   *
   * <h4>Selective use of blocks in MatrixFree</h4>
   *
   * MatrixFree allows to use several DoFHandler/AffineConstraints combinations
   * by passing a std::vector with pointers to the respective objects into
   * the MatrixFree::reinit function. This class supports to select only some
   * of the blocks in the underlying MatrixFree object by optional integer
   * lists that specify the chosen blocks.
   *
   * One application of constructing a matrix-free operator only on selected
   * blocks would be the setting of the step-32 tutorial program: This
   * problem has three <i>blocks</i>, one for the velocity, one for the
   * pressure, and one for temperature. The time lag scheme used for temporal
   * evolution splits the temperature equation away from the Stokes system in
   * velocity and pressure. However, there are cross terms like the velocity
   * that enters the temperature advection-diffusion equation or the
   * temperature that enters the right hand side of the velocity. In order to
   * be sure that MatrixFree uses the same integer indexing to the different
   * blocks, one needs to put all the three blocks into the same MatrixFree
   * object. However, when solving a linear system the operators involved
   * either address the first two in the Stokes solver, or the last one for
   * the temperature solver. In the former case, a BlockVector of two
   * components would be selected with a vector selecting the blocks {0, 1} in
   * MatrixFree, whereas in the latter, a non-block vector selecting the block
   * {2} would be used.
   *
   * A second application of selection is in problems with a Newton-type
   * iteration or problems with inhomogeneous boundary conditions. In such a
   * case, one has to deal with two different sets of constraints: One set of
   * constraints applies to the solution vector which might include hanging
   * node constraints or periodicity constraints but no constraints on
   * inhomogeneous Dirichlet boundaries. Before the nonlinear iteration, the
   * boundary values are set to the expected value in the vector, representing
   * the initial guess. In each iteration of the Newton method, a linear
   * system subject to zero Dirichlet boundary conditions is solved that is
   * then added to the initial guess. This setup can be realized by using a
   * vector of two pointers pointing to the same DoFHandler object and a
   * vector of two pointers to the two AffineConstraints objects. If the first
   * AffineConstraints object is the one including the zero Dirichlet
   * constraints, one would give a std::vector<unsigned int>(1, 0) to the
   * initialize() function, i.e., a vector of length 1 that selects exactly the
   * first AffineConstraints object with index 0.
   *
   * For systems of PDEs where the different blocks of MatrixFree are
   * associated with different physical components of the equations, adding
   * another block with a different AffineConstraints argument solely for the
   * purpose of boundary conditions might lead to cumbersome index
   * handling. Instead, one could set up a second MatrixFree instance with the
   * different AffineConstraints object but the same interpretation of blocks,
   * and use that for interpolating inhomogeneous boundary conditions (see also
   * the discussion in the results section of the step-37 tutorial program):
   *
   * @code
   * matrix_free_inhomogeneous.reinit(dof_handler, constraints_no_dirichlet,
   *                                  quadrature, additional_data);
   * operator_inhomogeneous.initialize(matrix_free_inhomogeneous,
   *                                   selected_blocks);
   * LinearAlgebra::distributed::Vector<double> inhomogeneity;
   * matrix_free_inhomogeneous.initialize_dof_vector(inhomogeneity);
   * constraints_with_dirichlet.distribute(inhomogeneity);
   * operator_inhomogeneous.vmult(system_rhs, inhomogeneity);
   * system_rhs *= -1.;
   * // proceed with other terms from right hand side...
   * @endcode
   *
   * @author Denis Davydov, Daniel Arndt, Martin Kronbichler, 2016, 2017
   */
  template <int dim,
            typename VectorType = LinearAlgebra::distributed::Vector<double>,
            typename VectorizedArrayType =
              VectorizedArray<typename VectorType::value_type>>
  class Base : public Subscriptor
  {
  public:
    /**
     * Number alias.
     */
    using value_type = typename VectorType::value_type;

    /**
     * size_type needed for preconditioner classes.
     */
    using size_type = typename VectorType::size_type;

    /**
     * Default constructor.
     */
    Base();

    /**
     * Virtual destructor.
     */
    virtual ~Base() override = default;

    /**
     * Release all memory and return to a state just like after having called
     * the default constructor.
     */
    virtual void
    clear();

    /**
     * Initialize operator on fine scale.
     *
     * The optional selection vector allows to choose only some components
     * from the underlying MatrixFree object, e.g. just a single one. The
     * entry @p selected_row_blocks[i] in the vector chooses the DoFHandler
     * and AffineConstraints object that was given as the
     * @p selected_row_blocks[i]-th argument to the MatrixFree::reinit() call.
     * Different arguments for rows and columns also make it possible to
     * select non-diagonal blocks or rectangular blocks. If the row vector is
     * empty, all components are selected, otherwise its size must be smaller
     * or equal to MatrixFree::n_components() and all indices need to be
     * unique and within the range of 0 and MatrixFree::n_components(). If the
     * column selection vector is empty, it is taken the same as the row
     * selection, defining a diagonal block.
     */
    void
    initialize(std::shared_ptr<
                 const MatrixFree<dim, value_type, VectorizedArrayType>> data,
               const std::vector<unsigned int> &selected_row_blocks =
                 std::vector<unsigned int>(),
               const std::vector<unsigned int> &selected_column_blocks =
                 std::vector<unsigned int>());

    /**
     * Initialize operator on a level @p level for a single FiniteElement.
     *
     * The optional selection vector allows to choose only some components
     * from the underlying MatrixFree object, e.g. just a single one. The
     * entry @p selected_row_blocks[i] in the vector chooses the DoFHandler
     * and AffineConstraints object that was given as the
     * @p selected_row_blocks[i]-th argument to the MatrixFree::reinit() call.
     * Since a multigrid operator is always associated to inverting a matrix
     * and thus represents a diagonal block, the same vector for rows and
     * columns is used as opposed to the non-level initialization function. If
     * empty, all components are selected.
     */
    void
    initialize(std::shared_ptr<
                 const MatrixFree<dim, value_type, VectorizedArrayType>> data,
               const MGConstrainedDoFs &        mg_constrained_dofs,
               const unsigned int               level,
               const std::vector<unsigned int> &selected_row_blocks =
                 std::vector<unsigned int>());

    /**
     * Initialize operator on a level @p level for multiple FiniteElement
     * objects.
     *
     * The optional selection vector allows to choose only some components
     * from the underlying MatrixFree object, e.g. just a single one. The
     * entry @p selected_row_blocks[i] in the vector chooses the DoFHandler
     * and AffineConstraints object that was given as the
     * @p selected_row_blocks[i]-th argument to the MatrixFree::reinit() call.
     * Since a multigrid operator is always associated to inverting a matrix
     * and thus represents a diagonal block, the same vector for rows and
     * columns is used as opposed to the non-level initialization function. If
     * empty, all components are selected.
     */
    void
    initialize(std::shared_ptr<
                 const MatrixFree<dim, value_type, VectorizedArrayType>> data_,
               const std::vector<MGConstrainedDoFs> &mg_constrained_dofs,
               const unsigned int                    level,
               const std::vector<unsigned int> &     selected_row_blocks =
                 std::vector<unsigned int>());

    /**
     * Return the dimension of the codomain (or range) space.
     */
    size_type
    m() const;

    /**
     * Return the dimension of the domain space.
     */
    size_type
    n() const;

    /**
     * vmult operator for interface.
     */
    void
    vmult_interface_down(VectorType &dst, const VectorType &src) const;

    /**
     * vmult operator for interface.
     */
    void
    vmult_interface_up(VectorType &dst, const VectorType &src) const;

    /**
     * Matrix-vector multiplication.
     */
    void
    vmult(VectorType &dst, const VectorType &src) const;

    /**
     * Transpose matrix-vector multiplication.
     */
    void
    Tvmult(VectorType &dst, const VectorType &src) const;

    /**
     * Adding Matrix-vector multiplication.
     */
    void
    vmult_add(VectorType &dst, const VectorType &src) const;

    /**
     * Adding transpose matrix-vector multiplication.
     */
    void
    Tvmult_add(VectorType &dst, const VectorType &src) const;

    /**
     * Return the value of the matrix entry (row,col). In matrix-free context
     * this function is valid only for row==col when diagonal is initialized.
     */
    value_type
    el(const unsigned int row, const unsigned int col) const;

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    virtual std::size_t
    memory_consumption() const;

    /**
     * A wrapper for initialize_dof_vector() of MatrixFree object.
     */
    void
    initialize_dof_vector(VectorType &vec) const;

    /**
     * Compute diagonal of this operator.
     *
     * A derived class needs to implement this function and resize and fill
     * the protected member inverse_diagonal_entries and/or diagonal_entries
     * accordingly.
     */
    virtual void
    compute_diagonal() = 0;

    /**
     * Get read access to the MatrixFree object stored with this operator.
     */
    std::shared_ptr<const MatrixFree<dim, value_type, VectorizedArrayType>>
    get_matrix_free() const;

    /**
     * Get read access to the inverse diagonal of this operator.
     */
    const std::shared_ptr<DiagonalMatrix<VectorType>> &
    get_matrix_diagonal_inverse() const;

    /**
     * Get read access to the diagonal of this operator.
     */
    const std::shared_ptr<DiagonalMatrix<VectorType>> &
    get_matrix_diagonal() const;


    /**
     * Apply the Jacobi preconditioner, which multiplies every element of the
     * <tt>src</tt> vector by the inverse of the respective diagonal element and
     * multiplies the result with the relaxation factor <tt>omega</tt>.
     */
    void
    precondition_Jacobi(VectorType &      dst,
                        const VectorType &src,
                        const value_type  omega) const;

  protected:
    /**
     * Perform necessary operations related to constraints before calling
     * apply_add() or Tapply_add() inside mult_add().
     */
    void
    preprocess_constraints(VectorType &dst, const VectorType &src) const;

    /**
     * Perform necessary operations related to constraints after calling
     * apply_add() or Tapply_add() inside mult_add().
     */
    void
    postprocess_constraints(VectorType &dst, const VectorType &src) const;

    /**
     * Set constrained entries (both from hanging nodes and edge constraints)
     * of @p dst to one.
     */
    void
    set_constrained_entries_to_one(VectorType &dst) const;

    /**
     * Apply operator to @p src and add result in @p dst.
     */
    virtual void
    apply_add(VectorType &dst, const VectorType &src) const = 0;

    /**
     * Apply transpose operator to @p src and add result in @p dst.
     *
     * Default implementation is to call apply_add().
     */
    virtual void
    Tapply_add(VectorType &dst, const VectorType &src) const;

    /**
     * MatrixFree object to be used with this operator.
     */
    std::shared_ptr<const MatrixFree<dim, value_type, VectorizedArrayType>>
      data;

    /**
     * A shared pointer to a diagonal matrix that stores the
     * diagonal elements as a vector.
     */
    std::shared_ptr<DiagonalMatrix<VectorType>> diagonal_entries;

    /**
     * A shared pointer to a diagonal matrix that stores the inverse of
     * diagonal elements as a vector.
     */
    std::shared_ptr<DiagonalMatrix<VectorType>> inverse_diagonal_entries;

    /**
     * A vector which defines the selection of sub-components of MatrixFree
     * for the rows of the matrix representation.
     */
    std::vector<unsigned int> selected_rows;

    /**
     * A vector which defines the selection of sub-components of MatrixFree
     * for the columns of the matrix representation.
     */
    std::vector<unsigned int> selected_columns;

  private:
    /**
     * Indices of DoFs on edge in case the operator is used in GMG context.
     */
    std::vector<std::vector<unsigned int>> edge_constrained_indices;

    /**
     * Auxiliary vector.
     */
    mutable std::vector<std::vector<std::pair<value_type, value_type>>>
      edge_constrained_values;

    /**
     * A flag which determines whether or not this operator has interface
     * matrices in GMG context.
     */
    bool have_interface_matrices;

    /**
     * Function which implements vmult_add (@p transpose = false) and
     * Tvmult_add (@p transpose = true).
     */
    void
    mult_add(VectorType &      dst,
             const VectorType &src,
             const bool        transpose) const;

    /**
     * Adjust the ghost range of the vectors to the storage requirements of
     * the underlying MatrixFree class. This is used inside the mult_add() as
     * well as vmult_interface_up() and vmult_interface_down() methods in
     * order to ensure that the cell loops will be able to access the ghost
     * indices with the correct local indices.
     */
    void
    adjust_ghost_range_if_necessary(const VectorType &vec,
                                    const bool        is_row) const;
  };



  /**
   * Auxiliary class to provide interface vmult/Tvmult methods required in
   * adaptive geometric multgrids. @p OperatorType class should be derived
   * from MatrixFreeOperators::Base class.
   *
   * The adaptive multigrid realization in deal.II implements an approach
   * called local smoothing. This means that the smoothing on the finest level
   * only covers the local part of the mesh defined by the fixed (finest) grid
   * level and ignores parts of the computational domain where the terminal
   * cells are coarser than this level. As the method progresses to coarser
   * levels, more and more of the global mesh will be covered. At some coarser
   * level, the whole mesh will be covered. Since all level matrices in the
   * multigrid method cover a single level in the mesh, no hanging nodes
   * appear on the level matrices. At the interface between multigrid levels,
   * homogeneous Dirichlet boundary conditions are set while smoothing. When
   * the residual is transferred to the next coarser level, however, the
   * coupling over the multigrid interface needs to be taken into account.
   * This is done by the so-called interface (or edge) matrices that compute
   * the part of the residual that is missed by the level matrix with
   * homogeneous Dirichlet conditions. We refer to the
   * @ref mg_paper "Multigrid paper by Janssen and Kanschat"
   * for more details.
   *
   * For the implementation of those interface matrices, most infrastructure
   * is already in place and provided by MatrixFreeOperators::Base through the
   * two multiplication routines vmult_interface_down() and
   * vmult_interface_up(). The only thing MGInterfaceOperator does is
   * wrapping those operations and make them accessible via
   * @p vmult() and @p Tvmult interface as expected by the multigrid routines
   * (that were originally written for matrices, hence expecting those names).
   * Note that the vmult_interface_down is used during the restriction phase of
   * the multigrid V-cycle, whereas vmult_interface_up is used during the
   * prolongation phase.
   *
   * @author Martin Kronbichler, 2016
   */
  template <typename OperatorType>
  class MGInterfaceOperator : public Subscriptor
  {
  public:
    /**
     * Number alias.
     */
    using value_type = typename OperatorType::value_type;

    /**
     * Size type.
     */
    using size_type = typename OperatorType::size_type;

    /**
     * Default constructor.
     */
    MGInterfaceOperator();

    /**
     * Clear the pointer to the OperatorType object.
     */
    void
    clear();

    /**
     * Initialize this class with an operator @p operator_in.
     */
    void
    initialize(const OperatorType &operator_in);

    /**
     * vmult operator, see class description for more info.
     */
    template <typename VectorType>
    void
    vmult(VectorType &dst, const VectorType &src) const;

    /**
     * Tvmult operator, see class description for more info.
     */
    template <typename VectorType>
    void
    Tvmult(VectorType &dst, const VectorType &src) const;

    /**
     * A wrapper for initialize_dof_vector() of OperatorType object.
     */
    template <typename VectorType>
    void
    initialize_dof_vector(VectorType &vec) const;


  private:
    /**
     * Const pointer to the operator class.
     */
    SmartPointer<const OperatorType> mf_base_operator;
  };



  /**
   * This class implements the operation of the action of the inverse of a
   * mass matrix on an element for the special case of an evaluation object
   * with as many quadrature points as there are cell degrees of freedom. It
   * uses algorithms from FEEvaluation and produces the exact mass matrix for
   * DGQ elements. This algorithm uses tensor products of inverse 1D shape
   * matrices over quadrature points, so the inverse operation is exactly as
   * expensive as applying the forward operator on each cell. Of course, for
   * continuous finite elements this operation does not produce the inverse of
   * a mass operation as the coupling between the elements cannot be
   * considered by this operation.
   *
   * The equation may contain variable coefficients, so the user is required
   * to provide an array for the inverse of the local coefficient (this class
   * provide a helper method 'fill_inverse_JxW_values' to get the inverse of a
   * constant-coefficient operator).
   *
   * @author Martin Kronbichler, 2014
   */
  template <int dim,
            int fe_degree,
            int n_components             = 1,
            typename Number              = double,
            typename VectorizedArrayType = VectorizedArray<Number>>
  class CellwiseInverseMassMatrix
  {
    static_assert(
      std::is_same<Number, typename VectorizedArrayType::value_type>::value,
      "Type of Number and of VectorizedArrayType do not match.");

  public:
    /**
     * Constructor. Initializes the shape information from the ShapeInfo field
     * in the class FEEval.
     */
    CellwiseInverseMassMatrix(
      const FEEvaluationBase<dim,
                             n_components,
                             Number,
                             false,
                             VectorizedArrayType> &fe_eval);

    /**
     * Applies the inverse mass matrix operation on an input array. It is
     * assumed that the passed input and output arrays are of correct size,
     * namely FEEvaluation::dofs_per_cell long. The inverse of the
     * local coefficient (also containing the inverse JxW values) must be
     * passed as first argument. Passing more than one component in the
     * coefficient is allowed.
     */
    void
    apply(const AlignedVector<VectorizedArrayType> &inverse_coefficient,
          const unsigned int                        n_actual_components,
          const VectorizedArrayType *               in_array,
          VectorizedArrayType *                     out_array) const;

    /**
     * Applies the inverse mass matrix operation on an input array, using the
     * inverse of the JxW values provided by the `fe_eval` argument passed to
     * the constructor of this class. Note that the user code must call
     * FEEvaluation::reinit() on the underlying evaluator to make the
     * FEEvaluationBase::JxW() method return the information of the correct
     * cell. It is assumed that the pointers of the input and output arrays
     * are valid over the length FEEvaluation::dofs_per_cell, which is the
     * number of entries processed by this function. The `in_array` and
     * `out_array` arguments may point to the same memory position.
     */
    void
    apply(const VectorizedArrayType *in_array,
          VectorizedArrayType *      out_array) const;

    /**
     * This operation performs a projection from the data given in quadrature
     * points to the actual basis underlying this object. This projection can
     * also be interpreted as a change of the basis from the Lagrange
     * interpolation polynomials in the quadrature points to the
     * basis underlying the current `fe_eval` object.
     *
     * Calling this function on an array as
     * @code
     * inverse_mass.transform_from_q_points_to_basis(1, array,
     *                                               phi.begin_dof_values());
     * @endcode
     * is equivalent to
     * @code
     * for (unsigned int q=0; q<phi.n_q_points; ++q)
     *   phi.submit_value(array[q], q);
     * phi.integrate(true, false);
     * inverse_mass.apply(coefficients, 1, phi.begin_dof_values(),
     *                    phi.begin_dof_values());
     * @endcode
     * provided that @p coefficients holds the inverse of the quadrature
     * weights and no additional coefficients. This setup highlights the
     * underlying projection, testing a right hand side and applying an
     * inverse mass matrix. This function works both for the scalar case as
     * described in the example or for multiple components that are laid out
     * component by component.
     *
     * Compared to the more verbose alternative, the given procedure is
     * considerably faster because it can bypass the @p integrate() step
     * and the first half of the transformation to the quadrature points,
     * reducing the number of tensor product calls from 3*dim*n_components to
     * dim*n_components.
     */
    void
    transform_from_q_points_to_basis(const unsigned int n_actual_components,
                                     const VectorizedArrayType *in_array,
                                     VectorizedArrayType *out_array) const;

    /**
     * Fills the given array with the inverse of the JxW values, i.e., a mass
     * matrix with coefficient 1. Non-unit coefficients must be multiplied (in
     * inverse form) to this array.
     */
    void
    fill_inverse_JxW_values(
      AlignedVector<VectorizedArrayType> &inverse_jxw) const;

  private:
    /**
     * A reference to the FEEvaluation object for getting the JxW_values.
     */
    const FEEvaluationBase<dim,
                           n_components,
                           Number,
                           false,
                           VectorizedArrayType> &fe_eval;
  };



  /**
   * This class implements the operation of the action of a mass matrix.
   *
   * Note that this class only supports the non-blocked vector variant of the
   * Base operator because only a single FEEvaluation object is used in the
   * apply function.
   *
   * @author Daniel Arndt, 2016
   */
  template <int dim,
            int fe_degree,
            int n_q_points_1d   = fe_degree + 1,
            int n_components    = 1,
            typename VectorType = LinearAlgebra::distributed::Vector<double>,
            typename VectorizedArrayType =
              VectorizedArray<typename VectorType::value_type>>
  class MassOperator : public Base<dim, VectorType, VectorizedArrayType>
  {
  public:
    /**
     * Number alias.
     */
    using value_type =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;

    /**
     * size_type needed for preconditioner classes.
     */
    using size_type =
      typename Base<dim, VectorType, VectorizedArrayType>::size_type;

    /**
     * Constructor.
     */
    MassOperator();

    /**
     * For preconditioning, we store a lumped mass matrix at the diagonal
     * entries.
     */
    virtual void
    compute_diagonal() override;

  private:
    /**
     * Applies the mass matrix operation on an input vector. It is
     * assumed that the passed input and output vector are correctly initialized
     * using initialize_dof_vector().
     */
    virtual void
    apply_add(VectorType &dst, const VectorType &src) const override;

    /**
     * For this operator, there is just a cell contribution.
     */
    void
    local_apply_cell(
      const MatrixFree<dim, value_type, VectorizedArrayType> &data,
      VectorType &                                            dst,
      const VectorType &                                      src,
      const std::pair<unsigned int, unsigned int> &           cell_range) const;
  };



  /**
   * This class implements the operation of the action of a Laplace matrix,
   * namely $ L_{ij} = \int_\Omega c(\mathbf x) \mathbf \nabla N_i(\mathbf x)
   * \cdot \mathbf \nabla N_j(\mathbf x)\,d \mathbf x$, where $c(\mathbf x)$ is
   * the scalar heterogeneity coefficient.
   *
   * Note that this class only supports the non-blocked vector variant of the
   * Base operator because only a single FEEvaluation object is used in the
   * apply function.
   *
   * @author Denis Davydov, 2016
   */
  template <int dim,
            int fe_degree,
            int n_q_points_1d   = fe_degree + 1,
            int n_components    = 1,
            typename VectorType = LinearAlgebra::distributed::Vector<double>,
            typename VectorizedArrayType =
              VectorizedArray<typename VectorType::value_type>>
  class LaplaceOperator : public Base<dim, VectorType, VectorizedArrayType>
  {
  public:
    /**
     * Number alias.
     */
    using value_type =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;

    /**
     * size_type needed for preconditioner classes.
     */
    using size_type =
      typename Base<dim, VectorType, VectorizedArrayType>::size_type;

    /**
     * Constructor.
     */
    LaplaceOperator();

    /**
     * The diagonal is approximated by computing a local diagonal matrix per
     * element and distributing it to the global diagonal. This will lead to
     * wrong results on element with hanging nodes but is still an acceptable
     * approximation to be used in preconditioners.
     */
    virtual void
    compute_diagonal();

    /**
     * Set the heterogeneous scalar coefficient @p scalar_coefficient to be used at
     * the quadrature points. The Table should be of correct size, consistent
     * with the total number of quadrature points in
     * <code>dim</code>-dimensions,
     * controlled by the @p n_q_points_1d template parameter. Here,
     * <code>(*scalar_coefficient)(cell,q)</code> corresponds to the value of
     * the coefficient, where <code>cell</code> is an index into a set of cell
     * batches as administered by the MatrixFree framework (which does not work
     * on individual cells, but instead of batches of cells at once), and
     * <code>q</code> is the number of the quadrature point within this batch.
     *
     * Such tables can be initialized by
     * @code
     * std::shared_ptr<Table<2, VectorizedArray<double> > > coefficient;
     * coefficient = std::make_shared<Table<2, VectorizedArray<double> > >();
     * {
     *   FEEvaluation<dim,fe_degree,n_q_points_1d,1,double> fe_eval(mf_data);
     *   const unsigned int n_cells = mf_data.n_macro_cells();
     *   const unsigned int n_q_points = fe_eval.n_q_points;
     *   coefficient->reinit(n_cells, n_q_points);
     *   for (unsigned int cell=0; cell<n_cells; ++cell)
     *     {
     *       fe_eval.reinit(cell);
     *       for (unsigned int q=0; q<n_q_points; ++q)
     *         (*coefficient)(cell,q) =
     *           function.value(fe_eval.quadrature_point(q));
     *     }
     * }
     * @endcode
     * where <code>mf_data</code> is a MatrixFree object and
     * <code>function</code> is a function which provides the following method
     * <code>VectorizedArray<double> value(const Point<dim,
     * VectorizedArray<double> > &p_vec)</code>.
     *
     * If this function is not called, the coefficient is assumed to be unity.
     *
     * The argument to this function is a shared pointer to such a table. The
     * class stores the shared pointer to this table, not a deep copy
     * and uses it to form the Laplace matrix. Consequently, you can update the
     * table and re-use the current object to obtain the action of a Laplace
     * matrix with this updated coefficient. Alternatively, if the table values
     * are only to be filled once, the original shared pointer can also go out
     * of scope in user code and the clear() command or destructor of this class
     * will delete the table.
     */
    void
    set_coefficient(
      const std::shared_ptr<Table<2, VectorizedArrayType>> &scalar_coefficient);

    virtual void
    clear();

    /**
     * Read/Write access to coefficients to be used in Laplace operator.
     *
     * The function will throw an error if coefficients are not previously set
     * by set_coefficient() function.
     */
    std::shared_ptr<Table<2, VectorizedArrayType>>
    get_coefficient();

  private:
    /**
     * Applies the laplace matrix operation on an input vector. It is
     * assumed that the passed input and output vector are correctly initialized
     * using initialize_dof_vector().
     */
    virtual void
    apply_add(VectorType &dst, const VectorType &src) const;

    /**
     * Applies the Laplace operator on a cell.
     */
    void
    local_apply_cell(
      const MatrixFree<dim, value_type, VectorizedArrayType> &data,
      VectorType &                                            dst,
      const VectorType &                                      src,
      const std::pair<unsigned int, unsigned int> &           cell_range) const;

    /**
     * Apply diagonal part of the Laplace operator on a cell.
     */
    void
    local_diagonal_cell(
      const MatrixFree<dim, value_type, VectorizedArrayType> &data,
      VectorType &                                            dst,
      const VectorType &,
      const std::pair<unsigned int, unsigned int> &cell_range) const;

    /**
     * Apply Laplace operator on a cell @p cell.
     */
    void
    do_operation_on_cell(
      FEEvaluation<dim, fe_degree, n_q_points_1d, n_components, value_type>
        &                phi,
      const unsigned int cell) const;

    /**
     * User-provided heterogeneity coefficient.
     */
    std::shared_ptr<Table<2, VectorizedArrayType>> scalar_coefficient;
  };



  // ------------------------------------ inline functions ---------------------

  template <int dim,
            int fe_degree,
            int n_components,
            typename Number,
            typename VectorizedArrayType>
  inline CellwiseInverseMassMatrix<dim,
                                   fe_degree,
                                   n_components,
                                   Number,
                                   VectorizedArrayType>::
    CellwiseInverseMassMatrix(
      const FEEvaluationBase<dim,
                             n_components,
                             Number,
                             false,
                             VectorizedArrayType> &fe_eval)
    : fe_eval(fe_eval)
  {}



  template <int dim,
            int fe_degree,
            int n_components,
            typename Number,
            typename VectorizedArrayType>
  inline void
  CellwiseInverseMassMatrix<dim,
                            fe_degree,
                            n_components,
                            Number,
                            VectorizedArrayType>::
    fill_inverse_JxW_values(
      AlignedVector<VectorizedArrayType> &inverse_jxw) const
  {
    constexpr unsigned int dofs_per_component_on_cell =
      Utilities::pow(fe_degree + 1, dim);
    Assert(inverse_jxw.size() > 0 &&
             inverse_jxw.size() % dofs_per_component_on_cell == 0,
           ExcMessage(
             "Expected diagonal to be a multiple of scalar dof per cells"));

    // temporarily reduce size of inverse_jxw to dofs_per_cell to get JxW values
    // from fe_eval (will not reallocate any memory)
    for (unsigned int q = 0; q < dofs_per_component_on_cell; ++q)
      inverse_jxw[q] = 1. / fe_eval.JxW(q);
    // copy values to rest of vector
    for (unsigned int q = dofs_per_component_on_cell; q < inverse_jxw.size();)
      for (unsigned int i = 0; i < dofs_per_component_on_cell; ++i, ++q)
        inverse_jxw[q] = inverse_jxw[i];
  }



  template <int dim,
            int fe_degree,
            int n_components,
            typename Number,
            typename VectorizedArrayType>
  inline void
  CellwiseInverseMassMatrix<
    dim,
    fe_degree,
    n_components,
    Number,
    VectorizedArrayType>::apply(const VectorizedArrayType *in_array,
                                VectorizedArrayType *      out_array) const
  {
    internal::CellwiseInverseMassMatrixImpl<
      dim,
      fe_degree,
      n_components,
      VectorizedArrayType>::apply(fe_eval, in_array, out_array);
  }



  template <int dim,
            int fe_degree,
            int n_components,
            typename Number,
            typename VectorizedArrayType>
  inline void
  CellwiseInverseMassMatrix<dim,
                            fe_degree,
                            n_components,
                            Number,
                            VectorizedArrayType>::
    apply(const AlignedVector<VectorizedArrayType> &inverse_coefficients,
          const unsigned int                        n_actual_components,
          const VectorizedArrayType *               in_array,
          VectorizedArrayType *                     out_array) const
  {
    internal::CellwiseInverseMassMatrixImpl<dim,
                                            fe_degree,
                                            n_components,
                                            VectorizedArrayType>::
      apply(fe_eval.get_shape_info().data.front().inverse_shape_values_eo,
            inverse_coefficients,
            n_actual_components,
            in_array,
            out_array);
  }



  template <int dim,
            int fe_degree,
            int n_components,
            typename Number,
            typename VectorizedArrayType>
  inline void
  CellwiseInverseMassMatrix<dim,
                            fe_degree,
                            n_components,
                            Number,
                            VectorizedArrayType>::
    transform_from_q_points_to_basis(const unsigned int n_actual_components,
                                     const VectorizedArrayType *in_array,
                                     VectorizedArrayType *      out_array) const
  {
    internal::CellwiseInverseMassMatrixImpl<dim,
                                            fe_degree,
                                            n_components,
                                            VectorizedArrayType>::
      transform_from_q_points_to_basis(
        fe_eval.get_shape_info().data.front().inverse_shape_values_eo,
        n_actual_components,
        in_array,
        out_array);
  }



  //----------------- Base operator -----------------------------
  template <int dim, typename VectorType, typename VectorizedArrayType>
  Base<dim, VectorType, VectorizedArrayType>::Base()
    : Subscriptor()
    , have_interface_matrices(false)
  {}



  template <int dim, typename VectorType, typename VectorizedArrayType>
  typename Base<dim, VectorType, VectorizedArrayType>::size_type
  Base<dim, VectorType, VectorizedArrayType>::m() const
  {
    Assert(data.get() != nullptr, ExcNotInitialized());
    typename Base<dim, VectorType, VectorizedArrayType>::size_type total_size =
      0;
    for (const unsigned int selected_row : selected_rows)
      total_size += data->get_vector_partitioner(selected_row)->size();
    return total_size;
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  typename Base<dim, VectorType, VectorizedArrayType>::size_type
  Base<dim, VectorType, VectorizedArrayType>::n() const
  {
    Assert(data.get() != nullptr, ExcNotInitialized());
    typename Base<dim, VectorType, VectorizedArrayType>::size_type total_size =
      0;
    for (const unsigned int selected_column : selected_columns)
      total_size += data->get_vector_partitioner(selected_column)->size();
    return total_size;
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::clear()
  {
    data.reset();
    inverse_diagonal_entries.reset();
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  typename Base<dim, VectorType, VectorizedArrayType>::value_type
  Base<dim, VectorType, VectorizedArrayType>::el(const unsigned int row,
                                                 const unsigned int col) const
  {
    (void)col;
    Assert(row == col, ExcNotImplemented());
    Assert(inverse_diagonal_entries.get() != nullptr &&
             inverse_diagonal_entries->m() > 0,
           ExcNotInitialized());
    return 1.0 / (*inverse_diagonal_entries)(row, row);
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::initialize_dof_vector(
    VectorType &vec) const
  {
    Assert(data.get() != nullptr, ExcNotInitialized());
    AssertDimension(BlockHelper::n_blocks(vec), selected_rows.size());
    for (unsigned int i = 0; i < BlockHelper::n_blocks(vec); ++i)
      {
        const unsigned int index = selected_rows[i];
        if (!BlockHelper::subblock(vec, index)
               .partitioners_are_compatible(
                 *data->get_dof_info(index).vector_partitioner))
          data->initialize_dof_vector(BlockHelper::subblock(vec, index), index);

        Assert(BlockHelper::subblock(vec, index)
                 .partitioners_are_globally_compatible(
                   *data->get_dof_info(index).vector_partitioner),
               ExcInternalError());
      }
    BlockHelper::collect_sizes(vec);
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::initialize(
    std::shared_ptr<const MatrixFree<dim, value_type, VectorizedArrayType>>
                                     data_,
    const std::vector<unsigned int> &given_row_selection,
    const std::vector<unsigned int> &given_column_selection)
  {
    data = data_;

    selected_rows.clear();
    selected_columns.clear();
    if (given_row_selection.empty())
      for (unsigned int i = 0; i < data_->n_components(); ++i)
        selected_rows.push_back(i);
    else
      {
        for (unsigned int i = 0; i < given_row_selection.size(); ++i)
          {
            AssertIndexRange(given_row_selection[i], data_->n_components());
            for (unsigned int j = 0; j < given_row_selection.size(); ++j)
              if (j != i)
                Assert(given_row_selection[j] != given_row_selection[i],
                       ExcMessage("Given row indices must be unique"));

            selected_rows.push_back(given_row_selection[i]);
          }
      }
    if (given_column_selection.size() == 0)
      selected_columns = selected_rows;
    else
      {
        for (unsigned int i = 0; i < given_column_selection.size(); ++i)
          {
            AssertIndexRange(given_column_selection[i], data_->n_components());
            for (unsigned int j = 0; j < given_column_selection.size(); ++j)
              if (j != i)
                Assert(given_column_selection[j] != given_column_selection[i],
                       ExcMessage("Given column indices must be unique"));

            selected_columns.push_back(given_column_selection[i]);
          }
      }

    edge_constrained_indices.clear();
    edge_constrained_indices.resize(selected_rows.size());
    edge_constrained_values.clear();
    edge_constrained_values.resize(selected_rows.size());
    have_interface_matrices = false;
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::initialize(
    std::shared_ptr<const MatrixFree<dim, value_type, VectorizedArrayType>>
                                     data_,
    const MGConstrainedDoFs &        mg_constrained_dofs,
    const unsigned int               level,
    const std::vector<unsigned int> &given_row_selection)
  {
    std::vector<MGConstrainedDoFs> mg_constrained_dofs_vector(
      1, mg_constrained_dofs);
    initialize(data_, mg_constrained_dofs_vector, level, given_row_selection);
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::initialize(
    std::shared_ptr<const MatrixFree<dim, value_type, VectorizedArrayType>>
                                          data_,
    const std::vector<MGConstrainedDoFs> &mg_constrained_dofs,
    const unsigned int                    level,
    const std::vector<unsigned int> &     given_row_selection)
  {
    AssertThrow(level != numbers::invalid_unsigned_int,
                ExcMessage("level is not set"));

    selected_rows.clear();
    selected_columns.clear();
    if (given_row_selection.empty())
      for (unsigned int i = 0; i < data_->n_components(); ++i)
        selected_rows.push_back(i);
    else
      {
        for (unsigned int i = 0; i < given_row_selection.size(); ++i)
          {
            AssertIndexRange(given_row_selection[i], data_->n_components());
            for (unsigned int j = 0; j < given_row_selection.size(); ++j)
              if (j != i)
                Assert(given_row_selection[j] != given_row_selection[i],
                       ExcMessage("Given row indices must be unique"));

            selected_rows.push_back(given_row_selection[i]);
          }
      }
    selected_columns = selected_rows;

    AssertDimension(mg_constrained_dofs.size(), selected_rows.size());
    edge_constrained_indices.clear();
    edge_constrained_indices.resize(selected_rows.size());
    edge_constrained_values.clear();
    edge_constrained_values.resize(selected_rows.size());

    data = data_;

    for (unsigned int j = 0; j < selected_rows.size(); ++j)
      {
        if (data_->n_macro_cells() > 0)
          {
            AssertDimension(level, data_->get_cell_iterator(0, 0, j)->level());
          }

        // setup edge_constrained indices
        std::vector<types::global_dof_index> interface_indices;
        mg_constrained_dofs[j]
          .get_refinement_edge_indices(level)
          .fill_index_vector(interface_indices);
        edge_constrained_indices[j].clear();
        edge_constrained_indices[j].reserve(interface_indices.size());
        edge_constrained_values[j].resize(interface_indices.size());
        const IndexSet &locally_owned =
          data->get_dof_handler(selected_rows[j]).locally_owned_mg_dofs(level);
        for (const auto interface_index : interface_indices)
          if (locally_owned.is_element(interface_index))
            edge_constrained_indices[j].push_back(
              locally_owned.index_within_set(interface_index));
        have_interface_matrices |=
          Utilities::MPI::max(
            static_cast<unsigned int>(edge_constrained_indices[j].size()),
            data->get_vector_partitioner()->get_mpi_communicator()) > 0;
      }
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::set_constrained_entries_to_one(
    VectorType &dst) const
  {
    for (unsigned int j = 0; j < BlockHelper::n_blocks(dst); ++j)
      {
        const std::vector<unsigned int> &constrained_dofs =
          data->get_constrained_dofs(selected_rows[j]);
        for (const auto constrained_dof : constrained_dofs)
          BlockHelper::subblock(dst, j).local_element(constrained_dof) = 1.;
        for (unsigned int i = 0; i < edge_constrained_indices[j].size(); ++i)
          BlockHelper::subblock(dst, j).local_element(
            edge_constrained_indices[j][i]) = 1.;
      }
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::vmult(VectorType &      dst,
                                                    const VectorType &src) const
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;
    dst = Number(0.);
    vmult_add(dst, src);
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::vmult_add(
    VectorType &      dst,
    const VectorType &src) const
  {
    mult_add(dst, src, false);
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::Tvmult_add(
    VectorType &      dst,
    const VectorType &src) const
  {
    mult_add(dst, src, true);
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::adjust_ghost_range_if_necessary(
    const VectorType &src,
    const bool        is_row) const
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;
    for (unsigned int i = 0; i < BlockHelper::n_blocks(src); ++i)
      {
        const unsigned int mf_component =
          is_row ? selected_rows[i] : selected_columns[i];
        // If both vectors use the same partitioner -> done
        if (BlockHelper::subblock(src, i).get_partitioner().get() ==
            data->get_dof_info(mf_component).vector_partitioner.get())
          continue;

        // If not, assert that the local ranges are the same and reset to the
        // current partitioner
        Assert(
          BlockHelper::subblock(src, i).get_partitioner()->local_size() ==
            data->get_dof_info(mf_component).vector_partitioner->local_size(),
          ExcMessage("The vector passed to the vmult() function does not have "
                     "the correct size for compatibility with MatrixFree."));

        // copy the vector content to a temporary vector so that it does not get
        // lost
        LinearAlgebra::distributed::Vector<Number> copy_vec(
          BlockHelper::subblock(src, i));
        BlockHelper::subblock(const_cast<VectorType &>(src), i)
          .reinit(data->get_dof_info(mf_component).vector_partitioner);
        BlockHelper::subblock(const_cast<VectorType &>(src), i)
          .copy_locally_owned_data_from(copy_vec);
      }
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::preprocess_constraints(
    VectorType &      dst,
    const VectorType &src) const
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;
    adjust_ghost_range_if_necessary(src, false);
    adjust_ghost_range_if_necessary(dst, true);

    // set zero Dirichlet values on the input vector (and remember the src and
    // dst values because we need to reset them at the end)
    for (unsigned int j = 0; j < BlockHelper::n_blocks(dst); ++j)
      {
        for (unsigned int i = 0; i < edge_constrained_indices[j].size(); ++i)
          {
            edge_constrained_values[j][i] = std::pair<Number, Number>(
              BlockHelper::subblock(src, j).local_element(
                edge_constrained_indices[j][i]),
              BlockHelper::subblock(dst, j).local_element(
                edge_constrained_indices[j][i]));
            BlockHelper::subblock(const_cast<VectorType &>(src), j)
              .local_element(edge_constrained_indices[j][i]) = 0.;
          }
      }
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::mult_add(
    VectorType &      dst,
    const VectorType &src,
    const bool        transpose) const
  {
    AssertDimension(dst.size(), src.size());
    AssertDimension(BlockHelper::n_blocks(dst), BlockHelper::n_blocks(src));
    AssertDimension(BlockHelper::n_blocks(dst), selected_rows.size());
    preprocess_constraints(dst, src);
    if (transpose)
      Tapply_add(dst, src);
    else
      apply_add(dst, src);
    postprocess_constraints(dst, src);
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::postprocess_constraints(
    VectorType &      dst,
    const VectorType &src) const
  {
    for (unsigned int j = 0; j < BlockHelper::n_blocks(dst); ++j)
      {
        const std::vector<unsigned int> &constrained_dofs =
          data->get_constrained_dofs(selected_rows[j]);
        for (const auto constrained_dof : constrained_dofs)
          BlockHelper::subblock(dst, j).local_element(constrained_dof) +=
            BlockHelper::subblock(src, j).local_element(constrained_dof);
      }

    // reset edge constrained values, multiply by unit matrix and add into
    // destination
    for (unsigned int j = 0; j < BlockHelper::n_blocks(dst); ++j)
      {
        for (unsigned int i = 0; i < edge_constrained_indices[j].size(); ++i)
          {
            BlockHelper::subblock(const_cast<VectorType &>(src), j)
              .local_element(edge_constrained_indices[j][i]) =
              edge_constrained_values[j][i].first;
            BlockHelper::subblock(dst, j).local_element(
              edge_constrained_indices[j][i]) =
              edge_constrained_values[j][i].second +
              edge_constrained_values[j][i].first;
          }
      }
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::vmult_interface_down(
    VectorType &      dst,
    const VectorType &src) const
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;
    AssertDimension(dst.size(), src.size());
    adjust_ghost_range_if_necessary(src, false);
    adjust_ghost_range_if_necessary(dst, true);

    dst = Number(0.);

    if (!have_interface_matrices)
      return;

    // set zero Dirichlet values on the input vector (and remember the src and
    // dst values because we need to reset them at the end)
    for (unsigned int j = 0; j < BlockHelper::n_blocks(dst); ++j)
      for (unsigned int i = 0; i < edge_constrained_indices[j].size(); ++i)
        {
          edge_constrained_values[j][i] = std::pair<Number, Number>(
            BlockHelper::subblock(src, j).local_element(
              edge_constrained_indices[j][i]),
            BlockHelper::subblock(dst, j).local_element(
              edge_constrained_indices[j][i]));
          BlockHelper::subblock(const_cast<VectorType &>(src), j)
            .local_element(edge_constrained_indices[j][i]) = 0.;
        }

    apply_add(dst, src);

    for (unsigned int j = 0; j < BlockHelper::n_blocks(dst); ++j)
      {
        unsigned int c = 0;
        for (unsigned int i = 0; i < edge_constrained_indices[j].size(); ++i)
          {
            for (; c < edge_constrained_indices[j][i]; ++c)
              BlockHelper::subblock(dst, j).local_element(c) = 0.;
            ++c;

            // reset the src values
            BlockHelper::subblock(const_cast<VectorType &>(src), j)
              .local_element(edge_constrained_indices[j][i]) =
              edge_constrained_values[j][i].first;
          }
        for (; c < BlockHelper::subblock(dst, j).local_size(); ++c)
          BlockHelper::subblock(dst, j).local_element(c) = 0.;
      }
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::vmult_interface_up(
    VectorType &      dst,
    const VectorType &src) const
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;
    AssertDimension(dst.size(), src.size());
    adjust_ghost_range_if_necessary(src, false);
    adjust_ghost_range_if_necessary(dst, true);

    dst = Number(0.);

    if (!have_interface_matrices)
      return;

    VectorType src_cpy(src);
    for (unsigned int j = 0; j < BlockHelper::n_blocks(dst); ++j)
      {
        unsigned int c = 0;
        for (unsigned int i = 0; i < edge_constrained_indices[j].size(); ++i)
          {
            for (; c < edge_constrained_indices[j][i]; ++c)
              BlockHelper::subblock(src_cpy, j).local_element(c) = 0.;
            ++c;
          }
        for (; c < BlockHelper::subblock(src_cpy, j).local_size(); ++c)
          BlockHelper::subblock(src_cpy, j).local_element(c) = 0.;
      }

    apply_add(dst, src_cpy);

    for (unsigned int j = 0; j < BlockHelper::n_blocks(dst); ++j)
      for (unsigned int i = 0; i < edge_constrained_indices[j].size(); ++i)
        BlockHelper::subblock(dst, j).local_element(
          edge_constrained_indices[j][i]) = 0.;
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::Tvmult(
    VectorType &      dst,
    const VectorType &src) const
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;
    dst = Number(0.);
    Tvmult_add(dst, src);
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  std::size_t
  Base<dim, VectorType, VectorizedArrayType>::memory_consumption() const
  {
    return inverse_diagonal_entries.get() != nullptr ?
             inverse_diagonal_entries->memory_consumption() :
             sizeof(*this);
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  std::shared_ptr<const MatrixFree<
    dim,
    typename Base<dim, VectorType, VectorizedArrayType>::value_type,
    VectorizedArrayType>>
  Base<dim, VectorType, VectorizedArrayType>::get_matrix_free() const
  {
    return data;
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  const std::shared_ptr<DiagonalMatrix<VectorType>> &
  Base<dim, VectorType, VectorizedArrayType>::get_matrix_diagonal_inverse()
    const
  {
    Assert(inverse_diagonal_entries.get() != nullptr &&
             inverse_diagonal_entries->m() > 0,
           ExcNotInitialized());
    return inverse_diagonal_entries;
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  const std::shared_ptr<DiagonalMatrix<VectorType>> &
  Base<dim, VectorType, VectorizedArrayType>::get_matrix_diagonal() const
  {
    Assert(diagonal_entries.get() != nullptr && diagonal_entries->m() > 0,
           ExcNotInitialized());
    return diagonal_entries;
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::Tapply_add(
    VectorType &      dst,
    const VectorType &src) const
  {
    apply_add(dst, src);
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::precondition_Jacobi(
    VectorType &                                                          dst,
    const VectorType &                                                    src,
    const typename Base<dim, VectorType, VectorizedArrayType>::value_type omega)
    const
  {
    Assert(inverse_diagonal_entries.get() && inverse_diagonal_entries->m() > 0,
           ExcNotInitialized());
    inverse_diagonal_entries->vmult(dst, src);
    dst *= omega;
  }



  //------------------------- MGInterfaceOperator ------------------------------

  template <typename OperatorType>
  MGInterfaceOperator<OperatorType>::MGInterfaceOperator()
    : Subscriptor()
    , mf_base_operator(nullptr)
  {}



  template <typename OperatorType>
  void
  MGInterfaceOperator<OperatorType>::clear()
  {
    mf_base_operator = nullptr;
  }



  template <typename OperatorType>
  void
  MGInterfaceOperator<OperatorType>::initialize(const OperatorType &operator_in)
  {
    mf_base_operator = &operator_in;
  }



  template <typename OperatorType>
  template <typename VectorType>
  void
  MGInterfaceOperator<OperatorType>::vmult(VectorType &      dst,
                                           const VectorType &src) const
  {
#ifndef DEAL_II_MSVC
    static_assert(
      std::is_same<typename VectorType::value_type, value_type>::value,
      "The vector type must be based on the same value type as this "
      "operator");
#endif

    Assert(mf_base_operator != nullptr, ExcNotInitialized());

    mf_base_operator->vmult_interface_down(dst, src);
  }



  template <typename OperatorType>
  template <typename VectorType>
  void
  MGInterfaceOperator<OperatorType>::Tvmult(VectorType &      dst,
                                            const VectorType &src) const
  {
#ifndef DEAL_II_MSVC
    static_assert(
      std::is_same<typename VectorType::value_type, value_type>::value,
      "The vector type must be based on the same value type as this "
      "operator");
#endif

    Assert(mf_base_operator != nullptr, ExcNotInitialized());

    mf_base_operator->vmult_interface_up(dst, src);
  }



  template <typename OperatorType>
  template <typename VectorType>
  void
  MGInterfaceOperator<OperatorType>::initialize_dof_vector(
    VectorType &vec) const
  {
    Assert(mf_base_operator != nullptr, ExcNotInitialized());

    mf_base_operator->initialize_dof_vector(vec);
  }



  //-----------------------------MassOperator----------------------------------

  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  MassOperator<dim,
               fe_degree,
               n_q_points_1d,
               n_components,
               VectorType,
               VectorizedArrayType>::MassOperator()
    : Base<dim, VectorType, VectorizedArrayType>()
  {}



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  void
  MassOperator<dim,
               fe_degree,
               n_q_points_1d,
               n_components,
               VectorType,
               VectorizedArrayType>::compute_diagonal()
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;
    Assert((Base<dim, VectorType, VectorizedArrayType>::data.get() != nullptr),
           ExcNotInitialized());

    this->inverse_diagonal_entries =
      std::make_shared<DiagonalMatrix<VectorType>>();
    this->diagonal_entries = std::make_shared<DiagonalMatrix<VectorType>>();
    VectorType &inverse_diagonal_vector =
      this->inverse_diagonal_entries->get_vector();
    VectorType &diagonal_vector = this->diagonal_entries->get_vector();
    this->initialize_dof_vector(inverse_diagonal_vector);
    this->initialize_dof_vector(diagonal_vector);
    inverse_diagonal_vector = Number(1.);
    apply_add(diagonal_vector, inverse_diagonal_vector);

    this->set_constrained_entries_to_one(diagonal_vector);
    inverse_diagonal_vector = diagonal_vector;

    const unsigned int local_size = inverse_diagonal_vector.local_size();
    for (unsigned int i = 0; i < local_size; ++i)
      inverse_diagonal_vector.local_element(i) =
        Number(1.) / inverse_diagonal_vector.local_element(i);

    inverse_diagonal_vector.update_ghost_values();
    diagonal_vector.update_ghost_values();
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  void
  MassOperator<dim,
               fe_degree,
               n_q_points_1d,
               n_components,
               VectorType,
               VectorizedArrayType>::apply_add(VectorType &      dst,
                                               const VectorType &src) const
  {
    Base<dim, VectorType, VectorizedArrayType>::data->cell_loop(
      &MassOperator::local_apply_cell, this, dst, src);
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  void
  MassOperator<dim,
               fe_degree,
               n_q_points_1d,
               n_components,
               VectorType,
               VectorizedArrayType>::
    local_apply_cell(
      const MatrixFree<
        dim,
        typename Base<dim, VectorType, VectorizedArrayType>::value_type,
        VectorizedArrayType> &                     data,
      VectorType &                                 dst,
      const VectorType &                           src,
      const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;
    FEEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components,
                 Number,
                 VectorizedArrayType>
      phi(data, this->selected_rows[0]);
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values(src);
        phi.evaluate(true, false, false);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_value(phi.get_value(q), q);
        phi.integrate(true, false);
        phi.distribute_local_to_global(dst);
      }
  }


  //-----------------------------LaplaceOperator----------------------------------

  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  LaplaceOperator<dim,
                  fe_degree,
                  n_q_points_1d,
                  n_components,
                  VectorType,
                  VectorizedArrayType>::LaplaceOperator()
    : Base<dim, VectorType, VectorizedArrayType>()
  {}



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  void
  LaplaceOperator<dim,
                  fe_degree,
                  n_q_points_1d,
                  n_components,
                  VectorType,
                  VectorizedArrayType>::clear()
  {
    Base<dim, VectorType, VectorizedArrayType>::clear();
    scalar_coefficient.reset();
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  void
  LaplaceOperator<dim,
                  fe_degree,
                  n_q_points_1d,
                  n_components,
                  VectorType,
                  VectorizedArrayType>::
    set_coefficient(
      const std::shared_ptr<Table<2, VectorizedArrayType>> &scalar_coefficient_)
  {
    scalar_coefficient = scalar_coefficient_;
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  std::shared_ptr<Table<2, VectorizedArrayType>>
  LaplaceOperator<dim,
                  fe_degree,
                  n_q_points_1d,
                  n_components,
                  VectorType,
                  VectorizedArrayType>::get_coefficient()
  {
    Assert(scalar_coefficient.get(), ExcNotInitialized());
    return scalar_coefficient;
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  void
  LaplaceOperator<dim,
                  fe_degree,
                  n_q_points_1d,
                  n_components,
                  VectorType,
                  VectorizedArrayType>::compute_diagonal()
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;
    Assert((Base<dim, VectorType, VectorizedArrayType>::data.get() != nullptr),
           ExcNotInitialized());

    this->inverse_diagonal_entries =
      std::make_shared<DiagonalMatrix<VectorType>>();
    this->diagonal_entries = std::make_shared<DiagonalMatrix<VectorType>>();
    VectorType &inverse_diagonal_vector =
      this->inverse_diagonal_entries->get_vector();
    VectorType &diagonal_vector = this->diagonal_entries->get_vector();
    this->initialize_dof_vector(inverse_diagonal_vector);
    this->initialize_dof_vector(diagonal_vector);

    this->data->cell_loop(&LaplaceOperator::local_diagonal_cell,
                          this,
                          diagonal_vector,
                          /*unused*/ diagonal_vector);
    this->set_constrained_entries_to_one(diagonal_vector);

    inverse_diagonal_vector = diagonal_vector;

    for (unsigned int i = 0; i < inverse_diagonal_vector.local_size(); ++i)
      if (std::abs(inverse_diagonal_vector.local_element(i)) >
          std::sqrt(std::numeric_limits<Number>::epsilon()))
        inverse_diagonal_vector.local_element(i) =
          1. / inverse_diagonal_vector.local_element(i);
      else
        inverse_diagonal_vector.local_element(i) = 1.;

    inverse_diagonal_vector.update_ghost_values();
    diagonal_vector.update_ghost_values();
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  void
  LaplaceOperator<dim,
                  fe_degree,
                  n_q_points_1d,
                  n_components,
                  VectorType,
                  VectorizedArrayType>::apply_add(VectorType &      dst,
                                                  const VectorType &src) const
  {
    Base<dim, VectorType, VectorizedArrayType>::data->cell_loop(
      &LaplaceOperator::local_apply_cell, this, dst, src);
  }

  namespace Implementation
  {
    template <typename VectorizedArrayType>
    bool
    non_negative(const VectorizedArrayType &n)
    {
      for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
        if (n[v] < 0.)
          return false;

      return true;
    }
  } // namespace Implementation



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  void
  LaplaceOperator<dim,
                  fe_degree,
                  n_q_points_1d,
                  n_components,
                  VectorType,
                  VectorizedArrayType>::
    do_operation_on_cell(
      FEEvaluation<
        dim,
        fe_degree,
        n_q_points_1d,
        n_components,
        typename Base<dim, VectorType, VectorizedArrayType>::value_type> &phi,
      const unsigned int cell) const
  {
    phi.evaluate(false, true, false);
    if (scalar_coefficient.get())
      {
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          {
            Assert(Implementation::non_negative((*scalar_coefficient)(cell, q)),
                   ExcMessage("Coefficient must be non-negative"));
            phi.submit_gradient((*scalar_coefficient)(cell, q) *
                                  phi.get_gradient(q),
                                q);
          }
      }
    else
      {
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          {
            phi.submit_gradient(phi.get_gradient(q), q);
          }
      }
    phi.integrate(false, true);
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  void
  LaplaceOperator<dim,
                  fe_degree,
                  n_q_points_1d,
                  n_components,
                  VectorType,
                  VectorizedArrayType>::
    local_apply_cell(
      const MatrixFree<
        dim,
        typename Base<dim, VectorType, VectorizedArrayType>::value_type,
        VectorizedArrayType> &                     data,
      VectorType &                                 dst,
      const VectorType &                           src,
      const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;
    FEEvaluation<dim, fe_degree, n_q_points_1d, n_components, Number> phi(
      data, this->selected_rows[0]);
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values(src);
        do_operation_on_cell(phi, cell);
        phi.distribute_local_to_global(dst);
      }
  }


  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  void
  LaplaceOperator<dim,
                  fe_degree,
                  n_q_points_1d,
                  n_components,
                  VectorType,
                  VectorizedArrayType>::
    local_diagonal_cell(
      const MatrixFree<
        dim,
        typename Base<dim, VectorType, VectorizedArrayType>::value_type,
        VectorizedArrayType> &data,
      VectorType &            dst,
      const VectorType &,
      const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;

    FEEvaluation<dim, fe_degree, n_q_points_1d, n_components, Number> phi(
      data, this->selected_rows[0]);
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        VectorizedArrayType local_diagonal_vector[phi.static_dofs_per_cell];
        for (unsigned int i = 0; i < phi.dofs_per_component; ++i)
          {
            for (unsigned int j = 0; j < phi.dofs_per_component; ++j)
              phi.begin_dof_values()[j] = VectorizedArrayType();
            phi.begin_dof_values()[i] = 1.;
            do_operation_on_cell(phi, cell);
            local_diagonal_vector[i] = phi.begin_dof_values()[i];
          }
        for (unsigned int i = 0; i < phi.dofs_per_component; ++i)
          for (unsigned int c = 0; c < phi.n_components; ++c)
            phi.begin_dof_values()[i + c * phi.dofs_per_component] =
              local_diagonal_vector[i];
        phi.distribute_local_to_global(dst);
      }
  }


} // end of namespace MatrixFreeOperators


DEAL_II_NAMESPACE_CLOSE

#endif

// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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

#ifndef dealii_non_matching_fe_field_function_h
#define dealii_non_matching_fe_field_function_h

#include <deal.II/base/function.h>
#include <deal.II/base/polynomial.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/tria.h>

DEAL_II_NAMESPACE_OPEN

namespace NonMatching
{
  /**
   * Interface for a scalar Function which has a
   * set_active_cell(..)-function. That is, a function which we in some way
   * need to associate with a given cell in order to evaluate.
   */
  template <int dim>
  class CellWiseFunction : public Function<dim>
  {
  public:
    /**
     * Destructor. Declared to make it virtual.
     */
    virtual ~CellWiseFunction() = default;

    /**
     * Set the cell that the function should be evaluated on.
     */
    virtual void
    set_active_cell(
      const typename Triangulation<dim>::active_cell_iterator &cell) = 0;
  };

  /**
   * This class evaluates a function defined by a solution vector and a
   * DoFHandler transformed to reference space. To be precise, if we let
   * $\hat{x}$ be a point on the reference cell, this class implements the
   * function
   *
   * $\hat{f}(\hat{x}) = \sum_{j=0}^{n-1} f_j \hat{\phi}_j(\hat{x})$,
   *
   * where $f_j$ are the local solution values and $\hat{\phi}_j(\hat(x))$
   * are the local reference space shape functions. The gradient and Hessian
   * of this function are thus derivatives with respect to the reference
   * space coordinates, $\hat{x}_0, \hat{x}_1, \ldots$.
   *
   * Note that this class is similar to FEFieldFunction, but that
   * FEFieldFunction implements the following function on a given cell, $K$,
   *
   * $f(x) = \sum_{j=0}^{n-1} f_j \hat{\phi}_j(F_K^{-1}(x))$,
   *
   * which has the same coefficients but uses real space basis functions.
   * Here, $F_K$ is the mapping from the reference cell to the real cell.
   *
   * Before calling the value/gradient/hessian function, the set_active_cell
   * function must be called to specify which cell the function should be
   * evaluated on.
   */
  template <int dim, class VectorType = Vector<double>>
  class RefSpaceFEFieldFunction : public CellWiseFunction<dim>
  {
  public:
    /**
     * Constructor. Takes the solution vector and the associated DoFHandler.
     *
     * Pointers to the input arguments are stored internally, so they must
     * have a longer lifetime than the created RefSpaceFEFieldFunction
     * object.
     */
    RefSpaceFEFieldFunction(const DoFHandler<dim> &dof_handler,
                            const VectorType &     dof_values);

    /**
     * @copydoc CellWiseFunction::set_active_cell()
     */
    void
    set_active_cell(
      const typename Triangulation<dim>::active_cell_iterator &cell) override;

    /**
     * @copydoc Function::value()
     *
     * @note The set_active_cell function must be called before this function.
     * The incoming point should be on the reference cell, but this is not
     * checked.
     */
    double
    value(const Point<dim> & point,
          const unsigned int component = 0) const override;

    /**
     * @copydoc Function::gradient()
     *
     * @note The set_active_cell function must be called before this function.
     * The incoming point should be on the reference cell, but this is not
     * checked.
     */
    Tensor<1, dim>
    gradient(const Point<dim> & point,
             const unsigned int component = 0) const override;

    /**
     * @copydoc Function::hessian()
     *
     * @note The set_active_cell function must be called before this function.
     * The incoming point should be on the reference cell, but this is not
     * checked.
     */
    SymmetricTensor<2, dim>
    hessian(const Point<dim> & point,
            const unsigned int component = 0) const override;

  private:
    /**
     * Return whether the set_active_cell function has been called.
     */
    bool
    cell_is_set() const;

    /**
     * Pointer to the DoFHandler passed to the constructor.
     */
    const SmartPointer<const DoFHandler<dim>> dof_handler;

    /**
     * Pointer to the vector of solution coefficients passed to the
     * constructor.
     */
    const SmartPointer<const VectorType> global_dof_values;

    /**
     * Pointer to the element associated with the cell in the last call to
     * set_active_cell().
     */
    SmartPointer<const FiniteElement<dim>> element;

    /**
     * DOF-indices of the cell in the last call to set_active_cell()
     */
    std::vector<types::global_dof_index> local_dof_indices;

    /**
     * Local solution values of the cell in the last call to
     * set_active_cell()
     */
    std::vector<typename VectorType::value_type> local_dof_values;

    /**
     * Description of the 1D polynomial basis for tensor product elements
     * used for the fast path of this class using tensor product
     * evaluators.
     */
    std::vector<Polynomials::Polynomial<double>> poly;

    /**
     * Renumbering for the tensor-product evaluator in the fast path.
     */
    std::vector<unsigned int> renumber;

    /**
     * Check whether the shape functions are linear.
     */
    bool polynomials_are_hat_functions;
  };
} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE

#endif // dealii_non_matching_fe_field_function_h

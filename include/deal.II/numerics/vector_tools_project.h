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

#ifndef dealii_vector_tools_project_h
#define dealii_vector_tools_project_h


#include <deal.II/base/config.h>

#include <deal.II/base/template_constraints.h>
#include <deal.II/base/vectorization.h>

#include <functional>
#include <memory>

DEAL_II_NAMESPACE_OPEN

template <typename number>
class AffineConstraints;

template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class DoFHandler;

template <int dim, typename Number>
class Function;
template <int dim, int spacedim>
class Mapping;
template <int dim, typename number, typename VectorizedArrayType>
class MatrixFree;
template <int dim>
class Quadrature;
template <int dim>
class QGauss;
template <typename Number, std::size_t width>
class VectorizedArray;
namespace hp
{
  template <int dim, int spacedim>
  class MappingCollection;
  template <int dim>
  class QCollection;
} // namespace hp


namespace VectorTools
{
  /**
   * @name Interpolation and projection
   * @{
   */

  /**
   * Compute the projection of @p function to the finite element space. In other
   * words, given a function $f(\mathbf x)$, the current function computes a
   * finite element function $f_h(\mathbf x)=\sum_j F_j \varphi_j(\mathbf x)$
   * characterized by the (output) vector of nodal values $F$ that satisfies
   * the equation
   * @f{align*}{
   *   (\varphi_i, f_h)_\Omega = (\varphi_i,f)_\Omega
   * @f}
   * for all test functions $\varphi_i$. This requires solving a linear system
   * involving the @ref GlossMassMatrix "mass matrix" since the equation above is equivalent to
   * the linear system
   * @f{align*}{
   *   \sum_j (\varphi_i, \varphi_j)_\Omega F_j = (\varphi_i,f)_\Omega
   * @f}
   * which can also be written as $MF = \Phi$ with
   * $M_{ij} = (\varphi_i, \varphi_j)_\Omega$ and
   * $\Phi_i = (\varphi_i,f)_\Omega$.
   *
   * By default, no boundary values for $f_h$ are needed nor
   * imposed, but there are optional parameters to this function that allow
   * imposing either zero boundary values or, in a first step, to project
   * the boundary values of $f$ onto the finite element space on the boundary
   * of the mesh in a similar way to above, and then using these values as the
   * imposed boundary values for $f_h$. The ordering of arguments to this
   * function is such that you need not give a second quadrature formula (of
   * type `Quadrature<dim-1>` and used for the computation of the matrix and
   * right hand side for the projection of boundary values) if you
   * don't want to project to the boundary first, but that you must if you want
   * to do so.
   *
   * A MatrixFree implementation is used if the following conditions are met:
   * - @p enforce_zero_boundary is false,
   * - @p project_to_boundary_first is false,
   * - the FiniteElement is supported by the MatrixFree class,
   * - the FiniteElement has less than five components
   * - the degree of the FiniteElement is less than nine.
   * - dim==spacedim
   *
   * In this case, this function performs numerical quadrature using the given
   * quadrature formula for integration of the right hand side $\Phi_i$ and
   * for the mass operator. In the case of hypercube cells, a
   * QGauss(fe_degree+2) object is used for the mass operator. You should
   * therefore make sure that the given quadrature formula is sufficiently
   * accurate for creating the right-hand side.
   *
   * Otherwise, only serial Triangulations are supported and the @ref GlossMassMatrix "mass matrix"
   * is assembled using MatrixTools::create_mass_matrix. The given
   * quadrature rule is then used for both the matrix and the right-hand side.
   * You should therefore make sure that the given quadrature formula is also
   * sufficient for creating the mass matrix. In particular, the degree of the
   * quadrature formula must be sufficiently high to ensure that the mass
   * matrix is invertible. For example, if you are using a FE_Q(k) element,
   * then the integrand of the matrix entries $M_{ij}$ is of polynomial
   * degree $2k$ in each variable, and you need a Gauss quadrature formula
   * with $k+1$ points in each coordinate direction to ensure that $M$
   * is invertible.
   *
   * See the general documentation of this namespace for further information.
   *
   * In 1d, the default value of the boundary quadrature formula is an invalid
   * object since integration on the boundary doesn't happen in 1d.
   *
   * @param[in] mapping The mapping object to use.
   * @param[in] dof The DoFHandler the describes the finite element space to
   * project into and that corresponds to @p vec.
   * @param[in] constraints Constraints to be used when assembling the mass
   * matrix, typically needed when you have hanging nodes.
   * @param[in] quadrature The quadrature formula to be used for assembling the
   * mass matrix.
   * @param[in] function The function to project into the finite element space.
   * @param[out] vec The output vector where the projected function will be
   * stored in. This vector is required to be already initialized and must not
   * have ghost elements.
   * @param[in] enforce_zero_boundary If true, @p vec will have zero boundary
   * conditions.
   * @param[in] q_boundary Quadrature rule to be used if @p project_to_boundary_first
   * is true.
   * @param[in] project_to_boundary_first If true, perform a projection on the
   * boundary before projecting the interior of the function.
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, typename VectorType, int spacedim>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void project(
    const Mapping<dim, spacedim>                              &mapping,
    const DoFHandler<dim, spacedim>                           &dof,
    const AffineConstraints<typename VectorType::value_type>  &constraints,
    const Quadrature<dim>                                     &quadrature,
    const Function<spacedim, typename VectorType::value_type> &function,
    VectorType                                                &vec,
    const bool                 enforce_zero_boundary = false,
    const Quadrature<dim - 1> &q_boundary = (dim > 1 ? QGauss<dim - 1>(2) :
                                                       Quadrature<dim - 1>()),
    const bool                 project_to_boundary_first = false);

  /**
   * Call the project() function above, with
   * <tt>mapping=MappingQ@<dim@>(1)</tt>.
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, typename VectorType, int spacedim>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void project(
    const DoFHandler<dim, spacedim>                           &dof,
    const AffineConstraints<typename VectorType::value_type>  &constraints,
    const Quadrature<dim>                                     &quadrature,
    const Function<spacedim, typename VectorType::value_type> &function,
    VectorType                                                &vec,
    const bool                 enforce_zero_boundary = false,
    const Quadrature<dim - 1> &q_boundary = (dim > 1 ? QGauss<dim - 1>(2) :
                                                       Quadrature<dim - 1>()),
    const bool                 project_to_boundary_first = false);

  /**
   * Same as above, but with hp-capabilities.
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, typename VectorType, int spacedim>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void project(
    const hp::MappingCollection<dim, spacedim>                &mapping,
    const DoFHandler<dim, spacedim>                           &dof,
    const AffineConstraints<typename VectorType::value_type>  &constraints,
    const hp::QCollection<dim>                                &quadrature,
    const Function<spacedim, typename VectorType::value_type> &function,
    VectorType                                                &vec,
    const bool                      enforce_zero_boundary = false,
    const hp::QCollection<dim - 1> &q_boundary = hp::QCollection<dim - 1>(
      dim > 1 ? QGauss<dim - 1>(2) : Quadrature<dim - 1>()),
    const bool project_to_boundary_first = false);

  /**
   * Call the project() function above, with a collection of $Q_1$ mapping
   * objects, i.e., with hp::StaticMappingQ1::mapping_collection.
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, typename VectorType, int spacedim>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void project(
    const DoFHandler<dim, spacedim>                           &dof,
    const AffineConstraints<typename VectorType::value_type>  &constraints,
    const hp::QCollection<dim>                                &quadrature,
    const Function<spacedim, typename VectorType::value_type> &function,
    VectorType                                                &vec,
    const bool                      enforce_zero_boundary = false,
    const hp::QCollection<dim - 1> &q_boundary = hp::QCollection<dim - 1>(
      dim > 1 ? QGauss<dim - 1>(2) : Quadrature<dim - 1>()),
    const bool project_to_boundary_first = false);

  /**
   * The same as above for projection of scalar-valued quadrature data.
   * The user provided function should return a value at the quadrature point
   * based on the cell iterator and quadrature number and of course should be
   * consistent with the provided @p quadrature object, which will be used
   * to assemble the right-hand-side.
   *
   * This function can be used with lambdas:
   * @code
   * VectorTools::project
   * (mapping,
   *  dof_handler,
   *  constraints,
   *  quadrature_formula,
   *  [&] (const typename DoFHandler<dim>::active_cell_iterator & cell,
   *       const unsigned int q) -> double
   *  {
   *    return qp_data.get_data(cell)[q]->density;
   *  },
   *  field);
   * @endcode
   * where <code>qp_data</code> is a CellDataStorage object, which stores
   * quadrature point data.
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, typename VectorType, int spacedim>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void project(
    const Mapping<dim, spacedim>                             &mapping,
    const DoFHandler<dim, spacedim>                          &dof,
    const AffineConstraints<typename VectorType::value_type> &constraints,
    const Quadrature<dim>                                    &quadrature,
    const std::function<typename VectorType::value_type(
      const typename DoFHandler<dim, spacedim>::active_cell_iterator &,
      const unsigned int)>                                   &func,
    VectorType                                               &vec_result);

  /**
   * The same as above for projection of scalar-valued MatrixFree quadrature
   * data.
   * The user provided function @p func should return a VectorizedArray value
   * at the quadrature point based on the cell number and quadrature number and
   * should be consistent with the @p n_q_points_1d.
   *
   * This function can be used with lambdas:
   * @code
   * VectorTools::project
   * (matrix_free_data,
   *  constraints,
   *  3,
   *  [&] (const unsigned int cell,
   *       const unsigned int q) -> VectorizedArray<double>
   *  {
   *    return qp_data(cell,q);
   *  },
   *  field);
   * @endcode
   * where <code>qp_data</code> is a an object of type Table<2,
   * VectorizedArray<double> >, which stores quadrature point data.
   *
   * @p fe_component allow to additionally specify which component of @p data
   * to use in case it was constructed with an <code>std::vector<const
   * DoFHandler<dim>*></code>. It will be used internally in constructor of
   * FEEvaluation object.
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void project(
    std::shared_ptr<
      const MatrixFree<dim,
                       typename VectorType::value_type,
                       VectorizedArray<typename VectorType::value_type>>> data,
    const AffineConstraints<typename VectorType::value_type> &constraints,
    const unsigned int                                        n_q_points_1d,
    const std::function<VectorizedArray<typename VectorType::value_type>(
      const unsigned int,
      const unsigned int)>                                   &func,
    VectorType                                               &vec_result,
    const unsigned int                                        fe_component = 0);

  /**
   * Same as above but for <code>n_q_points_1d =
   * matrix_free.get_dof_handler().get_fe().degree+1</code>.
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void project(
    std::shared_ptr<
      const MatrixFree<dim,
                       typename VectorType::value_type,
                       VectorizedArray<typename VectorType::value_type>>> data,
    const AffineConstraints<typename VectorType::value_type> &constraints,
    const std::function<VectorizedArray<typename VectorType::value_type>(
      const unsigned int,
      const unsigned int)>                                   &func,
    VectorType                                               &vec_result,
    const unsigned int                                        fe_component = 0);

  /** @} */

} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_project_h

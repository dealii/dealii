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

#ifndef dealii_vector_tools_mean_value_h
#define dealii_vector_tools_mean_value_h


#include <deal.II/base/config.h>

#include <deal.II/base/template_constraints.h>

#include <deal.II/lac/read_vector.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN
// forward declarations
template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class DoFHandler;
template <int dim, int spacedim>
class Mapping;
template <int dim>
class Quadrature;
namespace hp
{
  template <int dim, int spacedim>
  class MappingCollection;
  template <int dim>
  class QCollection;
} // namespace hp
#endif

namespace VectorTools
{
  /**
   * Mean value operations
   */
  /** @{ */

  /**
   * Subtract the (algebraic) mean value from a vector.
   *
   * This function is most frequently used as a mean-value filter for Stokes:
   * The pressure in Stokes' equations with only Dirichlet boundaries for the
   * velocities is only determined up to a constant. This function allows to
   * subtract the mean value of the pressure. It is usually called in a
   * preconditioner and generates updates with mean value zero. The mean value
   * is computed as the mean value of the degrees of freedom values as given
   * by the input vector; they are not weighted by the area of cells, i.e. the
   * mean is computed as $\sum_i v_i$, rather than as $\int_\Omega v(x) =
   * \int_\Omega \sum_i v_i \phi_i(x)$. The latter can be obtained from the
   * VectorTools::compute_mean_function, however.
   *
   * Apart from the vector @p v to operate on, this function takes a boolean
   * mask @p p_select that has a true entry for every element of the vector
   * for which the mean value shall be computed and later subtracted. The
   * argument is used to denote which components of the solution vector
   * correspond to the pressure, and avoid touching all other components of
   * the vector, such as the velocity components. (Note, however, that the
   * mask is not a
   * @ref GlossComponentMask
   * operating on the vector components of the finite element the solution
   * vector @p v may be associated with; rather, it is a mask on the entire
   * vector, without reference to what the vector elements mean.)
   *
   * The boolean mask @p p_select has an empty vector as default value, which
   * will be interpreted as selecting all vector elements, hence, subtracting
   * the algebraic mean value on the whole vector. This allows to call this
   * function without a boolean mask if the whole vector should be processed.
   *
   * @note In the context of using this function to filter out the kernel of
   * an operator (such as the null space of the Stokes operator that consists
   * of the constant pressures), this function only makes sense for finite
   * elements for which the null space indeed consists of the vector
   * $(1,1,\ldots,1)^T$. This is the case for example for the usual Lagrange
   * elements where the sum of all shape functions equals the function that is
   * constant one. However, it is not true for some other functions: for
   * example, for the FE_DGP element (another valid choice for the pressure in
   * Stokes discretizations), the first shape function on each cell is
   * constant while further elements are $L_2$ orthogonal to it (on the
   * reference cell); consequently, the sum of all shape functions is not
   * equal to one, and the vector that is associated with the constant mode is
   * not equal to $(1,1,\ldots,1)^T$. For such elements, a different procedure
   * has to be used when subtracting the mean value.
   *
   * @warning This function can only be used for distributed vector classes
   * provided the boolean mask is empty, i.e. selecting the whole vector.
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void subtract_mean_value(VectorType              &v,
                           const std::vector<bool> &p_select = {});

  /**
   * Add the constant @p constant_adjustment to the specified
   * component of the finite element function given by the coefficient
   * vector @p solution defined by the given DoFHandler.
   *
   * This operation is a common operation to compute a solution with
   * mean pressure zero for a Stokes flow problem. Here, one can
   * use VectorTools::compute_mean_value() to compute the value to
   * subtract.
   *
   * For a nodal finite element like FE_Q, this function will simply add
   * the value @p constant_adjustment to each coefficient of the corresponding
   * component. If you have the component in a separate block @p b, you
   * could directly use
   * <code>solution.block(b) += constant_adjustment</code>
   * instead of calling this function (and this would be more efficient).
   * For other finite element spaces like FE_DGP, the logic is more
   * complicated and handled correctly by this function.
   *
   * @note Not all kinds of finite elements are supported and the selected
   * component must not be part of a non-primitive element for this
   * implementation to work.
   *
   * @note In contrast to subtract_mean_value(), this function can
   * adjust a single component of a distributed vector.
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <typename VectorType, int dim, int spacedim = dim>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void add_constant(VectorType                           &solution,
                    const DoFHandler<dim, spacedim>      &dof_handler,
                    const unsigned int                    component,
                    const typename VectorType::value_type constant_adjustment);

  /**
   * Compute the mean value of one component of the solution.
   *
   * This function integrates the chosen component over the whole domain and
   * returns the result, i.e. it computes $\frac{1}{|\Omega|}\int_\Omega
   * [u_h(x)]_c \; dx$ where $c$ is the vector component and $u_h$ is the
   * function representation of the nodal vector given as fourth argument. The
   * integral is evaluated numerically using the quadrature formula given as
   * third argument.
   *
   * This function is used in the "Possibilities for extensions" part of the
   * results section of
   * @ref step_3 "step-3".
   *
   * @note The function is most often used when solving a problem whose
   * solution is only defined up to a constant, for example a pure Neumann
   * problem or the pressure in a Stokes or Navier-Stokes problem. In both
   * cases, subtracting the mean value as computed by the current function,
   * from the nodal vector does not generally yield the desired result of a
   * finite element function with mean value zero. In fact, it only works for
   * Lagrangian elements. For all other elements, you will need to compute the
   * mean value and subtract it right inside the evaluation routine.
   */
  template <int dim, typename Number, int spacedim>
  Number
  compute_mean_value(
    const hp::MappingCollection<dim, spacedim> &mapping_collection,
    const DoFHandler<dim, spacedim>            &dof,
    const hp::QCollection<dim>                 &q_collection,
    const ReadVector<Number>                   &v,
    const unsigned int                          component);

  /**
   * Calls the other compute_mean_value() function, see above, for the non-hp
   * case. That means, it requires a single FiniteElement, a single Quadrature,
   * and a single Mapping object.
   */
  template <int dim, typename Number, int spacedim>
  Number
  compute_mean_value(const Mapping<dim, spacedim>    &mapping,
                     const DoFHandler<dim, spacedim> &dof,
                     const Quadrature<dim>           &quadrature,
                     const ReadVector<Number>        &v,
                     const unsigned int               component);

  /**
   * Call the other compute_mean_value() function, see above, with
   * <tt>mapping=MappingQ@<dim@>(1)</tt>.
   */
  template <int dim, typename Number, int spacedim>
  Number
  compute_mean_value(const DoFHandler<dim, spacedim> &dof,
                     const Quadrature<dim>           &quadrature,
                     const ReadVector<Number>        &v,
                     const unsigned int               component);
  /** @} */
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_mean_value_h

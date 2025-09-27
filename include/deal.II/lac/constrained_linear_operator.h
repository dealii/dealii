// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_constrained_linear_operator_h
#define dealii_constrained_linear_operator_h

#include <deal.II/base/config.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/packaged_operation.h>


DEAL_II_NAMESPACE_OPEN


/**
 * @name Indirectly applying constraints to LinearOperator
 */
/** @{ */


/**
 * This function takes an AffineConstraints object @p constraints and
 * an operator exemplar @p exemplar (this exemplar is usually a linear
 * operator that describes the system matrix - it is only used to create
 * domain and range vectors of appropriate sizes, its action <tt>vmult</tt>
 * is never used). A LinearOperator object associated with the "homogeneous
 * action" of the underlying AffineConstraints object is returned:
 *
 * Applying the LinearOperator object on a vector <code>u</code> results in a
 * vector <code>v</code> that stores the result of calling
 * AffineConstraints::distribute() on <code>u</code> - with one important
 * difference: inhomogeneities are not applied, but always treated as 0
 * instead.
 *
 * The LinearOperator object created by this function is primarily used
 * internally in constrained_linear_operator() to build up a modified system
 * of linear equations. How to solve a linear system of equations with this
 * approach is explained in detail in the
 * @ref constraints
 * topic.
 *
 *
 * @note Currently, this function may not work correctly for distributed data
 * structures.
 *
 * @relatesalso LinearOperator
 * @ingroup constraints
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
distribute_constraints_linear_operator(
  const AffineConstraints<typename Range::value_type> &constraints,
  const LinearOperator<Range, Domain, Payload>        &exemplar)
{
  LinearOperator<Range, Domain, Payload> return_op = exemplar;

  return_op.vmult_add = [&constraints](Range &v, const Domain &u) {
    Assert(!dealii::PointerComparison::equal(&v, &u),
           dealii::ExcMessage("The domain and range vectors must be different "
                              "storage locations"));

    // First, add vector u to v unconditionally and clean up constrained
    // degrees of freedom later.
    v += u;

    const auto &locally_owned_elements = v.locally_owned_elements();
    for (const auto &line : constraints.get_lines())
      {
        const auto i = line.index;
        if (locally_owned_elements.is_element(i))
          {
            v(i) -= u(i);
            const auto &entries = line.entries;
            for (types::global_dof_index j = 0; j < entries.size(); ++j)
              {
                const auto pos = entries[j].first;
                v(i) += u(pos) * entries[j].second;
              }
          }
      }

    v.compress(VectorOperation::add);
  };

  return_op.Tvmult_add = [&constraints](Domain &v, const Range &u) {
    Assert(!dealii::PointerComparison::equal(&v, &u),
           dealii::ExcMessage("The domain and range vectors must be different "
                              "storage locations"));

    // First, add vector u to v unconditionally and clean up constrained
    // degrees of freedom later.
    v += u;

    const auto &locally_owned_elements = v.locally_owned_elements();
    for (const auto &line : constraints.get_lines())
      {
        const auto i = line.index;

        if (locally_owned_elements.is_element(i))
          {
            v(i) -= u(i);
          }

        const auto &entries = line.entries;
        for (types::global_dof_index j = 0; j < entries.size(); ++j)
          {
            const auto pos = entries[j].first;
            if (locally_owned_elements.is_element(pos))
              v(pos) += u(i) * entries[j].second;
          }
      }

    v.compress(VectorOperation::add);
  };

  return_op.vmult = [vmult_add = return_op.vmult_add](Range        &v,
                                                      const Domain &u) {
    v = 0.;
    vmult_add(v, u);
  };

  return_op.Tvmult = [Tvmult_add = return_op.Tvmult_add](Domain      &v,
                                                         const Range &u) {
    v = 0.;
    Tvmult_add(v, u);
  };

  return return_op;
}


/**
 * Given a AffineConstraints @p constraints and an operator exemplar @p
 * exemplar, return a LinearOperator that is the projection to the subspace of
 * constrained degrees of freedom, i.e. all entries of the result vector that
 * correspond to unconstrained degrees of freedom are set to zero.
 *
 *
 * @relatesalso LinearOperator
 * @ingroup constraints
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
project_to_constrained_linear_operator(
  const AffineConstraints<typename Range::value_type> &constraints,
  const LinearOperator<Range, Domain, Payload>        &exemplar)
{
  LinearOperator<Range, Domain, Payload> return_op = exemplar;

  return_op.vmult_add = [&constraints](Range &v, const Domain &u) {
    const auto &locally_owned_elements = v.locally_owned_elements();
    for (const auto &line : constraints.get_lines())
      {
        const auto i = line.index;
        if (locally_owned_elements.is_element(i))
          {
            v(i) += u(i);
          }
      }

    v.compress(VectorOperation::add);
  };

  return_op.Tvmult_add = [&constraints](Domain &v, const Range &u) {
    const auto &locally_owned_elements = v.locally_owned_elements();
    for (const auto &line : constraints.get_lines())
      {
        const auto i = line.index;
        if (locally_owned_elements.is_element(i))
          {
            v(i) += u(i);
          }
      }

    v.compress(VectorOperation::add);
  };

  return_op.vmult = [vmult_add = return_op.vmult_add](Range        &v,
                                                      const Domain &u) {
    v = 0.;
    vmult_add(v, u);
  };

  return_op.Tvmult = [Tvmult_add = return_op.Tvmult_add](Domain      &v,
                                                         const Range &u) {
    v = 0.;
    Tvmult_add(v, u);
  };

  return return_op;
}


/**
 * Given a AffineConstraints object @p constraints and a LinearOperator
 * @p linop, this function creates a LinearOperator object consisting of the
 * composition of three operations and a regularization:
 * @code
 *   Ct * linop * C + Id_c;
 * @endcode
 * with
 * @code
 *   C = distribute_constraints_linear_operator(constraints, linop);
 *   Ct = transpose_operator(C);
 *   Id_c = project_to_constrained_linear_operator(constraints, linop);
 * @endcode
 * and <code>Id_c</code> is the projection to the subspace consisting of all
 * vector entries associated with constrained degrees of freedom.
 *
 * This LinearOperator object is used together with
 * constrained_right_hand_side() to build up the following modified system of
 * linear equations:
 * @f[
 *   (C^T A C + Id_c) x = C^T (b - A\,k)
 * @f]
 * with a given (unconstrained) system matrix $A$, right hand side $b$, and
 * linear constraints $C$ with inhomogeneities $k$.
 *
 * A detailed explanation of this approach is given in the
 * @ref constraints
 * topic.
 *
 *
 * @note Currently, this function may not work correctly for distributed data
 * structures.
 *
 * @relatesalso LinearOperator
 * @ingroup constraints
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
constrained_linear_operator(
  const AffineConstraints<typename Range::value_type> &constraints,
  const LinearOperator<Range, Domain, Payload>        &linop)
{
  const auto C    = distribute_constraints_linear_operator(constraints, linop);
  const auto Ct   = transpose_operator(C);
  const auto Id_c = project_to_constrained_linear_operator(constraints, linop);
  return Ct * linop * C + Id_c;
}


/**
 * Given a AffineConstraints object @p constraints, a LinearOperator @p
 * linop and a right-hand side @p right_hand_side, this function creates a
 * PackagedOperation that stores the following computation:
 * @code
 *   Ct * (right_hand_side - linop * k)
 * @endcode
 * with
 * @code
 *   C = distribute_constraints_linear_operator(constraints, linop);
 *   Ct = transpose_operator(C);
 * @endcode
 *
 * This LinearOperator object is used together with
 * constrained_right_hand_side() to build up the following modified system of
 * linear equations:
 * @f[
 *   (C^T A C + Id_c) x = C^T (b - A\,k)
 * @f]
 * with a given (unconstrained) system matrix $A$, right hand side $b$, and
 * linear constraints $C$ with inhomogeneities $k$.
 *
 * A detailed explanation of this approach is given in the
 * @ref constraints
 * topic.
 *
 *
 * @note Currently, this function may not work correctly for distributed data
 * structures.
 *
 * @relatesalso LinearOperator
 * @ingroup constraints
 */
template <typename Range, typename Domain, typename Payload>
PackagedOperation<Range>
constrained_right_hand_side(
  const AffineConstraints<typename Range::value_type> &constraints,
  const LinearOperator<Range, Domain, Payload>        &linop,
  const Range                                         &right_hand_side)
{
  PackagedOperation<Range> return_comp;

  return_comp.reinit_vector = linop.reinit_range_vector;

  return_comp.apply_add = [&constraints, &linop, &right_hand_side](Range &v) {
    const auto C  = distribute_constraints_linear_operator(constraints, linop);
    const auto Ct = transpose_operator(C);

    GrowingVectorMemory<Domain>            vector_memory;
    typename VectorMemory<Domain>::Pointer k(vector_memory);
    linop.reinit_domain_vector(*k, /*bool fast=*/false);
    constraints.distribute(*k);

    v += Ct * (right_hand_side - linop * *k);
  };

  return_comp.apply = [apply_add = return_comp.apply_add](Range &v) {
    v = 0.;
    apply_add(v);
  };

  return return_comp;
}

/** @} */

DEAL_II_NAMESPACE_CLOSE

#endif

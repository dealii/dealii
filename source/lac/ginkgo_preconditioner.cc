// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2023 by the deal.II authors
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


#include <deal.II/lac/ginkgo_preconditioner.h>


#ifdef DEAL_II_WITH_GINKGO

#  include <ginkgo/core/base/types.hpp>
#  include <ginkgo/core/preconditioner/jacobi.hpp>
#  include <ginkgo/core/solver/ir.hpp>

DEAL_II_NAMESPACE_OPEN


namespace GinkgoWrappers
{
  template <typename ValueType, typename IndexType>
  PreconditionBase<ValueType, IndexType>::PreconditionBase(
    const AbstractMatrix<ValueType, IndexType> &A,
    std::shared_ptr<gko::LinOpFactory>          factory)
    : preconditioner_(factory->generate(A.get_shared_gko_object()))
  {}

  template <typename ValueType, typename IndexType>
  std::shared_ptr<const gko::LinOp>
  PreconditionBase<ValueType, IndexType>::get_shared_gko_object() const
  {
    return preconditioner_;
  }

  template <typename ValueType, typename IndexType>
  void
  PreconditionBase<ValueType, IndexType>::vmult(Vector<ValueType> &      u,
                                                const Vector<ValueType> &v)
  {
    if (preconditioner_)
      {
        preconditioner_->apply(v.get_gko_object(), u.get_gko_object());
      }
    else
      {
        u = v;
      }
  }


  template <typename ValueType, typename IndexType>
  PreconditionIdentity<ValueType, IndexType>::PreconditionIdentity(
    const PreconditionIdentity::AdditionalData &)
    : PreconditionBase<ValueType, IndexType>()
  {}


  template <typename ValueType, typename IndexType>
  PreconditionJacobi<ValueType, IndexType>::PreconditionJacobi(
    const Csr<ValueType, IndexType> &         A,
    const PreconditionJacobi::AdditionalData &data)
    : PreconditionBase<ValueType, IndexType>(
        A,
        gko::solver::Ir<ValueType>::build()
          .with_solver(
            gko::preconditioner::Jacobi<ValueType, IndexType>::build()
              .with_max_block_size(
                static_cast<gko::uint32>(data.max_block_size))
              .on(A.get_gko_object()->get_executor()))
          .with_relaxation_factor(data.relaxation)
          .with_criteria(gko::stop::Iteration::build()
                           .with_max_iters(data.n_iterations)
                           .on(A.get_gko_object()->get_executor()))
          .with_default_initial_guess(gko::solver::initial_guess_mode::zero)
          .on(A.get_gko_object()->get_executor()))
  {}

  using gko::int32;
  using gko::int64;

#  define DECLARE_GINKGO_PRECONDITIONER_BASE(_value_type, _index_type) \
    class PreconditionBase<_value_type, _index_type>
  GKO_INSTANTIATE_FOR_EACH_VALUE_AND_INDEX_TYPE(
    DECLARE_GINKGO_PRECONDITIONER_BASE);


#  define DECLARE_GINKGO_PRECONDITIONER_IDENTITY(_value_type, _index_type) \
    class PreconditionIdentity<_value_type, _index_type>
  GKO_INSTANTIATE_FOR_EACH_VALUE_AND_INDEX_TYPE(
    DECLARE_GINKGO_PRECONDITIONER_IDENTITY);


#  define DECLARE_GINKGO_PRECONDITIONER_JACOBI(_value_type, _index_type) \
    class PreconditionJacobi<_value_type, _index_type>
  GKO_INSTANTIATE_FOR_EACH_VALUE_AND_INDEX_TYPE(
    DECLARE_GINKGO_PRECONDITIONER_JACOBI);

} // namespace GinkgoWrappers

DEAL_II_NAMESPACE_CLOSE

#endif

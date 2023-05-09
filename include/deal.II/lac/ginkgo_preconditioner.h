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

#ifndef dealii_ginkgo_preconditioner_h
#define dealii_ginkgo_preconditioner_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_GINKGO

#  include <deal.II/base/subscriptor.h>

#  include <deal.II/lac/ginkgo_sparse_matrix.h>
#  include <deal.II/lac/ginkgo_vector.h>


DEAL_II_NAMESPACE_OPEN

namespace GinkgoWrappers
{
  template <typename ValueType, typename IndexType = int32_t>
  class PreconditionBase : public Subscriptor
  {
  public:
    PreconditionBase(const AbstractMatrix<ValueType, IndexType> &A,
                     std::shared_ptr<gko::LinOpFactory>          factory);

    std::shared_ptr<const gko::LinOp>
    get_shared_gko_object() const;

    void
    vmult(Vector<ValueType> &u, const Vector<ValueType> &v);

  protected:
    PreconditionBase() = default;

  private:
    std::shared_ptr<gko::LinOp> preconditioner_;
  };


  namespace detail
  {
    struct AdditionalRelaxationData
    {
      double       relaxation   = 1.0;
      unsigned int n_iterations = 1;
    };
  } // namespace detail


  template <typename ValueType, typename IndexType = int32_t>
  class PreconditionIdentity : public PreconditionBase<ValueType, IndexType>
  {
  public:
    struct AdditionalData
    {};

    explicit PreconditionIdentity(const AdditionalData & = {});
  };


  template <typename ValueType, typename IndexType = int32_t>
  class PreconditionJacobi : public PreconditionBase<ValueType, IndexType>
  {
  public:
    struct AdditionalData : detail::AdditionalRelaxationData
    {
      AdditionalData() = default;
      AdditionalData(double       relaxation,
                     unsigned int n_iterations,
                     unsigned int max_block_size)
        : detail::AdditionalRelaxationData{relaxation, n_iterations}
        , max_block_size(max_block_size)
      {}

      unsigned int max_block_size = 32;
    };

    explicit PreconditionJacobi(const Csr<ValueType, IndexType> &A,
                                const AdditionalData &           data = {});
  };

} // namespace GinkgoWrappers

DEAL_II_NAMESPACE_CLOSE

#endif

#endif // dealii_ginkgo_preconditioner_h

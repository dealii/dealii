// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Check the intermediate vector maps of rectangular Trilinos payload sums and
// products, including their transposed and inverse actions.

#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <Epetra_MultiVector.h>

#include <cmath>
#include <initializer_list>
#include <string>

#include "../tests.h"

using Payload =
  TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload;

void
initialize(TrilinosWrappers::SparseMatrix &matrix,
           const double                    a,
           const double                    b)
{
  matrix.set(0, 0, a);
  matrix.set(1, 1, b);
  matrix.compress(VectorOperation::insert);
}

void
set_inverse_actions(Payload &payload, const double a, const double b)
{
  const auto domain = payload.OperatorDomainMap();
  const auto range  = payload.OperatorRangeMap();
  // These diagonal embeddings have explicit Moore-Penrose inverses. Supply
  // their actions directly to check the payload's vector-map plumbing without
  // requiring a particular solver or preconditioner.
  payload.inv_vmult = [domain, range, a, b](Epetra_MultiVector       &dst,
                                            const Epetra_MultiVector &src) {
    AssertThrow(dst.Map().SameAs(domain) && src.Map().SameAs(range),
                ExcMessage("Incorrect inverse-action maps"));
    dst.PutScalar(0.0);
    dst[0][0] = src[0][0] / a;
    dst[0][1] = src[0][1] / b;
  };
  payload.inv_Tvmult = [domain, range, a, b](Epetra_MultiVector       &dst,
                                             const Epetra_MultiVector &src) {
    AssertThrow(dst.Map().SameAs(range) && src.Map().SameAs(domain),
                ExcMessage("Incorrect transposed inverse-action maps"));
    dst.PutScalar(0.0);
    dst[0][0] = src[0][0] / a;
    dst[0][1] = src[0][1] / b;
  };
}

void
check(Payload                            &payload,
      const bool                          inverse,
      const std::initializer_list<double> input,
      const std::initializer_list<double> expected)
{
  Epetra_MultiVector src(inverse ? payload.OperatorRangeMap() :
                                   payload.OperatorDomainMap(),
                         1);
  Epetra_MultiVector dst(inverse ? payload.OperatorDomainMap() :
                                   payload.OperatorRangeMap(),
                         1);
  AssertThrow(src.MyLength() == static_cast<int>(input.size()) &&
                dst.MyLength() == static_cast<int>(expected.size()),
              ExcInternalError());
  unsigned int i = 0;
  for (const double value : input)
    src[0][i++] = value;
  AssertThrow((inverse ? payload.ApplyInverse(src, dst) :
                         payload.Apply(src, dst)) == 0,
              ExcInternalError());
  i = 0;
  for (const double value : expected)
    AssertThrow(std::abs(dst[0][i++] - value) < 1e-12, ExcInternalError());
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  initlog();

  const auto                     two   = complete_index_set(2);
  const auto                     three = complete_index_set(3);
  const auto                     four  = complete_index_set(4);
  TrilinosWrappers::SparseMatrix a(three, two, MPI_COMM_SELF, 2);
  TrilinosWrappers::SparseMatrix b(three, two, MPI_COMM_SELF, 2);
  TrilinosWrappers::SparseMatrix c(two, four, MPI_COMM_SELF, 2);
  initialize(a, 2.0, 3.0);
  initialize(b, 4.0, 6.0);
  initialize(c, 5.0, 7.0);
  Payload pa(a, a), pb(b, b), pc(c, c);
  set_inverse_actions(pa, 2.0, 3.0);
  set_inverse_actions(pb, 4.0, 6.0);
  set_inverse_actions(pc, 5.0, 7.0);
  auto sum     = pa + pb;
  auto product = pa * pc;

  const std::string selected = argc > 1 ? argv[1] : "all";
  bool              executed = false;
  auto              run = [&](const std::string &name, const auto &operation) {
    if (selected == "all" || selected == name)
      {
        operation();
        deallog << name << " OK" << std::endl;
        executed = true;
      }
  };
  run("sum", [&]() { check(sum, false, {1, 2}, {6, 18, 0}); });
  // A sum payload adds the supplied inverse callbacks; this is not the
  // Moore-Penrose inverse of the sum matrix.
  run("sum_inverse", [&]() { check(sum, true, {4, 6, 8}, {3, 3}); });
  sum.SetUseTranspose(true);
  run("sum_transpose", [&]() { check(sum, false, {1, 2, 3}, {6, 18}); });
  run("sum_inverse_transpose", [&]() { check(sum, true, {4, 6}, {3, 3, 0}); });
  sum.SetUseTranspose(false);
  run("product", [&]() { check(product, false, {1, 2, 3, 4}, {10, 42, 0}); });
  run("product_inverse", [&]() {
    check(product, true, {10, 42, 5}, {1, 2, 0, 0});
  });
  product.SetUseTranspose(true);
  run("product_transpose", [&]() {
    check(product, false, {1, 2, 3}, {10, 42, 0, 0});
  });
  run("product_inverse_transpose", [&]() {
    check(product, true, {10, 42, 8, 9}, {1, 2, 0});
  });
  AssertThrow(executed, ExcMessage("Unknown test case"));
}

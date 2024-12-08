// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2023 by the deal.II Authors
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

// Tests the GinkgoWrappers interface

#include "../tests.h"

#include "test_macros.h"

// all include files you need here
#include <deal.II/base/utilities.h>

#include <deal.II/lac/ginkgo_interface.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/sparse_matrix.templates.h>

#include <Kokkos_Random.hpp>


using namespace dealii;


template <typename View, typename ValueType>
bool
check_equality(View view, gko::matrix::Dense<ValueType> *vector)
{
  using exec_space = typename View::execution_space;
  auto data        = vector->get_values();
  bool success     = true;
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<exec_space>(0, view.size()),
    KOKKOS_LAMBDA(const int i, bool &lsuccess) {
      lsuccess = lsuccess && (view(i) == data[i]);
    },
    Kokkos::LAnd<bool, exec_space>{success});
  return success;
}


template <typename View>
void
fill_view(View view)
{
  using exec_space = typename View::execution_space;
  Kokkos::Random_XorShift1024_Pool<exec_space> pool(/*seed=*/12345);
  Kokkos::parallel_for(
    Kokkos::RangePolicy<exec_space>(0, view.size()),
    KOKKOS_LAMBDA(const int i) {
      auto generator = pool.get_state();
      view(i)        = generator.drand(0., 1.);
      pool.free_state(generator);
    });
}



TEST(can_create_vector_host)
{
  Kokkos::View<double[5], MemorySpace::Host::kokkos_space> view{"name"};
  ArrayView<double, MemorySpace::Host> array_view(view.data(), view.size());
  fill_view(view);

  auto vector = GinkgoWrappers::create_vector(array_view);

  TEST_ASSERT(vector->get_values() == view.data());
}


TEST(can_create_vector_host_const)
{
  using array_view_type = ArrayView<double, MemorySpace::Host>;
  Kokkos::View<double[5], MemorySpace::Host::kokkos_space> view{"name"};
  array_view_type array_view(view.data(), view.size());
  fill_view(view);

  auto vector = GinkgoWrappers::create_vector(
    const_cast<const array_view_type &>(array_view));

  using gko_type = typename decltype(vector)::element_type;
  TEST_ASSERT(std::is_const_v<gko_type>);
}


TEST(can_create_vector_default)
{
  Kokkos::View<double[5], MemorySpace::Default::kokkos_space> view{"name"};
  ArrayView<double, MemorySpace::Default> array_view(view.data(), view.size());
  fill_view(view);

  auto vector = GinkgoWrappers::create_vector(array_view);

  TEST_ASSERT(vector->get_values() == array_view.begin());
}


TEST(can_create_vector_with_executor)
{
  if (!std::is_same_v<MemorySpace::Host::kokkos_space,
                      MemorySpace::Default::kokkos_space>)
    {
      Kokkos::View<double[5], MemorySpace::Default::kokkos_space> view{"name"};
      ArrayView<double, MemorySpace::Default> array_view(view.data(),
                                                         view.size());
      fill_view(view);

      auto ref    = gko::ReferenceExecutor::create();
      auto vector = GinkgoWrappers::create_vector(ref, array_view);

      auto mirror = Kokkos::create_mirror(view);
      Kokkos::deep_copy(mirror, view);
      TEST_ASSERT(vector->get_values() != array_view.begin());
      TEST_ASSERT(check_equality(mirror, vector.get()));
    }
}


TEST(can_create_csr)
{
  std::vector<std::vector<unsigned int>> col_idxs{{0, 3}, {1}, {0, 2}};
  dealii::SparsityPattern                pattern;
  pattern.copy_from(3, 4, col_idxs.begin(), col_idxs.end());
  dealii::SparseMatrix<double> dealii_obj(pattern);
  dealii_obj.set(0, 0, 1);
  dealii_obj.set(0, 3, 2);
  dealii_obj.set(1, 1, 3);
  dealii_obj.set(2, 0, 4);
  dealii_obj.set(2, 2, 5);
  auto exec = gko::ext::kokkos::get_default_executor();

  auto obj = GinkgoWrappers::create_csr_matrix(exec, dealii_obj);
  TEST_ASSERT(obj->get_num_stored_elements() == 5);
  TEST_ASSERT(obj->get_const_col_idxs()[0] == 0);
  TEST_ASSERT(obj->get_const_col_idxs()[1] == 3);
  TEST_ASSERT(obj->get_const_col_idxs()[2] == 1);
  TEST_ASSERT(obj->get_const_col_idxs()[3] == 0);
  TEST_ASSERT(obj->get_const_col_idxs()[4] == 2);
  TEST_ASSERT(obj->get_const_values()[0] == 1);
  TEST_ASSERT(obj->get_const_values()[1] == 2);
  TEST_ASSERT(obj->get_const_values()[2] == 3);
  TEST_ASSERT(obj->get_const_values()[3] == 4);
  TEST_ASSERT(obj->get_const_values()[4] == 5);
  TEST_ASSERT(obj->get_const_row_ptrs()[0] == 0);
  TEST_ASSERT(obj->get_const_row_ptrs()[1] == 2);
  TEST_ASSERT(obj->get_const_row_ptrs()[2] == 3);
  TEST_ASSERT(obj->get_const_row_ptrs()[3] == 5);
}


TEST(can_create_csr_with_strategy)
{
  std::vector<std::vector<unsigned int>> col_idxs{{0, 3}, {1}, {0, 2}};
  dealii::SparsityPattern                pattern;
  pattern.copy_from(3, 4, col_idxs.begin(), col_idxs.end());
  dealii::SparseMatrix<double> dealii_obj(pattern);
  dealii_obj.set(0, 0, 1);
  dealii_obj.set(0, 3, 2);
  dealii_obj.set(1, 1, 3);
  dealii_obj.set(2, 0, 4);
  dealii_obj.set(2, 2, 5);
  auto exec = gko::ext::kokkos::get_default_executor();

  auto obj =
    GinkgoWrappers::create_csr_matrix(exec,
                                      dealii_obj,
                                      GinkgoWrappers::csr_strategy::merge_path);

  using csr_type =
    std::remove_cv_t<std::remove_reference_t<decltype(obj)>>::element_type;
  TEST_ASSERT(std::dynamic_pointer_cast<typename csr_type::merge_path>(
    obj->get_strategy()));
}


TEST(can_create_csr_with_index_type)
{
  std::vector<std::vector<unsigned int>> col_idxs{{0, 3}, {1}, {0, 2}};
  dealii::SparsityPattern                pattern;
  pattern.copy_from(3, 4, col_idxs.begin(), col_idxs.end());
  dealii::SparseMatrix<double> dealii_obj(pattern);
  dealii_obj.set(0, 0, 1);
  dealii_obj.set(0, 3, 2);
  dealii_obj.set(1, 1, 3);
  dealii_obj.set(2, 0, 4);
  dealii_obj.set(2, 2, 5);
  auto exec = gko::ext::kokkos::get_default_executor();

  auto obj = GinkgoWrappers::create_csr_matrix<gko::int64>(
    exec, dealii_obj, GinkgoWrappers::csr_strategy::merge_path);

  using csr_type =
    std::remove_cv_t<std::remove_reference_t<decltype(obj)>>::element_type;
  auto same_index_type =
    std::is_same_v<typename csr_type::index_type, gko::int64>;
  TEST_ASSERT(same_index_type);
}


TEST(can_create_csr_with_value_index_type)
{
  std::vector<std::vector<unsigned int>> col_idxs{{0, 3}, {1}, {0, 2}};
  dealii::SparsityPattern                pattern;
  pattern.copy_from(3, 4, col_idxs.begin(), col_idxs.end());
  dealii::SparseMatrix<double> dealii_obj(pattern);
  dealii_obj.set(0, 0, 1);
  dealii_obj.set(0, 3, 2);
  dealii_obj.set(1, 1, 3);
  dealii_obj.set(2, 0, 4);
  dealii_obj.set(2, 2, 5);
  auto exec = gko::ext::kokkos::get_default_executor();

  auto obj = GinkgoWrappers::create_csr_matrix<float, gko::int64>(
    exec, dealii_obj, GinkgoWrappers::csr_strategy::merge_path);

  using csr_type =
    std::remove_cv_t<std::remove_reference_t<decltype(obj)>>::element_type;
  auto same_value_type = std::is_same_v<typename csr_type::value_type, float>;
  auto same_index_type =
    std::is_same_v<typename csr_type::index_type, gko::int64>;
  TEST_ASSERT(same_value_type);
  TEST_ASSERT(same_index_type);
}


template <typename MemorySpace>
struct linop_data
{
  linop_data()
  {
    std::vector<std::vector<unsigned int>> col_idxs{{0, 3}, {1}, {0, 2}};
    pattern.copy_from(3, 4, col_idxs.begin(), col_idxs.end());

    dealii_obj = dealii::SparseMatrix<double>(pattern);
    dealii_obj.set(0, 0, 1);
    dealii_obj.set(0, 3, 2);
    dealii_obj.set(1, 1, 3);
    dealii_obj.set(2, 0, 4);
    dealii_obj.set(2, 2, 5);

    exec = gko::ext::kokkos::get_executor(
      typename MemorySpace::kokkos_space::execution_space{});
    gko_obj = gko::share(GinkgoWrappers::create_csr_matrix(exec, dealii_obj));
  }

  dealii::SparsityPattern      pattern;
  dealii::SparseMatrix<double> dealii_obj;

  std::shared_ptr<gko::Executor> exec;
  std::shared_ptr<gko::LinOp>    gko_obj;
};

template <typename MemorySpace>
struct linop_vectors
{
  using Vector = LinearAlgebra::distributed::Vector<double, MemorySpace>;
  using VectorHost =
    LinearAlgebra::distributed::Vector<double, ::dealii::MemorySpace::Host>;

  linop_vectors(size_t n, size_t m)
  {
    auto init_vectors = [](auto &x, auto &x_ref, auto size) {
      auto get_view = [](auto &vector) {
        return Kokkos::View<double *,
                            typename MemorySpace::kokkos_space,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>>{
          vector.begin(), vector.size()};
      };

      x.reinit(size);
      auto view_x = get_view(x);
      fill_view(view_x);

      x_ref.reinit(size);
      auto view_x_ref = Kokkos::create_mirror(view_x);
      Kokkos::deep_copy(view_x, view_x_ref);
    };

    init_vectors(u, u_ref, n);
    init_vectors(v, v_ref, m);
  }

  Vector     u;
  Vector     v;
  VectorHost u_ref;
  VectorHost v_ref;
};



TEST(can_create_linop_host)
{
  linop_data<MemorySpace::Host> data{};

  auto linop = GinkgoWrappers::linear_operator(data.gko_obj);

  {
    linop_vectors<MemorySpace::Host> vectors(data.dealii_obj.n(),
                                             data.dealii_obj.m());

    linop.vmult(vectors.v, vectors.u);
    data.dealii_obj.vmult(vectors.v_ref, vectors.u_ref);

    vectors.v_ref -= vectors.v;
    TEST_ASSERT(vectors.v_ref.linfty_norm() < 1e-12);
  }
  {
    linop_vectors<MemorySpace::Host> vectors(data.dealii_obj.n(),
                                             data.dealii_obj.m());

    linop.vmult_add(vectors.v, vectors.u);
    data.dealii_obj.vmult_add(vectors.v_ref, vectors.u_ref);

    vectors.v_ref -= vectors.v;
    TEST_ASSERT(vectors.v_ref.linfty_norm() < 1e-12);
  }
  {
    linop_vectors<MemorySpace::Host> vectors(data.dealii_obj.n(),
                                             data.dealii_obj.m());

    linop.Tvmult(vectors.u, vectors.v);
    data.dealii_obj.Tvmult(vectors.u_ref, vectors.v_ref);

    vectors.u_ref -= vectors.u;
    TEST_ASSERT(vectors.u_ref.linfty_norm() < 1e-12);
  }
  {
    linop_vectors<MemorySpace::Host> vectors(data.dealii_obj.n(),
                                             data.dealii_obj.m());

    linop.Tvmult_add(vectors.u, vectors.v);
    data.dealii_obj.Tvmult_add(vectors.u_ref, vectors.v_ref);

    vectors.u_ref -= vectors.u;
    TEST_ASSERT(vectors.u_ref.linfty_norm() < 1e-12);
  }
}

TEST(can_create_linop_default)
{
  linop_data<MemorySpace::Default> data{};
  linop_data<MemorySpace::Host>    data_host{};

  auto linop = GinkgoWrappers::
    linear_operator<double, double, MemorySpace::Default, MemorySpace::Default>(
      data.gko_obj);
  auto host_linop = GinkgoWrappers::linear_operator(data_host.gko_obj);

  auto compare = [&](const auto &u, const auto &u_ref) {
    LinearAlgebra::distributed::Vector<double, MemorySpace::Default> mirror(
      u_ref.size());

    auto default_exec = gko::ext::kokkos::get_default_executor();
    auto host_exec    = gko::ext::kokkos::get_default_host_executor();

    default_exec->copy_from(host_exec,
                            mirror.size(),
                            u_ref.begin(),
                            mirror.begin());

    mirror -= u;
    TEST_ASSERT(mirror.linfty_norm() < 1e-12);
  };
  {
    linop_vectors<MemorySpace::Default> vectors(data.dealii_obj.n(),
                                                data.dealii_obj.m());

    linop.vmult(vectors.v, vectors.u);
    host_linop.vmult(vectors.v_ref, vectors.u_ref);

    compare(vectors.v, vectors.v_ref);
  }
  {
    linop_vectors<MemorySpace::Default> vectors(data.dealii_obj.n(),
                                                data.dealii_obj.m());

    linop.vmult_add(vectors.v, vectors.u);
    host_linop.vmult_add(vectors.v_ref, vectors.u_ref);

    compare(vectors.v, vectors.v_ref);
  }
  {
    linop_vectors<MemorySpace::Default> vectors(data.dealii_obj.n(),
                                                data.dealii_obj.m());

    linop.Tvmult(vectors.u, vectors.v);
    host_linop.Tvmult(vectors.u_ref, vectors.v_ref);

    compare(vectors.v, vectors.v_ref);
  }
  {
    linop_vectors<MemorySpace::Default> vectors(data.dealii_obj.n(),
                                                data.dealii_obj.m());

    linop.Tvmult_add(vectors.u, vectors.v);
    host_linop.Tvmult_add(vectors.u_ref, vectors.v_ref);

    compare(vectors.v, vectors.v_ref);
  }
}


TEST(can_create_inverse_linop)
{
  auto exec            = gko::ext::kokkos::get_default_executor();
  auto gko_obj         = gko::share(gko::initialize<gko::matrix::Csr<>>(
    {{2, -1, 0, 0}, {-1, 2, -1, 0}, {0, -1, 2, -1}, {0, 0, -1, 2}}, exec));
  auto linear_operator = GinkgoWrappers::
    linear_operator<double, double, MemorySpace::Default, MemorySpace::Default>(
      gko_obj);
  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> u(
    gko_obj->get_size()[0]);
  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> v(
    gko_obj->get_size()[0]);
  u             = 1.0;
  v             = 0.0;
  linear_operator.vmult(v, u);

  auto inverse_operator =
    GinkgoWrappers::inverse_operator<double,
                                     double,
                                     MemorySpace::Default,
                                     MemorySpace::Default>(
      gko_obj,
      gko::solver::Cg<>::build()
        .with_criteria(
          gko::stop::Iteration::build().with_max_iters(4u).on(exec))
        .on(exec));
  auto solution = u;
  u = 0.0;
  inverse_operator.vmult(u, v);

  solution -= u;
  TEST_ASSERT(solution.linfty_norm() < 1e-12);
}


int
main(int argc, char **argv)
{
  Kokkos::ScopeGuard guard(argc, argv);

  // Initialize deallog for test output.
  // This also reroutes deallog output to a file "output".
  initlog();

  can_create_vector_host();
  can_create_vector_host_const();
  can_create_vector_default();
  can_create_vector_with_executor();
  can_create_csr();
  can_create_csr_with_index_type();
  can_create_csr_with_value_index_type();
  can_create_linop_host();
  can_create_linop_default();
  can_create_inverse_linop();

  return 0;
}

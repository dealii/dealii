// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2022 by the deal.II authors
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


#ifndef dealii_ginkgo_interface_h
#define dealii_ginkgo_interface_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_GINKGO

#  include <deal.II/base/array_view.h>

#  include <deal.II/lac/sparse_matrix.h>

#  include <ginkgo/core/matrix/csr.hpp>
#  include <ginkgo/core/matrix/dense.hpp>
#  include <ginkgo/extensions/kokkos.hpp>



DEAL_II_NAMESPACE_OPEN

namespace GinkgoWrappers
{

  template <typename ValueType, typename MemorySpace>
  std::unique_ptr<gko::matrix::Dense<ValueType>>
  create_vector(std::shared_ptr<const gko::Executor> exec,
                ArrayView<ValueType, MemorySpace>   &array)
  {
    std::shared_ptr<const gko::Executor> array_exec =
      gko::ext::kokkos::get_executor(
        typename MemorySpace::kokkos_space::execution_space{});

    gko::dim<2> size = {static_cast<gko::size_type>(array.size()), 1};

    return gko::matrix::Dense<ValueType>::create(
      exec, size, gko::make_array_view(array_exec, size[0], array.data()), 1);
  }

  template <typename ValueType, typename MemorySpace>
  std::unique_ptr<gko::matrix::Dense<ValueType>>
  create_vector(ArrayView<ValueType, MemorySpace> &array)
  {
    std::shared_ptr<const gko::Executor> exec = gko::ext::kokkos::get_executor(
      typename MemorySpace::kokkos_space::execution_space{});

    return create_vector(std::move(exec), array);
  }

  template <typename ValueType, typename MemorySpace, typename = std::enable_if_t<!std::is_const_v<ValueType>>>
  std::unique_ptr<gko::matrix::Dense<ValueType>>
  create_vector(ArrayView<ValueType, MemorySpace> &&array)
  {
    std::shared_ptr<const gko::Executor> exec = gko::ext::kokkos::get_executor(
      typename MemorySpace::kokkos_space::execution_space{});

    return create_vector(std::move(exec), array);
  }

  template <typename ValueType, typename MemorySpace>
  std::unique_ptr<const gko::matrix::Dense<std::decay_t<ValueType>>>
  create_vector(std::shared_ptr<const gko::Executor>           exec,
                const ArrayView<ValueType, MemorySpace> &array)
  {
    std::shared_ptr<const gko::Executor> array_exec =
      gko::ext::kokkos::get_executor(
        typename MemorySpace::kokkos_space::execution_space{});

    gko::dim<2> size = {static_cast<gko::size_type>(array.size()), 1};

    return gko::matrix::Dense<std::decay_t<ValueType>>::create_const(
      exec,
      size,
      gko::make_const_array_view(array_exec, size[0], array.data()),
      1);
  }

  template <typename ValueType, typename MemorySpace>
  std::unique_ptr<const gko::matrix::Dense<std::decay_t<ValueType>>>
  create_vector(const ArrayView<ValueType, MemorySpace> &array)
  {
    std::shared_ptr<const gko::Executor> exec = gko::ext::kokkos::get_executor(
      typename MemorySpace::kokkos_space::execution_space{});

    return create_vector(std::move(exec), array);
  }

} // namespace GinkgoWrappers


DEAL_II_NAMESPACE_CLOSE

#endif

#endif // dealii_ginkgo_interface_h

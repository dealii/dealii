//-----------------------------------------------------------
//
//    Copyright (C) 2020 - 2021 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//-----------------------------------------------------------

// Test SUNDIALS' vector operations on N_Vector implementation. The N_Vectors
// are created by calling NVectorView on one of the internal vector types.

#include "../../include/deal.II/sundials/n_vector.h"

#include <deal.II/base/logstream.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include "../../include/deal.II/sundials/n_vector.templates.h"
#include <deal.II/sundials/n_vector.h>

#include "../tests.h"

using namespace SUNDIALS::internal;

// anonymous namespace groups helper functions for testing
namespace
{
  DeclExceptionMsg(NVectorTestError,
                   "The internal N_Vector implementation didn't pass a test.");


  IndexSet
  create_parallel_index_set()
  {
    const unsigned int current_process =
      Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
    const unsigned int n_processes =
      Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
    // give the zeroth process 2 dofs, the first 4, etc
    const types::global_dof_index n_local_dofs = 2 * (1 + current_process);
    const types::global_dof_index begin_dof =
      (current_process * (current_process + 1));
    const types::global_dof_index end_dof = begin_dof + n_local_dofs;

    const types::global_dof_index n_global_dofs =
      ((n_processes - 1) * (n_processes)) + 2 * n_processes;
    IndexSet local_dofs(n_global_dofs);
    local_dofs.add_range(begin_dof, end_dof);
    local_dofs.compress();
    AssertDimension(local_dofs.n_elements(), n_local_dofs);
    return local_dofs;
  }

  /**
   * Create a vector for testing and intialize all entries to @p value.
   */
  template <typename VectorType>
  VectorType
  create_test_vector(double value = 0.0);



  template <>
  Vector<double>
  create_test_vector(double value)
  {
    Vector<double> vector(3 /*size*/);
    vector = value;
    return vector;
  }



  template <>
  BlockVector<double>
  create_test_vector(double value)
  {
    const int           num_blocks = 2;
    const int           size_block = 3;
    BlockVector<double> vector(num_blocks, size_block);
    vector = value;
    return vector;
  }

  template <>
  LinearAlgebra::distributed::Vector<double>
  create_test_vector(double value)
  {
    IndexSet local_dofs = create_parallel_index_set();
    LinearAlgebra::distributed::Vector<double> vector(local_dofs,
                                                      MPI_COMM_WORLD);
    vector = value;
    return vector;
  }

  template <>
  LinearAlgebra::distributed::BlockVector<double>
  create_test_vector(double value)
  {
    const unsigned n_processes =
      Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
    const unsigned block_size =
      ((n_processes - 1) * (n_processes)) / 2 + n_processes;
    const auto partitioning =
      create_parallel_index_set().split_by_block({block_size, block_size});
    LinearAlgebra::distributed::BlockVector<double> vector(partitioning,
                                                           MPI_COMM_WORLD);
    vector = value;
    return vector;
  }


  template <>
  TrilinosWrappers::MPI::Vector
  create_test_vector(double value)
  {
    IndexSet local_dofs = create_parallel_index_set();

    TrilinosWrappers::MPI::Vector vector(local_dofs, MPI_COMM_WORLD);
    vector = value;
    return vector;
  }

  template <>
  TrilinosWrappers::MPI::BlockVector
  create_test_vector(double value)
  {
    const unsigned n_processes =
      Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
    const unsigned block_size =
      ((n_processes - 1) * (n_processes)) / 2 + n_processes;
    const auto partitioning =
      create_parallel_index_set().split_by_block({block_size, block_size});
    TrilinosWrappers::MPI::BlockVector vector(partitioning, MPI_COMM_WORLD);
    vector = value;
    return vector;
  }

  template <>
  PETScWrappers::MPI::Vector
  create_test_vector(double value)
  {
    IndexSet local_dofs = create_parallel_index_set();

    PETScWrappers::MPI::Vector vector(local_dofs, MPI_COMM_WORLD);
    vector = value;
    return vector;
  }

  template <>
  PETScWrappers::MPI::BlockVector
  create_test_vector(double value)
  {
    const unsigned n_processes =
      Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
    const unsigned block_size =
      ((n_processes - 1) * (n_processes)) / 2 + n_processes;
    const auto partitioning =
      create_parallel_index_set().split_by_block({block_size, block_size});
    PETScWrappers::MPI::BlockVector vector(partitioning, MPI_COMM_WORLD);
    vector = value;
    return vector;
  }

  template <typename VectorType>
  bool
  vector_equal(const VectorType &a, const VectorType &b)
  {
    auto elements_a = a.locally_owned_elements();
    auto elements_b = b.locally_owned_elements();
    if (elements_a.size() != elements_b.size())
      return false;

    const auto is_equal = [&](const IndexSet::size_type idxa,
                              const IndexSet::size_type idxb) {
      return std::fabs(a[idxa] - b[idxb]) < 1e-15;
    };

    return std::equal(elements_a.begin(),
                      elements_a.end(),
                      elements_b.begin(),
                      is_equal);
  }
} // namespace


template <typename VectorType>
void
test_nvector_view_unwrap()
{
  auto       vector         = create_test_vector<VectorType>();
  const auto const_vector   = create_test_vector<VectorType>();
  auto       n_vector       = make_nvector_view(vector);
  auto       const_n_vector = make_nvector_view(const_vector);

  Assert(n_vector != nullptr, NVectorTestError());
  Assert((n_vector)->content != nullptr, NVectorTestError());

  Assert(const_n_vector != nullptr, NVectorTestError());
  Assert((const_n_vector)->content != nullptr, NVectorTestError());

  // unwrap non-const as non-const
  auto *vector_unwrapped = unwrap_nvector<VectorType>(n_vector);
  Assert(vector_unwrapped == &vector, NVectorTestError());

  // unwrap non-const as const
  const auto *const_vector_unwrapped =
    unwrap_nvector_const<VectorType>(n_vector);
  Assert(const_vector_unwrapped == &vector, NVectorTestError());

  // unwrap const vector as const
  const_vector_unwrapped = unwrap_nvector_const<VectorType>(const_n_vector);
  Assert(const_vector_unwrapped == &const_vector, NVectorTestError());

  // unwrap const vector as non-const throws
  try
    {
      unwrap_nvector<VectorType>(const_n_vector);
      Assert(false, NVectorTestError());
    }
  catch (const ExcMessage &e)
    {
      const std::string msg(e.what());
      Assert(msg.find(
               "Tried to access a constant vector content as non-const.") !=
               std::string::npos,
             NVectorTestError());
    }

  deallog << "test_nvector_view_unwrap OK" << std::endl;
}



template <typename VectorType>
void
test_get_vector_id()
{
  auto vector   = create_test_vector<VectorType>();
  auto n_vector = make_nvector_view(vector);
  auto id       = N_VGetVectorID(n_vector);
  Assert(id == SUNDIALS_NVEC_CUSTOM, NVectorTestError());
  deallog << "test_get_vector_id OK" << std::endl;
}



template <typename VectorType>
void
test_clone()
{
  auto vector   = create_test_vector<VectorType>();
  auto n_vector = make_nvector_view(vector);
  auto cloned   = N_VClone(n_vector);

  Assert(cloned != nullptr, NVectorTestError());
  AssertDimension(unwrap_nvector<VectorType>(cloned)->size(), vector.size());

  N_VDestroy(cloned);
  deallog << "test_clone OK" << std::endl;
}



template <typename VectorType>
void
test_constraint_mask()
{
  // deal.II vectors
  auto constraint_0  = create_test_vector<VectorType>(0);
  auto constraint_1  = create_test_vector<VectorType>(1);
  auto constraint_2  = create_test_vector<VectorType>(2);
  auto constraint_m1 = create_test_vector<VectorType>(-1);
  auto constraint_m2 = create_test_vector<VectorType>(-2);

  auto vector0        = create_test_vector<VectorType>(0);
  auto vectorpositive = create_test_vector<VectorType>(3);
  auto vectornegative = create_test_vector<VectorType>(-3);
  auto mask           = create_test_vector<VectorType>(0);

  // N_Vectors
  auto view_constraint_0  = make_nvector_view(constraint_0);
  auto view_constraint_1  = make_nvector_view(constraint_1);
  auto view_constraint_2  = make_nvector_view(constraint_2);
  auto view_constraint_m1 = make_nvector_view(constraint_m1);
  auto view_constraint_m2 = make_nvector_view(constraint_m2);

  auto view_vector0        = make_nvector_view(vector0);
  auto view_vectorpositive = make_nvector_view(vectorpositive);
  auto view_vectornegative = make_nvector_view(vectornegative);

  auto view_mask = make_nvector_view(mask);

  auto check_true = [](bool in) {
    if (in == false)
      {
        deallog << "test_failed" << std::endl;
        AssertThrow(false, ExcInternalError());
      }
  };

  auto check_false = [](bool in) {
    if (in == true)
      {
        deallog << "test_failed" << std::endl;
        AssertThrow(false, ExcInternalError());
      }
  };

  check_true(N_VConstrMask(view_constraint_0, view_vector0, view_mask));
  check_true(N_VConstrMask(view_constraint_0, view_vectorpositive, view_mask));
  check_true(N_VConstrMask(view_constraint_0, view_vectornegative, view_mask));

  check_true(N_VConstrMask(view_constraint_1, view_vector0, view_mask));
  check_true(N_VConstrMask(view_constraint_1, view_vectorpositive, view_mask));
  check_false(N_VConstrMask(view_constraint_1, view_vectornegative, view_mask));

  check_false(N_VConstrMask(view_constraint_2, view_vector0, view_mask));
  check_true(N_VConstrMask(view_constraint_2, view_vectorpositive, view_mask));
  check_false(N_VConstrMask(view_constraint_2, view_vectornegative, view_mask));

  check_true(N_VConstrMask(view_constraint_m1, view_vector0, view_mask));
  check_false(
    N_VConstrMask(view_constraint_m1, view_vectorpositive, view_mask));
  check_true(N_VConstrMask(view_constraint_m1, view_vectornegative, view_mask));

  check_false(N_VConstrMask(view_constraint_m2, view_vector0, view_mask));
  check_false(
    N_VConstrMask(view_constraint_m2, view_vectorpositive, view_mask));
  check_true(N_VConstrMask(view_constraint_m2, view_vectornegative, view_mask));

  deallog << "test_constraint_mask OK" << std::endl;
}



template <typename VectorType>
void
test_min_quotient()
{
  auto vector1   = create_test_vector<VectorType>(12.);
  auto vector2   = create_test_vector<VectorType>(3.);
  auto n_vector1 = make_nvector_view(vector1);
  auto n_vector2 = make_nvector_view(vector2);

  auto m = N_VMinQuotient(n_vector1, n_vector2);
  if (m == 4)
    deallog << "test_min_quotient OK" << std::endl;
  else
    deallog << "test_min_quotient failed" << std::endl;
}



template <typename VectorType>
void
test_destroy()
{
  GrowingVectorMemory<VectorType>                   mem;
  typename GrowingVectorMemory<VectorType>::Pointer vector(mem);
  auto n_vector = make_nvector_view(*vector);

  Assert(n_vector != nullptr, NVectorTestError());
  auto cloned = N_VClone(n_vector);
  Assert(cloned != nullptr, NVectorTestError());

  // destroy memory-owning vector
  N_VDestroy(cloned);

  // clone empty(=without an actual vector attached) and destroy
  cloned = N_VCloneEmpty(n_vector);
  N_VDestroy(cloned);

  // NOTE: destroying a view vector is not possible and would lead to a double
  // delete at the end of scope when the destructor of NVectorView calls destroy
  // again
  // N_VDestroy(n_vector);

  // destroying a nullptr does nothing but is possible
  N_VDestroy(nullptr);

  // the destructor of NVectorView will do a destroy, so this is also
  // checked here
  deallog << "test_destroy OK" << std::endl;
}



template <typename VectorType,
          std::enable_if_t<is_serial_vector<VectorType>::value, int> = 0>
void
test_get_communicator()
{
  auto vector   = create_test_vector<VectorType>();
  auto n_vector = make_nvector_view(vector);
  // required by SUNDIALS: MPI-unaware vectors should return the nullptr
  Assert(N_VGetCommunicator(n_vector) == nullptr, NVectorTestError());

  deallog << "test_get_communicator OK" << std::endl;
}



template <typename VectorType,
          std::enable_if_t<!is_serial_vector<VectorType>::value, int> = 0>
void
test_get_communicator()
{
  auto vector   = create_test_vector<VectorType>();
  auto n_vector = make_nvector_view(vector);
  Assert(N_VGetCommunicator(n_vector) == MPI_COMM_WORLD, NVectorTestError());

  deallog << "test_get_communicator OK" << std::endl;
}



template <typename VectorType>
void
test_length()
{
  auto vector   = create_test_vector<VectorType>();
  auto n_vector = make_nvector_view(vector);
  Assert(N_VGetLength(n_vector) == static_cast<int>(vector.size()),
         NVectorTestError());

  deallog << "test_length OK" << std::endl;
}



template <typename VectorType>
void
test_linear_sum()
{
  auto va       = create_test_vector<VectorType>();
  auto vb       = create_test_vector<VectorType>();
  auto vc       = create_test_vector<VectorType>();
  auto expected = create_test_vector<VectorType>();

  auto nv_a = make_nvector_view(va);
  auto nv_b = make_nvector_view(vb);
  auto nv_c = make_nvector_view(vc);

  va = 1.0;
  vb = 2.0;

  expected = -3.0;

  // test sum into third vector
  N_VLinearSum(1.0, nv_a, -2.0, nv_b, nv_c);
  Assert(vector_equal(vc, expected), NVectorTestError());
  // repeat to test that sum overwrites initial content
  N_VLinearSum(1.0, nv_a, -2.0, nv_b, nv_c);
  Assert(vector_equal(vc, expected), NVectorTestError());

  // test store sum into one of the summands
  N_VLinearSum(1.0, nv_a, -2.0, nv_b, nv_a);
  Assert(vector_equal(va, expected), NVectorTestError());
  va = 1.0;

  N_VLinearSum(1.0, nv_a, -2.0, nv_b, nv_b);
  Assert(vector_equal(vb, expected), NVectorTestError());

  deallog << "test_linear_sum OK" << std::endl;
}



template <typename VectorType>
void
test_dot_product()
{
  const auto va   = create_test_vector<VectorType>(2.0);
  const auto vb   = create_test_vector<VectorType>(3.0);
  const auto size = va.size();

  auto nv_a = make_nvector_view(va);
  auto nv_b = make_nvector_view(vb);

  auto result   = N_VDotProd(nv_a, nv_b);
  auto expected = 6.0 * size;
  Assert(std::fabs(result - expected) < 1e-12, NVectorTestError());

  // test same vector
  result   = N_VDotProd(nv_a, nv_a);
  expected = 4.0 * size;
  Assert(std::fabs(result - expected) < 1e-12, NVectorTestError());

  deallog << "test_dot_product OK" << std::endl;
}



template <typename VectorType>
void
test_set_constant()
{
  auto vector   = create_test_vector<VectorType>();
  auto expected = create_test_vector<VectorType>();
  expected      = 1.0;

  auto n_vector = make_nvector_view(vector);

  N_VConst(1.0, n_vector);
  Assert(vector_equal(vector, expected), NVectorTestError());

  N_VConst(2.0, n_vector);
  expected = 2.0;
  Assert(vector_equal(vector, expected), NVectorTestError());

  deallog << "test_set_constant OK" << std::endl;
}



template <typename VectorType>
void
test_add_constant()
{
  const auto vector_a = create_test_vector<VectorType>(2.0);
  auto       result   = create_test_vector<VectorType>(-1.0);
  auto       expected = create_test_vector<VectorType>(3.0);

  auto nv_a      = make_nvector_view(vector_a);
  auto nv_result = make_nvector_view(result);

  N_VAddConst(nv_a, 1.0, nv_result);
  Assert(vector_equal(result, expected), NVectorTestError());

  // test overwrite self
  result = 2.0;
  N_VAddConst(nv_result, 1.0, nv_result);
  Assert(vector_equal(result, expected), NVectorTestError());

  deallog << "test_add_constant OK" << std::endl;
}



template <typename VectorType>
void
test_elementwise_product()
{
  const auto vector_a      = create_test_vector<VectorType>(2.0);
  const auto vector_b      = create_test_vector<VectorType>(3.0);
  auto       vector_result = create_test_vector<VectorType>();
  auto       expected      = create_test_vector<VectorType>(6.0);

  auto nv_a      = make_nvector_view(vector_a);
  auto nv_b      = make_nvector_view(vector_b);
  auto nv_result = make_nvector_view(vector_result);

  // result = a.*b
  N_VProd(nv_a, nv_b, nv_result);
  Assert(vector_equal(vector_result, expected), NVectorTestError());

  auto vector_c = create_test_vector<VectorType>(2.0);

  auto nv_c = make_nvector_view(vector_c);

  // c .*= b
  N_VProd(nv_c, nv_b, nv_c);
  Assert(vector_equal(vector_c, expected), NVectorTestError());

  // same as above with flipped arguments
  vector_c = 2.0;
  N_VProd(nv_b, nv_c, nv_c);
  Assert(vector_equal(vector_c, expected), NVectorTestError());

  // all operands the same
  vector_c = 2.0;
  expected = 4.0;
  N_VProd(nv_c, nv_c, nv_c);
  Assert(vector_equal(vector_c, expected), NVectorTestError());

  deallog << "test_elementwise_product OK" << std::endl;
}

template <typename VectorType>
void
test_elementwise_div()
{
  const auto vector_a      = create_test_vector<VectorType>(6.0);
  const auto vector_b      = create_test_vector<VectorType>(3.0);
  auto       vector_result = create_test_vector<VectorType>();
  auto       expected      = create_test_vector<VectorType>(2.0);

  auto nv_a      = make_nvector_view(vector_a);
  auto nv_b      = make_nvector_view(vector_b);
  auto nv_result = make_nvector_view(vector_result);

  N_VDiv(nv_a, nv_b, nv_result);
  Assert(vector_equal(vector_result, expected), NVectorTestError());

  deallog << "test_elementwise_div OK" << std::endl;
}



template <typename VectorType>
void
test_elementwise_inv()
{
  const auto vector_a      = create_test_vector<VectorType>(6.0);
  auto       vector_result = create_test_vector<VectorType>();
  auto       expected      = create_test_vector<VectorType>(1.0 / 6.0);

  auto nv_a      = make_nvector_view(vector_a);
  auto nv_result = make_nvector_view(vector_result);

  N_VInv(nv_a, nv_result);
  Assert(vector_equal(vector_result, expected), NVectorTestError());

  // test overwrite self
  vector_result = 6.0;
  N_VInv(nv_result, nv_result);
  Assert(vector_equal(vector_result, expected), NVectorTestError());

  deallog << "test_elementwise_inv OK" << std::endl;
}



template <typename VectorType>
void
test_elementwise_abs()
{
  const auto vector_a      = create_test_vector<VectorType>(-2.0);
  auto       vector_result = create_test_vector<VectorType>();
  auto       expected      = create_test_vector<VectorType>(2.0);

  auto nv_a      = make_nvector_view(vector_a);
  auto nv_result = make_nvector_view(vector_result);

  N_VAbs(nv_a, nv_result);
  Assert(vector_equal(vector_result, expected), NVectorTestError());

  // test overwrite self
  vector_result = -2.0;
  N_VAbs(nv_result, nv_result);
  Assert(vector_equal(vector_result, expected), NVectorTestError());

  deallog << "test_elementwise_abs OK" << std::endl;
}



template <typename VectorType>
void
test_weighted_rms_norm()
{
  const auto vector_a = create_test_vector<VectorType>(2.0);
  const auto vector_b = create_test_vector<VectorType>(3.0);

  auto nv_a = make_nvector_view(vector_a);
  auto nv_b = make_nvector_view(vector_b);

  const auto result = N_VWrmsNorm(nv_a, nv_b);
  Assert(std::fabs(result - 6.0) < 1e-12, NVectorTestError());

  deallog << "test_weighted_rms_norm OK" << std::endl;
}



template <typename VectorType>
void
test_weighted_rms_norm_mask()
{
  const auto vector_a = create_test_vector<VectorType>(2.0);
  const auto vector_b = create_test_vector<VectorType>(3.0);
  const auto vector_c = create_test_vector<VectorType>(1.0);
  const auto vector_d = create_test_vector<VectorType>(0.0);

  auto nv_a = make_nvector_view(vector_a);
  auto nv_b = make_nvector_view(vector_b);
  auto nv_c = make_nvector_view(vector_c);
  auto nv_d = make_nvector_view(vector_d);
  {
    const auto result = N_VWrmsNormMask(nv_a, nv_b, nv_c);
    Assert(std::fabs(result - 6.0) < 1e-12, NVectorTestError());
    deallog << "test_weighted_rms_norm_mask 1 OK" << std::endl;
  }
  {
    const auto result = N_VWrmsNormMask(nv_a, nv_b, nv_d);
    Assert(std::fabs(result) < 1e-12, NVectorTestError());
    deallog << "test_weighted_rms_norm_mask 0 OK" << std::endl;
  }
}



template <typename VectorType>
void
test_weighted_l2_norm()
{
  const auto vector_a = create_test_vector<VectorType>(2.0);
  const auto vector_b = create_test_vector<VectorType>(3.0);

  auto nv_a = make_nvector_view(vector_a);
  auto nv_b = make_nvector_view(vector_b);

  const auto result = N_VWL2Norm(nv_a, nv_b);

  auto vector_a_reference = create_test_vector<VectorType>(2.0);
  auto vector_b_reference = create_test_vector<VectorType>(3.0);
  vector_a_reference.scale(vector_b_reference);
  const auto result_reference = vector_a_reference.l2_norm();

  Assert(std::fabs(result - result_reference) < 1e-12, NVectorTestError());

  deallog << "test_weighted_l2_norm OK" << std::endl;
}



template <typename VectorType>
void
test_max_norm()
{
  const auto vector_a = create_test_vector<VectorType>(2.0);
  const auto vector_b = create_test_vector<VectorType>(-3.0);

  auto nv_a = make_nvector_view(vector_a);
  auto nv_b = make_nvector_view(vector_b);

  auto result = N_VMaxNorm(nv_a);
  Assert(std::fabs(result - 2.0) < 1e-12, NVectorTestError());

  result = N_VMaxNorm(nv_b);
  Assert(std::fabs(result - 3.0) < 1e-12, NVectorTestError());

  deallog << "test_max_norm OK" << std::endl;
}



template <typename VectorType>
void
test_l1_norm()
{
  const auto vector_a = create_test_vector<VectorType>(2.0);
  const auto vector_b = create_test_vector<VectorType>(-3.0);

  auto nv_a = make_nvector_view(vector_a);
  auto nv_b = make_nvector_view(vector_b);

  auto result = N_VL1Norm(nv_a);
  Assert(std::fabs(result - vector_a.l1_norm()) < 1e-12, NVectorTestError());

  result = N_VL1Norm(nv_b);
  Assert(std::fabs(result - vector_b.l1_norm()) < 1e-12, NVectorTestError());

  deallog << "test_l1_norm OK" << std::endl;
}



template <typename VectorType>
void
test_min_element()
{
  const auto vector_a = create_test_vector<VectorType>(2.0);
  const auto vector_b = create_test_vector<VectorType>(-3.0);

  auto nv_a = make_nvector_view(vector_a);
  auto nv_b = make_nvector_view(vector_b);

  auto result = N_VMin(nv_a);
  Assert(std::fabs(result - 2.0) < 1e-12, NVectorTestError());

  result = N_VMin(nv_b);
  Assert(std::fabs(result - (-3.0)) < 1e-12, NVectorTestError());

  deallog << "test_min_element OK" << std::endl;
}



template <typename VectorType>
void
test_scale()
{
  const auto vector_a = create_test_vector<VectorType>(2.0);
  auto       vector_b = create_test_vector<VectorType>(2.0);
  const auto expected = create_test_vector<VectorType>(4.0);

  auto nv_a = make_nvector_view(vector_a);
  auto nv_b = make_nvector_view(vector_b);

  // b *= 2
  N_VScale(2.0, nv_b, nv_b);
  Assert(vector_equal(vector_b, expected), NVectorTestError());

  // b = 2*a
  N_VScale(2.0, nv_a, nv_b);
  Assert(vector_equal(vector_b, expected), NVectorTestError());

  deallog << "test_scale OK" << std::endl;
}



template <typename VectorType>
void
run_all_tests(const std::string &prefix)
{
  LogStream::Prefix p(prefix);
  // test conversion between vectors
  test_nvector_view_unwrap<VectorType>();
  test_get_vector_id<VectorType>();

  // test vector operations
  test_add_constant<VectorType>();
  test_clone<VectorType>();
  test_constraint_mask<VectorType>();
  test_destroy<VectorType>();
  test_dot_product<VectorType>();
  test_elementwise_abs<VectorType>();
  test_elementwise_div<VectorType>();
  test_elementwise_inv<VectorType>();
  test_elementwise_product<VectorType>();
  test_get_communicator<VectorType>();
  test_l1_norm<VectorType>();
  test_length<VectorType>();
  test_linear_sum<VectorType>();
  test_max_norm<VectorType>();
  test_min_element<VectorType>();
  test_min_quotient<VectorType>();
  test_scale<VectorType>();
  test_set_constant<VectorType>();
  test_weighted_l2_norm<VectorType>();
  test_weighted_rms_norm_mask<VectorType>();
  test_weighted_rms_norm<VectorType>();
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log_all;

  using VectorType = Vector<double>;

  run_all_tests<Vector<double>>("Vector<double>");
  run_all_tests<BlockVector<double>>("BlockVector<double>");
  run_all_tests<LinearAlgebra::distributed::Vector<double>>(
    "LinearAlgebra::distributed::Vector<double>");
  run_all_tests<LinearAlgebra::distributed::BlockVector<double>>(
    "LinearAlgebra::distributed::BlockVector<double>");
  run_all_tests<TrilinosWrappers::MPI::Vector>("TrilinosWrappers::MPI::Vector");
  run_all_tests<TrilinosWrappers::MPI::BlockVector>(
    "TrilinosWrappers::MPI::BlockVector");
  run_all_tests<PETScWrappers::MPI::Vector>("PETScWrappers::MPI::Vector");
  run_all_tests<PETScWrappers::MPI::BlockVector>(
    "PETScWrappers::MPI::BlockVector");

  // although the memory would be cleared in ~MPI_InitFinalize it needs to be
  // done manually to satisfy the PETSc memory check inside ~MPILogInitAll,
  // which is invoked first
  GrowingVectorMemory<PETScWrappers::MPI::Vector>::release_unused_memory();
  GrowingVectorMemory<PETScWrappers::MPI::BlockVector>::release_unused_memory();
}

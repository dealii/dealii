// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check Utilities::MPI::mpi_type_id_for_type


#include <deal.II/base/mpi.h>
#include <deal.II/base/template_constraints.h>

#include <cstdint>

#include "../tests.h"



template <typename T>
using mpi_type_id_for_type_t =
  decltype(Utilities::MPI::internal::MPIDataTypes::mpi_type_id(
    &std::declval<const T &>()));


template <typename T>
void
test(const char *type, MPI_Datatype mpi_type)
{
  deallog << std::boolalpha;
  deallog << "--------- Checking <" << type << "> ---------------" << std::endl;
  deallog
    << "T is supported: "
    << internal::is_supported_operation<mpi_type_id_for_type_t, T> << std::endl;
  deallog << "T[] is supported: "
          << internal::is_supported_operation<mpi_type_id_for_type_t,
                                              T[]> << std::endl;
  deallog << "T* is supported: "
          << internal::is_supported_operation<mpi_type_id_for_type_t,
                                              T *> << std::endl;

  Assert(Utilities::MPI::mpi_type_id_for_type<T> == mpi_type,
         ExcInternalError());
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  initlog();

  // Go through the entire list of MPI-supported data types:
#define TEST(a, b) test<a>(#a, b)
  TEST(signed char, MPI_SIGNED_CHAR);
  TEST(unsigned char, MPI_UNSIGNED_CHAR);
  TEST(wchar_t, MPI_WCHAR);
  TEST(short, MPI_SHORT);
  TEST(unsigned short, MPI_UNSIGNED_SHORT);
  TEST(int, MPI_INT);
  TEST(unsigned int, MPI_UNSIGNED);
  TEST(long, MPI_LONG);
  TEST(unsigned long, MPI_UNSIGNED_LONG);
  TEST(long long, MPI_LONG_LONG_INT);
  TEST(long long int, MPI_LONG_LONG);
  TEST(unsigned long long, MPI_UNSIGNED_LONG_LONG);
  TEST(float, MPI_FLOAT);
  TEST(double, MPI_DOUBLE);
  TEST(long double, MPI_LONG_DOUBLE);

  // The following types have their own MPI tags, but these tags are
  // aliases for some of the types above and consequently are
  // indistinguishable even though the tags are distinguishable. That
  // is: while the type system cannot distinguish between 'unsigned
  // char' and 'std::uint8_t', the tags 'MPI_SIGNED_CHAR' and 'MPI_INT8_T'
  // are different. We cannot test, then, that the tag we get for
  // these types equals their tags...

  // TEST(int8_t, MPI_INT8_T);
  // TEST(int16_t, MPI_INT16_T);
  // TEST(int32_t, MPI_INT32_T);
  // TEST(int64_t, MPI_INT64_T);
  // TEST(std::uint8_t, MPI_UINT8_T);
  // TEST(std::uint16_t, MPI_UINT16_T);
  // TEST(std::uint32_t, MPI_UINT32_T);
  // TEST(std::uint64_t, MPI_UINT64_T);

  // Now make sure that the following type is not supported:
  struct X
  {};

  using T = X;
  deallog << "--------- Checking <X> ---------------" << std::endl;
  deallog
    << "T is supported: "
    << internal::is_supported_operation<mpi_type_id_for_type_t, T> << std::endl;
}

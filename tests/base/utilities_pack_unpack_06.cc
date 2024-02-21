// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// tests the 'allow for compression' feature of
// Utilities::pack/unpack
// (based upon "utilities_pack_unpack_04" and "*_05")


#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

// We need this header to define the exception when we compile with zlib
#include <boost/iostreams/filter/gzip.hpp>

#include "../tests.h"



template <int N, int dim>
void
check(const double (&array)[N], const Point<dim>(&point))
{
  // ----- PACK -----
  std::vector<char> array_compressed, array_uncompressed;
  std::vector<char> point_compressed, point_uncompressed;

  // default option should work for compression
  Utilities::pack(array, array_compressed);
  Utilities::pack(point, point_compressed);

  Utilities::pack(array, array_uncompressed, false);
  Utilities::pack(point, point_uncompressed, false);

  // check if compression has been invoked by comparing sizes
  deallog << "check packed array with compression: "
          << (sizeof(array) > sizeof(char) * array_compressed.size() ? "OK" :
                                                                       "Failed")
          << std::endl;
  deallog << "check packed point with compression: "
          << (sizeof(point) > sizeof(char) * point_compressed.size() ? "OK" :
                                                                       "Failed")
          << std::endl;

  deallog << "check packed array without compression: "
          << (sizeof(array) <= sizeof(char) * array_uncompressed.size() ?
                "OK" :
                "Failed")
          << std::endl;
  deallog << "check packed point without compression: "
          << (sizeof(point) <= sizeof(char) * point_uncompressed.size() ?
                "OK" :
                "Failed")
          << std::endl;



  // ----- UNPACK -----
  // default option should work for compression
  double unpacked_array_compressed[N];
  Utilities::unpack(array_compressed.cbegin(),
                    array_compressed.cend(),
                    unpacked_array_compressed);
  Point<dim> unpacked_point_compressed =
    Utilities::unpack<Point<dim>>(point_compressed.cbegin(),
                                  point_compressed.cend());

  double unpacked_array_uncompressed[N];
  Utilities::unpack(array_uncompressed.cbegin(),
                    array_uncompressed.cend(),
                    unpacked_array_uncompressed,
                    false);
  Point<dim> unpacked_point_uncompressed =
    Utilities::unpack<Point<dim>>(point_uncompressed.cbegin(),
                                  point_uncompressed.cend(),
                                  false);

  // check if unpacked results are okay
  bool equal_array = true;
  for (unsigned int i = 0; i < N; ++i)
    if (array[i] != unpacked_array_compressed[i])
      {
        equal_array = false;
        break;
      }
  deallog << "check unpacked array with compression: "
          << (equal_array ? "OK" : "Failed") << std::endl;
  deallog << "check unpacked point with compression: "
          << (point.distance(unpacked_point_compressed) < 1e-12 ? "OK" :
                                                                  "Failed")
          << std::endl;

  equal_array = true;
  for (unsigned int i = 0; i < N; ++i)
    if (array[i] != unpacked_array_uncompressed[i])
      {
        equal_array = false;
        break;
      }
  deallog << "check unpacked array without compression: "
          << (equal_array ? "OK" : "Failed") << std::endl;
  deallog << "check unpacked point without compression: "
          << (point.distance(unpacked_point_uncompressed) < 1e-12 ? "OK" :
                                                                    "Failed")
          << std::endl;



  // try something nasty
  try
    {
      Point<dim> forbidden =
        Utilities::unpack<Point<dim>>(point_compressed.cbegin(),
                                      point_compressed.cend(),
                                      false);
    }
  catch (const boost::archive::archive_exception &)
    {
      deallog << "unpacking compressed point without decompression failed!"
              << std::endl;
    }
  try
    {
      double forbidden[N];
      Utilities::unpack(array_compressed.cbegin(),
                        array_compressed.cend(),
                        forbidden,
                        false);
    }
  catch (const boost::archive::archive_exception &)
    {
      deallog << "unpacking compressed array without decompression failed!"
              << std::endl;
    }

  try
    {
      Point<dim> forbidden =
        Utilities::unpack<Point<dim>>(point_uncompressed.cbegin(),
                                      point_uncompressed.cend(),
                                      true);
    }
  catch (const boost::iostreams::gzip_error &)
    {
      deallog << "unpacking uncompressed point with decompression failed!"
              << std::endl;
    }
  try
    {
      double forbidden[N];
      Utilities::unpack(array_uncompressed.cbegin(),
                        array_uncompressed.cend(),
                        forbidden,
                        true);
    }
  catch (const boost::iostreams::gzip_error &)
    {
      deallog << "unpacking uncompressed array with decompression failed!"
              << std::endl;
    }
}


void
test()
{
  // pick large data types and arrays that could be compressed,
  // and check for both compression options
  const unsigned int N  = 10000;
  Point<N>           p2 = random_point<N>();
  double             x2[N];
  for (unsigned int i = 0; i < N; ++i)
    x2[i] = i;

  check(x2, p2);

  deallog << "OK!" << std::endl;
}

int
main()
{
  initlog();

  test();
}

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test that we can use BOOST iostreams

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "../tests.h"


void
test()
{
  // create a stream where we read from gzipped data
  boost::iostreams::filtering_istream in;
  in.push(boost::iostreams::basic_gzip_decompressor<>());
  in.push(
    boost::iostreams::file_source(SOURCE_DIR "/boost_iostreams_01.data.gz"));
  AssertThrow(in, ExcIO());

  // read the two numbers that are in this file
  int    i;
  double d;
  in >> i >> d;

  deallog << i << ' ' << d << std::endl;

  // make sure that we got here just fine but that there is nothing else
  AssertThrow(in, ExcIO());
  in >> i;
  AssertThrow(!in, ExcIO());
}



int
main()
{
  initlog();

  test();
}

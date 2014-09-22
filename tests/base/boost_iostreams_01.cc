// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// test that we can use BOOST iostreams

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>


void test ()
{
  // create a stream where we read from gzipped data
  boost::iostreams::filtering_istream in;
  in.push(boost::iostreams::basic_gzip_decompressor<>());
  in.push(boost::iostreams::file_source(SOURCE_DIR "/boost_iostreams_01.data.gz"));
  AssertThrow (in, ExcIO());

  // read the two numbers that are in this file
  int i;
  double d;
  in >> i >> d;

  deallog << i << ' ' << d << std::endl;

  // make sure that we got here just fine but that there is nothing else
  AssertThrow (in, ExcIO());
  in >> i;
  AssertThrow (!in, ExcIO());  
}




int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}

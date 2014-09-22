// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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

// Author: Guido Kanschat

// A little program reading a grid *.inp and writing it to *.eps.
// Some more functionality should be added som time.

#include <grid/tria.h>
#include <grid/grid_in.h>
#include <grid/grid_out.h>

#include <iostream>
#include <fstream>
#include <string>

#include <unistd.h>

using namespace std;

template <int dim>
void convert(const char* infile,
	     const char* outfile,
	     GridOut::OutputFormat oformat)
{
  Triangulation<dim> tr;
  GridIn<dim> gin;
  gin.attach_triangulation(tr);
  gin.read(infile);
  
  GridOut gout;
  GridOutFlags::DX dx_flags(true, true, true, false, true);
  gout.set_flags(dx_flags);
  
  ofstream out(outfile);
  gout.write(tr, out, oformat);
}


int main(int argc, char** argv)
{
  if (argc<4)
    {
      cerr << "Usage: " << argv[0]
	   << " dim infile outfile" << endl;
      exit(1);
    }

  const unsigned int d = atoi(argv[1]);
  
  std::string outstring(argv[3]);
  GridOut::OutputFormat oformat = GridOut::eps;
  
  const unsigned int dotpos = outstring.find_last_of('.');
  if (dotpos < outstring.length())
    {
      std::string ext = outstring.substr(dotpos+1);
      if (ext == "inp")
	ext = "ucd";
      
      oformat = GridOut::parse_output_format(ext);
    }

  if (d == 2)
    convert<2>(argv[2], argv[3], oformat);
  else if (d == 3)
    convert<3>(argv[2], argv[3], oformat);
}

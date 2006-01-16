//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

#include "../tests.h"
#include <base/data_out_base.h>
#include <base/logstream.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include "patches.h"

// Output data on repetitions of the unit hypercube

template <int dim, int spacedim>
void check(DataOutBase::EpsFlags flags,
	   std::ostream& out)
{
  const unsigned int np = 4;
  
  std::vector<DataOutBase::Patch<dim, spacedim> > patches(np);
  
  create_patches(patches);
  
  std::vector<std::string> names(5);
  names[0] = "x1";
  names[1] = "x2";
  names[2] = "x3";
  names[3] = "x4";
  names[4] = "i";
  DataOutBase::write_eps(patches, names, flags, out);
}


template<int dim, int spacedim>
void check_all()
{
  DataOutBase::EpsFlags flags;
  if (true) {
    std::ofstream out("data_out_base_dx/fff3ff22.dx");
    check<dim,spacedim>(flags, out);
  }
  flags.int_binary = true;
  if (true) {
    std::ofstream out("data_out_base_dx/tff3ff22.dx");
    check<dim,spacedim>(flags, out);
  }
  flags.coordinates_binary = true;
  if (true) {
    std::ofstream out("data_out_base_dx/ttf3ff22.dx");
    check<dim,spacedim>(flags, out);
  }
  flags.data_binary = true;
  if (true) {
    std::ofstream out("data_out_base_dx/ttt3ff22.dx");
    check<dim,spacedim>(flags, out);
  }
}

int main()
{
  check_all<1,1>();
  check_all<1,2>();
  check_all<2,2>();
  check_all<2,3>();
  check_all<3,3>();  
}

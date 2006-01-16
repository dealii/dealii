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
#include <stdio.h>

#include "patches.h"

// Output data on repetitions of the unit hypercube

template <int dim, int spacedim>
void check(DataOutBase::VtkFlags flags,
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
  DataOutBase::write_vtk(patches, names, flags, out);
}


template<int dim, int spacedim>
void check_all()
{
  char name[100];
  DataOutBase::VtkFlags flags;
  if (true) {
    sprintf(name, "data_out_base_vtk/%d%d.vtk", dim, spacedim);
    std::ofstream out(name);
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

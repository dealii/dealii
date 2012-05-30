//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006, 2007, 2010, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// write the pvd master record for parallel visualization through the
// vtu file format

#include "../tests.h"
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/logstream.h>

#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>

#include "patches.h"



std::vector<DataOutBase::Patch<2,2> > patches;

class DataOutX : public DataOutInterface<2,2>
{
    virtual
    const std::vector< ::DataOutBase::Patch<2,2> > &
    get_patches () const
      {
	return patches;
      }

    virtual
    std::vector<std::string>
    get_dataset_names () const
      {
	return std::vector<std::string>();
      }
};


template <int dim, int spacedim>
void check(std::ostream& out)
{
  std::vector<std::pair<double,std::string> > names(5);
  names[0] = std::make_pair(0,"x1");
  names[1] = std::make_pair(1,"x2");
  names[2] = std::make_pair(1e1,"x3");
  names[3] = std::make_pair(3.141,"d");
  names[4] = std::make_pair(42e19,"i");

  DataOutX x;
  x.write_pvd_record (out, names);
}



int main()
{
  std::ofstream logfile("data_out_base_pvd/output");
  check<2,2>(logfile);
}

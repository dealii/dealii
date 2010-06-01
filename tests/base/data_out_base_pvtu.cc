//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006, 2007, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// write the pvtu master record for parallel visualization through the
// vtu file format

#include "../tests.h"
#include <base/data_out_base.h>
#include <base/logstream.h>

#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>

#include "patches.h"



std::vector<DataOutBase::Patch<2,2> > patches;
std::vector<std::string> names;

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
	return names;
      }
};


template <int dim, int spacedim>
void check(std::ostream& out)
{
  const unsigned int np = 4;
  patches.resize (np);

  create_patches(patches);

  names.resize (5);
  names[0] = "x1";
  names[1] = "x2";
  names[2] = "x3";
  names[3] = "x4";
  names[4] = "i";

  std::vector<std::string> filenames = names;

  DataOutX x;
  x.write_pvtu_record (out, filenames);
}



int main()
{
  std::ofstream logfile("data_out_base_pvtu/output");
  check<2,2>(logfile);
}

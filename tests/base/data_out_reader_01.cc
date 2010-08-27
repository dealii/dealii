//-----------------------------------------------------------------------------
//    $Id: data_out_base_gnuplot.cc 20952 2010-04-06 15:02:46Z bangerth $
//    Version: $Name$ 
//
//    Copyright (C) 2006, 2007 by the deal.II authors
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
#include <iomanip>
#include <fstream>
#include <string>

#include "patches.h"

// test DataOutReader::merge

void cat_file(const char * filename)
{
  std::ifstream in(filename);
  while (in)
    {
      std::string s;
      std::getline(in, s);
      deallog.get_file_stream() << s << "\n";
    }
}

template <int dim, int spacedim>
void check(const char * name, DataOutBase::Deal_II_IntermediateFlags flags,
	   std::ostream& out)
{
  const unsigned int np = 1;
  
  std::vector<DataOutBase::Patch<dim, spacedim> > patches(np);
  
  create_patches(patches);
  
  std::vector<std::string> names(5);
  names[0] = "x1";
  names[1] = "x2";
  names[2] = "x3";
  names[3] = "x4";
  names[4] = "i";
  std::vector<std_cxx1x::tuple<unsigned int, unsigned int, std::string> > vectors;
  
  DataOutBase::write_deal_II_intermediate(patches, names, vectors, flags, out);

  DataOutReader<dim,spacedim> data;
  {
    std::ifstream input(name);    
    data.read (input);
  }  
  DataOutReader<dim,spacedim> additional_data;
  {
    std::ifstream input(name);
    additional_data.read (input);
  }
   
  data.merge (additional_data);

  {
    std::ofstream out2( "data_out_reader_01/outfile");
    data.write_deal_II_intermediate (out2);
  }  

  cat_file("data_out_reader_01/outfile");
}


template<int dim, int spacedim>
void check_all(std::ostream& log)
{
  char name[100];
  const char* format = "data_out_reader_01/%d%d.d2";
  DataOutBase::Deal_II_IntermediateFlags flags;
    {
      sprintf(name, format, dim, spacedim, "");

      std::ofstream out(name);

      check<dim,spacedim>(name, flags, out);
    }
}

int main()
{
  std::ofstream logfile("data_out_reader_01/output");
  deallog.attach(logfile);

  //  check_all<1,1>(logfile);
  //  check_all<1,2>(logfile);
  check_all<2,2>(logfile);
  //  check_all<2,3>(logfile);
  //  check_all<3,3>(logfile);  
}

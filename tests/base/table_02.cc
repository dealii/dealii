// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2014 by the deal.II authors
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


// check TableBase::fill


#include "../tests.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>

#include <deal.II/base/logstream.h>
#include <deal.II/base/table.h>


int
main ()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(0);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  // rank=1
  {
    deallog << "rank=1" << std::endl;
    
    const double entries[] = { 1, 2, 3 };

    Table<1,double> t (3);
    t.fill (entries, true);

    for (unsigned int i=0; i<t.size()[0]; ++i)
      deallog << t[i] << ' ';
    deallog << std::endl;

    // passing false as second argument shouldn't
    // make a difference for rank-1 tables
    t.fill (entries, false);
    for (unsigned int i=0; i<t.size()[0]; ++i)
      deallog << t[i] << ' ';
    deallog << std::endl;
  }


  // rank=2
  {
    deallog << "rank=2" << std::endl;
    
    const double entries[] = { 1, 2, 3, 4, 5, 6 };

    // create a 2x3 table from this
    Table<2,double> t (2,3);
    t.fill (entries, true);

    for (unsigned int i=0; i<t.size()[0]; ++i)
      {
	for (unsigned int j=0; j<t.size()[1]; ++j)
	  deallog << t[i][j] << ' ';
	deallog << std::endl;
      }

    // same data, same table, but filled in transpose ordering
    t.fill (entries, false);

    for (unsigned int i=0; i<t.size()[0]; ++i)
      {
	for (unsigned int j=0; j<t.size()[1]; ++j)
	  deallog << t[i][j] << ' ';
	deallog << std::endl;
      }
  }


  // rank=3
  {
    deallog << "rank=3" << std::endl;
    
    const double entries[] = { 1, 2, 3, 4, 5, 6,
			       7, 8, 9, 10, 11, 12 };

    // create a 2x3x2 table from this
    Table<3,double> t (2,3,2);
    t.fill (entries, true);

    for (unsigned int i=0; i<t.size()[0]; ++i)
      {
	for (unsigned int j=0; j<t.size()[1]; ++j)
	  {
	    deallog << '(';
	    for (unsigned int k=0; k<t.size()[2]; ++k)
	      deallog << t[i][j][k] << ' ';
	    deallog << ')';
	  }
	deallog << std::endl;
      }

    // same data, same table, but filled in transpose ordering
    t.fill (entries, false);

    for (unsigned int i=0; i<t.size()[0]; ++i)
      {
	for (unsigned int j=0; j<t.size()[1]; ++j)
	  {
	    deallog << '(';
	    for (unsigned int k=0; k<t.size()[2]; ++k)
	      deallog << t[i][j][k] << ' ';
	    deallog << ')';
	  }
	deallog << std::endl;
      }
  }
}




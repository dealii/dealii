//----------------------------  $RCSfile$  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  $RCSfile$  ---------------------------


#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include <base/logstream.h>
#include <base/vector2d.h>

const int entries[9] = { 11,12,13,21,22,23,31,32,33 };

int
main ()
{
  std::ofstream logfile("vector2d.output");
  logfile.setf(std::ios::fixed);
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);

  vector2d<double> Td(3,3);
  vector2d<int> Ti(3,3);
  
  for (unsigned int i=0; i<9; ++i)
    {
      Td[i/3][i%3] = entries[i];
      Ti[i/3][i%3] = entries[i];
    };

  for (unsigned int i=0; i<3; ++i)
    for (unsigned int j=0; j<3; ++j)
      {
	Assert (Td[i][j] == Td(i,j), ExcInternalError());
	Assert (Ti[i][j] == Ti(i,j), ExcInternalError());
	Assert (Ti[i][j] == Td(i,j), ExcInternalError());

	logfile << i << " " << j << " " << Td[i][j] << " ok" << std::endl;
      };
};


      

//----------------------------  fe_tools_01c.cc  ---------------------------
//    fe_tools_01c.cc,v 1.1 2003/02/24 15:48:02 wolf Exp
//    Version: 
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_tools_01c.cc  ---------------------------


// common framework for the various dof_tools_*.cc tests

#include "../tests.h"
#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_handler.h>
#include <fe/fe_tools.h>
#include <fe/fe_q.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>


// check invertability of the map from 
//   DoFTools::hierarchic_to_lexicographic_numbering
// to  
//   DoFTools::lexicographic_to_hierarchic_numbering


template <int dim>
void
check (const FE_Q<dim>   &fe,
       const std::string &name)
{
  deallog << "Checking " << name
          << " in " << dim << "d:"
          << std::endl;

  std::vector<unsigned int> n1(fe.dofs_per_cell);
  FETools::hierarchic_to_lexicographic_numbering (fe, n1);

  std::vector<unsigned int> n2(fe.dofs_per_cell);
  FETools::lexicographic_to_hierarchic_numbering (fe, n2);

  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
    {
      Assert (n1[i] < fe.dofs_per_cell,
              ExcInternalError());
      Assert (n2[i] < fe.dofs_per_cell,
              ExcInternalError());
      Assert (n1[n2[i]] == i,
              ExcInternalError());
      Assert (n2[n1[i]] == i,
              ExcInternalError());      
      
      deallog << n1[n2[i]] << " ";
    }
  deallog << std::endl;
}





#define CHECK(EL,deg,dim)\
 { FE_ ## EL<dim> EL(deg);   \
   check(EL, #EL #deg); }

#define CHECK_ALL(EL,deg)\
 { CHECK(EL,deg,1); \
   CHECK(EL,deg,2); \
   CHECK(EL,deg,3); \
 }


int
main()
{
  try
    {
      std::ofstream logfile("fe_tools_01c.output");
      logfile.precision (2);
      deallog.attach(logfile);
      deallog.depth_console(0);

      CHECK_ALL(Q,1);
      CHECK_ALL(Q,2);
      CHECK_ALL(Q,3);
      CHECK_ALL(Q,4);
  
      return 0;
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    }
  catch (...) 
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };
}


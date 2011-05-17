//----------------------------  fe_tools_01c.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_tools_01c.cc  ---------------------------


// common framework for the various fe_tools_*.cc tests

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <iomanip>
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
      std::ofstream logfile("fe_tools_01c/output");
      deallog << std::setprecision (2);
      deallog.attach(logfile);
      deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

      CHECK_ALL(Q,1);
      CHECK_ALL(Q,2);
      CHECK_ALL(Q,3);
      CHECK_ALL(Q,4);
  
      return 0;
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      deallog << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    }
  catch (...) 
    {
      deallog << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      deallog << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };
}


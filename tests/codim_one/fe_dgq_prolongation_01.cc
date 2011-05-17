//----------------------------  fe_dgq_prolongation_01.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_dgq_prolongation_01.cc  ---------------------------

// this is the root cause for solution_tranfer_01: the prolongation
// matrices for FE_DGQ<dim-1,dim> were not computed at all

#include "../tests.h"
#include <deal.II/fe/fe_dgq.h>


#include <fstream>

using namespace dealii;


int main ()
{
  std::ofstream logfile("fe_dgq_prolongation_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  const unsigned int spacedim = 2;
  const unsigned int dim = spacedim-1;

  for (unsigned int degree=0; degree<3; ++degree)
    {
      deallog << "Degree=" << degree << std::endl;
      FE_DGQ<dim,spacedim> fe(degree);
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
	  deallog << fe.get_prolongation_matrix(0)(i,j) << std::endl;
    }

  return 0;
}

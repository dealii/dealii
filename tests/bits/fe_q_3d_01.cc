//----------------------------  fe_q_3d_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_q_3d_01.cc  ---------------------------


// we used to have constraint matrices for the FE_Q elements precomputed and
// stored, but later moved to compute them on the fly. make sure these
// correspond to the previously available precomputed ones for 3d and q=1,2


#include "../tests.h"
#include <fe/fe_q.h>

#include <fstream>
#include <iostream>


// matrices taken from the old file deal.II/source/fe/fe_q_3d.cc
namespace FE_Q_3d 
{
  static const double constraint_q1[] =
  {
	.25,.25,.25,.25,
	.5,.5,0.,0.,
	0.,.5,.5,0.,
	0.,0.,.5,.5,
	.5,0.,0.,.5
  };

  static const double constraint_q2[] =
  {
	0,0,0,0,0,0,0,0,1,
	0,0,0,0,1,0,0,0,0,
	0,0,0,0,0,1,0,0,0,
	0,0,0,0,0,0,1,0,0,
	0,0,0,0,0,0,0,1,0,
	0,0,0,0,.375,0,-.125,0,.75,
	0,0,0,0,0,.375,0,-.125,.75,
	0,0,0,0,-.125,0,.375,0,.75,
	0,0,0,0,0,-.125,0,.375,.75,
	.375,-.125,0,0,.75,0,0,0,0,
	-.125,.375,0,0,.75,0,0,0,0,
	0,.375,-.125,0,0,.75,0,0,0,
	0,-.125,.375,0,0,.75,0,0,0,
	0,0,-.125,.375,0,0,.75,0,0,
	0,0,.375,-.125,0,0,.75,0,0,
	.375,0,0,-.125,0,0,0,.75,0,
	-.125,0,0,.375,0,0,0,.75,0,
	.140625,-.046875,.015625,-.046875,.28125,-.09375,-.09375,.28125,.5625,
	-.046875,.140625,-.046875,.015625,.28125,.28125,-.09375,-.09375,.5625,
	.015625,-.046875,.140625,-.046875,-.09375,.28125,.28125,-.09375,.5625,
	-.046875,.015625,-.046875,.140625,-.09375,-.09375,.28125,.28125,.5625
  };
}


namespace Matrices
{
  const double * const 
  constraint_matrices[] =
  {
        FE_Q_3d::constraint_q1,
        FE_Q_3d::constraint_q2
  };

  const unsigned int 
  n_constraint_matrices
  = sizeof(constraint_matrices) /
    sizeof(constraint_matrices[0]);
}



void check ()
{
                                   // check for q=1,2
  for (unsigned int q=1; q<=2; ++q)
    {
      deallog << "q=" << q << std::endl;
      
      FE_Q<3> fe(q);
    
      FullMatrix<double> x(fe.constraints().m(),
                           fe.constraints().n());

      Assert (q<=Matrices::n_constraint_matrices,
              ExcInternalError());
      x.fill (Matrices::constraint_matrices[q-1]);

      for (unsigned int i=0; i<x.m(); ++i)
        for (unsigned int j=0; j<x.n(); ++j)
          {
            deallog << i << ' ' << j << ' '
                    << x(i,j) << ' '
                    << fe.constraints()(i,j)
                    << std::endl;
            Assert (std::fabs (x(i,j) - fe.constraints()(i,j))
                    <
                    1e-14,
                    ExcInternalError());
          }
    }
  deallog << "OK" << std::endl;
}


int main () 
{
  std::ofstream logfile("fe_q_3d_01.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  check ();
}


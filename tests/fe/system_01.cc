// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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



// previously, the FETools::get_interpolation_matrix function would
// compute its result itself by interpolation. now, the different
// finite elements do that themselves, if they can. make sure the
// result doesn't change, in this case for this function in FESystem
// that does it's job by calling the respective function in the base
// elements and then munging the results

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>

#include <fstream>
#include <string>

#define PRECISION 5



template<int dim>
void
check(const FiniteElement<dim> &fe1,
      const FiniteElement<dim> &fe2)
{
  deallog << fe1.get_name () << " to " << fe2.get_name ()
          << std::endl;

  FullMatrix<float> m (fe2.dofs_per_cell,
                       fe1.dofs_per_cell);
  FETools::get_interpolation_matrix (fe1, fe2, m);

  for (unsigned int i=0; i<m.m(); ++i)
    {
      for (unsigned int j=0; j<m.n(); ++j)
        deallog << m(i,j) << ' ';

      deallog << std::endl;
    }

  deallog << std::endl;
}



#define CHECK_SYS1(sub1_1,N1_1,sub2_1,N2_1,dim) \
  { FESystem<dim> fe1(sub1_1, N1_1);   \
    FESystem<dim> fe2(sub2_1, N2_1);   \
    check(fe1, fe2); \
    check(fe2, fe1); \
  }

#define CHECK_SYS3(sub1,N1,sub2,N2,sub3,N3,dim)   \
  { FESystem<dim> q(sub1, N1, sub2, N2, sub3, N3); \
    check(q, q); }


int
main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  CHECK_SYS1(FE_Q<1>(1),  3,
             FE_Q<1>(2),  3,
             1);
  CHECK_SYS1(FE_DGQ<1>(2),2,
             FE_DGQ<1>(3),2,
             1);

  CHECK_SYS1(FE_Q<2>(1),  3,
             FE_Q<2>(2),  3,
             2);
  CHECK_SYS1(FE_DGQ<2>(2),2,
             FE_DGQ<2>(3),2,
             2);

  CHECK_SYS1(FE_Q<3>(1),  3,
             FE_Q<3>(2),  3,
             3);
  CHECK_SYS1(FE_DGQ<3>(2),2,
             FE_DGQ<3>(3),2,
             3);

  // systems in 1d
  CHECK_SYS3(FE_Q<1>(1),  3,FE_DGQ<1>(3),1,FE_Q<1>(1),3,1);
  CHECK_SYS3(FE_DGQ<1>(2),2,FE_DGQ<1>(2),2,FE_Q<1>(3),3,1);
  CHECK_SYS3(FE_DGQ<1>(3),1,FE_DGQ<1>(3),1,FE_Q<1>(2),3,1);

  // systems in 2d
  CHECK_SYS3(FE_Q<2>(1),  3,FE_DGQ<2>(3),1,FE_Q<2>(1),3,2);
  CHECK_SYS3(FE_DGQ<2>(2),2,FE_DGQ<2>(2),2,FE_Q<2>(3),3,2);

  // systems in 3d
  CHECK_SYS3(FE_Q<3>(1),  3,FE_DGQ<3>(3),1,FE_Q<3>(1),3,3);

  // systems of systems
  CHECK_SYS3((FESystem<2>(FE_Q<2>(1),3)), 3,
             FE_DGQ<2>(3), 1,
             FE_Q<2>(1), 3,
             2);
  CHECK_SYS3(FE_DGQ<2>(3), 1,
             FESystem<2>(FE_DGQ<2>(3),3), 1,
             FESystem<2>(FE_Q<2>(2),3,
                         FE_DGQ<2>(0),1),2,
             2);

  return 0;
}




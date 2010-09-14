//----------------------------  n_dofs.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  n_dofs.cc  ---------------------------


// check that even for wickedly complicated finite elements the
// hp::DoFHandler allocates the same number of degrees of freedom as a
// ::DoFHandler, as long as the same element is used on all cells


#include "../tests.h"
#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_out.h>
#include <grid/tria_iterator.h>
#include <hp/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_dgp.h>
#include <fe/fe_nedelec.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_system.h>

#include <fstream>



template <int dim>
void test (const FiniteElement<dim> &fe)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global (2);
  tria.begin_active()->set_refine_flag ();
  tria.execute_coarsening_and_refinement ();
  tria.refine_global (1);  

  const hp::FECollection<dim> fe_collection (fe);

  hp::DoFHandler<dim> hp_dof_handler(tria);
  hp_dof_handler.distribute_dofs(fe_collection);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs (fe);

  deallog << fe.get_name() << ' ' << dof_handler.n_dofs()
	  << ' ' << hp_dof_handler.n_dofs() << std::endl;

  Assert (dof_handler.n_dofs() == hp_dof_handler.n_dofs(),
	  ExcInternalError());
}



#define CHECK(EL,deg,dim)\
 { FE_ ## EL<dim> EL(deg);   \
   test(EL); }

#define CHECK_SYS1(sub1,N1,dim) \
 { FESystem<dim> q(sub1, N1);   \
   test(q); }

#define CHECK_SYS2(sub1,N1,sub2,N2,dim) \
 { FESystem<dim> q(sub1, N1, sub2, N2); \
   test(q); }

#define CHECK_SYS3(sub1,N1,sub2,N2,sub3,N3,dim)   \
 { FESystem<dim> q(sub1, N1, sub2, N2, sub3, N3); \
   test(q); }


#define CHECK_ALL(EL,deg)\
 { CHECK(EL,deg,1); \
   CHECK(EL,deg,2); \
   CHECK(EL,deg,3); \
 }



int main ()
{
  std::ofstream logfile("n_dofs/output");
  logfile.precision(2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);  

  CHECK_ALL(Q,1);
  CHECK_ALL(Q,2);
  CHECK_ALL(Q,3);

  CHECK_ALL(DGQ,0);
  CHECK_ALL(DGQ,1);
  CHECK_ALL(DGQ,3);

  CHECK_ALL(DGP,0);
  CHECK_ALL(DGP,1);
  CHECK_ALL(DGP,3);

  CHECK(Nedelec, 0, 2);
  CHECK(Nedelec, 0, 3);

  CHECK(RaviartThomas, 0, 2);
  CHECK(RaviartThomas, 1, 2);
  CHECK(RaviartThomas, 2, 2);
  CHECK(RaviartThomas, 0, 3);
  CHECK(RaviartThomas, 1, 3);
  CHECK(RaviartThomas, 2, 3);

  CHECK_SYS1(FE_Q<1>(1),  3,1);
  CHECK_SYS1(FE_DGQ<1>(2),2,1);
  CHECK_SYS1(FE_DGP<1>(3),1,1);

  CHECK_SYS1(FE_Q<2>(1),  3,2);
  CHECK_SYS1(FE_DGQ<2>(2),2,2);
  CHECK_SYS1(FE_DGP<2>(3),1,2);

  CHECK_SYS1(FE_Q<3>(1),  3,3);
  CHECK_SYS1(FE_DGQ<3>(2),2,3);
  CHECK_SYS1(FE_DGP<3>(3),1,3);


  CHECK_SYS2(FE_Q<1>(1),  3,FE_DGQ<1>(2),2,1);
  CHECK_SYS2(FE_DGQ<1>(2),2,FE_DGP<1>(3),1,1);
  CHECK_SYS2(FE_DGP<1>(3),1,FE_DGQ<1>(2),2,1);

  CHECK_SYS2(FE_Q<2>(1),  3,FE_DGQ<2>(2),2,2);
  CHECK_SYS2(FE_DGQ<2>(2),2,FE_DGP<2>(3),1,2);
  CHECK_SYS2(FE_DGP<2>(3),1,FE_DGQ<2>(2),2,2);

  CHECK_SYS3(FE_Q<1>(1),  3,FE_DGP<1>(3),1,FE_Q<1>(1),3,1);
  CHECK_SYS3(FE_DGQ<1>(2),2,FE_DGQ<1>(2),2,FE_Q<1>(3),3,1);
  CHECK_SYS3(FE_DGP<1>(3),1,FE_DGP<1>(3),1,FE_Q<1>(2),3,1);

  CHECK_SYS3(FE_Q<2>(1),  3,FE_DGP<2>(3),1,FE_Q<2>(1),3,2);
  CHECK_SYS3(FE_DGQ<2>(2),2,FE_DGQ<2>(2),2,FE_Q<2>(3),3,2);
  CHECK_SYS3(FE_DGP<2>(3),1,FE_DGP<2>(3),1,FE_Q<2>(2),3,2);

  CHECK_SYS3(FE_DGQ<3>(1),  3,FE_DGP<3>(3),1,FE_Q<3>(1),3,3);

				   // systems of systems  
  CHECK_SYS3((FESystem<2>(FE_Q<2>(1),3)), 3,
	     FE_DGQ<2>(0), 1,
	     FE_Q<2>(1), 3,
	     2);
  CHECK_SYS3(FE_DGQ<2>(3), 1,
	     FESystem<2>(FE_DGQ<2>(0),3), 1,
	     FESystem<2>(FE_Q<2>(2),1,
			 FE_DGQ<2>(0),1),2,
	     2);

				   // systems with Nedelec elements
  CHECK_SYS2 (FE_DGQ<2>(3), 1,
	      FE_Nedelec<2>(0), 2,
	      2);
  CHECK_SYS3(FE_Nedelec<2>(0), 1,
	     FESystem<2>(FE_DGQ<2>(1),2), 1,
	     FESystem<2>(FE_Q<2>(2),1,
			 FE_Nedelec<2>(0),2),2,
	     2);
  
  deallog << "OK" << std::endl;
}

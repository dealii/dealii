//----------------------------  dof_tools_common.cc  ---------------------------
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
//----------------------------  dof_tools_common.cc  ---------------------------


// common framework for the various dof_tools_*.cc tests

#include "../tests.h"
#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_dgp.h>
#include <fe/fe_nedelec.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_system.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>


// forward declaration of the function that must be provided in the
// .cc files
template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler);

// forward declaration of a variable with the name of the output file
extern std::string output_file_name;


void
output_bool_vector (std::vector<bool> &v)
{
for (unsigned int i=0; i<v.size(); ++i)
    deallog << (v[i] ? '1' : '0');
  deallog << std::endl;
}



template <int dim>
void
set_boundary_ids (Triangulation<dim> &tria)
{  
  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
    tria.begin_active()->face(f)->set_boundary_indicator (f);
}


void
set_boundary_ids (Triangulation<1> &)
{}
  


template <int dim>
void
check (const FiniteElement<dim> &fe,
       const std::string        &name)
{
  deallog << "Checking " << name
          << " in " << dim << "d:"
          << std::endl;
  
                                   // create tria and dofhandler
                                   // objects. set different boundary
                                   // and sub-domain ids
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);
  set_boundary_ids (tria);
  tria.refine_global (1);
  for (int i=0; i<2; ++i)
    {
      tria.begin_active()->set_refine_flag();
      tria.execute_coarsening_and_refinement ();
    }
  for (typename Triangulation<dim>::active_cell_iterator
         cell=tria.begin_active();
       cell!=tria.end(); ++cell)
    cell->set_subdomain_id (cell->level());
  DoFHandler<dim> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

                                   // call main function in .cc files
  check_this (dof_handler);
}





#define CHECK(EL,deg,dim)\
 { FE_ ## EL<dim> EL(deg);   \
   check(EL, #EL #deg); }

#define CHECK_SYS1(sub1,N1,dim) \
 { FESystem<dim> q(sub1, N1);   \
   check(q, #sub1 #N1); }

#define CHECK_SYS2(sub1,N1,sub2,N2,dim) \
 { FESystem<dim> q(sub1, N1, sub2, N2); \
   check(q, #sub1 #N1 #sub2 #N2); }

#define CHECK_SYS3(sub1,N1,sub2,N2,sub3,N3,dim)   \
 { FESystem<dim> q(sub1, N1, sub2, N2, sub3, N3); \
   check(q, #sub1 #N1 #sub2 #N2 #sub3 #N3); }


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
      std::ofstream logfile(output_file_name.c_str());
      logfile.precision (2);
      deallog.attach(logfile);
      deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

      CHECK_ALL(Q,1);
      CHECK_ALL(Q,2);
      CHECK_ALL(Q,3);

      CHECK_ALL(DGQ,0);
      CHECK_ALL(DGQ,1);
      CHECK_ALL(DGQ,2);
      CHECK_ALL(DGQ,3);
      CHECK_ALL(DGQ,4);

      CHECK_ALL(DGP,0);
      CHECK_ALL(DGP,1);
      CHECK_ALL(DGP,2);
      CHECK_ALL(DGP,3);
      CHECK_ALL(DGP,4);

      CHECK(Nedelec, 1, 2);
      CHECK(Nedelec, 1, 3);

      CHECK(RaviartThomas, 0, 2);
      CHECK(RaviartThomas, 1, 2);
      CHECK(RaviartThomas, 2, 2);

      CHECK(RaviartThomasNodal, 0, 2);
      CHECK(RaviartThomasNodal, 1, 2);
      CHECK(RaviartThomasNodal, 2, 2);
      CHECK(RaviartThomasNodal, 0, 3);
      CHECK(RaviartThomasNodal, 1, 3);
      CHECK(RaviartThomasNodal, 2, 3);
      
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

      CHECK_SYS2(FE_Q<3>(1)  ,3,FE_DGQ<3>(2),2,3);
      CHECK_SYS2(FE_DGQ<3>(2),2,FE_DGP<3>(3),1,3);
      CHECK_SYS2(FE_DGP<3>(3),1,FE_DGQ<3>(2),2,3);


      CHECK_SYS3(FE_Q<1>(1),  3,FE_DGP<1>(3),1,FE_Q<1>(1),3,1);
      CHECK_SYS3(FE_DGQ<1>(2),2,FE_DGQ<1>(2),2,FE_Q<1>(3),3,1);
      CHECK_SYS3(FE_DGP<1>(3),1,FE_DGP<1>(3),1,FE_Q<1>(2),3,1);

      CHECK_SYS3(FE_Q<2>(1),  3,FE_DGP<2>(3),1,FE_Q<2>(1),3,2);
      CHECK_SYS3(FE_DGQ<2>(2),2,FE_DGQ<2>(2),2,FE_Q<2>(3),3,2);
      CHECK_SYS3(FE_DGP<2>(3),1,FE_DGP<2>(3),1,FE_Q<2>(2),3,2);

      CHECK_SYS3(FE_Q<3>(1),  3,FE_DGP<3>(3),1,FE_Q<3>(1),3,3);
      CHECK_SYS3(FE_DGQ<3>(2),2,FE_DGQ<3>(2),2,FE_Q<3>(3),3,3);
      CHECK_SYS3(FE_DGP<3>(3),1,FE_DGP<3>(3),1,FE_Q<3>(2),3,3);

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

                                       // systems with Nedelec elements
      CHECK_SYS2 (FE_DGQ<2>(3), 1,
                  FE_Nedelec<2>(1), 2,
                  2);
      CHECK_SYS3(FE_Nedelec<2>(1), 1,
                 FESystem<2>(FE_DGQ<2>(3),3), 1,
                 FESystem<2>(FE_Q<2>(2),3,
                             FE_Nedelec<2>(1),2),2,
                 2);
  
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


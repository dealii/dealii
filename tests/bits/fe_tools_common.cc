//----------------------------  fe_tools_common.cc  ---------------------------
//    $Id$
//    Version: 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_tools_common.cc  ---------------------------


// common framework for the various fe_tools_*.cc tests

#include "../tests.h"
#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_constraints.h>
#include <lac/vector.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_dgp.h>
#include <fe/fe_nedelec.h>
#include <fe/fe_system.h>
#include <fe/fe_tools.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <memory>


// forward declaration of the function that must be provided in the
// .cc files
template <int dim>
void
check_this (const FiniteElement<dim> &fe1,
            const FiniteElement<dim> &fe2);

// forward declaration of a variable with the name of the output file
extern std::string output_file_name;

  

// output some indicators for a given matrix. we don't write out the
// entire matrix since this would blow up our output files beyond
// reasonable limits
void
output_matrix (const FullMatrix<double> &m)
{
  deallog << m.l1_norm() << ' ' << m.linfty_norm()
          << std::endl;
  if (m.m() == m.n())
    deallog << m.norm2() << std::endl;
  
  for (unsigned int i=0; i<std::min(m.m(),m.n()); ++i)
    deallog << m(i,i) << ' ' << m(i,std::min(m.m(),m.n())-i-1) << ' ';
  deallog << std::endl;
}


// output some indicators for a given vector
template <class VECTOR>
void
output_vector (const VECTOR &v)
{
  deallog << v.l1_norm() << ' ' << v.l2_norm() << ' ' << v.linfty_norm()
          << std::endl;

                                   // write out at most 20 equispaced
                                   // elements of the vector
  for (unsigned int i=0; i<v.size(); i+=std::max(1U,v.size()/20))
    deallog << v(i) << ' ';
  deallog << std::endl;
}




template <int dim>
Triangulation<dim> * make_tria ()
{
  Triangulation<dim> *tria = new Triangulation<dim>();
  GridGenerator::hyper_cube(*tria, 0., 1.);
  tria->refine_global (1);
  for (int i=0; i<2; ++i)
    {
      tria->begin_active()->set_refine_flag();
      tria->execute_coarsening_and_refinement ();
    }
  return tria;
}



template <int dim>
DoFHandler<dim> * make_dof_handler (const Triangulation<dim> &tria,
                                    const FiniteElement<dim> &fe)
{
  DoFHandler<dim> *dof_handler = new DoFHandler<dim>(tria);
  dof_handler->distribute_dofs (fe);
  return dof_handler;
}



template <int dim>
void
check (const FiniteElement<dim> &fe1,
       const FiniteElement<dim> &fe2,
       const std::string        &name)
{
  deallog << "Checking " << name
          << " in " << dim << "d:"
          << std::endl;
  
                                   // call main function in .cc files
  check_this (fe1, fe2);
}





#define CHECK(EL1,deg1,EL2,deg2,dim)\
 { FE_ ## EL1<dim> fe1(deg1);   \
   FE_ ## EL2<dim> fe2(deg2);   \
   check(fe1, fe2, #EL1 #deg1 " against " #EL2 #deg2); \
   check(fe2, fe1, #EL2 #deg2 " against " #EL1 #deg1); \
 }

#define CHECK_SYS1(sub1_1,N1_1,sub2_1,N2_1,dim) \
 { FESystem<dim> fe1(sub1_1, N1_1);   \
   FESystem<dim> fe2(sub2_1, N2_1);   \
   check(fe1, fe2, #sub1_1 #N1_1 " against " #sub2_1 #N2_1); \
   check(fe2, fe1, #sub2_1 #N2_1 " against " #sub1_1 #N1_1); \
 }

/*
#define CHECK_SYS2(sub1,N1,sub2,N2,dim) \
 { FESystem<dim> q(sub1, N1, sub2, N2); \
   check(q, #sub1 #N1 #sub2 #N2); }

#define CHECK_SYS3(sub1,N1,sub2,N2,sub3,N3,dim)   \
 { FESystem<dim> q(sub1, N1, sub2, N2, sub3, N3); \
   check(q, #sub1 #N1 #sub2 #N2 #sub3 #N3); }

*/
#define CHECK_ALL(EL1,deg1,EL2,deg2)\
 { CHECK(EL1,deg1,EL2,deg2,1); \
   CHECK(EL1,deg1,EL2,deg2,2); \
   CHECK(EL1,deg1,EL2,deg2,3); \
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

      CHECK_ALL(Q,1,Q,1);
      CHECK_ALL(Q,1,Q,2);
      CHECK_ALL(Q,1,Q,3);
      CHECK_ALL(Q,2,Q,2);
      CHECK_ALL(Q,2,Q,3);
      CHECK_ALL(Q,3,Q,3);

      CHECK_ALL(DGQ,0,DGQ,0);
      CHECK_ALL(DGQ,0,DGQ,1);
      CHECK_ALL(DGQ,0,DGQ,2);
      CHECK_ALL(DGQ,0,DGQ,4);
      CHECK_ALL(DGQ,1,DGQ,1);
      CHECK_ALL(DGQ,1,DGQ,3);
      CHECK_ALL(DGQ,1,DGQ,5);
      CHECK_ALL(DGQ,2,DGQ,2);
      CHECK_ALL(DGQ,2,DGQ,2);
      CHECK_ALL(DGQ,2,DGQ,4);
      CHECK_ALL(DGQ,3,DGQ,3);
      CHECK_ALL(DGQ,3,DGQ,5);
      CHECK_ALL(DGQ,4,DGQ,5);

      CHECK_ALL(DGP,0,DGP,0);
      CHECK_ALL(DGP,0,DGP,1);
      CHECK_ALL(DGP,0,DGP,2);
      CHECK_ALL(DGP,0,DGP,4);
      CHECK_ALL(DGP,1,DGP,1);
      CHECK_ALL(DGP,1,DGP,3);
      CHECK_ALL(DGP,1,DGP,5);
      CHECK_ALL(DGP,2,DGP,2);
      CHECK_ALL(DGP,2,DGP,2);
      CHECK_ALL(DGP,2,DGP,4);
      CHECK_ALL(DGP,3,DGP,3);
      CHECK_ALL(DGP,3,DGP,5);
      CHECK_ALL(DGP,4,DGP,5);
      
      CHECK(Nedelec, 1, Nedelec, 1, 2);
      CHECK(Nedelec, 1, Nedelec, 1, 3);
  

      CHECK_SYS1(FE_Q<1>(1),  3,
                 FE_Q<1>(2),  3,
                 1);
      CHECK_SYS1(FE_DGQ<1>(2),2,
                 FE_DGQ<1>(3),2,
                 1);
      CHECK_SYS1(FE_DGP<1>(3),1,
                 FE_DGP<1>(1),1,
                 1);

      CHECK_SYS1(FE_Q<2>(1),  3,
                 FE_Q<2>(2),  3,
                 2);
      CHECK_SYS1(FE_DGQ<2>(2),2,
                 FE_DGQ<2>(3),2,
                 2);
      CHECK_SYS1(FE_DGP<2>(3),1,
                 FE_DGP<2>(1),1,
                 2);

      CHECK_SYS1(FE_Q<3>(1),  3,
                 FE_Q<3>(2),  3,
                 3);
      CHECK_SYS1(FE_DGQ<3>(2),2,
                 FE_DGQ<3>(3),2,
                 3);
      CHECK_SYS1(FE_DGP<3>(3),1,
                 FE_DGP<3>(1),1,
                 3);


/*
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
*/  
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


//----------------------------  dof_tools_common.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_tools_common.cc  ---------------------------


// common framework for the various dof_tools_*.cc tests

#include <base/logstream.h>
#include <lac/vector.h>
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
#include <fe/fe_system.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>


// forward declaration of the function that must be provided in the
// .cc files
template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler,
            const Vector<double>  &v_node,
            const Vector<double>  &v_cell);

// forward declaration of a variable with the name of the output file
extern std::string output_file_name;

  


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
  tria.refine_global (1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement ();

  DoFHandler<dim> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  Vector<double> v_node (dof_handler.n_dofs());
  for (unsigned int i=0; i<v_node.size(); ++i) v_node(i) = i;

  Vector<double> v_cell (dof_handler.get_tria().n_active_cells());
  for (unsigned int i=0; i<v_cell.size(); ++i) v_cell(i) = i;
  
                                   // call main function in .cc files
  check_this (dof_handler, v_node, v_cell);
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
      std::ofstream logfile(output_file_name.c_str());
      logfile.precision (2);
      deallog.attach(logfile);
      deallog.depth_console(0);

      CHECK_ALL(Q,1);
      CHECK_ALL(Q,2);
      CHECK_ALL(Q,3);

      CHECK_ALL(DGQ,0);
      CHECK_ALL(DGQ,1);
      CHECK_ALL(DGQ,2);

      CHECK(Nedelec, 1, 2);
      CHECK(Nedelec, 1, 3);
  
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


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


// common framework for the various dof_tools_*.cc tests

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_system.h>

#include <fstream>
#include <iomanip>
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


// take a vector, and make a block vector out of it
void
make_block_vector (const Vector<double> &in,
                   BlockVector<double>  &out)
{
  std::vector<types::global_dof_index> block_sizes(2);
  block_sizes[0] = in.size() / 2;
  block_sizes[1] = in.size() - block_sizes[0];

  out.reinit (block_sizes);
  std::copy (in.begin(), in.end(), out.begin());
}




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
      deallog << std::setprecision (2);
      logfile << std::setprecision (2);
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      CHECK_ALL(Q,1);
      CHECK_ALL(Q,2);
      CHECK_ALL(Q,3);

      CHECK_ALL(DGQ,0);
      CHECK_ALL(DGQ,1);
      CHECK_ALL(DGQ,2);

      CHECK(Nedelec, 0, 2);
      CHECK(Nedelec, 0, 3);

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


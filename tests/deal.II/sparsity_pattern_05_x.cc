// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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



// check that generating a sparsity pattern without constraints and later
// calling constraints.condense() on it results in the same pattern as when
// creating the condensed pattern right away


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/compressed_set_sparsity_pattern.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <fstream>



bool operator == (const BlockSparsityPattern &sp1,
                  const BlockSparsityPattern &sp2)
{
  if (sp1.n_block_rows() != sp2.n_block_rows())
    return false;

  if (sp1.n_block_cols() != sp2.n_block_cols())
    return false;

  for (unsigned int i=0; i<sp1.n_block_rows(); ++i)
    for (unsigned int j=0; j<sp1.n_block_cols(); ++j)
      if (!(sp1.block(i,j) == sp2.block(i,j)))
        {
          Assert (false, ExcInternalError());
          return false;
        }

  return true;
}




template <int dim>
void
check ()
{
  Triangulation<dim> tr;
  if (dim==2)
    GridGenerator::hyper_ball(tr, Point<dim>(), 1);
  else
    GridGenerator::hyper_cube(tr, -1,1);
  tr.refine_global (1);
  tr.begin_active()->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();
  tr.begin_active(2)->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();
  if (dim==1)
    tr.refine_global(2);

  // create a system element composed
  // of one Q1 and one Q2 element
  FESystem<dim> element(FE_Q<dim>(1), 1,
                        FE_Q<dim>(2), 1);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof, constraints);
  constraints.close ();


//--------------- Regular sparsity pattern checks -----------------

  // first way: directly
  SparsityPattern sparsity_1 (dof.n_dofs(), dof.n_dofs());
  DoFTools::make_sparsity_pattern (dof, sparsity_1);
  constraints.condense (sparsity_1);
  sparsity_1.compress ();

  // second way: via direct elimination of
  // constraints
  SparsityPattern sparsity_2;
  CompressedSetSparsityPattern csp_2 (dof.n_dofs());
  DoFTools::make_sparsity_pattern (dof, csp_2, constraints);
  sparsity_2.copy_from (csp_2);


  // the exact content of sparsity
  // patterns is checked in other
  // tests, so only make sure that
  // sparsity_[12] are equal
  deallog << "Check 1:"
          << " -- "
          << (sparsity_1 == sparsity_2 ? "ok" : "failed")
          << std::endl;



//--------------- Block sparsity pattern checks -----------------

  const unsigned int n  = dof.n_dofs();
  const unsigned int n1 = n/3;
  const unsigned int n2 = n - n1;

  BlockSparsityPattern sparsity_3(2,2);
  sparsity_3.block(0,0).reinit (n1,n1,n);
  sparsity_3.block(1,0).reinit (n2,n1,n);
  sparsity_3.block(0,1).reinit (n1,n2,n);
  sparsity_3.block(1,1).reinit (n2,n2,n);
  sparsity_3.collect_sizes ();

  DoFTools::make_sparsity_pattern (dof, sparsity_3);
  constraints.condense (sparsity_3);
  sparsity_3.compress ();

  BlockSparsityPattern sparsity_4;
  BlockCompressedSetSparsityPattern csp_4(2,2);
  csp_4.block(0,0).reinit (n1,n1);
  csp_4.block(1,0).reinit (n2,n1);
  csp_4.block(0,1).reinit (n1,n2);
  csp_4.block(1,1).reinit (n2,n2);
  csp_4.collect_sizes ();

  DoFTools::make_sparsity_pattern (dof, csp_4, constraints);
  csp_4.compress ();

  sparsity_4.copy_from (csp_4);

  deallog << "Check 2:"
          << " -- "
          << (sparsity_3 == sparsity_4 ? "ok" : "failed")
          << std::endl;
}



int main ()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);

  deallog.push ("1d");
  check<1> ();
  deallog.pop ();
  deallog.push ("2d");
  check<2> ();
  deallog.pop ();
  deallog.push ("3d");
  check<3> ();
  deallog.pop ();
}

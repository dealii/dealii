//----------------------------  sparsity_pattern.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparsity_pattern.cc  ---------------------------


/* Author: Wolfgang Bangerth, University of Heidelberg, 2001 */

// check that the direct generation of the sparsity pattern and that
// via the CompressedSparsityPattern result in the same


#include <base/logstream.h>
#include <lac/sparsity_pattern.h>
#include <lac/compressed_sparsity_pattern.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_constraints.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>

#include <fstream>



bool operator == (const SparsityPattern &sp1,
		  const SparsityPattern &sp2)
{
  if (sp1.n_nonzero_elements() != sp2.n_nonzero_elements())
    return false;
  
  for (unsigned int i=0; i<sp1.n_nonzero_elements(); ++i)
    if (sp1.get_column_numbers()[i] !=
	sp2.get_column_numbers()[i])
      return false;

  for (unsigned int i=0; i<sp1.n_rows(); ++i)
    if (sp1.get_rowstart_indices()[i] !=
	sp2.get_rowstart_indices()[i])
      return false;
  
  return true;
};



template <int dim>
void
check_boundary (const DoFHandler<dim> &dof)
{
  std::vector<unsigned int> dof_to_boundary_mapping;
  DoFTools::map_dof_to_boundary_indices (dof,
					 dof_to_boundary_mapping);

				   // first way: direct generation
  SparsityPattern sparsity_1(dof.n_boundary_dofs(),
			     dof.max_couplings_between_boundary_dofs());
  DoFTools::make_boundary_sparsity_pattern (dof,
					    dof_to_boundary_mapping,
					    sparsity_1);
  sparsity_1.compress ();

				   // second way: via a CompressedSparsityPattern
  SparsityPattern sparsity_2;
  CompressedSparsityPattern csp(dof.n_boundary_dofs());
  DoFTools::make_boundary_sparsity_pattern (dof,
					    dof_to_boundary_mapping,
					    csp);
  sparsity_2.copy_from (csp);

				   // the exact content of sparsity
				   // patterns is checked in other
				   // tests, so only make sure that
				   // sparsity_[12] are equal
  deallog << __PRETTY_FUNCTION__
	  << " -- "
	  << (sparsity_1 == sparsity_2 ? "ok" : "failed")
	  << std::endl;
};




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

				   // create sparsity pattern. note
				   // that different components should
				   // not couple, so use pattern
  std::vector<std::vector<bool> > mask (2, std::vector<bool>(2, false));
  mask[0][0] = mask[1][1] = true;

				   // first way: directly
  SparsityPattern sparsity_1 (dof.n_dofs(), dof.n_dofs());
  DoFTools::make_sparsity_pattern (dof, mask, sparsity_1);
  constraints.condense (sparsity_1);
  sparsity_1.compress ();

				   // second way: via CompressedSparsityPattern
  SparsityPattern sparsity_2;
  CompressedSparsityPattern csp (dof.n_dofs());
  DoFTools::make_sparsity_pattern (dof, mask, csp);
  constraints.condense (csp);
  sparsity_2.copy_from (csp);


				   // the exact content of sparsity
				   // patterns is checked in other
				   // tests, so only make sure that
				   // sparsity_[12] are equal
  deallog << __PRETTY_FUNCTION__
	  << " -- "
	  << (sparsity_1 == sparsity_2 ? "ok" : "failed")
	  << std::endl;
  
  check_boundary (dof);
};



int main ()
{
  std::ofstream logfile ("sparsity_pattern.output");
  logfile.precision (2);
  logfile.setf(std::ios::fixed);  
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

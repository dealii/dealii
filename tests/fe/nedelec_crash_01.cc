//--------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2005, 2006, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//--------------------------------------------------------------------

// Test, whether differential operators produce a cochain complex on
// the standard Hilbert space sequence

#include "../tests.h"
#include "../lib/test_grids.h"

#include <base/logstream.h>

#include <lac/matrix_block.h>
#include <lac/sparse_matrix.h>
#include <lac/block_sparsity_pattern.h>
#include <lac/compressed_sparsity_pattern.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_tools.h>

#include <fe/fe_q.h>
#include <fe/fe_nedelec.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>
#include <fe/mapping_q1.h>

#include <numerics/mesh_worker_info.h>
#include <numerics/mesh_worker_assembler.h>
#include <numerics/mesh_worker_loop.h>


int main()
{
  const std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  struct MyFE : FE_Nedelec<3> {
      MyFE () : FE_Nedelec<3>(0) {};
      virtual bool hp_constraints_are_implemented () const
	{
	  return false;
	}
  } fe;

  Triangulation<3> tr(Triangulation<3>::limit_level_difference_at_vertices);
  TestGrids::hypercube(tr, 2, true);

  DoFHandler<3> dof(tr);
  dof.distribute_dofs(fe);

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof, constraints);
  constraints.close();
}

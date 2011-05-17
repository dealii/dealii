//----------------------------  face_orientation_crash.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  face_orientation_crash.cc  ---------------------------

// trip up the new code handling hanging node constraints with a face in 3d
// that has face_orientation==false



#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/vectors.h>
#include <deal.II/numerics/error_estimator.h>

#include <fstream>


template <int dim>
void
check ()
{
				   // create a mesh with at least one cell
				   // that has a face with
				   // face_orientation==false. refine each of
				   // the 7 cells in turn, to make sure we
				   // have a face with hanging nodes that has
				   // face_orientation==false at least once
  for (unsigned int i=0; i<7; ++i)
    {
      deallog << "Check " << i << std::endl;
      
      Triangulation<dim> tria;
      GridGenerator::hyper_ball(tria);

      typename Triangulation<dim>::active_cell_iterator
	cell = tria.begin_active();
      std::advance(cell,i);
      cell->set_refine_flag();
      tria.execute_coarsening_and_refinement ();

				       // attach a DoFHandler
      FE_Q<dim> element(1);
      DoFHandler<dim> dof(tria);
      dof.distribute_dofs(element);

				       // then build hanging node
				       // constraints. this should trip the
				       // new code using the hp constraints,
				       // added in late July 2006
      ConstraintMatrix constraints;
      DoFTools::make_hanging_node_constraints (dof,
					       constraints);

      for (unsigned int j=0; j<dof.n_dofs(); ++j)
	if (constraints.is_constrained (j))
	  deallog << j << std::endl;
    }
  
  deallog << "OK" << std::endl;
}


int main ()
{
  std::ofstream logfile ("face_orientation_crash/output");
  deallog << std::setprecision (2);
  deallog << std::fixed;  
  deallog.attach(logfile);
  deallog.depth_console (0);

  check<3> ();
}

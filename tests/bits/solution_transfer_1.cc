//----------------------------  solution_transfer_1.cc  ---------------------------
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
//----------------------------  solution_transfer_1.cc  ---------------------------


// something went wrong with the SolutionTransfer class: when we had a
// vector or vectors as input and an empty vector of vectors as
// output, then the output was a zero vector even if the input vector
// was nonzero
//
// this was due to us computing the output vector size before we
// resized the output vector
//
// reported by Brent Bailey (bailey@utias.utoronto.ca) on Wed Oct 22
// 17:02:59 2003

#include <base/logstream.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <numerics/solution_transfer.h>

#include <fstream>
    

int main () 
{
  std::ofstream logfile("solution_transfer_1.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  Triangulation<2> tria;
  GridGenerator::hyper_cube (tria);

  FE_Q<2> fe(1);
  
  DoFHandler<2> dof_handler (tria);
  dof_handler.distribute_dofs (fe);
  
  SolutionTransfer<2, double> soltrans (dof_handler);
  std::vector<Vector<double> > solution_in, solution_out;
  Vector<double> solution (dof_handler.n_dofs());
  solution(0) = 1;
  solution_in.push_back (solution);

  tria.begin_active()->set_refine_flag();
  tria.prepare_coarsening_and_refinement ();
  soltrans.prepare_for_coarsening_and_refinement(solution_in);

  tria.execute_coarsening_and_refinement ();
  dof_handler.distribute_dofs (fe);

  soltrans.interpolate(solution_in, solution_out);

  Assert (solution_out.size() == solution_in.size(),
	  ExcInternalError());
  Assert (solution_out[0].size() == dof_handler.n_dofs(),
	  ExcInternalError());

				   // this is the assertion that used
				   // to fail: the input vector is
				   // non-zero, so the output vector
				   // should not be either!
  Assert (solution_out[0].l2_norm() > 0,
	  ExcInternalError());

  deallog << "OK" << std::endl;
}

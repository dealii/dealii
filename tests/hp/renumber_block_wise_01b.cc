//----------------------------  renumber_block_wise_01a.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2003, 2004, 2009, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  renumber_block_wise_01b.cc  ---------------------------


// A further redux of the _01a test
//
// we distribute Q2 and Q1 dofs on a 2x1x1 mesh in 3d. this leads to
// 27+8-4=31 total dofs (-4 because of the unification of vertex dofs)
// but it turns out that at the time of writing this test not all dof
// were renumbered after unification, leading some to remain at values
// above n_dofs(). tsk tsk tsk...


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/multigrid/mg_dof_accessor.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/hp/fe_collection.h>

#include <fstream>



template <int dim>
void
check ()
{
  Triangulation<dim> tr;
  std::vector<unsigned int> sub(3,1U);
  sub[0] = 2;
  GridGenerator::subdivided_hyper_rectangle(tr, sub, Point<dim>(), Point<dim>(2,1,1));

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back (FE_Q<dim>(2));
  fe_collection.push_back (FE_Q<dim>(1));

  hp::DoFHandler<dim> dof(tr);
  {
    bool coin = false;
    for (typename hp::DoFHandler<dim>::active_cell_iterator cell=dof.begin_active();
	 cell != dof.end(); ++cell)
      {
	cell->set_active_fe_index (coin ? 0 : 1);
	coin = !coin;
      }
  }
  dof.distribute_dofs(fe_collection);

  std::vector<bool> touched(dof.n_dofs(), false);
  std::vector<unsigned int> local_dof_indices;
  for (typename hp::DoFHandler<dim>::active_cell_iterator cell=dof.begin_active();
	 cell != dof.end(); ++cell)
    {
      const unsigned int fe_index = cell->active_fe_index();
      const unsigned int dofs_per_cell = fe_collection[fe_index].dofs_per_cell;
      local_dof_indices.resize (dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);

      deallog << "cell=" << cell << std::endl;
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  Assert (local_dof_indices[i] < dof.n_dofs(),
		  ExcInternalError());
	  touched[local_dof_indices[i]] = true;
	  deallog << local_dof_indices[i] << ' ';
	}
      deallog << std::endl;
    }

				   // make sure all dof indices have
				   // actually been used
  for (unsigned int i=0; i<touched.size(); ++i)
    Assert (touched[i] == true, ExcInternalError());

  deallog << "OK" << std::endl;
}


int main ()
{
  std::ofstream logfile ("renumber_block_wise_01b/output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);

  check<3> ();
}

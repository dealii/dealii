/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 2000 */
/* 
   Purpose: check some things with the intergrid map
*/

#include <base/logstream.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_constraints.h>
#include <dofs/dof_tools.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_handler.h>
#include <grid/intergrid_map.h>
#include <fe/fe_lib.lagrange.h>
#include <fe/fe_lib.dg.h>
#include <fe/fe_system.h>

#include <fstream>


template <int dim>
void check ()
{
  deallog << "Checking in " << dim << " space dimensions"
	  << endl
	  << "---------------------------------------" << endl;
  
				   // create two grids
  Triangulation<dim> tria_1, tria_2;
  GridGenerator::hyper_cube (tria_1, -1, 1);
  tria_1.refine_global (5-dim);
  tria_2.copy_triangulation (tria_1);

				   // create two different fe's
  const FEQ1<dim> fe_linear;
  const FEQ2<dim> fe_quadratic;
  const FEDG_Q2<dim> fe_dq_quadratic;
  const FESystem<dim> fe_1 (fe_dq_quadratic, 3,
			    fe_quadratic, 2,
			    fe_linear, 1);
  const FESystem<dim> fe_2 (fe_quadratic, 1,
			    fe_linear, 2,
			    fe_dq_quadratic, 1);
  
  
				   // make several loops to refine the
				   // two grids
  for (unsigned int i=0; i<3; ++i)
    {
      deallog << "Refinement step " << i << endl;
      
      DoFHandler<dim> dof_1 (tria_1);
      DoFHandler<dim> dof_2 (tria_2);

      dof_1.distribute_dofs (fe_1);
      dof_2.distribute_dofs (fe_2);

				       // now compute intergrid
				       // constraints
      InterGridMap<DoFHandler,dim> intergrid_map;
      intergrid_map.make_mapping (dof_1, dof_2);
      ConstraintMatrix intergrid_constraints;
      DoFTools::compute_intergrid_constraints (dof_1, 4, dof_2, 0, 
					       intergrid_map,
					       intergrid_constraints);
      


				       // now refine grids a little,
				       // but differently. this
				       // produces quite random grids
      DoFHandler<dim>::cell_iterator cell, endc;
      cell = dof_1.begin();
      endc = dof_1.end();
      for (unsigned int index=0; cell!=endc; ++cell)
	if (cell->active())
	  {
	    ++index;
	    if (index % 3 == 0)
	      cell->set_refine_flag ();
	  };

      cell = dof_2.begin();
      endc = dof_2.end();
      for (unsigned int index=0; cell!=endc; ++cell)
	if (cell->active())
	  {
	    ++index;
	    if (index % 3 == 1)
	      cell->set_refine_flag ();
	  };
      
      tria_1.execute_coarsening_and_refinement ();
      tria_2.execute_coarsening_and_refinement ();
    };
};



int main ()
{
  ofstream logfile("intergrid_map.output");
  logfile.precision(4);
  
  deallog.attach(logfile);
  deallog.depth_console(0);

  check<1> ();
  check<2> ();
  check<3> ();
};


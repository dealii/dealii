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
#include <numerics/dof_renumbering.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_handler.h>
#include <grid/intergrid_map.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/mapping_q1.h>
#include <fe/fe_system.h>

#include <fstream>



ofstream logfile("intergrid_constraints.output");


template <int dim>
void check ()
{
  deallog << "Checking in " << dim << " space dimensions"
	  << endl
	  << "---------------------------------------" << endl;
  
				   // create two grids
  Triangulation<dim> tria_1, tria_2;
  GridGenerator::hyper_cube (tria_1, -1, 1);
  tria_1.refine_global (4-dim);
  tria_2.copy_triangulation (tria_1);
  
				   // create two different fe's:
				   // create two really perverse
				   // finite elements to check for
				   // obscure effects when mixing
				   // different FEs and creating
				   // constraints in between them
				   //
				   // note that presently dq2 elements
				   // are not implemented in 3d, so
				   // take special measures
  const FE_DGQ<dim> fe_constant(0);
  const FE_Q<dim>    fe_quadratic(2);
  const FE_DGQ<dim> fe_dq_linear(1);
  const FE_DGQ<dim> *fe_dq_quadratic = (dim != 3
					 ?
					 new FE_DGQ<dim>(2)
					 :
					 0);

  const FESystem<dim>
    *fe_1 = (dim != 3
	     ?
	     new FESystem<dim>(fe_quadratic,     4,
			       *fe_dq_quadratic, 2,
			       fe_constant,      12)
	     :
	     new FESystem<dim>(fe_quadratic,     4,
			       fe_dq_linear,     2,
			       fe_constant,      12)
	     );

  const FESystem<dim>
    *fe_2 = (dim != 3
	     ?
	     new FESystem<dim>(fe_constant,      1,
			       *fe_dq_quadratic, 2,
			       fe_quadratic,     5)
	     :
	     new FESystem<dim>(fe_constant,      1,
			       fe_dq_linear,     2,
			       fe_quadratic,     5)
	     );
  
				   // make several loops to refine the
				   // two grids
  for (unsigned int i=0; i<3; ++i)
    {
      deallog << "Refinement step " << i << endl;

      DoFHandler<dim> dof_1 (tria_1);
      DoFHandler<dim> dof_2 (tria_2);

      dof_1.distribute_dofs (*fe_1);
      dof_2.distribute_dofs (*fe_2);
      DoFRenumbering::Cuthill_McKee (dof_1);
      DoFRenumbering::Cuthill_McKee (dof_2);

      deallog << "  Grid 1: " << tria_1.n_active_cells() << " cells, "
	      << dof_1.n_dofs() << " dofs" << endl;
      deallog << "  Grid 2: " << tria_2.n_active_cells() << " cells, "
	      << dof_2.n_dofs() << " dofs" << endl;
      
				       // now compute intergrid
				       // constraints
      InterGridMap<DoFHandler,dim> intergrid_map;
      intergrid_map.make_mapping (dof_1, dof_2);
      ConstraintMatrix intergrid_constraints;
				       // dq quadratic
      DoFTools::compute_intergrid_constraints (dof_1, 5, dof_2, 2, 
					       intergrid_map,
					       intergrid_constraints);
				       // dq constant
      DoFTools::compute_intergrid_constraints (dof_1, 8, dof_2, 0, 
					       intergrid_map,
					       intergrid_constraints);
				       // continuous
				       // quadratic. continuous
				       // elements only work for 1d at
				       // present!
      if (dim == 1)
	DoFTools::compute_intergrid_constraints (dof_1, 3, dof_2, 5, 
						 intergrid_map,
						 intergrid_constraints);
      
      intergrid_constraints.print (logfile);
 
      

				       // now refine grids a little,
				       // but make sure that grid 2 is
				       // always more refined than
				       // grid 1
      typename DoFHandler<dim>::cell_iterator cell, endc;
      cell = dof_1.begin();
      endc = dof_1.end();
      for (unsigned int index=0; cell!=endc; ++cell)
	if (cell->active())
	  {
	    ++index;
	    if (index % 3 == 0) 
	      {
		cell->set_refine_flag ();

						 // make sure that the
						 // corrsponding cell
						 // on grid 2 is also
						 // refined if that
						 // has not yet
						 // happened
		typename DoFHandler<dim>::cell_iterator cell2 (&tria_2,
							       cell->level(),
							       cell->index(),
							       &dof_2);
		if (!cell2->has_children())
		  cell2->set_refine_flag();
	      };
	  };

      tria_1.execute_coarsening_and_refinement ();
      tria_2.execute_coarsening_and_refinement ();

				       // next refine grid 2 a little more
      cell = dof_2.begin();
      endc = dof_2.end();
      for (unsigned int index=0; cell!=endc; ++cell)
	if (cell->active())
	  {
	    ++index;
	    if (index % 3 == 1)
	      cell->set_refine_flag ();
	  };
      
      tria_2.execute_coarsening_and_refinement ();

				       // use the following instead of
				       // the previous block when
				       // testing continuous elements
				       // in 2d
//      tria_2.refine_global(1);
    };

  delete fe_1;
  delete fe_2;
  if (fe_dq_quadratic != 0)
    delete fe_dq_quadratic;
};



int main ()
{
  logfile.precision(4);
  
  deallog.attach(logfile);
  deallog.depth_console(0);

  check<1> ();
  check<2> ();
  check<3> ();
};


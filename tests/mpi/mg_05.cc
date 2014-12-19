// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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



// simplify test to show bug: missing level dofs:
/*
An error occurred in line <126> of file <mg_05.cc> in function
    void test() [with int dim = 2]
The violated condition was: 
    cell->neighbor(f)->level_subdomain_id()<100
The name and call sequence of the exception was:
    ExcInternalError()
Additional Information: 
(none)
*/

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/lac/trilinos_vector.h>

#include <fstream>

template<int dim>
void output(parallel::distributed::Triangulation<dim> &tr)
{
  const std::string filename = ("mesh." +
                                Utilities::int_to_string
                                (tr.locally_owned_subdomain(), 4) +
                                ".svg");
  std::ofstream stream(filename.c_str());
  /*
  GridOutFlags::XFig flags;
  flags.color_by = GridOutFlags::XFig::level_subdomain_id;
  GridOut out;
  out.set_flags(flags);

  out.write_xfig(tr, stream);
  */
  GridOut grid_out;
  GridOutFlags::Svg svg_flags;
  svg_flags.coloring = GridOutFlags::Svg::subdomain_id;
  svg_flags.label_material_id = false;
  svg_flags.label_level_number = false;
  svg_flags.label_cell_index = false;
  svg_flags.label_subdomain_id = false;
  svg_flags.label_level_subdomain_id = true;
  svg_flags.background = GridOutFlags::Svg::transparent;
  svg_flags.polar_angle = 60;
  svg_flags.convert_level_number_to_height = true;
  grid_out.set_flags(svg_flags);
  grid_out.write_svg(tr, stream);


}

template<int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD,
                                               Triangulation<dim>::none,
                                               parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(3);

  for (unsigned int ii=0;ii<5;++ii)
    {
      typename Triangulation<dim>::active_cell_iterator
        cell = tr.begin_active(),
        endc = tr.end();
      
      for (; cell!=endc; ++cell)
	if (Testing::rand()%42==1)
	  cell->set_refine_flag ();
      
      tr.execute_coarsening_and_refinement ();
      
      DoFHandler<dim> dofh(tr);
      
				       //output(tr);
      
      static const FE_Q<dim> fe(1);
      dofh.distribute_dofs (fe);
      dofh.distribute_mg_dofs (fe);
      
      {
	for (unsigned int lvl=0; lvl<tr.n_levels(); ++lvl)
	  {
//	    deallog << "level " << lvl << ": "<< std::endl;
	    typename DoFHandler<dim>::cell_iterator
	      cell = dofh.begin(lvl),
	      endc = dofh.end(lvl);

	    for (; cell!=endc; ++cell)
	      {
		if (cell->level_subdomain_id()!=tr.locally_owned_subdomain())
		  continue;
		for (unsigned int f=0;f<GeometryInfo<dim>::faces_per_cell;++f)
		  {
		    if (cell->at_boundary(f))
		      continue;

//		    deallog << cell->neighbor(f)->level_subdomain_id() << std::endl;
						     // is cell level-artificial?
		    Assert(cell->neighbor(f)->level_subdomain_id()<100, ExcInternalError());
		
		    std::vector<types::global_dof_index> dofs(fe.n_dofs_per_cell());
		    cell->neighbor(f)->get_mg_dof_indices(dofs);
		    for (unsigned int i=0;i<fe.n_dofs_per_cell();++i)
		      {
			Assert(dofs[i]!=numbers::invalid_dof_index, ExcInternalError());
		      }
		  
		  }
	      }
	  }
      }
    }
  
  if (myid==0)
    deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll log;

  deallog.push("2d");
  test<2>();
  deallog.pop();
}

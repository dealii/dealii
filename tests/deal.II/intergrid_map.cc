// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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

/*
/* Author: Wolfgang Bangerth, University of Heidelberg, 2000
/*
 * Purpose: check some things with the intergrid map
 */

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <fstream>


template <int dim>
void check ()
{
  deallog << "Checking in " << dim << " space dimensions"
          << std::endl
          << "---------------------------------------" << std::endl;

  // create two grids
  Triangulation<dim> tria_1, tria_2;
  GridGenerator::hyper_cube (tria_1, -1, 1);
  tria_1.refine_global (5-dim);
  tria_2.copy_triangulation (tria_1);

  FE_Q<dim> fe_1(1);
  FE_Q<dim> fe_2(2);

  // make several loops to refine the
  // two grids
  for (unsigned int i=0; i<3; ++i)
    {
      deallog << "Refinement step " << i << std::endl;

      DoFHandler<dim> dof_1 (tria_1);
      DoFHandler<dim> dof_2 (tria_2);

      dof_1.distribute_dofs (fe_1);
      dof_2.distribute_dofs (fe_2);

      // create some mapping
      InterGridMap<DoFHandler<dim> > intergrid_map_1;
      InterGridMap<DoFHandler<dim> > intergrid_map_2;
      intergrid_map_1.make_mapping (dof_1, dof_2);
      intergrid_map_2.make_mapping (dof_2, dof_1);

      // write out the mapping
      typename DoFHandler<dim>::cell_iterator cell=dof_1.begin(),
                                              endc=dof_1.end();
      for (; cell!=endc; ++cell)
        {
          deallog << cell
                  << "->"
                  << intergrid_map_1[cell]
                  << "->"
                  << intergrid_map_2[intergrid_map_1[cell]]
                  << std::endl;
// note that not necessarily intergrid_map_2[intergrid_map_1[cell]] ==
// cell, since the meshes have different refinement steps.
        };



      // now refine grids a little,
      // but differently. this
      // produces quite random grids
      cell = dof_1.begin();
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
}



int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(2);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<1> ();
  check<2> ();
  check<3> ();
}


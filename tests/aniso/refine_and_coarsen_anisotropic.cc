// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2017 by the deal.II authors
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


// simply perform some random refinement and coarsening operations. if we do
// this for several cycles chances are high, that any error in the refinement
// and coarsening procedure will show up in an Assert at one place or
// another. This is especially true because of the totally random refinement and
// coarsening flags which also produce unrealistic situations which were not
// correctly accounted for in the implementation phase...

#include "../tests.h"
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/timer.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_generator.h>

#include <iostream>
#include <sstream>

// flag to decide, whether
// 1) a single output_file should
// be generated for automated tests or
// 2) multiple eps-files and output to
// std::cout in order to check
// changes in the library
bool single_file=true;
// maximum number of refinement steps and cells
// to be created during refinement. whichever
// is reached first terminates the loop
const unsigned int max_cycle = 200;
const unsigned int max_cells = 50000;


template <int dim>
void test_isotropic (int type, std::ostream *logfile)
{
  const RefinementCase<dim> ref_cases[7] =
  {
    RefinementCase<dim>::cut_x,
    RefinementCase<dim>(dim > 1 ? RefinementCase<2>::cut_y : RefinementCase<2>::no_refinement),
    RefinementCase<dim>(dim > 1 ? RefinementCase<2>::cut_xy : RefinementCase<2>::no_refinement),
    RefinementCase<dim>(dim > 2 ? RefinementCase<3>::cut_z : RefinementCase<3>::no_refinement),
    RefinementCase<dim>(dim > 2 ? RefinementCase<3>::cut_xz : RefinementCase<3>::no_refinement),
    RefinementCase<dim>(dim > 2 ? RefinementCase<3>::cut_yz : RefinementCase<3>::no_refinement),
    RefinementCase<dim>(dim > 2 ? RefinementCase<3>::cut_xyz : RefinementCase<3>::no_refinement)
  };

  Triangulation<dim> tria(Triangulation<dim>::allow_anisotropic_smoothing);
  if (type == 0)
    GridGenerator::hyper_cube (tria);
  else
    GridGenerator::hyper_ball (tria);

  tria.reset_manifold(0);
  tria.refine_global(1);

  *logfile << "cycle: 0, number of cells: "<<tria.n_cells()
           << std::endl;

  typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active(),
                                                    endc = tria.end();
  for (unsigned int cycle=1; cycle<max_cycle+1; ++cycle)
    {
      cell = tria.begin_active();
      for (; cell!=endc; ++cell)
        if (Testing::rand()%5==0)
          {
            if (Testing::rand()%2==0)
              cell->set_refine_flag(ref_cases[Testing::rand()%
                                              RefinementCase<dim>::isotropic_refinement]);
          }
        else
          cell->set_coarsen_flag();

      tria.execute_coarsening_and_refinement();

      if (!single_file)
        {
          // graphical output
          GridOut grid_out;
          GridOutFlags::Eps<2> eps2(GridOutFlags::EpsFlagsBase::width,
                                    300, .01, false, 5, false);
          grid_out.set_flags (eps2);
          std::ostringstream filename;
          filename << "grid_" << type << "_" << cycle << ".eps";
          std::ofstream outfile(filename.str().c_str());
          grid_out.write_eps (tria, outfile);
        }

      *logfile << "cycle: "<<cycle
               <<", number of cells: "<<tria.n_cells()
               << std::endl;

      if (tria.n_cells() > max_cells || cycle == max_cycle)
        break;
    }
}



int main ()
{
  std::ostream *logfile;

  if (single_file)
    logfile = new std::ofstream("output");
  else
    logfile = &std::cout;

  *logfile<<std::endl
          <<"         2D"<<std::endl
          <<"-----------------------"<<std::endl;

  *logfile<<"HyperCube:"<<std::endl;
  test_isotropic<2> (0,logfile);

  *logfile<<"HyperBall:"<<std::endl;
  test_isotropic<2> (1,logfile);

  *logfile<<std::endl;
  *logfile<<"         3D"<<std::endl
          <<"-----------------------"<<std::endl;

  *logfile<<"HyperCube:"<<std::endl;
  test_isotropic<3> (0,logfile);

  *logfile<<"HyperBall:"<<std::endl;
  test_isotropic<3> (1,logfile);

  // clear whatever is in the buffer
  *logfile << std::flush;
  if (single_file)
    delete logfile;
}


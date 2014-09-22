// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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



// create and print a bunch of ConstrainMatrices

#include "../tests.h"
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/base/logstream.h>

#include <fstream>
#include <cmath>
#include <cstdlib>


std::ofstream logfile("output");


void make_tria (Triangulation<3> &tria, int step)
{
  switch (step)
    {
    case 0:
    case 1:
    {
      // two cells packed behind each
      // other. if step==0, refine back one,
      // otherwise the one in front
      const Point<3> vertices[12] = { Point<3>(0,0,0),
                                      Point<3>(1,0,0),
                                      Point<3>(1,0,1),
                                      Point<3>(0,0,1),

                                      Point<3>(0,1,0),
                                      Point<3>(1,1,0),
                                      Point<3>(1,1,1),
                                      Point<3>(0,1,1),

                                      Point<3>(0,2,0),
                                      Point<3>(1,2,0),
                                      Point<3>(1,2,1),
                                      Point<3>(0,2,1)
                                    };
      const int cell_vertices[2][8] = { { 0,1,2,3,4,5,6,7 }, { 4,5,6,7,8,9,10,11 } };
      std::vector<CellData<3> > cells (2, CellData<3>());
      for (unsigned int cell=0; cell<2; ++cell)
        for (unsigned int j=0; j<8; ++j)
          cells[cell].vertices[j] = cell_vertices[cell][j];
      cells[0].material_id = 0;
      cells[1].material_id = 0;

      tria.create_triangulation_compatibility (
        std::vector<Point<3> >(&vertices[0], &vertices[12]),
        cells,
        SubCellData());       // no boundary information

      if (step==0)
        tria.last_active()->set_refine_flag();
      else
        tria.begin_active()->set_refine_flag();
      tria.execute_coarsening_and_refinement ();

      break;
    };

    case 2:
    case 3:
    {
      // two cells packed next to each
      // other. if step==2, refine right one,
      // otherwise the left one
      const Point<3> vertices[12] = { Point<3>(0,0,0),
                                      Point<3>(1,0,0),
                                      Point<3>(1,0,1),
                                      Point<3>(0,0,1),

                                      Point<3>(0,1,0),
                                      Point<3>(1,1,0),
                                      Point<3>(1,1,1),
                                      Point<3>(0,1,1),

                                      Point<3>(2,0,0),
                                      Point<3>(2,0,1),
                                      Point<3>(2,1,0),
                                      Point<3>(2,1,1)
                                    };
      const int cell_vertices[2][8] = { { 0,1,2,3,4,5,6,7 }, { 1,8,9,2,5,10,11,6 } };
      std::vector<CellData<3> > cells (2, CellData<3>());
      for (unsigned int cell=0; cell<2; ++cell)
        for (unsigned int j=0; j<8; ++j)
          cells[cell].vertices[j] = cell_vertices[cell][j];
      cells[0].material_id = 0;
      cells[1].material_id = 0;

      tria.create_triangulation_compatibility (
        std::vector<Point<3> >(&vertices[0], &vertices[12]),
        cells,
        SubCellData());       // no boundary information

      if (step==2)
        tria.last_active()->set_refine_flag();
      else
        tria.begin_active()->set_refine_flag();
      tria.execute_coarsening_and_refinement ();

      break;
    };

    case 4:
    case 5:
    {
      // two cells packed on top of each
      // other. if step==4, refine top one,
      // otherwise the bottom one
      const Point<3> vertices[12] = { Point<3>(0,0,0),
                                      Point<3>(1,0,0),
                                      Point<3>(1,0,1),
                                      Point<3>(0,0,1),

                                      Point<3>(0,1,0),
                                      Point<3>(1,1,0),
                                      Point<3>(1,1,1),
                                      Point<3>(0,1,1),

                                      Point<3>(1,0,2),
                                      Point<3>(0,0,2),
                                      Point<3>(1,1,2),
                                      Point<3>(0,1,2)
                                    };
      const int cell_vertices[2][8] = { { 0,1,2,3,4,5,6,7 }, { 3, 2, 8, 9 , 7, 6, 10, 11} };
      std::vector<CellData<3> > cells (2, CellData<3>());
      for (unsigned int cell=0; cell<2; ++cell)
        for (unsigned int j=0; j<8; ++j)
          cells[cell].vertices[j] = cell_vertices[cell][j];
      cells[0].material_id = 0;
      cells[1].material_id = 0;

      tria.create_triangulation_compatibility (
        std::vector<Point<3> >(&vertices[0], &vertices[12]),
        cells,
        SubCellData());       // no boundary information

      if (step==4)
        tria.last_active()->set_refine_flag();
      else
        tria.begin_active()->set_refine_flag();
      tria.execute_coarsening_and_refinement ();

      break;
    };


    case 6:
    case 7:
    case 8:
    {
      // four cells, with several refined
      // (see below)
      const Point<3> vertices[18] = { Point<3>(0,0,0),
                                      Point<3>(1,0,0),
                                      Point<3>(1,0,1),
                                      Point<3>(0,0,1),

                                      Point<3>(0,1,0),
                                      Point<3>(1,1,0),
                                      Point<3>(1,1,1),
                                      Point<3>(0,1,1),

                                      Point<3>(2,0,0),
                                      Point<3>(2,0,1),
                                      Point<3>(2,1,0),
                                      Point<3>(2,1,1),

                                      Point<3>(0,2,0),
                                      Point<3>(1,2,0),
                                      Point<3>(1,2,1),
                                      Point<3>(0,2,1),

                                      Point<3>(2,2,0),
                                      Point<3>(2,2,1)
                                    };

      const int cell_vertices[4][8] = { { 0,1,2,3,4,5,6,7 },
        { 1,8,9,2,5,10,11,6 },
        { 4,5,6,7,12,13,14,15},
        { 5,10,11,6,13,16,17,14}
      };
      std::vector<CellData<3> > cells (4, CellData<3>());
      for (unsigned int cell=0; cell<4; ++cell)
        for (unsigned int j=0; j<8; ++j)
          cells[cell].vertices[j] = cell_vertices[cell][j];
      cells[0].material_id = 0;
      cells[1].material_id = 0;
      cells[2].material_id = 0;
      cells[3].material_id = 0;

      tria.create_triangulation_compatibility (
        std::vector<Point<3> >(&vertices[0], &vertices[18]),
        cells,
        SubCellData());       // no boundary information

      switch (step)
        {
        case 6:
          tria.begin_active()->set_refine_flag ();
          break;

        case 7:
          tria.begin_active()->set_refine_flag ();
          (++tria.begin_active())->set_refine_flag ();
          break;
        case 8:
          tria.begin_active()->set_refine_flag ();
          (++(++(++tria.begin_active())))->set_refine_flag ();
          break;
        };

      tria.execute_coarsening_and_refinement ();

      break;
    };
    };
}




int main ()
{
  deallog << std::setprecision (2);
  logfile << std::setprecision (2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  FiniteElement<3> *fe = 0;

  for (unsigned int element=0; element<2; ++element)
    {
      switch (element)
        {
        case 0:
          fe = new FE_Q<3>(1);
          break;
        case 1:
          fe = new FE_Q<3>(2);
          break;
        };

      for (int step=0; step<9; ++step)
        {
          deallog << "Element=" << element << ", Step=" << step << std::endl;

          Triangulation<3> tria;
          make_tria (tria, step);
          GridOut().write_gnuplot (tria, logfile);

          DoFHandler<3> dof (tria);
          dof.distribute_dofs (*fe);

          ConstraintMatrix constraints;
          DoFTools::make_hanging_node_constraints (dof, constraints);
          constraints.close ();

          constraints.print (logfile);

          // release fe
          dof.clear ();

          deallog << std::endl;
        };

      delete fe;
    };
}


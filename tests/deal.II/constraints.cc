//----------------------------  constraints.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  constraints.cc  ---------------------------


#include <dofs/dof_handler.h>
#include <grid/tria.h>
#include <fe/mapping_q1.h>
#include <fe/fe_q.h>
#include <grid/tria_boundary.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <lac/sparse_matrix.h>
#include <base/parameter_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_constraints.h>
#include <dofs/dof_tools.h>
#include <grid/grid_out.h>
#include <base/logstream.h>

#include <fstream>
#include <cmath>
#include <cstdlib>


std::ofstream logfile("constraints.output");


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
					Point<3>(0,2,1)    };
	const int cell_vertices[2][8] = { { 0,1,2,3,4,5,6,7 }, { 4,5,6,7,8,9,10,11 } };
	std::vector<CellData<3> > cells (2, CellData<3>());
	for (unsigned int cell=0; cell<2; ++cell)
	  for (unsigned int j=0; j<8; ++j)
	    cells[cell].vertices[j] = cell_vertices[cell][j];
	cells[0].material_id = 0;
	cells[1].material_id = 0;
  
	tria.create_triangulation (std::vector<Point<3> >(&vertices[0], &vertices[12]),
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
					Point<3>(2,1,1)    };
	const int cell_vertices[2][8] = { { 0,1,2,3,4,5,6,7 }, { 1,8,9,2,5,10,11,6 } };
	std::vector<CellData<3> > cells (2, CellData<3>());
	for (unsigned int cell=0; cell<2; ++cell)
	  for (unsigned int j=0; j<8; ++j)
	    cells[cell].vertices[j] = cell_vertices[cell][j];
	cells[0].material_id = 0;
	cells[1].material_id = 0;
  
	tria.create_triangulation (std::vector<Point<3> >(&vertices[0], &vertices[12]),
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
					Point<3>(0,1,2)    };
	const int cell_vertices[2][8] = { { 0,1,2,3,4,5,6,7 }, { 3, 2, 8, 9 , 7, 6, 10, 11} };
	std::vector<CellData<3> > cells (2, CellData<3>());
	for (unsigned int cell=0; cell<2; ++cell)
	  for (unsigned int j=0; j<8; ++j)
	    cells[cell].vertices[j] = cell_vertices[cell][j];
	cells[0].material_id = 0;
	cells[1].material_id = 0;
	
	tria.create_triangulation (std::vector<Point<3> >(&vertices[0], &vertices[12]),
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
					Point<3>(2,2,1)  };

	const int cell_vertices[4][8] = { { 0,1,2,3,4,5,6,7 },
					  { 1,8,9,2,5,10,11,6 },
					  { 4,5,6,7,12,13,14,15},
					  { 5,10,11,6,13,16,17,14} };
	std::vector<CellData<3> > cells (4, CellData<3>());
	for (unsigned int cell=0; cell<4; ++cell)
	  for (unsigned int j=0; j<8; ++j)
	    cells[cell].vertices[j] = cell_vertices[cell][j];
	cells[0].material_id = 0;
	cells[1].material_id = 0;
	cells[2].material_id = 0;
	cells[3].material_id = 0;
  
	tria.create_triangulation (std::vector<Point<3> >(&vertices[0], &vertices[18]),
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
};



void merge_check ()
{
  deallog << "Checking ConstraintMatrix::merge" << std::endl;

				   // check twice, once with closed
				   // objects, once with open ones
  for (unsigned int run=0; run<2; ++run)
    {
      deallog << "Checking with " << (run == 0 ? "open" : "closed")
	      << " objects" << std::endl;
      
				       // check that the `merge' function
				       // works correctly
      ConstraintMatrix c1, c2;

				       // enter simple line
      c1.add_line (0);
      c1.add_entry (0, 11, 1.);
				       // add more complex line
      c1.add_line (1);
      c1.add_entry (1, 3, 0.5);
      c1.add_entry (1, 4, 0.5);
  
				       // fill second constraints
				       // object with one trivial line
				       // and one which further
				       // constrains one of the
				       // entries in the first object
      c2.add_line (10);
      c2.add_entry (10, 11, 1.);
      c2.add_line (3);
      c2.add_entry (3, 12, 0.25);
      c2.add_entry (3, 13, 0.75);

				       // in one of the two runs,
				       // close the objects
      if (run == 1)
	{
	  c1.close ();
	  c2.close ();
	};

				       // now merge the two and print the
				       // results
      c1.merge (c2);
      c1.print (logfile);
    };
};


int main ()
{
  logfile.precision (2);
  deallog.attach(logfile);
  deallog.depth_console(0);

  FiniteElement<3> *fe;
  
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

  merge_check ();
};


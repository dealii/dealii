// check TriaCellAccessor::neighbor_child_on_subface

#include "../tests.h"
#include <base/quadrature_lib.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_reordering.h>
#include <grid/grid_in.h>
#include <grid/grid_out.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_q.h>

#include <fstream>
#include <iostream>

using namespace std;

void check_this (Triangulation<3> &tria)
{
  FE_Q<3> fe(1);
  
  DoFHandler<3> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  const unsigned int dim=3;
  DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (unsigned int cell_no=0; cell!=endc; ++cell, ++cell_no)
    for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
      {
	const DoFHandler<dim>::face_iterator
	  face=cell->face(face_no);
	    
	if (!face->at_boundary())
	  {
	    Assert (cell->neighbor(face_no).state() == IteratorState::valid,
		    ExcInternalError());
	    DoFHandler<dim>::cell_iterator neighbor=
	      cell->neighbor(face_no);

	    if (face->has_children())
	      {
		const unsigned int neighbor2=cell->neighbor_of_neighbor(face_no);
	    
		deallog << "cell_no=" << cell_no << ", face_no=" << face_no << endl;
		deallog << "  cell->face_orientation(face_no)="
			<< cell->face_orientation(face_no) << endl;
		deallog << "  neighbor2=" << neighbor2 << endl;
		deallog << "  neighbor->face_orientation(neighbor2)="
			<< neighbor->face_orientation(neighbor2) << endl;
		
		
		for (unsigned int subface_no=0;
		     subface_no<GeometryInfo<dim>::subfaces_per_face;
		     ++subface_no)
		  {
		    deallog << "subface_no=" << subface_no << endl;
		    DoFHandler<dim>::cell_iterator neighbor_child=
		      cell->neighbor_child_on_subface(face_no, subface_no);
		  }
	      }
	  }
      }
}



void check (Triangulation<3> &tria)
{
  tria.begin_active()->set_refine_flag ();
  tria.execute_coarsening_and_refinement ();
  
  if (false)
    {
      GridOut grid_out;
      ofstream grid_stream("grid.gnuplot");
      grid_out.write_gnuplot(tria, grid_stream);
    }
  
  check_this (tria);
}

int main () 
{
  ofstream logfile("mesh_3d_17.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  Triangulation<3> coarse_grid;
  GridIn<3> in;
  in.attach_triangulation(coarse_grid);
  ifstream ucd_file("two_cubes.inp");
  in.read_ucd(ucd_file);
  ucd_file.close();
  
  check (coarse_grid);
}


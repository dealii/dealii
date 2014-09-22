// ---------------------------------------------------------------------
//
// Copyright (C) 2014 by the deal.II authors
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

// This script reproduces an issue with merging a hyper ball within a 
// hyper shell using dealII-8.0.0.
//
// Test by Bamdad Hosseini <bamdad.hosseini@gmail.com> (see mail to mailing
// list on 6/20/2014)


#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>

#include <fstream>
#include <iostream>
#include <cmath>
#include <deal.II/base/logstream.h>

std::ofstream logfile("output");


template <int dim>
void flatten_triangulation( Triangulation<dim> &tria_in, Triangulation<dim> &tria_out )
// takes a possibly refined triangulation and returns a new coarse triangulation
// created from the vertices of the input
{

  // we start by extracting the number of vertices and cell is 
  // the original triangulation

  unsigned int n_vertices = tria_in.n_vertices();
  unsigned int n_cells = tria_in.n_cells();
  unsigned int n_active_cells = tria_in.n_active_cells();
  
  // also get all vertices 
  const std::vector< Point<dim> > vertices = tria_in.get_vertices();

  // now we need to loop over all cells and extract the cell data.
  typename Triangulation<dim>::active_cell_iterator
    cell=tria_in.begin_active(),
    endc=tria_in.end();

  std::vector< CellData<dim> > cells(n_active_cells, CellData<dim>());

  unsigned int cell_counter = 0;
  for(; cell!=endc; ++cell)
    {
      for (unsigned int j=0; 
	   j < GeometryInfo<dim>::vertices_per_cell;
	   ++j) {
	cells[cell_counter].vertices[j] = cell->vertex_index(j);
      }
      cells[cell_counter].material_id=0;
      ++cell_counter;
    }

  // pass the info to create_triangulation to create the flat mesh
  tria_out.create_triangulation(vertices, cells, SubCellData() );
}

// this is a function to provide some basic information about the 
// triangulations. The code inside can be ignored. It produces an 
// eps image of the mesh so we can see the nodes are at the right 
// place.

template<int dim>
void mesh_info(const Triangulation<dim> &tria, 
	       const std::string        &filename)
{
  deallog << "Mesh info for " << filename << ":" << std::endl
            << " dimension: " << dim << std::endl
            << " no. of cells: " << tria.n_active_cells() << std::endl;

  {
    std::map<unsigned int, unsigned int> boundary_count;
    typename Triangulation<dim>::active_cell_iterator
      cell = tria.begin_active(),
      endc = tria.end();
    for (; cell!=endc; ++cell)
      {
	for (unsigned int face=0;
	     face<GeometryInfo<dim>::faces_per_cell; 
	     ++face)
	  {
	    if (cell->face(face)->at_boundary())
	      boundary_count[cell->face(face)->boundary_indicator()]++;
	  }
      }
    deallog << " boundary indicators: ";
    for (std::map<unsigned int, 
	   unsigned int>::iterator it=boundary_count.begin();
	 it!=boundary_count.end();
	 ++it)
      {
	deallog << it->first << "(" << it->second << " times) ";
      }
    deallog << std::endl;
  }
}


void test ()
{  
  Point<2> center;
  
  Triangulation<2> ball;
  Triangulation<2> flat_ball;
  Triangulation<2> shell;
  Triangulation<2> tria_out;

  // create the two meshes
  GridGenerator::hyper_ball( ball, center, 0.5 );
  GridGenerator::hyper_shell(shell, center, 0.5, 1, 8, true); //colorize flag set to true so outer bnd is 1 and inner is 0

  // set boundaries 
  static const HyperBallBoundary<2> inner_bnd(center, 0.5); 
  static const HyperBallBoundary<2> outer_bnd(center, 1);

  ball.set_boundary(0, inner_bnd);
  shell.set_boundary(0, inner_bnd);
  shell.set_boundary(1, outer_bnd);

  // first we need to refine the ball once to get the right 
  // nodes to merge the meshes
  ball.refine_global(1);
  
  // flatten the refined ball so we have a coarse mesh to merge 
  flatten_triangulation( ball, flat_ball );

  mesh_info( shell, "shell.eps" );
  mesh_info( ball, "ball.eps" );
  mesh_info( flat_ball, "flat_ball.eps" );

  // now we merge (this throws an exception with no error messages 
  GridGenerator::merge_triangulations( flat_ball, shell, tria_out );
  mesh_info(tria_out, "tria_out.eps");
}


int main ()
{
  deallog << std::setprecision(2);
  logfile << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
  
  return 0;
}

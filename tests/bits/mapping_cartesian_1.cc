//----------------------------  mapping_cartesian_1.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mapping_cartesian_1.cc  ---------------------------


// there used to be a bug in MappingCartesian, where we used the size
// of the quadrature points array to initialize something else. this
// yielded a wrong result and after factoring out some code an abort

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_cartesian.h>

#include <fstream>



template <int dim>
void check (const Triangulation<dim> &tria)
{
  MappingCartesian<dim> mapping;
  FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  QGauss<dim-1> q_face(3);
  
  FEFaceValues<dim>    fe_face_values (mapping, fe, q_face,
                                       update_normal_vectors);
  fe_face_values.reinit (dof_handler.begin_active(), 0);

  deallog << "OK" << std::endl;
}


int main () 
{
  std::ofstream logfile("mapping_cartesian_1/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  {  
    Triangulation<2> coarse_grid;
    GridGenerator::hyper_cube (coarse_grid);
    check (coarse_grid);
  }

  {  
    Triangulation<3> coarse_grid;
    GridGenerator::hyper_cube (coarse_grid);
    check (coarse_grid);
  }  
}

  
  

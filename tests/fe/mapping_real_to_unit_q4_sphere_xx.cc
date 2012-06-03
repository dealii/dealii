//-----------------------------------------------------------------------------
//    $Id: data_out_base_pvd.cc 25569 2012-05-30 12:53:31Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------


// a redux of mapping_real_to_unit_q4_sphere_x.
// It seems that the point is outside the cell??


#include "../tests.h"

#include <deal.II/base/utilities.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/numerics/data_out.h>

void test_real_to_unit_cell()
{
  const unsigned int dim = 3;

				   // define a spherical cap boundary
				   // to be used for one of the faces
  const double radius = Point<dim>(1.43757e-10, 4.48023e+06, 4.48023e+06).norm();
  HyperBallBoundary<dim> boundary (Point<dim>(), radius);

				   // create the mesh: a single cell
				   // with the following coordinates:
  std::vector<Point<dim> > vertices;
  vertices.push_back (Point<dim>( 6.70384e-11, 3.17728e+06, 3.17728e+06));
  vertices.push_back (Point<dim>( -1.46060e+06, 3.99043e+06, 1.46060e+06));
  vertices.push_back (Point<dim>( -1.46060e+06, 1.46060e+06, 3.99043e+06));
  vertices.push_back (Point<dim>( -2.59424e+06, 2.59424e+06, 2.59424e+06));
  vertices.push_back (Point<dim>( 1.43757e-10, 4.48023e+06, 4.48023e+06));
  vertices.push_back (Point<dim>( -2.05956e+06, 5.62684e+06, 2.05956e+06));
  vertices.push_back (Point<dim>( -2.05956e+06, 2.05956e+06, 5.62684e+06));
  vertices.push_back (Point<dim>( -3.65809e+06, 3.65809e+06, 3.65809e+06));
				   // the points above don't show
				   // enough digits to have the same
				   // outer radius, so normalize the
				   // four outer ones
  for (unsigned int v=4; v<8; ++v)
    vertices[v] *= radius/vertices[v].norm();
  std::vector<CellData<dim> > cells;
  {
    CellData<dim> d = {{0,1,2,3,4,5,6,7},{0}};
    cells.push_back(d);
  }
  Triangulation<dim> triangulation;
  triangulation.create_triangulation (vertices, cells,
				      SubCellData());


  FE_Q<dim> fe (1);
  DoFHandler<dim> dh(triangulation);
  dh.distribute_dofs (fe);
  Vector<double> dummy;
  dummy.reinit(dh.n_dofs());
  DataOut<dim, DoFHandler<dim> > data_out;
  data_out.attach_dof_handler (dh);
  data_out.add_data_vector (dummy,"dummy");
  data_out.build_patches ();
  std::ofstream output ("mapping_real_to_unit_q4_sphere_xx/plot.vtk");
  data_out.write_vtk (output);
  
				   
  // const Point<dim> p (-3.56413e+06, 1.74215e+06, 2.14762e+06);
  // MappingQ1<dim> map;
  // typename Triangulation<dim >::active_cell_iterator
  //   cell = triangulation.begin_active();
  // deallog << map.transform_real_to_unit_cell(cell,p) << std::endl;
}


int
main()
{
  std::ofstream logfile ("mapping_real_to_unit_q4_sphere_xx/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test_real_to_unit_cell();

  return 0;
}




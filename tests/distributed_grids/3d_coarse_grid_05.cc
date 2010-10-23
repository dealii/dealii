//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// Test interaction with p4est with a 3d mesh. this is a redux of the
// _03 test where p4est errors out for a mesh that has multiple
// disconnected components
//
// because it happens to be simple to do that here, we also test the
// corresponding 2d case here

#include "../tests.h"
#include "coarse_grid_common.h"
#include <base/logstream.h>
#include <base/tensor.h>
#include <grid/tria.h>
#include <distributed/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>
#include <grid/grid_in.h>



// create a mesh that consists of two disconnected hypercubes
template <int dim>
void create_disconnected_mesh (Triangulation<dim> &tria)
{
  std::vector<Point<dim> > vertices (2*GeometryInfo<dim>::vertices_per_cell);
  std::vector<CellData<dim> > cells (2);

				   // cell 1
  {
    Point<dim> p1 = Point<dim>();
    Point<dim> p2;
    for (unsigned int i=0; i<dim; ++i)
      p2[i] = 1;

    switch (dim)
      {
	case 2:
	      vertices[0] = vertices[1] = p1;
	      vertices[2] = vertices[3] = p2;

	      vertices[1](0) = p2(0);
	      vertices[2](0) = p1(0);
	      break;
	case 3:
	      vertices[0] = vertices[1] = vertices[2] = vertices[3] = p1;
	      vertices[4] = vertices[5] = vertices[6] = vertices[7] = p2;

	      vertices[1](0) = p2(0);
	      vertices[2](1) = p2(1);
	      vertices[3](0) = p2(0);
	      vertices[3](1) = p2(1);

	      vertices[4](0) = p1(0);
	      vertices[4](1) = p1(1);
	      vertices[5](1) = p1(1);
	      vertices[6](0) = p1(0);

	      break;
	default:
	      Assert (false, ExcNotImplemented ());
      }

				     // Prepare cell data
    for (unsigned int i=0;i<GeometryInfo<dim>::vertices_per_cell;++i)
      cells[0].vertices[i] = i;
    cells[0].material_id = 0;
  }

				   // cell 2. shifted 2 units in x-direction
  {
    Point<dim> p1 = Point<dim>();
    Point<dim> p2;
    for (unsigned int i=0; i<dim; ++i)
      p2[i] = 1;

    p1[0] += 2;
    p2[0] += 2;

    switch (dim)
      {
	case 2:
	      vertices[GeometryInfo<dim>::vertices_per_cell+0] = vertices[GeometryInfo<dim>::vertices_per_cell+1] = p1;
	      vertices[GeometryInfo<dim>::vertices_per_cell+2] = vertices[GeometryInfo<dim>::vertices_per_cell+3] = p2;

	      vertices[GeometryInfo<dim>::vertices_per_cell+1](0) = p2(0);
	      vertices[GeometryInfo<dim>::vertices_per_cell+2](0) = p1(0);
	      break;
	case 3:
	      vertices[GeometryInfo<dim>::vertices_per_cell+0] = vertices[GeometryInfo<dim>::vertices_per_cell+1] = vertices[GeometryInfo<dim>::vertices_per_cell+2] = vertices[GeometryInfo<dim>::vertices_per_cell+3] = p1;
	      vertices[GeometryInfo<dim>::vertices_per_cell+4] = vertices[GeometryInfo<dim>::vertices_per_cell+5] = vertices[GeometryInfo<dim>::vertices_per_cell+6] = vertices[GeometryInfo<dim>::vertices_per_cell+7] = p2;

	      vertices[GeometryInfo<dim>::vertices_per_cell+1](0) = p2(0);
	      vertices[GeometryInfo<dim>::vertices_per_cell+2](1) = p2(1);
	      vertices[GeometryInfo<dim>::vertices_per_cell+3](0) = p2(0);
	      vertices[GeometryInfo<dim>::vertices_per_cell+3](1) = p2(1);

	      vertices[GeometryInfo<dim>::vertices_per_cell+4](0) = p1(0);
	      vertices[GeometryInfo<dim>::vertices_per_cell+4](1) = p1(1);
	      vertices[GeometryInfo<dim>::vertices_per_cell+5](1) = p1(1);
	      vertices[GeometryInfo<dim>::vertices_per_cell+6](0) = p1(0);

	      break;
	default:
	      Assert (false, ExcNotImplemented ());
      }

				     // Prepare cell data
    for (unsigned int i=0;i<GeometryInfo<dim>::vertices_per_cell;++i)
      cells[1].vertices[i] = GeometryInfo<dim>::vertices_per_cell+i;
    cells[1].material_id = 0;
  }

  tria.create_triangulation (vertices, cells, SubCellData());
}



template<int dim>
void test(std::ostream& /*out*/)
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  create_disconnected_mesh (tr);

  write_vtk (tr, "3d_coarse_grid_05", "1");
}


int main(int argc, char *argv[])
{
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Init (&argc,&argv);
#else
  (void)argc;
  (void)argv;
#endif

  std::ofstream logfile("3d_coarse_grid_05/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("2d");
  test<2>(logfile);
  deallog.pop();

  deallog.push("3d");
  test<3>(logfile);
  deallog.pop();

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Finalize();
#endif
}

//----------------------------  rt_2.cc  ---------------------------
//    rt_2.cc,v 1.3 2003/06/09 16:00:38 wolf Exp
//    Version: 
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  rt_2.cc  ---------------------------

// Show the shape functions of the Raviart-Thomas element on a grid
// with only one cell. This cell is rotated, stretched, scaled, etc,
// and on each of these cells each time we evaluate the shape
// functions.

#include "../tests.h"
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <grid/grid_generator.h>
#include <grid/grid_tools.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_values.h>

#include <vector>
#include <fstream>
#include <string>

#define PRECISION 5


Point<2> stretch_coordinates (const Point<2> p)
{
  return Point<2>(2*p(0), p(1));
}



Point<2> tilt_coordinates (const Point<2> p)
{
  return Point<2>(p(0)+p(1), p(1));
}


void
transform_grid (Triangulation<2> &tria,
		const unsigned int  transform)
{
  switch (transform)
    {
                                       // first round: take
                                       // original grid
      case 0:
	    break;
            
					     // second round: rotate
					     // triangulation
      case 1:
	    GridTools::rotate (3.14159265358/2, tria);
	    break;

					     // third round: inflate
					     // by a factor of 2
      case 2:
	    GridTools::scale (2, tria);
	    break;

                                             // third round: scale
                                             // back, rotate back,
                                             // stretch
      case 3:
            GridTools::scale (.5, tria);
	    GridTools::rotate (-3.14159265358/2, tria);
            GridTools::transform (&stretch_coordinates, tria);
            
	    break;
            
      default:
	    Assert (false, ExcNotImplemented());
    };
}

	    


template<int dim>
void
plot_shape_functions(const unsigned int degree)
{
  FE_RaviartThomas<dim> fe_rt(degree);
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, 0., 1.);

				   // check the following with a
				   // number of transformed
				   // triangulations
  for (unsigned int transform=0; transform<4; ++transform)
    {
      deallog << "GRID TRANSFORMATION " << transform << std::endl;
      
      transform_grid (tr, transform);

      DoFHandler<dim> dof(tr);
      typename DoFHandler<dim>::cell_iterator c = dof.begin();
      dof.distribute_dofs(fe_rt);
      
      QTrapez<1> q_trapez;
      const unsigned int div=2;
      QIterated<dim> q(q_trapez, div);
      FEValues<dim> fe(fe_rt, q, update_values|update_gradients|update_q_points);
      fe.reinit(c);
      
      unsigned int q_point=0;
      for (unsigned int mz=0;mz<=((dim>2) ? div : 0) ;++mz)
	for (unsigned int my=0;my<=((dim>1) ? div : 0) ;++my)
	  for (unsigned int mx=0;mx<=div;++mx)
	    {
	      deallog << "q_point(" << q_point << ")=" << fe.quadrature_point(q_point)
		      << std::endl;
	      
	      for (unsigned int i=0;i<fe_rt.dofs_per_cell;++i)
		{
		  deallog << "  phi(" << i << ") = [";
		  for (unsigned int c=0; c<fe.get_fe().n_components(); ++c)
		    deallog << " " << fe.shape_value_component(i,q_point,c);
		  deallog << " ]" << std::endl;
		};
	      for (unsigned int i=0;i<fe_rt.dofs_per_cell;++i)
		{
		  deallog << "  grad phi(" << i << ") = ";
		  for (unsigned int c=0; c<dim; ++c)
                    {
                      deallog << "[";
                      for (unsigned int d=0; d<dim; ++d)
                        deallog << " " << fe.shape_grad_component(i,q_point,c)[d];
                      deallog << " ]" << std::endl;
                      if (c != dim-1)
                        deallog << "                ";
                    };
		};

              q_point++;
	    }

      deallog << std::endl;
    }
}


int
main()
{
  std::ofstream logfile ("rt_2.output");
  logfile.precision (PRECISION);
  logfile.setf(std::ios::fixed);  
  deallog.attach(logfile);
  deallog.depth_console(0);

  for (unsigned int degree=0; degree<4; ++degree)
    plot_shape_functions<2>(degree);
  
  return 0;
}




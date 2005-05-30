//----------------------------------------------------------------------
//    rt_1.cc,v 1.3 2003/06/09 16:00:38 wolf Exp
//    Version: 
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

// FE_RaviartThomasNodal
// Show the shape functions in support points

#include "../tests.h"
#include <base/logstream.h>
#include <grid/grid_generator.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/mapping_cartesian.h>
#include <fe/fe_values.h>

#include <vector>
#include <fstream>
#include <string>

#include <cmath>
extern "C" 
{
#include <math.h>
};

template <int dim>
void
check_support_points (const FiniteElement<dim>& fe)
{
  deallog << fe.get_name() << std::endl;
  
  const std::vector< Point<dim > > & points = fe.get_unit_support_points();
  std::vector<double> weights(points.size());

  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);
  
  MappingCartesian<dim> mapping;
  
  Quadrature<dim> support_quadrature(points, weights);
  UpdateFlags flags = update_values | update_gradients;
  FEValues<dim> vals(mapping, fe, support_quadrature, flags);
  vals.reinit(dof.begin_active());
  
  for (unsigned int k=0;k<points.size();++k)
    {
      const Point<dim>& p = points[k];
      deallog << p;
      for (unsigned int i=0;i<fe.dofs_per_cell;++i)
	{
	  deallog << '\t';
	  for (unsigned int d=0;d<dim;++d)
	    {
	      const double sf = fe.shape_value_component(i,p,d);
	      deallog << ' ' << (int) round(60*sf);
	      
 	      const double diff = std::abs(sf - vals.shape_value_component(i,k,d));
 	      if (diff > 1.e-12)
 		deallog << "Error values" << std::endl;
	      Tensor<1,dim> grad = fe.shape_grad_component(i,p,d);
	      grad -= vals.shape_grad_component(i,k,d);
	      if (grad.norm() > 1.e-12)
 		deallog << "Error grads" << std::endl;
	    }
	}
      deallog << std::endl;
    }
  deallog << std::endl;
}


template <int dim>
void
check_face_support_points (const FiniteElement<dim>& fe)
{
  deallog << "Face " << fe.get_name() << std::endl;
  
  const std::vector< Point<dim-1> > & sub_points = fe.get_unit_face_support_points();
  std::vector<double> weights(sub_points.size());
  Quadrature<dim-1> sub_quadrature(sub_points, weights);
  std::vector< Point<dim> > points(sub_points.size());
  
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);
  
  MappingCartesian<dim> mapping;
  
  UpdateFlags flags = update_values | update_gradients;
  FEFaceValues<dim> vals(mapping, fe, sub_quadrature, flags);

  for (unsigned int face=0;face<GeometryInfo<dim>::faces_per_cell;++face)
    {
      QProjector<dim>::project_to_face(sub_quadrature, face, points);
      vals.reinit(dof.begin_active(), face);
      
      for (unsigned int k=0;k<points.size();++k)
	{
	  const Point<dim>& p = points[k];
	  deallog << p;
	  for (unsigned int i=0;i<fe.dofs_per_cell;++i)
	    {
	      deallog << '\t';
	      for (unsigned int d=0;d<dim;++d)
		{
		  const double sf = fe.shape_value_component(i,p,d);
		  deallog << ' ' << (int) round(60*sf);
		  
		  const double diff = std::abs(sf - vals.shape_value_component(i,k,d));
		  if (diff > 1.e-12)
		    deallog << "Error values" << std::endl;
		  Tensor<1,dim> grad = fe.shape_grad_component(i,p,d);
		  grad -= vals.shape_grad_component(i,k,d);
		  if (grad.norm() > 1.e-12)
		    deallog << "Error grads" << std::endl;
		}
	    }
	  deallog << std::endl;
	}
      deallog << std::endl;
    }
}


int
main()
{
  std::ofstream logfile ("rtn_1.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  FE_RaviartThomasNodal<2> e20(0);
  check_support_points(e20);
  FE_RaviartThomasNodal<3> e30(0);
  check_support_points(e30);
  check_face_support_points(e20);
  FE_RaviartThomasNodal<2> e21(1);
  check_support_points(e21);
  check_face_support_points(e21);
}

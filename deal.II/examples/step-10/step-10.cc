/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 1999 */

#include <base/quadrature_lib.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_q.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <grid/tria_boundary_lib.h>
#include <fe/mapping_q.h>




template <int dim>
void compute_pi_by_area ()
{
  std::cout << "Computation of Pi by the area:" << std::endl
	    << "==============================" << std::endl;
  
  for (unsigned int order=1; order<5; ++order)
    {
      cout << "Order = " << order << endl;
      Triangulation<dim> triangulation;
      GridGenerator::hyper_ball (triangulation);
  
      static const HyperBallBoundary<dim> boundary;
      triangulation.set_boundary (0, boundary);

      const MappingQ<dim> mapping (order);
      const FE_Q<dim>     fe (1);

      DoFHandler<dim> dof_handler (triangulation);
  
      for (unsigned int refinement=0; refinement<5; ++refinement)
	{
	  triangulation.refine_global (1);
	  dof_handler.distribute_dofs (fe);
	  
	  QGauss4<dim> quadrature;
	  FEValues<dim> fe_values (mapping, fe, quadrature, update_JxW_values);
	  
	  typename DoFHandler<dim>::active_cell_iterator
	    cell = dof_handler.begin_active(),
	    endc = dof_handler.end();
	  double area = 0;
	  for (; cell!=endc; ++cell)
	    {
	      fe_values.reinit (cell);
	      for (unsigned int i=0; i<fe_values.n_quadrature_points; ++i)
		area += fe_values.JxW (i);
	    };
	  std::cout << "Pi=" << area
		    << ", error=" << fabs(area-3.141592653589793238462643)
		    << std::endl;
	};
      std::cout << std::endl;
    };
};



template <int dim>
void compute_pi_by_perimeter ()
{
  std::cout << "Computation of Pi by the perimeter:" << std::endl
	    << "===================================" << std::endl;

  for (unsigned int order=1; order<5; ++order)
    {
      cout << "Order = " << order << endl;
      Triangulation<dim> triangulation;
      GridGenerator::hyper_ball (triangulation);
  
      static const HyperBallBoundary<dim> boundary;
      triangulation.set_boundary (0, boundary);

      const MappingQ<dim> mapping (order);
      const FE_Q<dim>     fe (1);

      DoFHandler<dim> dof_handler (triangulation);
  
      for (unsigned int refinement=0; refinement<5; ++refinement)
	{
	  triangulation.refine_global (1);
	  dof_handler.distribute_dofs (fe);
	  
	  QGauss4<dim-1> quadrature;
	  FEFaceValues<dim> fe_face_values (mapping, fe, quadrature, update_JxW_values);
	  
	  typename DoFHandler<dim>::active_cell_iterator
	    cell = dof_handler.begin_active(),
	    endc = dof_handler.end();
	  double perimeter = 0;
	  for (; cell!=endc; ++cell)
	    for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
	      if (cell->face(face_no)->at_boundary())
		{
		  fe_face_values.reinit (cell, face_no);
		  for (unsigned int i=0; i<fe_face_values.n_quadrature_points; ++i)
		    perimeter += fe_face_values.JxW (i);
		};
	  std::cout << "Pi=" << perimeter/2
		    << ", error=" << fabs(perimeter/2-3.141592653589793238462643)
		    << std::endl;
	};
      std::cout << std::endl;
    };
};




int main () 
{
  cout.precision (16);

  compute_pi_by_area<2> ();
  compute_pi_by_perimeter<2> ();
  
  return 0;
};

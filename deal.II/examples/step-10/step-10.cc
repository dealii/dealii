/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 1999 */

#include <base/quadrature_lib.h>
#include <base/convergence_table.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/tria.h>
#include <grid/grid_out.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <fe/mapping_q.h>

#include <fstream>


static const long double pi=3.141592653589793238462643;


template <int dim>
void compute_pi_by_area ()
{
  std::cout << "Computation of Pi by the area:" << std::endl
	    << "==============================" << std::endl;
  
  const QGauss4<dim> quadrature;
  for (unsigned int order=1; order<5; ++order)
    {
      cout << "Order = " << order << std::endl;
      Triangulation<dim> triangulation;
      GridGenerator::hyper_ball (triangulation);
  
      static const HyperBallBoundary<dim> boundary;
      triangulation.set_boundary (0, boundary);

      const MappingQ<dim> mapping (order);
      const FE_Q<dim>     fe (1);

      DoFHandler<dim> dof_handler (triangulation);
  
      FEValues<dim> fe_values (mapping, fe, quadrature, update_JxW_values);
      ConvergenceTable table;
      
      for (unsigned int refinement=0; refinement<6;
	   ++refinement, triangulation.refine_global (1))
	{
	  table.add_value("cells", triangulation.n_active_cells());
	    
	  dof_handler.distribute_dofs (fe);  
	  
	  typename DoFHandler<dim>::active_cell_iterator
	    cell = dof_handler.begin_active(),
	    endc = dof_handler.end();
	  long double area = 0;
	  for (; cell!=endc; ++cell)
	    {
	      fe_values.reinit (cell);
	      for (unsigned int i=0; i<fe_values.n_quadrature_points; ++i)
		area += fe_values.JxW (i);
	    };
	  table.add_value("eval.pi", static_cast<double> (area));
	  table.add_value("error", fabs(area-pi));
	};

      table.omit_column_from_convergence_rate_evaluation("cells");
      table.omit_column_from_convergence_rate_evaluation("eval.pi");
      table.evaluate_all_convergence_rates(
	ConvergenceTable::reduction_rate_log2);

      table.set_precision("eval.pi", 16);
      table.set_scientific("error", true);

      table.write_text(cout);
    };
};



template <int dim>
void compute_pi_by_perimeter ()
{
  std::cout << "Computation of Pi by the perimeter:" << std::endl
	    << "===================================" << std::endl;

  const QGauss4<dim-1> quadrature;
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
      
      FEFaceValues<dim> fe_face_values (mapping, fe, quadrature, update_JxW_values);
      ConvergenceTable table;

      for (unsigned int refinement=0; refinement<6;
	   ++refinement, triangulation.refine_global (1))
	{
	  table.add_value("cells", triangulation.n_active_cells());

	  dof_handler.distribute_dofs (fe);
	  
	  typename DoFHandler<dim>::active_cell_iterator
	    cell = dof_handler.begin_active(),
	    endc = dof_handler.end();
	  long double perimeter = 0;
	  for (; cell!=endc; ++cell)
	    for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
	      if (cell->face(face_no)->at_boundary())
		{
		  fe_face_values.reinit (cell, face_no);
		  for (unsigned int i=0; i<fe_face_values.n_quadrature_points; ++i)
		    perimeter += fe_face_values.JxW (i);
		};
	  table.add_value("eval.pi", static_cast<double> (perimeter/2.));
	  table.add_value("error", fabs(perimeter/2.-pi));
	};

      table.omit_column_from_convergence_rate_evaluation("cells");
      table.omit_column_from_convergence_rate_evaluation("eval.pi");
      table.evaluate_all_convergence_rates(
	ConvergenceTable::reduction_rate_log2);

      table.set_precision("eval.pi", 16);
      table.set_scientific("error", true);

      table.write_text(cout);
    };
};


template <int dim>
void gnuplot_output()
{
  cout << "Output of grids as gnuplot file:" << std::endl;
  
  static const HyperBallBoundary<dim> boundary;
  
  GridOut grid_out;
				   // on boundary faces plot 30
				   // additional points per face.
  GridOutFlags::Gnuplot gnuplot_flags(false, 30);
  grid_out.set_flags(gnuplot_flags);
  

  for (unsigned int order=1; order<4; ++order)
    {
      cout << "Order = " << order << std::endl;

      Triangulation<dim> triangulation;
      GridGenerator::hyper_ball (triangulation);
      triangulation.set_boundary (0, boundary);
      
      const MappingQ<dim> mapping (order);
      string filename_base="ball_mapping_q";
      filename_base += ('0'+order);

      
      for (unsigned int refinement=0; refinement<2;
	   ++refinement, triangulation.refine_global(1))
	{	  
	  string filename=filename_base+"_ref";
	  filename += ('0'+refinement);
	  filename += ".dat";
	  ofstream gnuplot_file(filename.c_str());

	  grid_out.write_gnuplot(triangulation, gnuplot_file, &mapping);
	}
    }
}



int main () 
{
  cout.precision (16);

  gnuplot_output<2>();

  compute_pi_by_area<2> ();
  compute_pi_by_perimeter<2> ();
  
  return 0;
};

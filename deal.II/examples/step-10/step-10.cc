/* $Id$ */
/* Author: Wolfgang Bangerth, Ralf Hartmann, University of Heidelberg, 2001 */

				 // The first of the following include
				 // files are probably well-known by
				 // now and need no further
				 // explanation.
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

				 // This is the only new one: in it,
				 // we declare the ``MappingQ'' class
				 // which we will use for polynomial
				 // mappings of arbitrary order:
#include <fe/mapping_q.h>

				 // And this again is C++:
#include <fstream>


				 // Now, as we want to compute the
				 // value of pi, we have to compare to
				 // somewhat. These are the first few
				 // digits of pi, which we define
				 // beforehand for later use. Since we
				 // would like to compute the
				 // difference between two numbers
				 // which are quite accurate, with the
				 // accuracy of the computed
				 // approximation to pi being in the
				 // range of the number of digits
				 // which a double variable can hold,
				 // we rather declare the reference
				 // value as a ``long double'' and
				 // give it a number of extra digits:
const long double pi = 3.141592653589793238462643;



				 // Then, the first task will be to
				 // generate some output. Since this
				 // program is so small, we do not
				 // employ object oriented techniques
				 // in it and do not declare classes
				 // (although, of course, we use the
				 // object oriented features of the
				 // library). Rather, we just pack the
				 // functionality into separate
				 // functions. We make these functions
				 // templates on the number of space
				 // dimensions to conform to usual
				 // practice when using deal.II,
				 // although we will only use them for
				 // two space dimensions.
				 //
				 // The first of these functions just
				 // generates a triangulation of a
				 // circle (hyperball) and outputs the
				 // Qp mapping of its cells for
				 // different values of ``p''. Then,
				 // we refine the grid once and do so
				 // again.
template <int dim>
void gnuplot_output()
{
  std::cout << "Output of grids into gnuplot files:" << std::endl
	    << "===================================" << std::endl;

				   // So first generate a coarse
				   // triangulation of the circle and
				   // associate a suitable boundary
				   // description to it:
  Triangulation<dim> triangulation;
  GridGenerator::hyper_ball (triangulation);
  static const HyperBallBoundary<dim> boundary;
  triangulation.set_boundary (0, boundary);

				   // Next generate output for this
				   // grid and for a once refined
				   // grid. Note that we have hidden
				   // the mesh refinement in the loop
				   // header, which might be uncommon
				   // but nevertheless works. Also it
				   // is strangly consistent with
				   // incrementing the loop index
				   // denoting the refinement level.
  for (unsigned int refinement=0; refinement<2;
       ++refinement, triangulation.refine_global(1))
    {
      std::cout << "Refinement level: " << refinement << std::endl;

				       // Then have a string which
				       // denotes the base part of the
				       // names of the files into
				       // which we write the
				       // output. Note that in the
				       // parentheses in the
				       // initializer we do arithmetic
				       // on characters, which assumes
				       // that first the characters
				       // denoting numbers are placed
				       // consecutively (which is
				       // probably true for all
				       // reasonable character sets
				       // nowadays), but also assumes
				       // that the increment
				       // ``refinement'' is less than
				       // ten. This is therefore more
				       // a quick hack if we know
				       // exactly the values which the
				       // increment can assume. A
				       // better implementation would
				       // use the
				       // ``std::istringstream''
				       // class to generate a name.
      std::string filename_base = std::string("ball");
      filename_base += '0'+refinement;

				       // Then output the present grid
				       // for Q1, Q2, and Q3 mappings:
      for (unsigned int order=1; order<4; ++order)
	{
	  std::cout << "Order = " << order << std::endl;

					   // For this, first set up
					   // an object describing the
					   // mapping. This is done
					   // using the ``MappingQ''
					   // class, which takes as
					   // argument to the
					   // constructor the
					   // polynomial order which
					   // it shall use.
	  const MappingQ<dim> mapping (order);
					   // We note one interesting
					   // fact: if you want a
					   // piecewise linear
					   // mapping, then you could
					   // give a value of ``1'' to
					   // the
					   // constructor. However,
					   // for linear mappings, so
					   // many things can be
					   // generated simpler that
					   // there is another class,
					   // called ``MappingQ1''
					   // which does exactly the
					   // same is if you gave an
					   // order of ``1'' to the
					   // ``MappingQ'' class, but
					   // does so significantly
					   // faster. ``MappingQ1'' is
					   // also the class that is
					   // implicitely used
					   // throughout the library
					   // in many functions and
					   // classes if you do not
					   // specify another mapping
					   // explicitly.


					   // In order to actually
					   // write out the present
					   // grid with this mapping,
					   // we set up an object
					   // which we will use for
					   // output. We will generate
					   // Gnuplot output, which
					   // consists of a set of
					   // lines describing the
					   // mapped triangulation. By
					   // default, only one line
					   // is drawn for each face
					   // of the triangulation,
					   // but since we want to
					   // explicitely see the
					   // effect of the mapping,
					   // we want to have teh
					   // faces in more
					   // detail. This can be done
					   // by passing the output
					   // object a structure which
					   // contains some flags. In
					   // the present case, since
					   // Gnuplot can only draw
					   // straight lines, we
					   // output a number of
					   // additional points on the
					   // faces so that each face
					   // is drawn by 30 small
					   // lines instead of only
					   // one. This is sufficient
					   // to give us the
					   // impression of seeing a
					   // curved line, rather than
					   // a set of straight lines.
	  GridOut grid_out;
	  GridOutFlags::Gnuplot gnuplot_flags(false, 30);
	  grid_out.set_flags(gnuplot_flags);
  
					   // Finally, generate a
					   // filename and a file for
					   // output using the same
					   // evil hack as above:
	  std::string filename = filename_base+"_mapping_q";
	  filename += ('0'+order);
	  filename += ".dat";
	  std::ofstream gnuplot_file (filename.c_str());

					   // Then write out the
					   // triangulation to this
					   // file. The last argument
					   // of the function is a
					   // pointer to a mapping
					   // object. This argument
					   // has a default value, and
					   // if no value is given a
					   // simple ``MappingQ1''
					   // object is taken, which
					   // we briefly described
					   // above. This would then
					   // result in a piecewise
					   // linear approximation of
					   // the true boundary in the
					   // output.
	  grid_out.write_gnuplot (triangulation, gnuplot_file, &mapping);
	}
      std::cout << std::endl;
    }
}



template <int dim>
void compute_pi_by_area ()
{
  std::cout << "Computation of Pi by the area:" << std::endl
	    << "==============================" << std::endl;
  
  const QGauss4<dim> quadrature;
  for (unsigned int order=1; order<5; ++order)
    {
      std::cout << "Order = " << order << std::endl;
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

      table.write_text(std::cout);

      std::cout << std::endl;
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
      std::cout << "Order = " << order << std::endl;
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

      table.write_text(std::cout);

      std::cout << std::endl;
    };
};


int main () 
{
  std::cout.precision (16);

  gnuplot_output<2>();

  compute_pi_by_area<2> ();
  compute_pi_by_perimeter<2> ();
  
  return 0;
};

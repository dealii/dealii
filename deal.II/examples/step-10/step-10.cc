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
				   // description to it. Note that the
				   // default values of the
				   // HyperBallBoundary constructor
				   // are a center at the origin and a
				   // radius equals one.
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
      for (unsigned int degree=1; degree<4; ++degree)
	{
	  std::cout << "Degree = " << degree << std::endl;

					   // For this, first set up
					   // an object describing the
					   // mapping. This is done
					   // using the ``MappingQ''
					   // class, which takes as
					   // argument to the
					   // constructor the
					   // polynomial degree which
					   // it shall use.
	  const MappingQ<dim> mapping (degree);
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
					   // degree of ``1'' to the
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


					   // In degree to actually
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
					   // we want to have the
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
	  filename += ('0'+degree);
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

				 // Now we proceed with the main part
				 // of the code, the approximation of
				 // pi. The area of a circle is given
				 // by pi*radius^2, so having a circle
				 // of radius 1, the area represents
				 // just the number that is searched
				 // for. The numerical computation of
				 // the area is performed by
				 // integrating the constant function
				 // of value 1 over the whole
				 // computational domain, i.e. by
				 // computing the areas $\int_K 1
				 // dx=\int_{\hat K} 1 J(\hat x) d\hat
				 // x\approx\sum J(\hat x)w(\hat x)$
				 // of all active cells of
				 // triangulation and summing up these
				 // contributions to gain the area of
				 // the overall domain. The integrals
				 // on each cell are approximated by
				 // numerical quadrature, hence the
				 // only additional ingredient we need
				 // is to set up a FEValues object
				 // that provides the corresponding
				 // `JxW' values of each cell. We note
				 // that here we won't use the
				 // FEValues object in its original
				 // purpose that is computing the
				 // values of basis functions of a
				 // specific finite element. But here
				 // we use it only to gain the `JxW'
				 // at the quadrature points,
				 // irrespective of the (dummy) finite
				 // element we will give to the
				 // constructor of the FEValues
				 // object.
template <int dim>
void compute_pi_by_area ()
{
  std::cout << "Computation of Pi by the area:" << std::endl
	    << "==============================" << std::endl;

				   // For the numerical quadrature on
				   // all cells we employ a quadrature
				   // rule of sufficiently high
				   // degree. We choose QGauss4 that is
				   // of order 8, to be sure that the
				   // errors due to numerical
				   // quadrature are of higher order
				   // than the order (maximal 6) that
				   // will occur due to the order of
				   // the approximation of the
				   // boundary, i.e. the order of the
				   // mappings employed.
  const QGauss4<dim> quadrature;

				   // Now start by looping over
				   // degrees=1..4
  for (unsigned int degree=1; degree<5; ++degree)
    {
      std::cout << "Degree = " << degree << std::endl;

				       // Then we generate the
				       // triangulation, the Boundary
				       // and the Mapping object as
				       // already seen.
      Triangulation<dim> triangulation;
      GridGenerator::hyper_ball (triangulation);
  
      static const HyperBallBoundary<dim> boundary;
      triangulation.set_boundary (0, boundary);

      const MappingQ<dim> mapping (degree);

				       // We now create a dummy finite
				       // element. Here we could
				       // choose a finite element no
				       // matter which, as we are only
				       // interested in the `JxW'
				       // values provided by the
				       // FEValues object below.
      const FE_Q<dim>     dummy_fe (1);

				       // Then we create a DofHandler
				       // object. This object will
				       // provide us with
				       // `active_cell_iterators' that
				       // are needed to reinit the
				       // FEValues object on each cell
				       // of the triangulation.
      DoFHandler<dim> dof_handler (triangulation);

				       // Now we set up the FEValues
				       // object, giving the Mapping,
				       // the dummy finite element and
				       // the quadrature object to the
				       // constructor, together with
				       // the UpdateFlag asking for
				       // the `JxW' values at the
				       // quadrature points only.
      FEValues<dim> fe_values (mapping, dummy_fe, quadrature, update_JxW_values);

				       // We employ an object of the
				       // ConvergenceTable class to
				       // store all important data
				       // like the approximative
				       // values for pi and the error
				       // wrt. the true value of
				       // pi. We will use functions
				       // provided by the
				       // ConvergenceTable class to
				       // compute convergence rates of
				       // the approximations to pi.
      ConvergenceTable table;

				       // Now we loop over several
				       // refinement steps of the
				       // triangulation.
      for (unsigned int refinement=0; refinement<6;
	   ++refinement, triangulation.refine_global (1))
	{
					   // In this loop we first
					   // add the number of active
					   // cells of the current
					   // triangulation to the
					   // table. This function
					   // automatically creates a
					   // table column with
					   // superscription `cells',
					   // for the case this column
					   // was not created before.
	  table.add_value("cells", triangulation.n_active_cells());

					   // Then we distribute the
					   // degrees of freedoms for
					   // the dummy finite
					   // element. Strictly
					   // speaking we do not need
					   // this function call in
					   // our special case but we
					   // call it to make the
					   // DoFHandler happy --
					   // otherwise it would throw
					   // an assertion in the
					   // FEValues::reinit
					   // function below.
	  dof_handler.distribute_dofs (dummy_fe);

					   // We define the variable
					   // area as `long double'
					   // like we did for the pi
					   // variable before.
	  long double area = 0;

					   // Now we loop over all
					   // cells, reinit the
					   // FEValues object for each
					   // cell, add all `JxW'
					   // values to `area'
	  typename DoFHandler<dim>::active_cell_iterator
	    cell = dof_handler.begin_active(),
	    endc = dof_handler.end();
	  for (; cell!=endc; ++cell)
	    {
	      fe_values.reinit (cell);
	      for (unsigned int i=0; i<fe_values.n_quadrature_points; ++i)
		area += fe_values.JxW (i);
	    };

					   // and store the resulting
					   // area values and the
					   // errors in the table. We
					   // need a static cast to
					   // double as there is no
					   // add_value(string, long
					   // double) function
					   // implemented.
	  table.add_value("eval.pi", static_cast<double> (area));
	  table.add_value("error", fabs(area-pi));
	};

				       // We want to compute
				       // the convergence rates of the
				       // `error' column. Therefore we
				       // need to omit the other
				       // columns from the convergence
				       // rate evaluation before
				       // calling
				       // `evaluate_all_convergence_rates'
      table.omit_column_from_convergence_rate_evaluation("cells");
      table.omit_column_from_convergence_rate_evaluation("eval.pi");
      table.evaluate_all_convergence_rates(
	ConvergenceTable::reduction_rate_log2);

				       // Finally we set the precision
				       // and the scientific mode
      table.set_precision("eval.pi", 16);
      table.set_scientific("error", true);

				       // and write the whole table to
				       // cout.
      table.write_text(std::cout);

      std::cout << std::endl;
    };
};


				 // The following function also
				 // computes an approximation of pi
				 // but this time via the diameter
				 // 2*pi*radius of the domain instead
				 // of the area. This function is only
				 // a variation of the previous
				 // function. So we will mainly give
				 // documentation for the differences.
template <int dim>
void compute_pi_by_perimeter ()
{
  std::cout << "Computation of Pi by the perimeter:" << std::endl
	    << "===================================" << std::endl;

				   // We take the same order of
				   // quadrature but this time a
				   // `dim-1' dimensional quadrature
				   // as we will integrate over
				   // (boundary) lines rather than
				   // over cells.
  const QGauss4<dim-1> quadrature;

				   // We loop over all degrees, create
				   // the Triangulation, the Boundary,
				   // the Mapping, the dummy
				   // FiniteElement and the DoFHandler
				   // object as seen before.
  for (unsigned int degree=1; degree<5; ++degree)
    {
      std::cout << "Degree = " << degree << std::endl;
      Triangulation<dim> triangulation;
      GridGenerator::hyper_ball (triangulation);
  
      static const HyperBallBoundary<dim> boundary;
      triangulation.set_boundary (0, boundary);

      const MappingQ<dim> mapping (degree);
      const FE_Q<dim>     fe (1);

      DoFHandler<dim> dof_handler (triangulation);

				       // Then we create a FEFaceValues
				       // object instead of a FEValues
				       // object as in the previous
				       // function.
      FEFaceValues<dim> fe_face_values (mapping, fe, quadrature, update_JxW_values);
      ConvergenceTable table;

      for (unsigned int refinement=0; refinement<6;
	   ++refinement, triangulation.refine_global (1))
	{
	  table.add_value("cells", triangulation.n_active_cells());

	  dof_handler.distribute_dofs (fe);

					   // Now we run over all
					   // cells and over all faces
					   // of each cell. Only the
					   // contributions of the
					   // `JxW' values on boundary
					   // faces are added to the
					   // long double variable
					   // `perimeter'.
	  typename DoFHandler<dim>::active_cell_iterator
	    cell = dof_handler.begin_active(),
	    endc = dof_handler.end();
	  long double perimeter = 0;
	  for (; cell!=endc; ++cell)
	    for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
	      if (cell->face(face_no)->at_boundary())
		{
						   // We reinit the
						   // FEFaceValues
						   // object with the
						   // cell iterator
						   // and the number
						   // of the face.
		  fe_face_values.reinit (cell, face_no);
		  for (unsigned int i=0; i<fe_face_values.n_quadrature_points; ++i)
		    perimeter += fe_face_values.JxW (i);
		};
					   // We store the evaluated
					   // values in the table
	  table.add_value("eval.pi", static_cast<double> (perimeter/2.));
	  table.add_value("error", fabs(perimeter/2.-pi));
	};

				       // and we end this function as
				       // we did in the previous
				       // function.
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


				 // The following main function just
				 // calles the above functions in the
				 // order of appearance.
int main () 
{
  std::cout.precision (16);

  gnuplot_output<2>();

  compute_pi_by_area<2> ();
  compute_pi_by_perimeter<2> ();
  
  return 0;
};

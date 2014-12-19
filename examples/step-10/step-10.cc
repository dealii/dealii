/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2001 - 2013 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Authors: Wolfgang Bangerth, Ralf Hartmann, University of Heidelberg, 2001
 */


// The first of the following include files are probably well-known by now and
// need no further explanation.
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

// This is the only new one: in it, we declare the <code>MappingQ</code> class
// which we will use for polynomial mappings of arbitrary order:
#include <deal.II/fe/mapping_q.h>

// And this again is C++:
#include <iostream>
#include <fstream>
#include <cmath>

// The last step is as in previous programs:
namespace Step10
{
  using namespace dealii;

  // Now, as we want to compute the value of $\pi$, we have to compare to
  // somewhat. These are the first few digits of $\pi$, which we define
  // beforehand for later use. Since we would like to compute the difference
  // between two numbers which are quite accurate, with the accuracy of the
  // computed approximation to $\pi$ being in the range of the number of
  // digits which a double variable can hold, we rather declare the reference
  // value as a <code>long double</code> and give it a number of extra digits:
  const long double pi = 3.141592653589793238462643;



  // Then, the first task will be to generate some output. Since this program
  // is so small, we do not employ object oriented techniques in it and do not
  // declare classes (although, of course, we use the object oriented features
  // of the library). Rather, we just pack the functionality into separate
  // functions. We make these functions templates on the number of space
  // dimensions to conform to usual practice when using deal.II, although we
  // will only use them for two space dimensions.
  //
  // The first of these functions just generates a triangulation of a circle
  // (hyperball) and outputs the Qp mapping of its cells for different values
  // of <code>p</code>. Then, we refine the grid once and do so again.
  template <int dim>
  void gnuplot_output()
  {
    std::cout << "Output of grids into gnuplot files:" << std::endl
              << "===================================" << std::endl;

    // So first generate a coarse triangulation of the circle and associate a
    // suitable boundary description to it. Note that the default values of
    // the HyperBallBoundary constructor are a center at the origin and a
    // radius equals one.
    Triangulation<dim> triangulation;
    GridGenerator::hyper_ball (triangulation);
    static const HyperBallBoundary<dim> boundary;
    triangulation.set_boundary (0, boundary);

    // Next generate output for this grid and for a once refined grid. Note
    // that we have hidden the mesh refinement in the loop header, which might
    // be uncommon but nevertheless works. Also it is strangely consistent
    // with incrementing the loop index denoting the refinement level.
    for (unsigned int refinement=0; refinement<2;
         ++refinement, triangulation.refine_global(1))
      {
        std::cout << "Refinement level: " << refinement << std::endl;

        // Then have a string which denotes the base part of the names of the
        // files into which we write the output. Note that in the parentheses
        // in the initializer we do arithmetic on characters, which assumes
        // that first the characters denoting numbers are placed consecutively
        // (which is probably true for all reasonable character sets
        // nowadays), but also assumes that the increment
        // <code>refinement</code> is less than ten. This is therefore more a
        // quick hack if we know exactly the values which the increment can
        // assume. A better implementation would use the
        // <code>std::istringstream</code> class to generate a name.
        std::string filename_base = "ball";
        filename_base += '0'+refinement;

        // Then output the present grid for Q1, Q2, and Q3 mappings:
        for (unsigned int degree=1; degree<4; ++degree)
          {
            std::cout << "Degree = " << degree << std::endl;

            // For this, first set up an object describing the mapping. This
            // is done using the <code>MappingQ</code> class, which takes as
            // argument to the constructor the polynomial degree which it
            // shall use.
            const MappingQ<dim> mapping (degree);
            // We note one interesting fact: if you want a piecewise linear
            // mapping, then you could give a value of <code>1</code> to the
            // constructor. However, for linear mappings, so many things can
            // be generated simpler that there is another class, called
            // <code>MappingQ1</code> which does exactly the same is if you
            // gave an degree of <code>1</code> to the <code>MappingQ</code>
            // class, but does so significantly faster. <code>MappingQ1</code>
            // is also the class that is implicitly used throughout the
            // library in many functions and classes if you do not specify
            // another mapping explicitly.


            // In degree to actually write out the present grid with this
            // mapping, we set up an object which we will use for output. We
            // will generate Gnuplot output, which consists of a set of lines
            // describing the mapped triangulation. By default, only one line
            // is drawn for each face of the triangulation, but since we want
            // to explicitly see the effect of the mapping, we want to have
            // the faces in more detail. This can be done by passing the
            // output object a structure which contains some flags. In the
            // present case, since Gnuplot can only draw straight lines, we
            // output a number of additional points on the faces so that each
            // face is drawn by 30 small lines instead of only one. This is
            // sufficient to give us the impression of seeing a curved line,
            // rather than a set of straight lines.
            GridOut grid_out;
            GridOutFlags::Gnuplot gnuplot_flags(false, 30);
            grid_out.set_flags(gnuplot_flags);

            // Finally, generate a filename and a file for output using the
            // same evil hack as above:
            std::string filename = filename_base+"_mapping_q";
            filename += ('0'+degree);
            filename += ".dat";
            std::ofstream gnuplot_file (filename.c_str());

            // Then write out the triangulation to this file. The last
            // argument of the function is a pointer to a mapping object. This
            // argument has a default value, and if no value is given a simple
            // <code>MappingQ1</code> object is taken, which we briefly
            // described above. This would then result in a piecewise linear
            // approximation of the true boundary in the output.
            grid_out.write_gnuplot (triangulation, gnuplot_file, &mapping);
          }
        std::cout << std::endl;
      }
  }

  // Now we proceed with the main part of the code, the approximation of
  // $\pi$. The area of a circle is of course given by $\pi r^2$, so having a
  // circle of radius 1, the area represents just the number that is searched
  // for. The numerical computation of the area is performed by integrating
  // the constant function of value 1 over the whole computational domain,
  // i.e. by computing the areas $\int_K 1 dx=\int_{\hat K} 1
  // \ \textrm{det}\ J(\hat x) d\hat x \approx \sum_i \textrm{det}
  // \ J(\hat x_i)w(\hat x_i)$,
  // where the sum extends over all quadrature points on all active cells in
  // the triangulation, with $w(x_i)$ being the weight of quadrature point
  // $x_i$. The integrals on each cell are approximated by numerical
  // quadrature, hence the only additional ingredient we need is to set up a
  // FEValues object that provides the corresponding `JxW' values of each
  // cell. (Note that `JxW' is meant to abbreviate <code>Jacobian determinant
  // times weight</code>; since in numerical quadrature the two factors always
  // occur at the same places, we only offer the combined quantity, rather
  // than two separate ones.) We note that here we won't use the FEValues
  // object in its original purpose, i.e. for the computation of values of
  // basis functions of a specific finite element at certain quadrature
  // points. Rather, we use it only to gain the `JxW' at the quadrature
  // points, irrespective of the (dummy) finite element we will give to the
  // constructor of the FEValues object. The actual finite element given to
  // the FEValues object is not used at all, so we could give any.
  template <int dim>
  void compute_pi_by_area ()
  {
    std::cout << "Computation of Pi by the area:" << std::endl
              << "==============================" << std::endl;

    // For the numerical quadrature on all cells we employ a quadrature rule
    // of sufficiently high degree. We choose QGauss that is of order 8 (4
    // points), to be sure that the errors due to numerical quadrature are of
    // higher order than the order (maximal 6) that will occur due to the
    // order of the approximation of the boundary, i.e. the order of the
    // mappings employed. Note that the integrand, the Jacobian determinant,
    // is not a polynomial function (rather, it is a rational one), so we do
    // not use Gauss quadrature in order to get the exact value of the
    // integral as done often in finite element computations, but could as
    // well have used any quadrature formula of like order instead.
    const QGauss<dim> quadrature(4);

    // Now start by looping over polynomial mapping degrees=1..4:
    for (unsigned int degree=1; degree<5; ++degree)
      {
        std::cout << "Degree = " << degree << std::endl;

        // First generate the triangulation, the boundary and the mapping
        // object as already seen.
        Triangulation<dim> triangulation;
        GridGenerator::hyper_ball (triangulation);

        static const HyperBallBoundary<dim> boundary;
        triangulation.set_boundary (0, boundary);

        const MappingQ<dim> mapping (degree);

        // We now create a dummy finite element. Here we could choose any
        // finite element, as we are only interested in the `JxW' values
        // provided by the FEValues object below. Nevertheless, we have to
        // provide a finite element since in this example we abuse the
        // FEValues class a little in that we only ask it to provide us with
        // the weights of certain quadrature points, in contrast to the usual
        // purpose (and name) of the FEValues class which is to provide the
        // values of finite elements at these points.
        const FE_Q<dim>     dummy_fe (1);

        // Likewise, we need to create a DoFHandler object. We do not actually
        // use it, but it will provide us with `active_cell_iterators' that
        // are needed to reinitialize the FEValues object on each cell of the
        // triangulation.
        DoFHandler<dim> dof_handler (triangulation);

        // Now we set up the FEValues object, giving the Mapping, the dummy
        // finite element and the quadrature object to the constructor,
        // together with the update flags asking for the `JxW' values at the
        // quadrature points only. This tells the FEValues object that it
        // needs not compute other quantities upon calling the
        // <code>reinit</code> function, thus saving computation time.
        //
        // The most important difference in the construction of the FEValues
        // object compared to previous example programs is that we pass a
        // mapping object as first argument, which is to be used in the
        // computation of the mapping from unit to real cell. In previous
        // examples, this argument was omitted, resulting in the implicit use
        // of an object of type MappingQ1.
        FEValues<dim> fe_values (mapping, dummy_fe, quadrature,
                                 update_JxW_values);

        // We employ an object of the ConvergenceTable class to store all
        // important data like the approximated values for $\pi$ and the error
        // with respect to the true value of $\pi$. We will also use functions
        // provided by the ConvergenceTable class to compute convergence rates
        // of the approximations to $\pi$.
        ConvergenceTable table;

        // Now we loop over several refinement steps of the triangulation.
        for (unsigned int refinement=0; refinement<6;
             ++refinement, triangulation.refine_global (1))
          {
            // In this loop we first add the number of active cells of the
            // current triangulation to the table. This function automatically
            // creates a table column with superscription `cells', in case
            // this column was not created before.
            table.add_value("cells", triangulation.n_active_cells());

            // Then we distribute the degrees of freedom for the dummy finite
            // element. Strictly speaking we do not need this function call in
            // our special case but we call it to make the DoFHandler happy --
            // otherwise it would throw an assertion in the FEValues::reinit
            // function below.
            dof_handler.distribute_dofs (dummy_fe);

            // We define the variable area as `long double' like we did for
            // the pi variable before.
            long double area = 0;

            // Now we loop over all cells, reinitialize the FEValues object
            // for each cell, and add up all the `JxW' values for this cell to
            // `area'...
            typename DoFHandler<dim>::active_cell_iterator
            cell = dof_handler.begin_active(),
            endc = dof_handler.end();
            for (; cell!=endc; ++cell)
              {
                fe_values.reinit (cell);
                for (unsigned int i=0; i<fe_values.n_quadrature_points; ++i)
                  area += fe_values.JxW (i);
              };

            // ...and store the resulting area values and the errors in the
            // table. We need a static cast to double as there is no
            // add_value(string, long double) function implemented. Note that
            // this also concerns the second call as the <code>fabs</code>
            // function in the <code>std</code> namespace is overloaded on its
            // argument types, so there exists a version taking and returning
            // a <code>long double</code>, in contrast to the global namespace
            // where only one such function is declared (which takes and
            // returns a double).
            table.add_value("eval.pi", static_cast<double> (area));
            table.add_value("error",   static_cast<double> (std::fabs(area-pi)));
          };

        // We want to compute the convergence rates of the `error'
        // column. Therefore we need to omit the other columns from the
        // convergence rate evaluation before calling
        // `evaluate_all_convergence_rates'
        table.omit_column_from_convergence_rate_evaluation("cells");
        table.omit_column_from_convergence_rate_evaluation("eval.pi");
        table.evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);

        // Finally we set the precision and scientific mode for output of some
        // of the quantities...
        table.set_precision("eval.pi", 16);
        table.set_scientific("error", true);

        // ...and write the whole table to std::cout.
        table.write_text(std::cout);

        std::cout << std::endl;
      };
  }


  // The following, second function also computes an approximation of $\pi$
  // but this time via the perimeter $2\pi r$ of the domain instead of the
  // area. This function is only a variation of the previous function. So we
  // will mainly give documentation for the differences.
  template <int dim>
  void compute_pi_by_perimeter ()
  {
    std::cout << "Computation of Pi by the perimeter:" << std::endl
              << "===================================" << std::endl;

    // We take the same order of quadrature but this time a `dim-1'
    // dimensional quadrature as we will integrate over (boundary) lines
    // rather than over cells.
    const QGauss<dim-1> quadrature(4);

    // We loop over all degrees, create the triangulation, the boundary, the
    // mapping, the dummy finite element and the DoFHandler object as seen
    // before.
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

        // Then we create a FEFaceValues object instead of a FEValues object
        // as in the previous function. Again, we pass a mapping as first
        // argument.
        FEFaceValues<dim> fe_face_values (mapping, fe, quadrature,
                                          update_JxW_values);
        ConvergenceTable table;

        for (unsigned int refinement=0; refinement<6;
             ++refinement, triangulation.refine_global (1))
          {
            table.add_value("cells", triangulation.n_active_cells());

            dof_handler.distribute_dofs (fe);

            // Now we run over all cells and over all faces of each cell. Only
            // the contributions of the `JxW' values on boundary faces are
            // added to the long double variable `perimeter'.
            typename DoFHandler<dim>::active_cell_iterator
            cell = dof_handler.begin_active(),
            endc = dof_handler.end();
            long double perimeter = 0;
            for (; cell!=endc; ++cell)
              for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
                if (cell->face(face_no)->at_boundary())
                  {
                    // We reinit the FEFaceValues object with the cell
                    // iterator and the number of the face.
                    fe_face_values.reinit (cell, face_no);
                    for (unsigned int i=0; i<fe_face_values.n_quadrature_points; ++i)
                      perimeter += fe_face_values.JxW (i);
                  };
            // Then store the evaluated values in the table...
            table.add_value("eval.pi", static_cast<double> (perimeter/2.));
            table.add_value("error",   static_cast<double> (std::fabs(perimeter/2.-pi)));
          };

        // ...and end this function as we did in the previous one:
        table.omit_column_from_convergence_rate_evaluation("cells");
        table.omit_column_from_convergence_rate_evaluation("eval.pi");
        table.evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);

        table.set_precision("eval.pi", 16);
        table.set_scientific("error", true);

        table.write_text(std::cout);

        std::cout << std::endl;
      };
  }
}


// The following main function just calls the above functions in the order of
// their appearance. Apart from this, it looks just like the main functions of
// previous tutorial programs.
int main ()
{
  try
    {
      std::cout.precision (16);

      Step10::gnuplot_output<2>();

      Step10::compute_pi_by_area<2> ();
      Step10::compute_pi_by_perimeter<2> ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}

/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 2001, 2002 */

/*    $Id$       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 2001, 2002, 2003, 2004 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */


				 // As in all programs, we start with
				 // a list of include files from the
				 // library, and as usual they are in
				 // the standard order which is
				 // ``base'' - ``lac'' - ``grid'' -
				 // ``dofs'' - ``fe'' - ``numerics''
				 // (as each of these categories
				 // roughly builds upon previous
				 // ones), then C++ standard headers:
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <base/table_handler.h>
#include <base/thread_management.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_refinement.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_constraints.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <numerics/error_estimator.h>

				 // Now for the C++ standard headers:
#include <iostream>
#include <fstream>
#include <list>

				 // Just as in the step-5 example
				 // program (see there for a lengthy
				 // discussion of the subject), we
				 // have to work around some
				 // historical confusion with the
				 // files declaring the stringstream
				 // classes:
#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif


				 // @sect3{Evaluation of the solution}

				 // As for the program itself, we
				 // first define classes that evaluate
				 // the solutions of a Laplace
				 // equation. In fact, they can
				 // evaluate every kind of solution,
				 // as long as it is described by a
				 // ``DoFHandler'' object, and a
				 // solution vector. We define them
				 // here first, even before the
				 // classes that actually generate the
				 // solution to be evaluated, since we
				 // need to declare an abstract base
				 // class that the solver classes can
				 // refer to.
				 //
				 // From an abstract point of view, we
				 // declare a pure base class
				 // that provides an evaluation
				 // operator ``operator()'' which will
				 // do the evaluation of the solution
				 // (whatever derived classes might
				 // consider an ``evaluation''). Since
				 // this is the only real function of
				 // this base class (except for some
				 // bookkeeping machinery), one
				 // usually terms such a class that
				 // only has an ``operator()'' a
				 // ``functor'' in C++ terminology,
				 // since it is used just like a
				 // function object.
				 //
				 // Objects of this functor type will
				 // then later be passed to the solver
				 // object, which applies it to the
				 // solution just computed. The
				 // evaluation objects may then
				 // extract any quantity they like
				 // from the solution. The advantage
				 // of putting these evaluation
				 // functions into a separate
				 // hierarchy of classes is that by
				 // design they cannot use the
				 // internals of the solver object and
				 // are therefore independent of
				 // changes to the way the solver
				 // works. Furthermore, it is trivial
				 // to write another evaluation class
				 // without modifying the solver
				 // class, which speeds up programming
				 // (not being able to use internals
				 // of another class also means that
				 // you do not have to worry about
				 // them -- programming evaluators is
				 // usually a rather quickly done
				 // task), as well as compilation (if
				 // solver and evaluation classes are
				 // put into different files: the
				 // solver only needs to see the
				 // declaration of the abstract base
				 // class, and therefore does not need
				 // to be recompiled upon addition of
				 // a new evaluation class, or
				 // modification of an old one).
				 // On a related note, you can reuse
				 // the evaluation classes for other
				 // projects, solving different
				 // equations.
				 //
				 // In order to improve separation of
				 // code into different modules, we
				 // put the evaluation classes into a
				 // namespace of their own. This makes
				 // it easier to actually solve
				 // different equations in the same
				 // program, by assembling it from
				 // existing building blocks. The
				 // reason for this is that classes
				 // for similar purposes tend to have
				 // the same name, although they were
				 // developed in different
				 // contexts. In order to be able to
				 // use them together in one program,
				 // it is necessary that they are
				 // placed in different
				 // namespaces. This we do here:
namespace Evaluation
{

				   // Now for the abstract base class
				   // of evaluation classes: its main
				   // purpose is to declare a pure
				   // virtual function ``operator()''
				   // taking a ``DoFHandler'' object,
				   // and the solution vector. In
				   // order to be able to use pointers
				   // to this base class only, it also
				   // has to declare a virtual
				   // destructor, which however does
				   // nothing. Besides this, it only
				   // provides for a little bit of
				   // bookkeeping: since we usually
				   // want to evaluate solutions on
				   // subsequent refinement levels, we
				   // store the number of the present
				   // refinement cycle, and provide a
				   // function to change this number.
  template <int dim>
  class EvaluationBase 
  {
    public:
      virtual ~EvaluationBase ();

      void set_refinement_cycle (const unsigned int refinement_cycle);
      
      virtual void operator () (const DoFHandler<dim> &dof_handler,
				const Vector<double>  &solution) const = 0;
    protected:
      unsigned int refinement_cycle;
  };


				   // After the declaration has been
				   // discussed above, the
				   // implementation is rather
				   // straightforward:
  template <int dim>
  EvaluationBase<dim>::~EvaluationBase ()
  {}
  

  
  template <int dim>
  void
  EvaluationBase<dim>::set_refinement_cycle (const unsigned int step)
  {
    refinement_cycle = step;
  }


				   // @sect4{Point evaluation}

				   // The next thing is to implement
				   // actual evaluation classes. As
				   // noted in the introduction, we'd
				   // like to extract a point value
				   // from the solution, so the first
				   // class does this in its
				   // ``operator()''. The actual point
				   // is given to this class through
				   // the constructor, as well as a
				   // table object into which it will
				   // put its findings.
				   //
				   // Finding out the value of a
				   // finite element field at an
				   // arbitrary point is rather
				   // difficult, if we cannot rely on
				   // knowing the actual finite
				   // element used, since then we
				   // cannot, for example, interpolate
				   // between nodes. For simplicity,
				   // we therefore assume here that
				   // the point at which we want to
				   // evaluate the field is actually a
				   // node. If, in the process of
				   // evaluating the solution, we find
				   // that we did not encounter this
				   // point upon looping over all
				   // vertices, we then have to throw
				   // an exception in order to signal
				   // to the calling functions that
				   // something has gone wrong, rather
				   // than silently ignore this error.
				   //
				   // In the step-9 example program,
				   // we have already seen how such an
				   // exception class can be declared,
				   // using the ``DeclExceptionN''
				   // macros. We use this mechanism
				   // here again.
				   //
				   // From this, the actual
				   // declaration of this class should
				   // be evident. Note that of course
				   // even if we do not list a
				   // destructor explicitely, an
				   // implicit destructor is generated
				   // from the compiler, and it is
				   // virtual just as the one of the
				   // base class.
  template <int dim>
  class PointValueEvaluation : public EvaluationBase<dim>
  {
    public:
      PointValueEvaluation (const Point<dim>   &evaluation_point,
			    TableHandler       &results_table);
      
      virtual void operator () (const DoFHandler<dim> &dof_handler,
				const Vector<double>  &solution) const;
      
      DeclException1 (ExcEvaluationPointNotFound,
		      Point<dim>,
		      << "The evaluation point " << arg1
		      << " was not found among the vertices of the present grid.");
    private:
      const Point<dim>  evaluation_point;
      TableHandler     &results_table;
  };


				   // As for the definition, the
				   // constructor is trivial, just
				   // taking data and storing it in
				   // object-local ones:
  template <int dim>
  PointValueEvaluation<dim>::
  PointValueEvaluation (const Point<dim>   &evaluation_point,
			TableHandler       &results_table)
		  :
		  evaluation_point (evaluation_point),
		  results_table (results_table)
  {}
  


				   // Now for the function that is
				   // mainly of interest in this
				   // class, the computation of the
				   // point value:
  template <int dim>
  void
  PointValueEvaluation<dim>::
  operator () (const DoFHandler<dim> &dof_handler,
	       const Vector<double>  &solution) const 
  {
				     // First allocate a variable that
				     // will hold the point
				     // value. Initialize it with a
				     // value that is clearly bogus,
				     // so that if we fail to set it
				     // to a reasonable value, we will
				     // note at once. This may not be
				     // necessary in a function as
				     // small as this one, since we
				     // can easily see all possible
				     // paths of execution here, but
				     // it proved to be helpful for
				     // more complex cases, and so we
				     // employ this strategy here as
				     // well.
    double point_value = 1e20;

				     // Then loop over all cells and
				     // all their vertices, and check
				     // whether a vertex matches the
				     // evaluation point. If this is
				     // the case, then extract the
				     // point value, set a flag that
				     // we have found the point of
				     // interest, and exit the loop.
    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
    bool evaluation_point_found = false;
    for (; (cell!=endc) && !evaluation_point_found; ++cell)
      for (unsigned int vertex=0;
	   vertex<GeometryInfo<dim>::vertices_per_cell;
	   ++vertex)
	if (cell->vertex(vertex) == evaluation_point)
	  {
					     // In order to extract
					     // the point value from
					     // the global solution
					     // vector, pick that
					     // component that belongs
					     // to the vertex of
					     // interest, and, in case
					     // the solution is
					     // vector-valued, take
					     // the first component of
					     // it:
	    point_value = solution(cell->vertex_dof_index(vertex,0));
					     // Note that by this we
					     // have made an
					     // assumption that is not
					     // valid always and
					     // should be documented
					     // in the class
					     // declaration if this
					     // were code for a real
					     // application rather
					     // than a tutorial
					     // program: we assume
					     // that the finite
					     // element used for the
					     // solution we try to
					     // evaluate actually has
					     // degrees of freedom
					     // associated with
					     // vertices. This, for
					     // example, does not hold
					     // for discontinuous
					     // elements, were the
					     // support points for the
					     // shape functions
					     // happen to be located
					     // at the vertices, but
					     // are not associated
					     // with the vertices but
					     // rather with the cell
					     // interior, since
					     // association with
					     // vertices would imply
					     // continuity there. It
					     // would also not hold
					     // for edge oriented
					     // elements, and the
					     // like.
					     //
					     // Ideally, we would
					     // check this at the
					     // beginning of the
					     // function, for example
					     // by a statement like
					     // ``Assert
					     // (dof_handler.get_fe().dofs_per_vertex
					     // > 0,
					     // ExcNotImplemented())'',
					     // which should make it
					     // quite clear what is
					     // going wrong when the
					     // exception is
					     // triggered. In this
					     // case, we omit it
					     // (which is indeed bad
					     // style), but knowing
					     // that that does not
					     // hurt here, since the
					     // statement
					     // ``cell->vertex_dof_index(vertex,0)''
					     // would fail if we asked
					     // it to give us the DoF
					     // index of a vertex if
					     // there were none.
					     //
					     // We stress again that
					     // this restriction on
					     // the allowed finite
					     // elements should be
					     // stated in the class
					     // documentation.

					     // Since we found the
					     // right point, we now
					     // set the respective
					     // flag and exit the
					     // innermost loop. The
					     // outer loop will the
					     // also be terminated due
					     // to the set flag.
	    evaluation_point_found = true;
	    break;
	  };

				     // Finally, we'd like to make
				     // sure that we have indeed found
				     // the evaluation point, since if
				     // that were not so we could not
				     // give a reasonable value of the
				     // solution there and the rest of
				     // the computations were useless
				     // anyway. So make sure through
				     // the ``AssertThrow'' macro
				     // already used in the step-9
				     // program that we have indeed
				     // found this point. If this is
				     // not so, the macro throws an
				     // exception of the type that is
				     // given to it as second
				     // argument, but compared to a
				     // straightforward ``throw''
				     // statement, it fills the
				     // exception object with a set of
				     // additional information, for
				     // example the source file and
				     // line number where the
				     // exception was generated, and
				     // the condition that failed. If
				     // you have a ``catch'' clause in
				     // your main function (as this
				     // program has), you will catch
				     // all exceptions that are not
				     // caught somewhere in between
				     // and thus already handled, and
				     // this additional information
				     // will help you find out what
				     // happened and where it went
				     // wrong.
    AssertThrow (evaluation_point_found,
		 ExcEvaluationPointNotFound(evaluation_point));
				     // Note that we have used the
				     // ``Assert'' macro in other
				     // example programs as well. It
				     // differed from the
				     // ``AssertThrow'' macro used
				     // here in that it simply aborts
				     // the program, rather than
				     // throwing an exception, and
				     // that it did so only in debug
				     // mode. It was the right macro
				     // to use to check about the size
				     // of vectors passed as arguments
				     // to functions, and the like.
				     //
				     // However, here the situation is
				     // different: whether we find the
				     // evaluation point or not may
				     // change from refinement to
				     // refinement (for example, if
				     // the four cells around point
				     // are coarsened away, then the
				     // point may vanish after
				     // refinement and
				     // coarsening). This is something
				     // that cannot be predicted from
				     // a few number of runs of the
				     // program in debug mode, but
				     // should be checked always, also
				     // in production runs. Thus the
				     // use of the ``AssertThrow''
				     // macro here.
    
				     // Now, if we are sure that we
				     // have found the evaluation
				     // point, we can add the results
				     // into the table of results:
    results_table.add_value ("DoFs", dof_handler.n_dofs());
    results_table.add_value ("u(x_0)", point_value);
  }




				   // @sect4{Generating output}

				   // A different, maybe slightly odd
				   // kind of ``evaluation'' of a
				   // solution is to output it to a
				   // file in a graphical
				   // format. Since in the evaluation
				   // functions we are given a
				   // ``DoFHandler'' object and the
				   // solution vector, we have all we
				   // need to do this, so we can do it
				   // in an evaluation class. The
				   // reason for actually doing so
				   // instead of putting it into the
				   // class that computed the solution
				   // is that this way we have more
				   // flexibility: if we choose to
				   // only output certain aspects of
				   // it, or not output it at all. In
				   // any case, we do not need to
				   // modify the solver class, we just
				   // have to modify one of the
				   // modules out of which we build
				   // this program. This form of
				   // encapsulation, as above, helps
				   // us to keep each part of the
				   // program rather simple as the
				   // interfaces are kept simple, and
				   // no access to hidden data is
				   // possible.
				   //
				   // Since this class which generates
				   // the output is derived from the
				   // common ``EvaluationBase'' base
				   // class, its main interface is the
				   // ``operator()''
				   // function. Furthermore, it has a
				   // constructor taking a string that
				   // will be used as the base part of
				   // the file name to which output
				   // will be sent (we will augment it
				   // by a number indicating the
				   // number of the refinement cycle
				   // -- the base class has this
				   // information at hand --, and a
				   // suffix), and the constructor
				   // also takes a value that
				   // indicates which format is
				   // requested, i.e. for which
				   // graphics program we shall
				   // generate output (from this we
				   // will then also generate the
				   // suffix of the filename to which
				   // we write).
				   //
				   // Regarding the output format, the
				   // ``DataOutInterface'' class
				   // (which is a base class of
				   // ``DataOut'' through which we
				   // will access its fields) provides
				   // an enumeration field
				   // ``OutputFormat'', which lists
				   // names for all supported output
				   // formats. At the time of writing
				   // of this program, the supported
				   // graphics formats are represented
				   // by the enum values ``ucd'',
				   // ``gnuplot'', ``povray'',
				   // ``eps'', ``gmv'', ``tecplot'',
				   // ``tecplot_binary'', ``dx'', and
				   // ``vtk'', but this list will
				   // certainly grow over time. Now,
				   // within various functions of that
				   // base class, you can use values
				   // of this type to get information
				   // about these graphics formats
				   // (for example the default suffix
				   // used for files of each format),
				   // and you can call a generic
				   // ``write'' function, which then
				   // branches to the
				   // ``write_gnuplot'',
				   // ``write_ucd'', etc functions
				   // which we have used in previous
				   // examples already, based on the
				   // value of a second argument given
				   // to it denoting the required
				   // output format. This mechanism
				   // makes it simple to write an
				   // extensible program that can
				   // decide which output format to
				   // use at runtime, and it also
				   // makes it rather simple to write
				   // the program in a way such that
				   // it takes advantage of newly
				   // implemented output formats,
				   // without the need to change the
				   // application program.
				   //
				   // Of these two fields, the base
				   // name and the output format
				   // descriptor, the constructor
				   // takes values and stores them for
				   // later use by the actual
				   // evaluation function.
  template <int dim>
  class SolutionOutput : public EvaluationBase<dim>
  {
    public:
      SolutionOutput (const std::string                         &output_name_base,
		      const typename DataOut<dim>::OutputFormat  output_format);
      
      virtual void operator () (const DoFHandler<dim> &dof_handler,
				const Vector<double>  &solution) const;
    private:
      const std::string                         output_name_base;
      const typename DataOut<dim>::OutputFormat output_format;
  };


  template <int dim>
  SolutionOutput<dim>::
  SolutionOutput (const std::string                         &output_name_base,
		  const typename DataOut<dim>::OutputFormat  output_format)
		  :
		  output_name_base (output_name_base),
		  output_format (output_format)
  {}
  

				   // After the description above, the
				   // function generating the actual
				   // output is now relatively
				   // straightforward. The only
				   // particularly interesting feature
				   // over previous example programs
				   // is the use of the
				   // ``DataOut::default_suffix''
				   // function, returning the usual
				   // suffix for files of a given
				   // format (e.g. ".eps" for
				   // encapsulated postscript files,
				   // ".gnuplot" for Gnuplot files),
				   // and of the generic
				   // ``DataOut::write'' function with
				   // a second argument, which
				   // branches to the actual output
				   // functions for the different
				   // graphics formats, based on the
				   // value of the format descriptor
				   // passed as second argument.
				   //
				   // The somewhat complicated use of
				   // the stringstream class,
				   // involving support from the
				   // preprocessor, is as already
				   // explained in the step-5 example
				   // program.
				   //
				   // Also note that we have to prefix
				   // ``this->'' to access a member
				   // variable of the template
				   // dependent base class. The reason
				   // here, and further down in the
				   // program is the same as the one
				   // described in the step-7 example
				   // program (look for ``two-stage
				   // name lookup'' there).
  template <int dim>
  void
  SolutionOutput<dim>::operator () (const DoFHandler<dim> &dof_handler,
				    const Vector<double>  &solution) const
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "solution");
    data_out.build_patches ();
  
#ifdef HAVE_STD_STRINGSTREAM
    std::ostringstream filename;
#else
    std::ostrstream filename;
#endif
    filename << output_name_base << "-"
	     << this->refinement_cycle
	     << data_out.default_suffix (output_format)
	     << std::ends;
#ifdef HAVE_STD_STRINGSTREAM
    std::ofstream out (filename.str().c_str());
#else
    std::ofstream out (filename.str());
#endif
    
    data_out.write (out, output_format);
  }



				   // @sect4{Other evaluations}
  
				   // In practical applications, one
				   // would add here a list of other
				   // possible evaluation classes,
				   // representing quantities that one
				   // may be interested in. For this
				   // example, that much shall be
				   // sufficient, so we close the
				   // namespace.
}

  
				 // @sect3{The Laplace solver classes}

				 // After defining what we want to
				 // know of the solution, we should
				 // now care how to get at it. We will
				 // pack everything we need into a
				 // namespace of its own, for much the
				 // same reasons as for the
				 // evaluations above.
				 //
				 // Since we have discussed Laplace
				 // solvers already in considerable
				 // detail in previous examples, there
				 // is not much new stuff
				 // following. Rather, we have to a
				 // great extent cannibalized previous
				 // examples and put them, in slightly
				 // different form, into this example
				 // program. We will therefore mostly
				 // be concerned with discussing the
				 // differences to previous examples.
				 //
				 // Basically, as already said in the
				 // introduction, the lack of new
				 // stuff in this example is
				 // deliberate, as it is more to
				 // demonstrate software design
				 // practices, rather than
				 // mathematics. The emphasis in
				 // explanations below will therefore
				 // be more on the actual
				 // implementation.
namespace LaplaceSolver
{
				   // @sect4{An abstract base class}

				   // In defining a Laplace solver, we
				   // start out by declaring an
				   // abstract base class, that has no
				   // functionality itself except for
				   // taking and storing a pointer to
				   // the triangulation to be used
				   // later.
				   //
				   // This base class is very general,
				   // and could as well be used for
				   // any other stationary problem. It
				   // provides declarations of
				   // functions that shall, in derived
				   // classes, solve a problem,
				   // postprocess the solution with a
				   // list of evaluation objects, and
				   // refine the grid,
				   // respectively. None of these
				   // functions actually does
				   // something itself in the base
				   // class.
				   //
				   // Due to the lack of actual
				   // functionality, the programming
				   // style of declaring very abstract
				   // base classes reminds of the
				   // style used in Smalltalk or Java
				   // programs, where all classes are
				   // derived from entirely abstract
				   // classes ``Object'', even number
				   // representations. The author
				   // admits that he does not
				   // particularly like the use of
				   // such a style in C++, as it puts
				   // style over reason. Furthermore,
				   // it promotes the use of virtual
				   // functions for everything (for
				   // example, in Java, all functions
				   // are virtual per se), which,
				   // however, has proven to be rather
				   // inefficient in many applications
				   // where functions are often only
				   // accessing data, not doing
				   // computations, and therefore
				   // quickly return; the overhead of
				   // virtual functions can then be
				   // significant. The opinion of the
				   // author is to have abstract base
				   // classes wherever at least some
				   // part of the code of actual
				   // implementations can be shared
				   // and thus separated into the base
				   // class.
				   //
				   // Besides all these theoretical
				   // questions, we here have a good
				   // reason, which will become
				   // clearer to the reader
				   // below. Basically, we want to be
				   // able to have a family of
				   // different Laplace solvers that
				   // differ so much that no larger
				   // common subset of functionality
				   // could be found. We therefore
				   // just declare such an abstract
				   // base class, taking a pointer to
				   // a triangulation in the
				   // constructor and storing it
				   // henceforth. Since this
				   // triangulation will be used
				   // throughout all computations, we
				   // have to make sure that the
				   // triangulation exists until the
				   // destructor exits. We do this by
				   // keeping a ``SmartPointer'' to
				   // this triangulation, which uses a
				   // counter in the triangulation
				   // class to denote the fact that
				   // there is still an object out
				   // there using this triangulation,
				   // thus leading to an abort in case
				   // the triangulation is attempted
				   // to be destructed while this
				   // object still uses it.
				   //
				   // Note that while the pointer
				   // itself is declared constant
				   // (i.e. throughout the lifetime of
				   // this object, the pointer points
				   // to the same object), it is not
				   // declared as a pointer to a
				   // constant triangulation. In fact,
				   // by this we allow that derived
				   // classes refine or coarsen the
				   // triangulation within the
				   // ``refine_grid'' function.
				   //
				   // Finally, we have a function
				   // ``n_dofs'' is only a tool for
				   // the driver functions to decide
				   // whether we want to go on with
				   // mesh refinement or not. It
				   // returns the number of degrees of
				   // freedom the present simulation
				   // has.
  template <int dim>
  class Base
  {
    public:
      Base (Triangulation<dim> &coarse_grid);
      virtual ~Base ();

      virtual void solve_problem () = 0;
      virtual void postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const = 0;
      virtual void refine_grid () = 0;
      virtual unsigned int n_dofs () const = 0;
      
    protected:
      const SmartPointer<Triangulation<dim> > triangulation;
  };


				   // The implementation of the only
				   // two non-abstract functions is
				   // then rather boring:
  template <int dim>
  Base<dim>::Base (Triangulation<dim> &coarse_grid)
		  :
		  triangulation (&coarse_grid)
  {}


  template <int dim>
  Base<dim>::~Base () 
  {}
  

				   // @sect4{A general solver class}

				   // Following now the main class
				   // that implements assembling the
				   // matrix of the linear system,
				   // solving it, and calling the
				   // postprocessor objects on the
				   // solution. It implements the
				   // ``solve_problem'' and
				   // ``postprocess'' functions
				   // declared in the base class. It
				   // does not, however, implement the
				   // ``refine_grid'' method, as mesh
				   // refinement will be implemented
				   // in a number of derived classes.
				   //
				   // It also declares a new abstract
				   // virtual function,
				   // ``assemble_rhs'', that needs to
				   // be overloaded in subclasses. The
				   // reason is that we will implement
				   // two different classes that will
				   // implement different methods to
				   // assemble the right hand side
				   // vector. This function might also
				   // be interesting in cases where
				   // the right hand side depends not
				   // simply on a continuous function,
				   // but on something else as well,
				   // for example the solution of
				   // another discretized problem,
				   // etc. The latter happens
				   // frequently in non-linear
				   // problems.
				   //
				   // As we mentioned previously, the
				   // actual content of this class is
				   // not new, but a mixture of
				   // various techniques already used
				   // in previous examples. We will
				   // therefore not discuss them in
				   // detail, but refer the reader to
				   // these programs.
				   //
				   // Basically, in a few words, the
				   // constructor of this class takes
				   // pointers to a triangulation, a
				   // finite element, and a function
				   // object representing the boundary
				   // values. These are either passed
				   // down to the base class's
				   // constructor, or are stored and
				   // used to generate a
				   // ``DoFHandler'' object
				   // later. Since finite elements and
				   // quadrature formula should match,
				   // it is also passed a quadrature
				   // object.
				   //
				   // The ``solve_problem'' sets up
				   // the data structures for the
				   // actual solution, calls the
				   // functions to assemble the linear
				   // system, and solves it.
				   //
				   // The ``postprocess'' function
				   // finally takes an evaluation
				   // object and applies it to the
				   // computed solution.
				   //
				   // The ``n_dofs'' function finally
				   // implements the pure virtual
				   // function of the base class.
  template <int dim>
  class Solver : public virtual Base<dim>
  {
    public:
      Solver (Triangulation<dim>       &triangulation,
	      const FiniteElement<dim> &fe,
	      const Quadrature<dim>    &quadrature,
	      const Function<dim>      &boundary_values);
      virtual
      ~Solver ();

      virtual
      void
      solve_problem ();

      virtual
      void
      postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const;

      virtual
      unsigned int
      n_dofs () const;
      
				       // In the protected section of
				       // this class, we first have a
				       // number of member variables,
				       // of which the use should be
				       // clear from the previous
				       // examples:
    protected:
      const SmartPointer<const FiniteElement<dim> >  fe;
      const SmartPointer<const Quadrature<dim> >     quadrature;
      DoFHandler<dim>                                dof_handler;
      Vector<double>                                 solution;
      const SmartPointer<const Function<dim> >       boundary_values;

				       // Then we declare an abstract
				       // function that will be used
				       // to assemble the right hand
				       // side. As explained above,
				       // there are various cases for
				       // which this action differs
				       // strongly in what is
				       // necessary, so we defer this
				       // to derived classes:
      virtual void assemble_rhs (Vector<double> &rhs) const = 0;
    
				       // Next, in the private
				       // section, we have a small
				       // class which represents an
				       // entire linear system, i.e. a
				       // matrix, a right hand side,
				       // and a solution vector, as
				       // well as the constraints that
				       // are applied to it, such as
				       // those due to hanging
				       // nodes. Its constructor
				       // initializes the various
				       // subobjects, and there is a
				       // function that implements a
				       // conjugate gradient method as
				       // solver.
    private:
      struct LinearSystem
      {
	  LinearSystem (const DoFHandler<dim> &dof_handler);

	  void solve (Vector<double> &solution) const;
	
	  ConstraintMatrix     hanging_node_constraints;
	  SparsityPattern      sparsity_pattern;
	  SparseMatrix<double> matrix;
	  Vector<double>       rhs;
      };

				       // Finally, there is a pair of
				       // functions which will be used
				       // to assemble the actual
				       // system matrix. It calls the
				       // virtual function assembling
				       // the right hand side, and
				       // installs a number threads
				       // each running the second
				       // function which assembles
				       // part of the system
				       // matrix. The mechanism for
				       // doing so is the same as in
				       // the step-9 example program.
      void
      assemble_linear_system (LinearSystem &linear_system);

      void
      assemble_matrix (LinearSystem                                         &linear_system,
		       const typename DoFHandler<dim>::active_cell_iterator &begin_cell,
		       const typename DoFHandler<dim>::active_cell_iterator &end_cell,
		       Threads::ThreadMutex                                 &mutex) const;
  };



				   // Now here comes the constructor
				   // of the class. It does not do
				   // much except store pointers to
				   // the objects given, and generate
				   // ``DoFHandler'' object
				   // initialized with the given
				   // pointer to a triangulation. This
				   // causes the DoF handler to store
				   // that pointer, but does not
				   // already generate a finite
				   // element numbering (we only ask
				   // for that in the
				   // ``solve_problem'' function).
  template <int dim>
  Solver<dim>::Solver (Triangulation<dim>       &triangulation,
		       const FiniteElement<dim> &fe,
		       const Quadrature<dim>    &quadrature,
		       const Function<dim>      &boundary_values)
		  :
		  Base<dim> (triangulation),
		  fe (&fe),
                  quadrature (&quadrature),
		  dof_handler (triangulation),
		  boundary_values (&boundary_values)
  {}


				   // The destructor is simple, it
				   // only clears the information
				   // stored in the DoF handler object
				   // to release the memory.
  template <int dim>
  Solver<dim>::~Solver () 
  {
    dof_handler.clear ();
  }


				   // The next function is the one
				   // which delegates the main work in
				   // solving the problem: it sets up
				   // the DoF handler object with the
				   // finite element given to the
				   // constructor of this object, the
				   // creates an object that denotes
				   // the linear system (i.e. the
				   // matrix, the right hand side
				   // vector, and the solution
				   // vector), calls the function to
				   // assemble it, and finally solves
				   // it:
  template <int dim>
  void
  Solver<dim>::solve_problem ()
  {
    dof_handler.distribute_dofs (*fe);
    solution.reinit (dof_handler.n_dofs());

    LinearSystem linear_system (dof_handler);
    assemble_linear_system (linear_system);
    linear_system.solve (solution);
  }


				   // As stated above, the
				   // ``postprocess'' function takes
				   // an evaluation object, and
				   // applies it to the computed
				   // solution. This function may be
				   // called multiply, once for each
				   // evaluation of the solution which
				   // the user required.
  template <int dim>
  void
  Solver<dim>::
  postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const
  {
    postprocessor (dof_handler, solution);
  }


				   // The ``n_dofs'' function should
				   // be self-explanatory:
  template <int dim>
  unsigned int
  Solver<dim>::n_dofs () const
  {
    return dof_handler.n_dofs();
  }
  

				   // The following function assembles
				   // matrix and right hand side of
				   // the linear system to be solved
				   // in each step. It goes along the
				   // same lines as used in previous
				   // examples, so we explain it only
				   // briefly:
  template <int dim>
  void
  Solver<dim>::assemble_linear_system (LinearSystem &linear_system)
  {
				     // First define a convenience
				     // abbreviation for these lengthy
				     // iterator names...
    typedef
      typename DoFHandler<dim>::active_cell_iterator
      active_cell_iterator;

				     // ... and use it to split up the
				     // set of cells into a number of
				     // pieces of equal size. The
				     // number of blocks is set to the
				     // default number of threads to
				     // be used, which by default is
				     // set to the number of
				     // processors found in your
				     // computer at startup of the
				     // program:
    const unsigned int n_threads = multithread_info.n_default_threads;
    std::vector<std::pair<active_cell_iterator,active_cell_iterator> >
      thread_ranges 
      = Threads::split_range<active_cell_iterator> (dof_handler.begin_active (),
						    dof_handler.end (),
						    n_threads);

				     // These ranges are then assigned
				     // to a number of threads which
				     // we create next. Each will
				     // assemble the local cell
				     // matrices on the assigned
				     // cells, and fill the matrix
				     // object with it. Since there is
				     // need for synchronization when
				     // filling the same matrix from
				     // different threads, we need a
				     // mutex here:
    Threads::ThreadMutex mutex;
    Threads::ThreadGroup<> threads;
    for (unsigned int thread=0; thread<n_threads; ++thread)
      threads += Threads::spawn (*this, &Solver<dim>::assemble_matrix)
                 (linear_system,
                  thread_ranges[thread].first,
                  thread_ranges[thread].second,
                  mutex);

				     // While the spawned threads
				     // assemble the system matrix, we
				     // can already compute the right
				     // hand side vector in the main
				     // thread, and condense away the
				     // constraints due to hanging
				     // nodes:
    assemble_rhs (linear_system.rhs);
    linear_system.hanging_node_constraints.condense (linear_system.rhs);

				     // And while we're already at it
				     // to compute things in parallel,
				     // interpolating boundary values
				     // is one more thing that can be
				     // done independently, so we do
				     // it here:
    std::map<unsigned int,double> boundary_value_map;
    VectorTools::interpolate_boundary_values (dof_handler,
					      0,
					      *boundary_values,
					      boundary_value_map);
    
    
				     // If this is done, wait for the
				     // matrix assembling threads, and
				     // condense the constraints in
				     // the matrix as well:
    threads.join_all ();
    linear_system.hanging_node_constraints.condense (linear_system.matrix);

				     // Now that we have the linear
				     // system, we can also treat
				     // boundary values, which need to
				     // be eliminated from both the
				     // matrix and the right hand
				     // side:
    MatrixTools::apply_boundary_values (boundary_value_map,
					linear_system.matrix,
					solution,
					linear_system.rhs);

  }


				   // The second of this pair of
				   // functions takes a range of cell
				   // iterators, and assembles the
				   // system matrix on this part of
				   // the domain. Since it's actions
				   // have all been explained in
				   // previous programs, we do not
				   // comment on it any more, except
				   // for one pointe below.
  template <int dim>
  void
  Solver<dim>::assemble_matrix (LinearSystem                                         &linear_system,
				const typename DoFHandler<dim>::active_cell_iterator &begin_cell,
				const typename DoFHandler<dim>::active_cell_iterator &end_cell,
				Threads::ThreadMutex                                 &mutex) const
  {
    FEValues<dim> fe_values (*fe, *quadrature, 
			     UpdateFlags(update_gradients |
					 update_JxW_values));

    const unsigned int   dofs_per_cell = fe->dofs_per_cell;
    const unsigned int   n_q_points    = quadrature->n_quadrature_points;

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    for (typename DoFHandler<dim>::active_cell_iterator cell=begin_cell;
	 cell!=end_cell; ++cell)
      {
	cell_matrix = 0;

	fe_values.reinit (cell);

	for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += (fe_values.shape_grad(i,q_point) *
				   fe_values.shape_grad(j,q_point) *
				   fe_values.JxW(q_point));


	cell->get_dof_indices (local_dof_indices);

                                         // In the step-9 program, we
                                         // have shown that you have
                                         // to use the mutex to lock
                                         // the matrix when copying
                                         // the elements from the
                                         // local to the global
                                         // matrix. This was necessary
                                         // to avoid that two threads
                                         // access it at the same
                                         // time, eventually
                                         // overwriting their
                                         // respective
                                         // work. Previously, we have
                                         // used the `acquire'' and
                                         // ``release'' functions of
                                         // the mutex to lock and
                                         // unlock the mutex,
                                         // respectively. While this
                                         // is valid, there is one
                                         // possible catch: if between
                                         // the locking operation and
                                         // the unlocking operation an
                                         // exception is thrown, the
                                         // mutex remains in the
                                         // locked state, and in some
                                         // cases this might lead to
                                         // deadlocks. A similar
                                         // situation arises, when one
                                         // changes the code to have a
                                         // return statement somewhere
                                         // in the middle of the
                                         // locked block, and forgets
                                         // that before we call
                                         // ``return'', we also have
                                         // to unlock the mutex. This
                                         // all is not be a problem
                                         // here, but we want to show
                                         // the general technique to
                                         // cope with these problems
                                         // nevertheless: have an
                                         // object that upon
                                         // initialization (i.e. in
                                         // its constructor) locks the
                                         // mutex, and on running the
                                         // destructor unlocks it
                                         // again. This is called the
                                         // ``scoped lock'' pattern
                                         // (apparently invented by
                                         // Doug Schmidt originally),
                                         // and it works because
                                         // destructors of local
                                         // objects are also run when
                                         // we exit the function
                                         // either through a
                                         // ``return'' statement, or
                                         // when an exception is
                                         // raised. Thus, it is
                                         // guaranteed that the mutex
                                         // will always be unlocked
                                         // when we exit this part of
                                         // the program, whether the
                                         // operation completed
                                         // successfully or not,
                                         // whether the exit path was
                                         // something we implemented
                                         // willfully or whether the
                                         // function was exited by an
                                         // exception that we did not
                                         // forsee.
                                         //
                                         // deal.II implements the
                                         // scoped locking pattern in
                                         // the
                                         // ThreadMutex::ScopedLock
                                         // class: it takes the mutex
                                         // in the constructor and
                                         // locks it; in its
                                         // destructor, it unlocks it
                                         // again. So here is how it
                                         // is used:
        Threads::ThreadMutex::ScopedLock lock (mutex);
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    linear_system.matrix.add (local_dof_indices[i],
				      local_dof_indices[j],
				      cell_matrix(i,j));
                                         // Here, at the brace, the
                                         // current scope ends, so the
                                         // ``lock'' variable goes out
                                         // of existence and its
                                         // destructor the mutex is
                                         // unlocked.
      };
  }


				   // Now for the functions that
				   // implement actions in the linear
				   // system class. First, the
				   // constructor initializes all data
				   // elements to their correct sizes,
				   // and sets up a number of
				   // additional data structures, such
				   // as constraints due to hanging
				   // nodes. Since setting up the
				   // hanging nodes and finding out
				   // about the nonzero elements of
				   // the matrix is independent, we do
				   // that in parallel (if the library
				   // was configured to use
				   // concurrency, at least;
				   // otherwise, the actions are
				   // performed sequentially). Note
				   // that we spawn only one thread,
				   // and do the second action in the
				   // main thread. Since only one
				   // thread is generated, we don't
				   // use the ``Threads::ThreadGroup''
				   // class here, but rather use the
				   // one created thread object
				   // directly to wait for this
				   // particular thread's exit.
				   //
				   // Note that taking up the address
				   // of the
				   // ``DoFTools::make_hanging_node_constraints''
				   // function is a little tricky,
				   // since there are actually three
				   // of them, one for each supported
				   // space dimension. Taking
				   // addresses of overloaded
				   // functions is somewhat
				   // complicated in C++, since the
				   // address-of operator ``&'' in
				   // that case returns more like a
				   // set of values (the addresses of
				   // all functions with that name),
				   // and selecting the right one is
				   // then the next step. If the
				   // context dictates which one to
				   // take (for example by assigning
				   // to a function pointer of known
				   // type), then the compiler can do
				   // that by itself, but if this set
				   // of pointers shall be given as
				   // the argument to a function that
				   // takes a template, the compiler
				   // could choose all without having
				   // a preference for one. We
				   // therefore have to make it clear
				   // to the compiler which one we
				   // would like to have; for this, we
				   // could use a cast, but for more
				   // clarity, we assign it to a
				   // temporary ``mhnc_p'' (short for
				   // ``pointer to
				   // make_hanging_node_constraints'')
				   // with the right type, and using
				   // this pointer instead.
  template <int dim>
  Solver<dim>::LinearSystem::
  LinearSystem (const DoFHandler<dim> &dof_handler)
  {
    hanging_node_constraints.clear ();

    void (*mhnc_p) (const DoFHandler<dim> &,
		    ConstraintMatrix      &)
      = &DoFTools::make_hanging_node_constraints;
    
    Threads::Thread<>
      mhnc_thread = Threads::spawn (mhnc_p)(dof_handler,
                                            hanging_node_constraints);

    sparsity_pattern.reinit (dof_handler.n_dofs(),
			     dof_handler.n_dofs(),
			     dof_handler.max_couplings_between_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

				     // Wait until the
				     // ``hanging_node_constraints''
				     // object is fully set up, then
				     // close it and use it to
				     // condense the sparsity pattern:
    mhnc_thread.join ();
    hanging_node_constraints.close ();
    hanging_node_constraints.condense (sparsity_pattern);

				     // Finally, close the sparsity
				     // pattern, initialize the
				     // matrix, and set the right hand
				     // side vector to the right size.
    sparsity_pattern.compress();
    matrix.reinit (sparsity_pattern);
    rhs.reinit (dof_handler.n_dofs());
  }



				   // The second function of this
				   // class simply solves the linear
				   // system by a preconditioned
				   // conjugate gradient method. This
				   // has been extensively discussed
				   // before, so we don't dwell into
				   // it any more.
  template <int dim>
  void
  Solver<dim>::LinearSystem::solve (Vector<double> &solution) const
  {
    SolverControl           solver_control (1000, 1e-12);
    PrimitiveVectorMemory<> vector_memory;
    SolverCG<>              cg (solver_control, vector_memory);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(matrix, 1.2);

    cg.solve (matrix, solution, rhs, preconditioner);

    hanging_node_constraints.distribute (solution);
  }




				   // @sect4{A primal solver}

				   // In the previous section, a base
				   // class for Laplace solvers was
				   // implemented, that lacked the
				   // functionality to assemble the
				   // right hand side vector, however,
				   // for reasons that were explained
				   // there. Now we implement a
				   // corresponding class that can do
				   // this for the case that the right
				   // hand side of a problem is given
				   // as a function object.
				   //
				   // The actions of the class are
				   // rather what you have seen
				   // already in previous examples
				   // already, so a brief explanation
				   // should suffice: the constructor
				   // takes the same data as does that
				   // of the underlying class (to
				   // which it passes all information)
				   // except for one function object
				   // that denotes the right hand side
				   // of the problem. A pointer to
				   // this object is stored (again as
				   // a ``SmartPointer'', in order to
				   // make sure that the function
				   // object is not deleted as long as
				   // it is still used by this class).
				   //
				   // The only functional part of this
				   // class is the ``assemble_rhs''
				   // method that does what its name
				   // suggests.
  template <int dim>
  class PrimalSolver : public Solver<dim>
  {
    public:
      PrimalSolver (Triangulation<dim>       &triangulation,
		    const FiniteElement<dim> &fe,
		    const Quadrature<dim>    &quadrature,
		    const Function<dim>      &rhs_function,
		    const Function<dim>      &boundary_values);
    protected:
      const SmartPointer<const Function<dim> > rhs_function;
      virtual void assemble_rhs (Vector<double> &rhs) const;
  };


				   // The constructor of this class
				   // basically does what it is
				   // announced to do above...
  template <int dim>
  PrimalSolver<dim>::
  PrimalSolver (Triangulation<dim>       &triangulation,
		const FiniteElement<dim> &fe,
		const Quadrature<dim>    &quadrature,
		const Function<dim>      &rhs_function,
		const Function<dim>      &boundary_values)
		  :
		  Base<dim> (triangulation),
		  Solver<dim> (triangulation, fe,
			       quadrature, boundary_values),
                  rhs_function (&rhs_function)
  {}



				   // ... as does the ``assemble_rhs''
				   // function. Since this is
				   // explained in several of the
				   // previous example programs, we
				   // leave it at that.
  template <int dim>
  void
  PrimalSolver<dim>::
  assemble_rhs (Vector<double> &rhs) const 
  {
    FEValues<dim> fe_values (*this->fe, *this->quadrature, 
			     UpdateFlags(update_values    |
					 update_q_points  |
					 update_JxW_values));

    const unsigned int   dofs_per_cell = this->fe->dofs_per_cell;
    const unsigned int   n_q_points    = this->quadrature->n_quadrature_points;

    Vector<double>       cell_rhs (dofs_per_cell);
    std::vector<double>  rhs_values (n_q_points);
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
      cell = this->dof_handler.begin_active(),
      endc = this->dof_handler.end();
    for (; cell!=endc; ++cell)
      {
	cell_rhs = 0;
	fe_values.reinit (cell);
	rhs_function->value_list (fe_values.get_quadrature_points(),
				  rhs_values);
      
	for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    cell_rhs(i) += (fe_values.shape_value(i,q_point) *
			    rhs_values[q_point] *
			    fe_values.JxW(q_point));

	cell->get_dof_indices (local_dof_indices);
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  rhs(local_dof_indices[i]) += cell_rhs(i);
      };
  }


				   // @sect4{Global refinement}

				   // By now, all functions of the
				   // abstract base class except for
				   // the ``refine_grid'' function
				   // have been implemented. We will
				   // now have two classes that
				   // implement this function for the
				   // ``PrimalSolver'' class, one
				   // doing global refinement, one a
				   // form of local refinement.
				   //
				   // The first, doing global
				   // refinement, is rather simple:
				   // its main function just calls
				   // ``triangulation->refine_global
				   // (1);'', which does all the work.
				   //
				   // Note that since the ``Base''
				   // base class of the ``Solver''
				   // class is virtual, we have to
				   // declare a constructor that
				   // initializes the immediate base
				   // class as well as the abstract
				   // virtual one.
				   //
				   // Apart from this technical
				   // complication, the class is
				   // probably simple enough to be
				   // left without further comments.
  template <int dim>
  class RefinementGlobal : public PrimalSolver<dim>
  {
    public:
      RefinementGlobal (Triangulation<dim>       &coarse_grid,
			const FiniteElement<dim> &fe,
			const Quadrature<dim>    &quadrature,
			const Function<dim>      &rhs_function,
			const Function<dim>      &boundary_values);

      virtual void refine_grid ();
  };



  template <int dim>
  RefinementGlobal<dim>::
  RefinementGlobal (Triangulation<dim>       &coarse_grid,
		    const FiniteElement<dim> &fe,
		    const Quadrature<dim>    &quadrature,
		    const Function<dim>      &rhs_function,
		    const Function<dim>      &boundary_values)
		  :
		  Base<dim> (coarse_grid),
                  PrimalSolver<dim> (coarse_grid, fe, quadrature,
				     rhs_function, boundary_values)
  {}



  template <int dim>
  void
  RefinementGlobal<dim>::refine_grid ()
  {
    this->triangulation->refine_global (1);
  }


				   // @sect4{Local refinement by the Kelly error indicator}

				   // The second class implementing
				   // refinement strategies uses the
				   // Kelly refinemet indicator used
				   // in various example programs
				   // before. Since this indicator is
				   // already implemented in a class
				   // of its own inside the deal.II
				   // library, there is not much t do
				   // here except cal the function
				   // computing the indicator, then
				   // using it to select a number of
				   // cells for refinement and
				   // coarsening, and refinement the
				   // mesh accordingly.
				   //
				   // Again, this should now be
				   // sufficiently standard to allow
				   // the omission of further
				   // comments.
  template <int dim>
  class RefinementKelly : public PrimalSolver<dim>
  {
    public:
      RefinementKelly (Triangulation<dim>       &coarse_grid,
		       const FiniteElement<dim> &fe,
		       const Quadrature<dim>    &quadrature,
		       const Function<dim>      &rhs_function,
		       const Function<dim>      &boundary_values);

      virtual void refine_grid ();
  };



  template <int dim>
  RefinementKelly<dim>::
  RefinementKelly (Triangulation<dim>       &coarse_grid,
		   const FiniteElement<dim> &fe,
		   const Quadrature<dim>    &quadrature,
		   const Function<dim>      &rhs_function,
		   const Function<dim>      &boundary_values)
		  :
		  Base<dim> (coarse_grid),
                  PrimalSolver<dim> (coarse_grid, fe, quadrature,
				     rhs_function, boundary_values)
  {}



  template <int dim>
  void
  RefinementKelly<dim>::refine_grid ()
  {
    Vector<float> estimated_error_per_cell (this->triangulation->n_active_cells());
    KellyErrorEstimator<dim>::estimate (this->dof_handler,
					QGauss<dim-1>(3),
					typename FunctionMap<dim>::type(),
					this->solution,
					estimated_error_per_cell);
    GridRefinement::refine_and_coarsen_fixed_number (*this->triangulation,
						     estimated_error_per_cell,
						     0.3, 0.03);
    this->triangulation->execute_coarsening_and_refinement ();
  }

}




				 // @sect3{Equation data}

				 // As this is one more academic
				 // example, we'd like to compare
				 // exact and computed solution
				 // against each other. For this, we
				 // need to declare function classes
				 // representing the exact solution
				 // (for comparison and for the
				 // Dirichlet boundary values), as
				 // well as a class that denotes the
				 // right hand side of the equation
				 // (this is simply the Laplace
				 // operator applied to the exact
				 // solution we'd like to recover).
				 //
				 // For this example, let us choose as
				 // exact solution the function
				 // u(x,y)=exp(x+sin(10y+5x^2)). In more
				 // than two dimensions, simply repeat
				 // the sine-factor with ``y''
				 // replaced by ''z'' and so on. Given
				 // this, the following two classes
				 // are probably straightforward from
				 // the previous examples.
				 //
				 // As in previous examples, the C++
				 // language forces us to declare and
				 // define a constructor to the
				 // following classes even though they
				 // are empty. This is due to the fact
				 // that the base class has no default
				 // constructor (i.e. one without
				 // arguments), even though it has a
				 // constructor which has default
				 // values for all arguments.
template <int dim>
class Solution : public Function<dim>
{
  public:
    Solution () : Function<dim> () {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component) const;
};


template <int dim>
double
Solution<dim>::value (const Point<dim>   &p,
		      const unsigned int  /*component*/) const
{
  double q = p(0);
  for (unsigned int i=1; i<dim; ++i)
    q += std::sin(10*p(i)+5*p(0)*p(0));
  const double exponential = std::exp(q);
  return exponential;
}



template <int dim>
class RightHandSide : public Function<dim>
{
  public:
    RightHandSide () : Function<dim> () {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component) const;
};


template <int dim>
double
RightHandSide<dim>::value (const Point<dim>   &p,
			   const unsigned int  /*component*/) const
{
  double q = p(0);
  for (unsigned int i=1; i<dim; ++i)
    q += std::sin(10*p(i)+5*p(0)*p(0));
  const double u = std::exp(q);
  double t1 = 1,
	 t2 = 0,
	 t3 = 0;
  for (unsigned int i=1; i<dim; ++i)
    {
      t1 += std::cos(10*p(i)+5*p(0)*p(0)) * 10 * p(0);
      t2 += 10*std::cos(10*p(i)+5*p(0)*p(0)) -
	    100*std::sin(10*p(i)+5*p(0)*p(0)) * p(0)*p(0);
      t3 += 100*std::cos(10*p(i)+5*p(0)*p(0))*std::cos(10*p(i)+5*p(0)*p(0)) -
	    100*std::sin(10*p(i)+5*p(0)*p(0));
    };
  t1 = t1*t1;
  
  return -u*(t1+t2+t3);
}



				 // @sect3{The driver routines}

				 // What is now missing are only the
				 // functions that actually select the
				 // various options, and run the
				 // simulation on successively finer
				 // grids to monitor the progress as
				 // the mesh is refined.
				 //
				 // This we do in the following
				 // function: it takes a solver
				 // object, and a list of
				 // postprocessing (evaluation)
				 // objects, and runs them with
				 // intermittent mesh refinement:
template <int dim>
void
run_simulation (LaplaceSolver::Base<dim>                     &solver,
		const std::list<Evaluation::EvaluationBase<dim> *> &postprocessor_list)
{
				   // We will give an indicator of the
				   // step we are presently computing,
				   // in order to keep the user
				   // informed that something is still
				   // happening, and that the program
				   // is not in an endless loop. This
				   // is the head of this status line:
  std::cout << "Refinement cycle: ";

				   // Then start a loop which only
				   // terminates once the number of
				   // degrees of freedom is larger
				   // than 20,000 (you may of course
				   // change this limit, if you need
				   // more -- or less -- accuracy from
				   // your program).
  for (unsigned int step=0; true; ++step)
    {
				       // Then give the ``alive''
				       // indication for this
				       // iteration. Note that the
				       // ``std::flush'' is needed to
				       // have the text actually
				       // appear on the screen, rather
				       // than only in some buffer
				       // that is only flushed the
				       // next time we issue an
				       // end-line.
      std::cout << step << " " << std::flush;

				       // Now solve the problem on the
				       // present grid, and run the
				       // evaluators on it. The long
				       // type name of iterators into
				       // the list is a little
				       // annoying, but could be
				       // shortened by a typedef, if
				       // so desired.
      solver.solve_problem ();

      for (typename std::list<Evaluation::EvaluationBase<dim> *>::const_iterator
	     i = postprocessor_list.begin();
	   i != postprocessor_list.end(); ++i)
	{
	  (*i)->set_refinement_cycle (step);
	  solver.postprocess (**i);
	};


				       // Now check whether more
				       // iterations are required, or
				       // whether the loop shall be
				       // ended:
      if (solver.n_dofs() < 20000)
	solver.refine_grid ();
      else
	break;
    };

				   // Finally end the line in which we
				   // displayed status reports:
  std::cout << std::endl;
}



				 // The final function is one which
				 // takes the name of a solver
				 // (presently "kelly" and "global"
				 // are allowed), creates a solver
				 // object out of it using a coarse
				 // grid (in this case the ubiquitous
				 // unit square) and a finite element
				 // object (here the likewise
				 // ubiquitous bilinear one), and uses
				 // that solver to ask for the
				 // solution of the problem on a
				 // sequence of successively refined
				 // grids.
				 //
				 // The function also sets up two of
				 // evaluation functions, one
				 // evaluating the solution at the
				 // point (0.5,0.5), the other writing
				 // out the solution to a file.
template <int dim>
void solve_problem (const std::string &solver_name) 
{
				   // First minor task: tell the user
				   // what is going to happen. Thus
				   // write a header line, and a line
				   // with all '-' characters of the
				   // same length as the first one
				   // right below.
  const std::string header = "Running tests with \"" + solver_name +
			     "\" refinement criterion:";
  std::cout << header << std::endl
	    << std::string (header.size(), '-') << std::endl;

				   // Then set up triangulation,
				   // finite element, etc.
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (2);
  const FE_Q<dim>          fe(1);
  const QGauss<dim>       quadrature(4);
  const RightHandSide<dim> rhs_function;
  const Solution<dim>      boundary_values;

				   // Create a solver object of the
				   // kind indicated by the argument
				   // to this function. If the name is
				   // not recognized, throw an
				   // exception!
  LaplaceSolver::Base<dim> * solver = 0;
  if (solver_name == "global")
    solver = new LaplaceSolver::RefinementGlobal<dim> (triangulation, fe,
						       quadrature,
						       rhs_function,
						       boundary_values);
  else if (solver_name == "kelly")
    solver = new LaplaceSolver::RefinementKelly<dim> (triangulation, fe,
						      quadrature,
						      rhs_function,
						      boundary_values);
  else
    AssertThrow (false, ExcNotImplemented());

				   // Next create a table object in
				   // which the values of the
				   // numerical solution at the point
				   // (0.5,0.5) will be stored, and
				   // create a respective evaluation
				   // object:
  TableHandler results_table;
  Evaluation::PointValueEvaluation<dim>
    postprocessor1 (Point<dim>(0.5,0.5), results_table);

				   // Also generate an evaluator which
				   // writes out the solution:
  Evaluation::SolutionOutput<dim>
    postprocessor2 (std::string("solution-")+solver_name,
		    DataOut<dim>::gnuplot);

				   // Take these two evaluation
				   // objects and put them in a
				   // list...
  std::list<Evaluation::EvaluationBase<dim> *> postprocessor_list;
  postprocessor_list.push_back (&postprocessor1);
  postprocessor_list.push_back (&postprocessor2);

				   // ... which we can then pass on to
				   // the function that actually runs
				   // the simulation on successively
				   // refined grids:
  run_simulation (*solver, postprocessor_list);

				   // When this all is done, write out
				   // the results of the point
				   // evaluations, and finally delete
				   // the solver object:
  results_table.write_text (std::cout);
  delete solver;

				   // And one blank line after all
				   // results:
  std::cout << std::endl;
}



				 // There is not much to say about the
				 // main function. It follows the same
				 // pattern as in all previous
				 // examples, with attempts to catch
				 // thrown exceptions, and displaying
				 // as much information as possible if
				 // we should get some. The rest is
				 // self-explanatory.
int main () 
{
  try
    {
      deallog.depth_console (0);

      solve_problem<2> ("global");
      solve_problem<2> ("kelly");      
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
    };

  return 0;
}

/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 2001, 2002 */

/*    $Id$       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 2001, 2002 by the deal.II authors */
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
				 // declare an abstract base class
				 // that provides and evaluation
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
				 // it easier to actually solver
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
  {};
  

  
  template <int dim>
  void
  EvaluationBase<dim>::set_refinement_cycle (const unsigned int step)
  {
    refinement_cycle = step;
  };


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
  {};
  


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
					     // happend to be located
					     // at the vertices, but
					     // are not associated
					     // with the vertices bur
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
					     // We briefly note that
					     // this restriction on
					     // the allowed finite
					     // elements should be
					     // stated in the class
					     // documentation.
	    
	    evaluation_point_found = true;
	    break;
	  };

				     // Finally, we'd like to make
				     // sure that we have indeed found
				     // the evaluation point, since if
				     // that were not so we could not
				     // give a reasonable value of the
				     // solution there and the rest of
				     // the computation were useless
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
				     // caught somewhere between and
				     // thus already handled, and this
				     // additional information will
				     // help you find out what
				     // happened and where it went
				     // wrong.
    AssertThrow (evaluation_point_found,
		 ExcEvaluationPointNotFound(evaluation_point));

				     // If we are sure that we have
				     // found the evaluation point, we
				     // can add the results into the
				     // table of results:
    results_table.add_value ("DoFs", dof_handler.n_dofs());
    results_table.add_value ("u(x_0)", point_value);
  };




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
				   // ``eps'', ``gmv'', and ``vtk'',
				   // but this list will certainly
				   // grow over time. Now, within
				   // various functions of that base
				   // class, you can use values of
				   // this type to get information
				   // about these graphics formats
				   // (for example the default suffix
				   // used for files of each format),
				   // and you can call a generic
				   // ``write'' function, which the
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
  {};
  

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
				   // preprocessor, as already
				   // explained in the step-5 example
				   // program.
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
	     << refinement_cycle
	     << data_out.default_suffix (output_format)
	     << std::ends;
#ifdef HAVE_STD_STRINGSTREAM
    std::ofstream out (filename.str().c_str());
#else
    std::ofstream out (filename.str());
#endif
    
    data_out.write (out, output_format);
  };


				   // In practical applications, one
				   // would add here a list of other
				   // possible evaluation classes,
				   // representing quantities of
				   // interest that one is interested
				   // in. For this examples, that much
				   // shall be sufficient, so we close
				   // the namespace.
};

  
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
				 // detail in previous examples, the
				 // is not much new stuff
				 // following. Rather, we have to a
				 // great extent cannibalized previous
				 // examples and put them, in slightly
				 // different form, into this examples
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
				   // classes, solver a problem,
				   // postprocess the solution with a
				   // list of evaluation objects, and
				   // refine the grid,
				   // respectively. None of these
				   // functions actually does
				   // something itself.
				   //
				   // Due to the lack of actual
				   // functionality, the programming
				   // style of declaring very abstract
				   // base classes reminds of the
				   // style used in Smalltalk or Java
				   // programs, where all classes are
				   // even derived from entirely
				   // abstract classes ``Object'',
				   // even number representations. The
				   // author admits that he does not
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
				   // virtual functions then can be
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
  template <int dim>
  class Base
  {
    public:
      Base (Triangulation<dim> &coarse_grid);
      virtual ~Base ();

      virtual void solve_problem () = 0;
      virtual void postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const = 0;
      virtual void refine_grid () = 0;

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
  {};


  template <int dim>
  Base<dim>::~Base () 
  {};
  

				   // @sect3{A general solver class}

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
  template <int dim>
  class Solver : public virtual Base<dim>
  {
    public:
      Solver (Triangulation<dim>       &triangulation,
	      const FiniteElement<dim> &fe,
	      const Function<dim>      &boundary_values);
      virtual ~Solver ();
      virtual void solve_problem ();
      virtual void postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const;

    protected:
      const SmartPointer<const FiniteElement<dim> >  fe;
      DoFHandler<dim>                                dof_handler;
      Vector<double>                                 solution;
      const SmartPointer<const Function<dim> >       boundary_values;
      
      virtual void assemble_rhs (Vector<double> &rhs) const = 0;
    
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

      void assemble_linear_system (LinearSystem &linear_system);

      void assemble_matrix (LinearSystem                                &linear_system,
			    const DoFHandler<dim>::active_cell_iterator &begin_cell,
			    const DoFHandler<dim>::active_cell_iterator &end_cell,
			    Threads::ThreadMutex                        &mutex) const      ;
  };




  template <int dim>
  Solver<dim>::Solver (Triangulation<dim>       &triangulation,
		       const FiniteElement<dim> &fe,
		       const Function<dim>      &boundary_values)
		  :
		  Base<dim> (triangulation),
		  fe (&fe),
		  dof_handler (triangulation),
		  boundary_values (&boundary_values)
  {};


  template <int dim>
  Solver<dim>::~Solver () 
  {
    dof_handler.clear ();
  };



  template <int dim>
  void
  Solver<dim>::solve_problem ()
  {
    dof_handler.distribute_dofs (*fe);
    solution.reinit (dof_handler.n_dofs());

    LinearSystem linear_system (dof_handler);
    assemble_linear_system (linear_system);
    linear_system.solve (solution);
  };



  template <int dim>
  Solver<dim>::LinearSystem::
  LinearSystem (const DoFHandler<dim> &dof_handler)
  {
    hanging_node_constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler,
					     hanging_node_constraints);
    hanging_node_constraints.close ();

    sparsity_pattern.reinit (dof_handler.n_dofs(),
			     dof_handler.n_dofs(),
			     dof_handler.max_couplings_between_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

    hanging_node_constraints.condense (sparsity_pattern);

    sparsity_pattern.compress();

    matrix.reinit (sparsity_pattern);
    rhs.reinit (dof_handler.n_dofs());
  };



  template <int dim>
  void
  Solver<dim>::assemble_linear_system (LinearSystem &linear_system)
  {
    typedef typename DoFHandler<dim>::active_cell_iterator active_cell_iterator;

    const unsigned int n_threads = multithread_info.n_default_threads;
    std::vector<std::pair<active_cell_iterator,active_cell_iterator> >
      thread_ranges 
      = Threads::split_range<active_cell_iterator> (dof_handler.begin_active (),
						    dof_handler.end (),
						    n_threads);
    Threads::ThreadMutex mutex;
    Threads::ThreadManager thread_manager;
    for (unsigned int thread=0; thread<n_threads; ++thread)
      Threads::spawn (thread_manager,
		      Threads::encapsulate(&Solver<dim>::assemble_matrix)
		      .collect_args (this,
				     linear_system,
				     thread_ranges[thread].first,
				     thread_ranges[thread].second,
				     mutex));
    assemble_rhs (linear_system.rhs);
    linear_system.hanging_node_constraints.condense (linear_system.rhs);

    thread_manager.wait ();
    linear_system.hanging_node_constraints.condense (linear_system.matrix);

    std::map<unsigned int,double> boundary_value_map;
    VectorTools::interpolate_boundary_values (dof_handler,
					      0,
					      *boundary_values,
					      boundary_value_map);
    MatrixTools::apply_boundary_values (boundary_value_map,
					linear_system.matrix,
					solution,
					linear_system.rhs);

  };

  
  template <int dim>
  void
  Solver<dim>::assemble_matrix (LinearSystem                                &linear_system,
				const DoFHandler<dim>::active_cell_iterator &begin_cell,
				const DoFHandler<dim>::active_cell_iterator &end_cell,
				Threads::ThreadMutex                        &mutex) const
  {
				     //TODO: adaptive
    QGauss4<dim>  quadrature_formula;

    FEValues<dim> fe_values (*fe, quadrature_formula, 
			     UpdateFlags(update_gradients |
					 update_JxW_values));

    const unsigned int   dofs_per_cell = fe->dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.n_quadrature_points;

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    for (typename DoFHandler<dim>::active_cell_iterator cell=begin_cell;
	 cell!=end_cell; ++cell)
      {
	cell_matrix.clear ();

	fe_values.reinit (cell);
	const std::vector<std::vector<Tensor<1,dim> > >
	  & shape_grads  = fe_values.get_shape_grads();
	const std::vector<double>
	  & JxW_values   = fe_values.get_JxW_values();

	for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += (shape_grads[i][q_point] *
				   shape_grads[j][q_point] *
				   JxW_values[q_point]);


	cell->get_dof_indices (local_dof_indices);
	mutex.acquire ();
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    linear_system.matrix.add (local_dof_indices[i],
				      local_dof_indices[j],
				      cell_matrix(i,j));
	mutex.release ();
      };
  };



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
  };



  template <int dim>
  void
  Solver<dim>::
  postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const
  {
    postprocessor (dof_handler, solution);
  };
  

//----------------------------------------------------------    

  template <int dim>
  class PrimalSolver : public Solver<dim>
  {
    public:
      PrimalSolver (Triangulation<dim>       &triangulation,
		    const FiniteElement<dim> &fe,
		    const Function<dim>      &rhs_function,
		    const Function<dim>      &boundary_values);
    protected:
      const SmartPointer<const Function<dim> > rhs_function;
      virtual void assemble_rhs (Vector<double> &rhs) const;
  };



  template <int dim>
  PrimalSolver<dim>::
  PrimalSolver (Triangulation<dim>       &triangulation,
		const FiniteElement<dim> &fe,
		const Function<dim>      &rhs_function,
		const Function<dim>      &boundary_values)
		  :
		  Base<dim> (triangulation),
		  Solver<dim> (triangulation, fe, boundary_values),
                  rhs_function (&rhs_function)
  {};



  template <int dim>
  void
  PrimalSolver<dim>::
  assemble_rhs (Vector<double> &rhs) const 
  {
				     //TODO: adaptive
    QGauss4<dim>  quadrature_formula;

    FEValues<dim> fe_values (*fe, quadrature_formula, 
			     UpdateFlags(update_values    |
					 update_q_points  |
					 update_JxW_values));

    const unsigned int   dofs_per_cell = fe->dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.n_quadrature_points;

    Vector<double>       cell_rhs (dofs_per_cell);
    std::vector<double>  rhs_values (n_q_points);
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
	cell_rhs.clear ();

	fe_values.reinit (cell);
	const FullMatrix<double> 
	  & shape_values = fe_values.get_shape_values();
	const std::vector<double>
	  & JxW_values   = fe_values.get_JxW_values();
	const std::vector<Point<dim> >
	  & q_points     = fe_values.get_quadrature_points();

	rhs_function->value_list (q_points, rhs_values);
      
	for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    cell_rhs(i) += (shape_values (i,q_point) *
			    rhs_values[q_point] *
			    JxW_values[q_point]);

	cell->get_dof_indices (local_dof_indices);
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  rhs(local_dof_indices[i]) += cell_rhs(i);
      };
  };


//----------------------------------------------------------    

  template <int dim>
  class RefinementKelly : public PrimalSolver<dim>
  {
    public:
      RefinementKelly (Triangulation<dim>       &coarse_grid,
		       const FiniteElement<dim> &fe,
		       const Function<dim>      &rhs_function,
		       const Function<dim>      &boundary_values);

      virtual void refine_grid ();
  };



  template <int dim>
  RefinementKelly<dim>::
  RefinementKelly (Triangulation<dim>       &coarse_grid,
		   const FiniteElement<dim> &fe,
		   const Function<dim>      &rhs_function,
		   const Function<dim>      &boundary_values)
		  :
		  Base<dim> (coarse_grid),
    PrimalSolver<dim> (coarse_grid, fe, rhs_function, boundary_values)
  {};



  template <int dim>
  void
  RefinementKelly<dim>::refine_grid ()
  {
    Vector<float> estimated_error_per_cell (triangulation->n_active_cells());
    KellyErrorEstimator<dim>::estimate (dof_handler,
					QGauss3<dim-1>(),
					typename FunctionMap<dim>::type(),
					solution,
					estimated_error_per_cell);
    GridRefinement::refine_and_coarsen_fixed_number (*triangulation,
						     estimated_error_per_cell,
						     0.3, 0.03);
    triangulation->execute_coarsening_and_refinement ();
  };



//----------------------------------------------------------    

  template <int dim>
  class RefinementGlobal : public PrimalSolver<dim>
  {
    public:
      RefinementGlobal (Triangulation<dim>       &coarse_grid,
			const FiniteElement<dim> &fe,
			const Function<dim>      &rhs_function,
			const Function<dim>      &boundary_values);

      virtual void refine_grid ();
  };



  template <int dim>
  RefinementGlobal<dim>::
  RefinementGlobal (Triangulation<dim>       &coarse_grid,
		    const FiniteElement<dim> &fe,
		    const Function<dim>      &rhs_function,
		    const Function<dim>      &boundary_values)
		  :
		  Base<dim> (coarse_grid),
    PrimalSolver<dim> (fe, rhs_function, boundary_values)
  {};



  template <int dim>
  void
  RefinementGlobal<dim>::refine_grid ()
  {
    triangulation->refine_global (1);
  };
};




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
template <int dim>
class Solution : public Function<dim>
{
  public:
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
    q += sin(10*p(i)+5*p(0)*p(0));
  const double exponential = exp(q);
  return exponential;
};



template <int dim>
class RightHandSide : public Function<dim>
{
  public:
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
    q += sin(10*p(i)+5*p(0)*p(0));
  const double u = exp(q);
  double t1 = 1,
	 t2 = 0,
	 t3 = 0;
  for (unsigned int i=1; i<dim; ++i)
    {
      t1 += cos(10*p(i)+5*p(0)*p(0)) * 10 * p(0);
      t2 += 10*cos(10*p(i)+5*p(0)*p(0)) -
	    100*sin(10*p(i)+5*p(0)*p(0)) * p(0)*p(0);
      t3 += 100*cos(10*p(i)+5*p(0)*p(0))*cos(10*p(i)+5*p(0)*p(0)) -
	    100*sin(10*p(i)+5*p(0)*p(0));
    };
  t1 = t1*t1;
  
  return -u*(t1+t2+t3);
};



				 // @sect3{The driver routines}


template <int dim>
void
run_simulation (LaplaceSolver::Base<dim>                     &solver,
		const std::list<Evaluation::EvaluationBase<dim> *> &postprocessor_list)
{
  const unsigned int max_steps = 10;
  for (unsigned int step=0; step<max_steps; ++step)
    {
      std::cout << "Refinement cycle " << step << std::endl;
      
      solver.solve_problem ();

      for (typename std::list<Evaluation::EvaluationBase<dim> *>::const_iterator
	     i = postprocessor_list.begin();
	   i != postprocessor_list.end(); ++i)
	{
	  (*i)->set_refinement_cycle (step);
	  solver.postprocess (**i);
	};

      if (step!=max_steps-1)
	solver.refine_grid ();
    };
};


template <int dim>
void solve_problem_kelly () 
{      
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (2);
  FE_Q<dim> fe(1);
  const RightHandSide<dim> rhs_function;
  const Solution<dim>      boundary_values;
      
  LaplaceSolver::RefinementKelly<dim> kelly (triangulation, fe,
					     rhs_function,
					     boundary_values);
  TableHandler results_table;
  
  Evaluation::PointValueEvaluation<dim>
    postprocessor1 (Point<dim>(.5,.5), results_table);
  Evaluation::SolutionOutput<dim>
    postprocessor2 ("solution-kelly", DataOut<dim>::gnuplot);
  std::list<Evaluation::EvaluationBase<dim> *> postprocessor_list;
  postprocessor_list.push_back (&postprocessor1);
  postprocessor_list.push_back (&postprocessor2);
  
  run_simulation (kelly, postprocessor_list);

  results_table.write_text (std::cout);
};


    
int main () 
{
  try
    {
      deallog.depth_console (0);

      solve_problem_kelly<2> ();
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
};

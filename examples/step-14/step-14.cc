/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2002 - 2014 by the deal.II authors
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
 * Author: Wolfgang Bangerth, ETH Zurich, 2002
 */


// Start out with well known things...
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>
#include <numeric>
#include <sstream>

// The last step is as in all previous programs:
namespace Step14
{
  using namespace dealii;

  // @sect3{Evaluating the solution}

  // As mentioned in the introduction, significant parts of the program have
  // simply been taken over from the step-13 example program. We therefore
  // only comment on those things that are new.
  //
  // First, the framework for evaluation of solutions is unchanged, i.e. the
  // base class is the same, and the class to evaluate the solution at a grid
  // point is unchanged:
  namespace Evaluation
  {
    // @sect4{The EvaluationBase class}
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


    template <int dim>
    EvaluationBase<dim>::~EvaluationBase ()
    {}



    template <int dim>
    void
    EvaluationBase<dim>::set_refinement_cycle (const unsigned int step)
    {
      refinement_cycle = step;
    }


    // @sect4{The PointValueEvaluation class}
    template <int dim>
    class PointValueEvaluation : public EvaluationBase<dim>
    {
    public:
      PointValueEvaluation (const Point<dim>   &evaluation_point);

      virtual void operator () (const DoFHandler<dim> &dof_handler,
                                const Vector<double>  &solution) const;

      DeclException1 (ExcEvaluationPointNotFound,
                      Point<dim>,
                      << "The evaluation point " << arg1
                      << " was not found among the vertices of the present grid.");
    private:
      const Point<dim>  evaluation_point;
    };


    template <int dim>
    PointValueEvaluation<dim>::
    PointValueEvaluation (const Point<dim>   &evaluation_point)
      :
      evaluation_point (evaluation_point)
    {}



    template <int dim>
    void
    PointValueEvaluation<dim>::
    operator () (const DoFHandler<dim> &dof_handler,
                 const Vector<double>  &solution) const
    {
      double point_value = 1e20;

      typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
      bool evaluation_point_found = false;
      for (; (cell!=endc) && !evaluation_point_found; ++cell)
        for (unsigned int vertex=0;
             vertex<GeometryInfo<dim>::vertices_per_cell;
             ++vertex)
          if (cell->vertex(vertex).distance (evaluation_point)
              <
              cell->diameter() * 1e-8)
            {
              point_value = solution(cell->vertex_dof_index(vertex,0));

              evaluation_point_found = true;
              break;
            }

      AssertThrow (evaluation_point_found,
                   ExcEvaluationPointNotFound(evaluation_point));

      std::cout << "   Point value=" << point_value
                << std::endl;
    }


    // @sect4{The PointXDerivativeEvaluation class}

    // Besides the class implementing the evaluation of the solution at one
    // point, we here provide one which evaluates the gradient at a grid
    // point. Since in general the gradient of a finite element function is
    // not continuous at a vertex, we have to be a little bit more careful
    // here. What we do is to loop over all cells, even if we have found the
    // point already on one cell, and use the mean value of the gradient at
    // the vertex taken from all adjacent cells.
    //
    // Given the interface of the <code>PointValueEvaluation</code> class, the
    // declaration of this class provides little surprise, and neither does
    // the constructor:
    template <int dim>
    class PointXDerivativeEvaluation : public EvaluationBase<dim>
    {
    public:
      PointXDerivativeEvaluation (const Point<dim>   &evaluation_point);

      virtual void operator () (const DoFHandler<dim> &dof_handler,
                                const Vector<double>  &solution) const;

      DeclException1 (ExcEvaluationPointNotFound,
                      Point<dim>,
                      << "The evaluation point " << arg1
                      << " was not found among the vertices of the present grid.");
    private:
      const Point<dim>  evaluation_point;
    };


    template <int dim>
    PointXDerivativeEvaluation<dim>::
    PointXDerivativeEvaluation (const Point<dim>   &evaluation_point)
      :
      evaluation_point (evaluation_point)
    {}


    // The more interesting things happen inside the function doing the actual
    // evaluation:
    template <int dim>
    void
    PointXDerivativeEvaluation<dim>::
    operator () (const DoFHandler<dim> &dof_handler,
                 const Vector<double>  &solution) const
    {
      // This time initialize the return value with something useful, since we
      // will have to add up a number of contributions and take the mean value
      // afterwards...
      double point_derivative = 0;

      // ...then have some objects of which the meaning will become clear
      // below...
      QTrapez<dim>  vertex_quadrature;
      FEValues<dim> fe_values (dof_handler.get_fe(),
                               vertex_quadrature,
                               update_gradients | update_quadrature_points);
      std::vector<Tensor<1,dim> >
      solution_gradients (vertex_quadrature.size());

      // ...and next loop over all cells and their vertices, and count how
      // often the vertex has been found:
      typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
      unsigned int evaluation_point_hits = 0;
      for (; cell!=endc; ++cell)
        for (unsigned int vertex=0;
             vertex<GeometryInfo<dim>::vertices_per_cell;
             ++vertex)
          if (cell->vertex(vertex) == evaluation_point)
            {
              // Things are now no more as simple, since we can't get the
              // gradient of the finite element field as before, where we
              // simply had to pick one degree of freedom at a vertex.
              //
              // Rather, we have to evaluate the finite element field on this
              // cell, and at a certain point. As you know, evaluating finite
              // element fields at certain points is done through the
              // <code>FEValues</code> class, so we use that. The question is:
              // the <code>FEValues</code> object needs to be a given a
              // quadrature formula and can then compute the values of finite
              // element quantities at the quadrature points. Here, we don't
              // want to do quadrature, we simply want to specify some points!
              //
              // Nevertheless, the same way is chosen: use a special
              // quadrature rule with points at the vertices, since these are
              // what we are interested in. The appropriate rule is the
              // trapezoidal rule, so that is the reason why we used that one
              // above.
              //
              // Thus: initialize the <code>FEValues</code> object on this
              // cell,
              fe_values.reinit (cell);
              // and extract the gradients of the solution vector at the
              // vertices:
              fe_values.get_function_grads (solution,
                                            solution_gradients);

              // Now we have the gradients at all vertices, so pick out that
              // one which belongs to the evaluation point (note that the
              // order of vertices is not necessarily the same as that of the
              // quadrature points):
              unsigned int q_point = 0;
              for (; q_point<solution_gradients.size(); ++q_point)
                if (fe_values.quadrature_point(q_point) ==
                    evaluation_point)
                  break;

              // Check that the evaluation point was indeed found,
              Assert (q_point < solution_gradients.size(),
                      ExcInternalError());
              // and if so take the x-derivative of the gradient there as the
              // value which we are interested in, and increase the counter
              // indicating how often we have added to that variable:
              point_derivative += solution_gradients[q_point][0];
              ++evaluation_point_hits;

              // Finally break out of the innermost loop iterating over the
              // vertices of the present cell, since if we have found the
              // evaluation point at one vertex it cannot be at a following
              // vertex as well:
              break;
            }

      // Now we have looped over all cells and vertices, so check whether the
      // point was found:
      AssertThrow (evaluation_point_hits > 0,
                   ExcEvaluationPointNotFound(evaluation_point));

      // We have simply summed up the contributions of all adjacent cells, so
      // we still have to compute the mean value. Once this is done, report
      // the status:
      point_derivative /= evaluation_point_hits;
      std::cout << "   Point x-derivative=" << point_derivative
                << std::endl;
    }



    // @sect4{The GridOutput class}

    // Since this program has a more difficult structure (it computed a dual
    // solution in addition to a primal one), writing out the solution is no
    // more done by an evaluation object since we want to write both solutions
    // at once into one file, and that requires some more information than
    // available to the evaluation classes.
    //
    // However, we also want to look at the grids generated. This again can be
    // done with one such class. Its structure is analog to the
    // <code>SolutionOutput</code> class of the previous example program, so
    // we do not discuss it here in more detail. Furthermore, everything that
    // is used here has already been used in previous example programs.
    template <int dim>
    class GridOutput : public EvaluationBase<dim>
    {
    public:
      GridOutput (const std::string &output_name_base);

      virtual void operator () (const DoFHandler<dim> &dof_handler,
                                const Vector<double>  &solution) const;
    private:
      const std::string output_name_base;
    };


    template <int dim>
    GridOutput<dim>::
    GridOutput (const std::string &output_name_base)
      :
      output_name_base (output_name_base)
    {}


    template <int dim>
    void
    GridOutput<dim>::operator () (const DoFHandler<dim> &dof_handler,
                                  const Vector<double>  &/*solution*/) const
    {
      std::ostringstream filename;
      filename << output_name_base << "-"
               << this->refinement_cycle
               << ".eps"
               << std::ends;

      std::ofstream out (filename.str().c_str());
      GridOut().write_eps (dof_handler.get_tria(), out);
    }
  }


  // @sect3{The Laplace solver classes}

  // Next are the actual solver classes. Again, we discuss only the
  // differences to the previous program.
  namespace LaplaceSolver
  {
    // Before everything else, forward-declare one class that we will have
    // later, since we will want to make it a friend of some of the classes
    // that follow, which requires the class to be known:
    template <int dim> class WeightedResidual;


    // @sect4{The Laplace solver base class}

    // This class is almost unchanged, with the exception that it declares two
    // more functions: <code>output_solution</code> will be used to generate
    // output files from the actual solutions computed by derived classes, and
    // the <code>set_refinement_cycle</code> function by which the testing
    // framework sets the number of the refinement cycle to a local variable
    // in this class; this number is later used to generate filenames for the
    // solution output.
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

      virtual void set_refinement_cycle (const unsigned int cycle);

      virtual void output_solution () const = 0;

    protected:
      const SmartPointer<Triangulation<dim> > triangulation;

      unsigned int refinement_cycle;
    };


    template <int dim>
    Base<dim>::Base (Triangulation<dim> &coarse_grid)
      :
      triangulation (&coarse_grid)
    {}


    template <int dim>
    Base<dim>::~Base ()
    {}



    template <int dim>
    void
    Base<dim>::set_refinement_cycle (const unsigned int cycle)
    {
      refinement_cycle = cycle;
    }


    // @sect4{The Laplace Solver class}

    // Likewise, the <code>Solver</code> class is entirely unchanged and will
    // thus not be discussed.
    template <int dim>
    class Solver : public virtual Base<dim>
    {
    public:
      Solver (Triangulation<dim>       &triangulation,
              const FiniteElement<dim> &fe,
              const Quadrature<dim>    &quadrature,
              const Quadrature<dim-1>  &face_quadrature,
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

    protected:
      const SmartPointer<const FiniteElement<dim> >  fe;
      const SmartPointer<const Quadrature<dim> >     quadrature;
      const SmartPointer<const Quadrature<dim-1> >   face_quadrature;
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


      // The remainder of the class is essentially a copy of step-13
      // as well, including the data structures and functions
      // necessary to compute the linear system in parallel using the
      // WorkStream framework:
      struct AssemblyScratchData
      {
        AssemblyScratchData (const FiniteElement<dim> &fe,
                             const Quadrature<dim>    &quadrature);
        AssemblyScratchData (const AssemblyScratchData &scratch_data);

        FEValues<dim>     fe_values;
      };

      struct AssemblyCopyData
      {
        FullMatrix<double> cell_matrix;
        std::vector<types::global_dof_index> local_dof_indices;
      };


      void
      assemble_linear_system (LinearSystem &linear_system);

      void
      local_assemble_matrix (const typename DoFHandler<dim>::active_cell_iterator &cell,
                             AssemblyScratchData                                  &scratch_data,
                             AssemblyCopyData                                     &copy_data) const;


      void
      copy_local_to_global(const AssemblyCopyData &copy_data,
                           LinearSystem           &linear_system) const;
    };



    template <int dim>
    Solver<dim>::Solver (Triangulation<dim>       &triangulation,
                         const FiniteElement<dim> &fe,
                         const Quadrature<dim>    &quadrature,
                         const Quadrature<dim-1>  &face_quadrature,
                         const Function<dim>      &boundary_values)
      :
      Base<dim> (triangulation),
      fe (&fe),
      quadrature (&quadrature),
      face_quadrature (&face_quadrature),
      dof_handler (triangulation),
      boundary_values (&boundary_values)
    {}


    template <int dim>
    Solver<dim>::~Solver ()
    {
      dof_handler.clear ();
    }


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


    template <int dim>
    void
    Solver<dim>::
    postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const
    {
      postprocessor (dof_handler, solution);
    }


    template <int dim>
    unsigned int
    Solver<dim>::n_dofs () const
    {
      return dof_handler.n_dofs();
    }


    // The following few functions and constructors are verbatim
    // copies taken from step-13:
    template <int dim>
    void
    Solver<dim>::assemble_linear_system (LinearSystem &linear_system)
    {
      Threads::Task<> rhs_task = Threads::new_task (&Solver<dim>::assemble_rhs,
                                                    *this,
                                                    linear_system.rhs);

      WorkStream::run(dof_handler.begin_active(),
                      dof_handler.end(),
                      std_cxx11::bind(&Solver<dim>::local_assemble_matrix,
                                      this,
                                      std_cxx11::_1,
                                      std_cxx11::_2,
                                      std_cxx11::_3),
                      std_cxx11::bind(&Solver<dim>::copy_local_to_global,
                                      this,
                                      std_cxx11::_1,
                                      std_cxx11::ref(linear_system)),
                      AssemblyScratchData(*fe, *quadrature),
                      AssemblyCopyData());

      rhs_task.join ();

      linear_system.hanging_node_constraints.condense (linear_system.rhs);

      std::map<types::global_dof_index,double> boundary_value_map;
      VectorTools::interpolate_boundary_values (dof_handler,
                                                0,
                                                *boundary_values,
                                                boundary_value_map);

      linear_system.hanging_node_constraints.condense (linear_system.matrix);

      MatrixTools::apply_boundary_values (boundary_value_map,
                                          linear_system.matrix,
                                          solution,
                                          linear_system.rhs);
    }


    template <int dim>
    Solver<dim>::AssemblyScratchData::
    AssemblyScratchData (const FiniteElement<dim> &fe,
                         const Quadrature<dim>    &quadrature)
      :
      fe_values (fe,
                 quadrature,
                 update_gradients | update_JxW_values)
    {}


    template <int dim>
    Solver<dim>::AssemblyScratchData::
    AssemblyScratchData (const AssemblyScratchData &scratch_data)
      :
      fe_values (scratch_data.fe_values.get_fe(),
                 scratch_data.fe_values.get_quadrature(),
                 update_gradients | update_JxW_values)
    {}


    template <int dim>
    void
    Solver<dim>::local_assemble_matrix (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                        AssemblyScratchData                                  &scratch_data,
                                        AssemblyCopyData                                     &copy_data) const
    {
      const unsigned int dofs_per_cell = fe->dofs_per_cell;
      const unsigned int n_q_points    = quadrature->size();

      copy_data.cell_matrix.reinit (dofs_per_cell, dofs_per_cell);

      copy_data.local_dof_indices.resize(dofs_per_cell);

      scratch_data.fe_values.reinit (cell);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            copy_data.cell_matrix(i,j) += (scratch_data.fe_values.shape_grad(i,q_point) *
                                           scratch_data.fe_values.shape_grad(j,q_point) *
                                           scratch_data.fe_values.JxW(q_point));

      cell->get_dof_indices (copy_data.local_dof_indices);
    }



    template <int dim>
    void
    Solver<dim>::copy_local_to_global(const AssemblyCopyData &copy_data,
                                      LinearSystem           &linear_system) const
    {
      for (unsigned int i=0; i<copy_data.local_dof_indices.size(); ++i)
        for (unsigned int j=0; j<copy_data.local_dof_indices.size(); ++j)
          linear_system.matrix.add (copy_data.local_dof_indices[i],
                                    copy_data.local_dof_indices[j],
                                    copy_data.cell_matrix(i,j));
    }


    // Now for the functions that implement actions in the linear
    // system class. First, the constructor initializes all data
    // elements to their correct sizes, and sets up a number of
    // additional data structures, such as constraints due to hanging
    // nodes. Since setting up the hanging nodes and finding out about
    // the nonzero elements of the matrix is independent, we do that
    // in parallel (if the library was configured to use concurrency,
    // at least; otherwise, the actions are performed
    // sequentially). Note that we start only one thread, and do the
    // second action in the main thread. Since only one thread is
    // generated, we don't use the <code>Threads::ThreadGroup</code>
    // class here, but rather use the one created thread object
    // directly to wait for this particular thread's exit. The
    // approach is generally the same as the one we have used in
    // <code>Solver::assemble_linear_system()</code> above.
    //
    // Note that taking the address of the
    // <code>DoFTools::make_hanging_node_constraints</code> function
    // is a little tricky, since there are actually three functions of
    // this name, one for each supported space dimension. Taking
    // addresses of overloaded functions is somewhat complicated in
    // C++, since the address-of operator <code>&</code> in that case
    // returns a set of values (the addresses of all
    // functions with that name), and selecting the right one is then
    // the next step. If the context dictates which one to take (for
    // example by assigning to a function pointer of known type), then
    // the compiler can do that by itself, but if this set of pointers
    // shall be given as the argument to a function that takes a
    // template, the compiler could choose all without having a
    // preference for one. We therefore have to make it clear to the
    // compiler which one we would like to have; for this, we could
    // use a cast, but for more clarity, we assign it to a temporary
    // <code>mhnc_p</code> (short for <code>pointer to
    // make_hanging_node_constraints</code>) with the right type, and
    // using this pointer instead.
    template <int dim>
    Solver<dim>::LinearSystem::
    LinearSystem (const DoFHandler<dim> &dof_handler)
    {
      hanging_node_constraints.clear ();

      void (*mhnc_p) (const DoFHandler<dim> &,
                      ConstraintMatrix &)
        = &DoFTools::make_hanging_node_constraints;

      Threads::Task<> side_task
        = Threads::new_task (mhnc_p,
                             dof_handler,
                             hanging_node_constraints);

      sparsity_pattern.reinit (dof_handler.n_dofs(),
                               dof_handler.n_dofs(),
                               dof_handler.max_couplings_between_dofs());
      DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

      side_task.join();

      hanging_node_constraints.close ();
      hanging_node_constraints.condense (sparsity_pattern);

      sparsity_pattern.compress();
      matrix.reinit (sparsity_pattern);
      rhs.reinit (dof_handler.n_dofs());
    }



    template <int dim>
    void
    Solver<dim>::LinearSystem::solve (Vector<double> &solution) const
    {
      SolverControl           solver_control (5000, 1e-12);
      SolverCG<>              cg (solver_control);

      PreconditionSSOR<> preconditioner;
      preconditioner.initialize(matrix, 1.2);

      cg.solve (matrix, solution, rhs, preconditioner);

      hanging_node_constraints.distribute (solution);
    }




    // @sect4{The PrimalSolver class}

    // The <code>PrimalSolver</code> class is also mostly unchanged except for
    // overloading the functions <code>solve_problem</code>,
    // <code>n_dofs</code>, and <code>postprocess</code> of the base class,
    // and implementing the <code>output_solution</code> function. These
    // overloaded functions do nothing particular besides calling the
    // functions of the base class -- that seems superfluous, but works around
    // a bug in a popular compiler which requires us to write such functions
    // for the following scenario: Besides the <code>PrimalSolver</code>
    // class, we will have a <code>DualSolver</code>, both derived from
    // <code>Solver</code>. We will then have a final classes which derived
    // from these two, which will then have two instances of the
    // <code>Solver</code> class as its base classes. If we want, for example,
    // the number of degrees of freedom of the primal solver, we would have to
    // indicate this like so: <code>PrimalSolver::n_dofs()</code>.  However,
    // the compiler does not accept this since the <code>n_dofs</code>
    // function is actually from a base class of the <code>PrimalSolver</code>
    // class, so we have to inject the name from the base to the derived class
    // using these additional functions.
    //
    // Regarding the implementation of the <code>output_solution</code>
    // function, we keep the <code>GlobalRefinement</code> and
    // <code>RefinementKelly</code> classes in this program, and they can then
    // rely on the default implementation of this function which simply
    // outputs the primal solution. The class implementing dual weighted error
    // estimators will overload this function itself, to also output the dual
    // solution.
    //
    // Except for this, the class is unchanged with respect to the previous
    // example.
    template <int dim>
    class PrimalSolver : public Solver<dim>
    {
    public:
      PrimalSolver (Triangulation<dim>       &triangulation,
                    const FiniteElement<dim> &fe,
                    const Quadrature<dim>    &quadrature,
                    const Quadrature<dim-1>  &face_quadrature,
                    const Function<dim>      &rhs_function,
                    const Function<dim>      &boundary_values);

      virtual
      void solve_problem ();

      virtual
      unsigned int n_dofs () const;

      virtual
      void postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const;

      virtual
      void output_solution () const;

    protected:
      const SmartPointer<const Function<dim> > rhs_function;
      virtual void assemble_rhs (Vector<double> &rhs) const;

      // Now, in order to work around some problems in one of the compilers
      // this library can be compiled with, we will have to declare a class
      // that is actually derived from the present one, as a friend (strange
      // as that seems). The full rationale will be explained below.
      friend class WeightedResidual<dim>;
    };


    template <int dim>
    PrimalSolver<dim>::
    PrimalSolver (Triangulation<dim>       &triangulation,
                  const FiniteElement<dim> &fe,
                  const Quadrature<dim>    &quadrature,
                  const Quadrature<dim-1>  &face_quadrature,
                  const Function<dim>      &rhs_function,
                  const Function<dim>      &boundary_values)
      :
      Base<dim> (triangulation),
      Solver<dim> (triangulation, fe,
                   quadrature, face_quadrature,
                   boundary_values),
      rhs_function (&rhs_function)
    {}


    template <int dim>
    void
    PrimalSolver<dim>::solve_problem ()
    {
      Solver<dim>::solve_problem ();
    }



    template <int dim>
    unsigned int
    PrimalSolver<dim>::n_dofs() const
    {
      return Solver<dim>::n_dofs();
    }


    template <int dim>
    void
    PrimalSolver<dim>::
    postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const
    {
      Solver<dim>::postprocess(postprocessor);
    }


    template <int dim>
    void
    PrimalSolver<dim>::output_solution () const
    {
      DataOut<dim> data_out;
      data_out.attach_dof_handler (this->dof_handler);
      data_out.add_data_vector (this->solution, "solution");
      data_out.build_patches ();

      std::ostringstream filename;
      filename << "solution-"
               << this->refinement_cycle
               << ".gnuplot"
               << std::ends;

      std::ofstream out (filename.str().c_str());
      data_out.write (out, DataOutBase::gnuplot);
    }



    template <int dim>
    void
    PrimalSolver<dim>::
    assemble_rhs (Vector<double> &rhs) const
    {
      FEValues<dim> fe_values (*this->fe, *this->quadrature,
                               update_values  | update_quadrature_points  |
                               update_JxW_values);

      const unsigned int   dofs_per_cell = this->fe->dofs_per_cell;
      const unsigned int   n_q_points    = this->quadrature->size();

      Vector<double>       cell_rhs (dofs_per_cell);
      std::vector<double>  rhs_values (n_q_points);
      std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

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
        }
    }


    // @sect4{The RefinementGlobal and RefinementKelly classes}

    // For the following two classes, the same applies as for most of the
    // above: the class is taken from the previous example as-is:
    template <int dim>
    class RefinementGlobal : public PrimalSolver<dim>
    {
    public:
      RefinementGlobal (Triangulation<dim>       &coarse_grid,
                        const FiniteElement<dim> &fe,
                        const Quadrature<dim>    &quadrature,
                        const Quadrature<dim-1>  &face_quadrature,
                        const Function<dim>      &rhs_function,
                        const Function<dim>      &boundary_values);

      virtual void refine_grid ();
    };



    template <int dim>
    RefinementGlobal<dim>::
    RefinementGlobal (Triangulation<dim>       &coarse_grid,
                      const FiniteElement<dim> &fe,
                      const Quadrature<dim>    &quadrature,
                      const Quadrature<dim-1>  &face_quadrature,
                      const Function<dim>      &rhs_function,
                      const Function<dim>      &boundary_values)
      :
      Base<dim> (coarse_grid),
      PrimalSolver<dim> (coarse_grid, fe, quadrature,
                         face_quadrature, rhs_function,
                         boundary_values)
    {}



    template <int dim>
    void
    RefinementGlobal<dim>::refine_grid ()
    {
      this->triangulation->refine_global (1);
    }



    template <int dim>
    class RefinementKelly : public PrimalSolver<dim>
    {
    public:
      RefinementKelly (Triangulation<dim>       &coarse_grid,
                       const FiniteElement<dim> &fe,
                       const Quadrature<dim>    &quadrature,
                       const Quadrature<dim-1>  &face_quadrature,
                       const Function<dim>      &rhs_function,
                       const Function<dim>      &boundary_values);

      virtual void refine_grid ();
    };



    template <int dim>
    RefinementKelly<dim>::
    RefinementKelly (Triangulation<dim>       &coarse_grid,
                     const FiniteElement<dim> &fe,
                     const Quadrature<dim>    &quadrature,
                     const Quadrature<dim-1>  &face_quadrature,
                     const Function<dim>      &rhs_function,
                     const Function<dim>      &boundary_values)
      :
      Base<dim> (coarse_grid),
      PrimalSolver<dim> (coarse_grid, fe, quadrature,
                         face_quadrature,
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



    // @sect4{The RefinementWeightedKelly class}

    // This class is a variant of the previous one, in that it allows to
    // weight the refinement indicators we get from the library's Kelly
    // indicator by some function. We include this class since the goal of
    // this example program is to demonstrate automatic refinement criteria
    // even for complex output quantities such as point values or stresses. If
    // we did not solve a dual problem and compute the weights thereof, we
    // would probably be tempted to give a hand-crafted weighting to the
    // indicators to account for the fact that we are going to evaluate these
    // quantities. This class accepts such a weighting function as argument to
    // its constructor:
    template <int dim>
    class RefinementWeightedKelly : public PrimalSolver<dim>
    {
    public:
      RefinementWeightedKelly (Triangulation<dim>       &coarse_grid,
                               const FiniteElement<dim> &fe,
                               const Quadrature<dim>    &quadrature,
                               const Quadrature<dim-1>  &face_quadrature,
                               const Function<dim>      &rhs_function,
                               const Function<dim>      &boundary_values,
                               const Function<dim>      &weighting_function);

      virtual void refine_grid ();

    private:
      const SmartPointer<const Function<dim> > weighting_function;
    };



    template <int dim>
    RefinementWeightedKelly<dim>::
    RefinementWeightedKelly (Triangulation<dim>       &coarse_grid,
                             const FiniteElement<dim> &fe,
                             const Quadrature<dim>    &quadrature,
                             const Quadrature<dim-1>  &face_quadrature,
                             const Function<dim>      &rhs_function,
                             const Function<dim>      &boundary_values,
                             const Function<dim>      &weighting_function)
      :
      Base<dim> (coarse_grid),
      PrimalSolver<dim> (coarse_grid, fe, quadrature,
                         face_quadrature,
                         rhs_function, boundary_values),
      weighting_function (&weighting_function)
    {}



    // Now, here comes the main function, including the weighting:
    template <int dim>
    void
    RefinementWeightedKelly<dim>::refine_grid ()
    {
      // First compute some residual based error indicators for all cells by a
      // method already implemented in the library. What exactly is computed
      // can be read in the documentation of that class.
      Vector<float> estimated_error (this->triangulation->n_active_cells());
      KellyErrorEstimator<dim>::estimate (this->dof_handler,
                                          *this->face_quadrature,
                                          typename FunctionMap<dim>::type(),
                                          this->solution,
                                          estimated_error);

      // Now we are going to weight these indicators by the value of the
      // function given to the constructor:
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->dof_handler.begin_active(),
      endc = this->dof_handler.end();
      for (unsigned int cell_index=0; cell!=endc; ++cell, ++cell_index)
        estimated_error(cell_index)
        *= weighting_function->value (cell->center());

      GridRefinement::refine_and_coarsen_fixed_number (*this->triangulation,
                                                       estimated_error,
                                                       0.3, 0.03);
      this->triangulation->execute_coarsening_and_refinement ();
    }

  }


  // @sect3{Equation data}
  //
  // In this example program, we work with the same data sets as in the
  // previous one, but as it may so happen that someone wants to run the
  // program with different boundary values and right hand side functions, or
  // on a different grid, we show a simple technique to do exactly that. For
  // more clarity, we furthermore pack everything that has to do with equation
  // data into a namespace of its own.
  //
  // The underlying assumption is that this is a research program, and that
  // there we often have a number of test cases that consist of a domain, a
  // right hand side, boundary values, possibly a specified coefficient, and a
  // number of other parameters. They often vary all at the same time when
  // shifting from one example to another. To make handling such sets of
  // problem description parameters simple is the goal of the following.
  //
  // Basically, the idea is this: let us have a structure for each set of
  // data, in which we pack everything that describes a test case: here, these
  // are two subclasses, one called <code>BoundaryValues</code> for the
  // boundary values of the exact solution, and one called
  // <code>RightHandSide</code>, and then a way to generate the coarse
  // grid. Since the solution of the previous example program looked like
  // curved ridges, we use this name here for the enclosing class. Note that
  // the names of the two inner classes have to be the same for all enclosing
  // test case classes, and also that we have attached the dimension template
  // argument to the enclosing class rather than to the inner ones, to make
  // further processing simpler.  (From a language viewpoint, a namespace
  // would be better to encapsulate these inner classes, rather than a
  // structure. However, namespaces cannot be given as template arguments, so
  // we use a structure to allow a second object to select from within its
  // given argument. The enclosing structure, of course, has no member
  // variables apart from the classes it declares, and a static function to
  // generate the coarse mesh; it will in general never be instantiated.)
  //
  // The idea is then the following (this is the right time to also take a
  // brief look at the code below): we can generate objects for boundary
  // values and right hand side by simply giving the name of the outer class
  // as a template argument to a class which we call here
  // <code>Data::SetUp</code>, and it then creates objects for the inner
  // classes. In this case, to get all that characterizes the curved ridge
  // solution, we would simply generate an instance of
  // <code>Data::SetUp@<Data::CurvedRidge@></code>, and everything we need to
  // know about the solution would be static member variables and functions of
  // that object.
  //
  // This approach might seem like overkill in this case, but will become very
  // handy once a certain set up is not only characterized by Dirichlet
  // boundary values and a right hand side function, but in addition by
  // material properties, Neumann values, different boundary descriptors,
  // etc. In that case, the <code>SetUp</code> class might consist of a dozen
  // or more objects, and each descriptor class (like the
  // <code>CurvedRidges</code> class below) would have to provide them. Then,
  // you will be happy to be able to change from one set of data to another by
  // only changing the template argument to the <code>SetUp</code> class at
  // one place, rather than at many.
  //
  // With this framework for different test cases, we are almost finished, but
  // one thing remains: by now we can select statically, by changing one
  // template argument, which data set to choose. In order to be able to do
  // that dynamically, i.e. at run time, we need a base class. This we provide
  // in the obvious way, see below, with virtual abstract functions. It forces
  // us to introduce a second template parameter <code>dim</code> which we
  // need for the base class (which could be avoided using some template
  // magic, but we omit that), but that's all.
  //
  // Adding new testcases is now simple, you don't have to touch the framework
  // classes, only a structure like the <code>CurvedRidges</code> one is
  // needed.
  namespace Data
  {
    // @sect4{The SetUpBase and SetUp classes}

    // Based on the above description, the <code>SetUpBase</code> class then
    // looks as follows. To allow using the <code>SmartPointer</code> class
    // with this class, we derived from the <code>Subscriptor</code> class.
    template <int dim>
    struct SetUpBase : public Subscriptor
    {
      virtual
      const Function<dim>   &get_boundary_values () const = 0;

      virtual
      const Function<dim>   &get_right_hand_side () const = 0;

      virtual
      void create_coarse_grid (Triangulation<dim> &coarse_grid) const = 0;
    };


    // And now for the derived class that takes the template argument as
    // explained above. For some reason, C++ requires us to define a
    // constructor (which maybe empty), as otherwise a warning is generated
    // that some data is not initialized.
    //
    // Here we pack the data elements into private variables, and allow access
    // to them through the methods of the base class.
    template <class Traits, int dim>
    struct SetUp : public SetUpBase<dim>
    {
      SetUp () {}

      virtual
      const Function<dim>   &get_boundary_values () const;

      virtual
      const Function<dim>   &get_right_hand_side () const;


      virtual
      void create_coarse_grid (Triangulation<dim> &coarse_grid) const;

    private:
      static const typename Traits::BoundaryValues boundary_values;
      static const typename Traits::RightHandSide  right_hand_side;
    };

    // We have to provide definitions for the static member variables of the
    // above class:
    template <class Traits, int dim>
    const typename Traits::BoundaryValues  SetUp<Traits,dim>::boundary_values;
    template <class Traits, int dim>
    const typename Traits::RightHandSide   SetUp<Traits,dim>::right_hand_side;

    // And definitions of the member functions:
    template <class Traits, int dim>
    const Function<dim> &
    SetUp<Traits,dim>::get_boundary_values () const
    {
      return boundary_values;
    }


    template <class Traits, int dim>
    const Function<dim> &
    SetUp<Traits,dim>::get_right_hand_side () const
    {
      return right_hand_side;
    }


    template <class Traits, int dim>
    void
    SetUp<Traits,dim>::
    create_coarse_grid (Triangulation<dim> &coarse_grid) const
    {
      Traits::create_coarse_grid (coarse_grid);
    }


    // @sect4{The CurvedRidges class}

    // The class that is used to describe the boundary values and right hand
    // side of the <code>curved ridge</code> problem already used in the
    // step-13 example program is then like so:
    template <int dim>
    struct CurvedRidges
    {
      class BoundaryValues : public Function<dim>
      {
      public:
        BoundaryValues () : Function<dim> () {}

        virtual double value (const Point<dim>   &p,
                              const unsigned int  component) const;
      };


      class RightHandSide : public Function<dim>
      {
      public:
        RightHandSide () : Function<dim> () {}

        virtual double value (const Point<dim>   &p,
                              const unsigned int  component) const;
      };

      static
      void
      create_coarse_grid (Triangulation<dim> &coarse_grid);
    };


    template <int dim>
    double
    CurvedRidges<dim>::BoundaryValues::
    value (const Point<dim>   &p,
           const unsigned int  /*component*/) const
    {
      double q = p(0);
      for (unsigned int i=1; i<dim; ++i)
        q += std::sin(10*p(i)+5*p(0)*p(0));
      const double exponential = std::exp(q);
      return exponential;
    }



    template <int dim>
    double
    CurvedRidges<dim>::RightHandSide::value (const Point<dim>   &p,
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
        }
      t1 = t1*t1;

      return -u*(t1+t2+t3);
    }


    template <int dim>
    void
    CurvedRidges<dim>::
    create_coarse_grid (Triangulation<dim> &coarse_grid)
    {
      GridGenerator::hyper_cube (coarse_grid, -1, 1);
      coarse_grid.refine_global (2);
    }


    // @sect4{The Exercise_2_3 class}

    // This example program was written while giving practical courses for a
    // lecture on adaptive finite element methods and duality based error
    // estimates. For these courses, we had one exercise, which required to
    // solve the Laplace equation with constant right hand side on a square
    // domain with a square hole in the center, and zero boundary
    // values. Since the implementation of the properties of this problem is
    // so particularly simple here, lets do it. As the number of the exercise
    // was 2.3, we take the liberty to retain this name for the class as well.
    template <int dim>
    struct Exercise_2_3
    {
      // We need a class to denote the boundary values of the problem. In this
      // case, this is simple: it's the zero function, so don't even declare a
      // class, just a typedef:
      typedef ZeroFunction<dim> BoundaryValues;

      // Second, a class that denotes the right hand side. Since they are
      // constant, just subclass the corresponding class of the library and be
      // done:
      class RightHandSide : public ConstantFunction<dim>
      {
      public:
        RightHandSide () : ConstantFunction<dim> (1.) {}
      };

      // Finally a function to generate the coarse grid. This is somewhat more
      // complicated here, see immediately below.
      static
      void
      create_coarse_grid (Triangulation<dim> &coarse_grid);
    };


    // As stated above, the grid for this example is the square [-1,1]^2 with
    // the square [-1/2,1/2]^2 as hole in it. We create the coarse grid as 4
    // times 4 cells with the middle four ones missing.
    //
    // Of course, the example has an extension to 3d, but since this function
    // cannot be written in a dimension independent way we choose not to
    // implement this here, but rather only specialize the template for
    // dim=2. If you compile the program for 3d, you'll get a message from the
    // linker that this function is not implemented for 3d, and needs to be
    // provided.
    //
    // For the creation of this geometry, the library has no predefined
    // method. In this case, the geometry is still simple enough to do the
    // creation by hand, rather than using a mesh generator.
    template <>
    void
    Exercise_2_3<2>::
    create_coarse_grid (Triangulation<2> &coarse_grid)
    {
      // First define the space dimension, to allow those parts of the
      // function that are actually dimension independent to use this
      // variable. That makes it simpler if you later takes this as a starting
      // point to implement the 3d version.
      const unsigned int dim = 2;

      // Then have a list of vertices. Here, they are 24 (5 times 5, with the
      // middle one omitted). It is probably best to draw a sketch here. Note
      // that we leave the number of vertices open at first, but then let the
      // compiler compute this number afterwards. This reduces the possibility
      // of having the dimension to large and leaving the last ones
      // uninitialized.
      static const Point<2> vertices_1[]
        = {  Point<2> (-1.,   -1.),
             Point<2> (-1./2, -1.),
             Point<2> (0.,    -1.),
             Point<2> (+1./2, -1.),
             Point<2> (+1,    -1.),

             Point<2> (-1.,   -1./2.),
             Point<2> (-1./2, -1./2.),
             Point<2> (0.,    -1./2.),
             Point<2> (+1./2, -1./2.),
             Point<2> (+1,    -1./2.),

             Point<2> (-1.,   0.),
             Point<2> (-1./2, 0.),
             Point<2> (+1./2, 0.),
             Point<2> (+1,    0.),

             Point<2> (-1.,   1./2.),
             Point<2> (-1./2, 1./2.),
             Point<2> (0.,    1./2.),
             Point<2> (+1./2, 1./2.),
             Point<2> (+1,    1./2.),

             Point<2> (-1.,   1.),
             Point<2> (-1./2, 1.),
             Point<2> (0.,    1.),
             Point<2> (+1./2, 1.),
             Point<2> (+1,    1.)
          };
      const unsigned int
      n_vertices = sizeof(vertices_1) / sizeof(vertices_1[0]);

      // From this static list of vertices, we generate an STL vector of the
      // vertices, as this is the data type the library wants to see.
      const std::vector<Point<dim> > vertices (&vertices_1[0],
                                               &vertices_1[n_vertices]);

      // Next, we have to define the cells and the vertices they
      // contain. Here, we have 8 vertices, but leave the number open and let
      // it be computed afterwards:
      static const int cell_vertices[][GeometryInfo<dim>::vertices_per_cell]
      = {{0, 1, 5, 6},
        {1, 2, 6, 7},
        {2, 3, 7, 8},
        {3, 4, 8, 9},
        {5, 6, 10, 11},
        {8, 9, 12, 13},
        {10, 11, 14, 15},
        {12, 13, 17, 18},
        {14, 15, 19, 20},
        {15, 16, 20, 21},
        {16, 17, 21, 22},
        {17, 18, 22, 23}
      };
      const unsigned int
      n_cells = sizeof(cell_vertices) / sizeof(cell_vertices[0]);

      // Again, we generate a C++ vector type from this, but this time by
      // looping over the cells (yes, this is boring). Additionally, we set
      // the material indicator to zero for all the cells:
      std::vector<CellData<dim> > cells (n_cells, CellData<dim>());
      for (unsigned int i=0; i<n_cells; ++i)
        {
          for (unsigned int j=0;
               j<GeometryInfo<dim>::vertices_per_cell;
               ++j)
            cells[i].vertices[j] = cell_vertices[i][j];
          cells[i].material_id = 0;
        }

      // Finally pass all this information to the library to generate a
      // triangulation. The last parameter may be used to pass information
      // about non-zero boundary indicators at certain faces of the
      // triangulation to the library, but we don't want that here, so we give
      // an empty object:
      coarse_grid.create_triangulation (vertices,
                                        cells,
                                        SubCellData());

      // And since we want that the evaluation point (3/4,3/4) in this example
      // is a grid point, we refine once globally:
      coarse_grid.refine_global (1);
    }
  }

  // @sect4{Discussion}
  //
  // As you have now read through this framework, you may be wondering why we
  // have not chosen to implement the classes implementing a certain setup
  // (like the <code>CurvedRidges</code> class) directly as classes derived
  // from <code>Data::SetUpBase</code>. Indeed, we could have done very well
  // so. The only reason is that then we would have to have member variables
  // for the solution and right hand side classes in the
  // <code>CurvedRidges</code> class, as well as member functions overloading
  // the abstract functions of the base class giving access to these member
  // variables. The <code>SetUp</code> class has the sole reason to relieve us
  // from the need to reiterate these member variables and functions that
  // would be necessary in all such classes. In some way, the template
  // mechanism here only provides a way to have default implementations for a
  // number of functions that depend on external quantities and can thus not
  // be provided using normal virtual functions, at least not without the help
  // of templates.
  //
  // However, there might be good reasons to actually implement classes
  // derived from <code>Data::SetUpBase</code>, for example if the solution or
  // right hand side classes require constructors that take arguments, which
  // the <code>Data::SetUpBase</code> class cannot provide. In that case,
  // subclassing is a worthwhile strategy. Other possibilities for special
  // cases are to derive from <code>Data::SetUp@<SomeSetUp@></code> where
  // <code>SomeSetUp</code> denotes a class, or even to explicitly specialize
  // <code>Data::SetUp@<SomeSetUp@></code>. The latter allows to transparently
  // use the way the <code>SetUp</code> class is used for other set-ups, but
  // with special actions taken for special arguments.
  //
  // A final observation favoring the approach taken here is the following: we
  // have found numerous times that when starting a project, the number of
  // parameters (usually boundary values, right hand side, coarse grid, just
  // as here) was small, and the number of test cases was small as well. One
  // then starts out by handcoding them into a number of <code>switch</code>
  // statements. Over time, projects grow, and so does the number of test
  // cases. The number of <code>switch</code> statements grows with that, and
  // their length as well, and one starts to find ways to consider impossible
  // examples where domains, boundary values, and right hand sides do not fit
  // together any more, and starts losing the overview over the whole
  // structure. Encapsulating everything belonging to a certain test case into
  // a structure of its own has proven worthwhile for this, as it keeps
  // everything that belongs to one test case in one place. Furthermore, it
  // allows to put these things all in one or more files that are only devoted
  // to test cases and their data, without having to bring their actual
  // implementation into contact with the rest of the program.


  // @sect3{Dual functionals}

  // As with the other components of the program, we put everything we need to
  // describe dual functionals into a namespace of its own, and define an
  // abstract base class that provides the interface the class solving the
  // dual problem needs for its work.
  //
  // We will then implement two such classes, for the evaluation of a point
  // value and of the derivative of the solution at that point. For these
  // functionals we already have the corresponding evaluation objects, so they
  // are complementary.
  namespace DualFunctional
  {
    // @sect4{The DualFunctionalBase class}

    // First start with the base class for dual functionals. Since for linear
    // problems the characteristics of the dual problem play a role only in
    // the right hand side, we only need to provide for a function that
    // assembles the right hand side for a given discretization:
    template <int dim>
    class DualFunctionalBase : public Subscriptor
    {
    public:
      virtual
      void
      assemble_rhs (const DoFHandler<dim> &dof_handler,
                    Vector<double>        &rhs) const = 0;
    };


    // @sect4{The PointValueEvaluation class}

    // As a first application, we consider the functional corresponding to the
    // evaluation of the solution's value at a given point which again we
    // assume to be a vertex. Apart from the constructor that takes and stores
    // the evaluation point, this class consists only of the function that
    // implements assembling the right hand side.
    template <int dim>
    class PointValueEvaluation : public DualFunctionalBase<dim>
    {
    public:
      PointValueEvaluation (const Point<dim> &evaluation_point);

      virtual
      void
      assemble_rhs (const DoFHandler<dim> &dof_handler,
                    Vector<double>        &rhs) const;

      DeclException1 (ExcEvaluationPointNotFound,
                      Point<dim>,
                      << "The evaluation point " << arg1
                      << " was not found among the vertices of the present grid.");

    protected:
      const Point<dim> evaluation_point;
    };


    template <int dim>
    PointValueEvaluation<dim>::
    PointValueEvaluation (const Point<dim> &evaluation_point)
      :
      evaluation_point (evaluation_point)
    {}


    // As for doing the main purpose of the class, assembling the right hand
    // side, let us first consider what is necessary: The right hand side of
    // the dual problem is a vector of values J(phi_i), where J is the error
    // functional, and phi_i is the i-th shape function. Here, J is the
    // evaluation at the point x0, i.e. J(phi_i)=phi_i(x0).
    //
    // Now, we have assumed that the evaluation point is a vertex. Thus, for
    // the usual finite elements we might be using in this program, we can
    // take for granted that at such a point exactly one shape function is
    // nonzero, and in particular has the value one. Thus, we set the right
    // hand side vector to all-zeros, then seek for the shape function
    // associated with that point and set the corresponding value of the right
    // hand side vector to one:
    template <int dim>
    void
    PointValueEvaluation<dim>::
    assemble_rhs (const DoFHandler<dim> &dof_handler,
                  Vector<double>        &rhs) const
    {
      // So, first set everything to zeros...
      rhs.reinit (dof_handler.n_dofs());

      // ...then loop over cells and find the evaluation point among the
      // vertices (or very close to a vertex, which may happen due to floating
      // point round-off):
      typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
      for (; cell!=endc; ++cell)
        for (unsigned int vertex=0;
             vertex<GeometryInfo<dim>::vertices_per_cell;
             ++vertex)
          if (cell->vertex(vertex).distance(evaluation_point)
              < cell->diameter()*1e-8)
            {
              // Ok, found, so set corresponding entry, and leave function
              // since we are finished:
              rhs(cell->vertex_dof_index(vertex,0)) = 1;
              return;
            }

      // Finally, a sanity check: if we somehow got here, then we must have
      // missed the evaluation point, so raise an exception unconditionally:
      AssertThrow (false, ExcEvaluationPointNotFound(evaluation_point));
    }


    // @sect4{The PointXDerivativeEvaluation class}

    // As second application, we again consider the evaluation of the
    // x-derivative of the solution at one point. Again, the declaration of
    // the class, and the implementation of its constructor is not too
    // interesting:
    template <int dim>
    class PointXDerivativeEvaluation : public DualFunctionalBase<dim>
    {
    public:
      PointXDerivativeEvaluation (const Point<dim> &evaluation_point);

      virtual
      void
      assemble_rhs (const DoFHandler<dim> &dof_handler,
                    Vector<double>        &rhs) const;

      DeclException1 (ExcEvaluationPointNotFound,
                      Point<dim>,
                      << "The evaluation point " << arg1
                      << " was not found among the vertices of the present grid.");

    protected:
      const Point<dim> evaluation_point;
    };


    template <int dim>
    PointXDerivativeEvaluation<dim>::
    PointXDerivativeEvaluation (const Point<dim> &evaluation_point)
      :
      evaluation_point (evaluation_point)
    {}


    // What is interesting is the implementation of this functional: here,
    // J(phi_i)=d/dx phi_i(x0).
    //
    // We could, as in the implementation of the respective evaluation object
    // take the average of the gradients of each shape function phi_i at this
    // evaluation point. However, we take a slightly different approach: we
    // simply take the average over all cells that surround this point. The
    // question which cells <code>surrounds</code> the evaluation point is
    // made dependent on the mesh width by including those cells for which the
    // distance of the cell's midpoint to the evaluation point is less than
    // the cell's diameter.
    //
    // Taking the average of the gradient over the area/volume of these cells
    // leads to a dual solution which is very close to the one which would
    // result from the point evaluation of the gradient. It is simple to
    // justify theoretically that this does not change the method
    // significantly.
    template <int dim>
    void
    PointXDerivativeEvaluation<dim>::
    assemble_rhs (const DoFHandler<dim> &dof_handler,
                  Vector<double>        &rhs) const
    {
      // Again, first set all entries to zero:
      rhs.reinit (dof_handler.n_dofs());

      // Initialize a <code>FEValues</code> object with a quadrature formula,
      // have abbreviations for the number of quadrature points and shape
      // functions...
      QGauss<dim> quadrature(4);
      FEValues<dim>  fe_values (dof_handler.get_fe(), quadrature,
                                update_gradients |
                                update_quadrature_points  |
                                update_JxW_values);
      const unsigned int n_q_points = fe_values.n_quadrature_points;
      const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

      // ...and have two objects that are used to store the global indices of
      // the degrees of freedom on a cell, and the values of the gradients of
      // the shape functions at the quadrature points:
      Vector<double> cell_rhs (dofs_per_cell);
      std::vector<unsigned int> local_dof_indices (dofs_per_cell);

      // Finally have a variable in which we will sum up the area/volume of
      // the cells over which we integrate, by integrating the unit functions
      // on these cells:
      double total_volume = 0;

      // Then start the loop over all cells, and select those cells which are
      // close enough to the evaluation point:
      typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
      for (; cell!=endc; ++cell)
        if (cell->center().distance(evaluation_point) <=
            cell->diameter())
          {
            // If we have found such a cell, then initialize the
            // <code>FEValues</code> object and integrate the x-component of
            // the gradient of each shape function, as well as the unit
            // function for the total area/volume.
            fe_values.reinit (cell);
            cell_rhs = 0;

            for (unsigned int q=0; q<n_q_points; ++q)
              {
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  cell_rhs(i) += fe_values.shape_grad(i,q)[0] *
                                 fe_values.JxW (q);
                total_volume += fe_values.JxW (q);
              }

            // If we have the local contributions, distribute them to the
            // global vector:
            cell->get_dof_indices (local_dof_indices);
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              rhs(local_dof_indices[i]) += cell_rhs(i);
          }

      // After we have looped over all cells, check whether we have found any
      // at all, by making sure that their volume is non-zero. If not, then
      // the results will be botched, as the right hand side should then still
      // be zero, so throw an exception:
      AssertThrow (total_volume > 0,
                   ExcEvaluationPointNotFound(evaluation_point));

      // Finally, we have by now only integrated the gradients of the shape
      // functions, not taking their mean value. We fix this by dividing by
      // the measure of the volume over which we have integrated:
      rhs /= total_volume;
    }


  }


  // @sect3{Extending the LaplaceSolver namespace}
  namespace LaplaceSolver
  {

    // @sect4{The DualSolver class}

    // In the same way as the <code>PrimalSolver</code> class above, we now
    // implement a <code>DualSolver</code>. It has all the same features, the
    // only difference is that it does not take a function object denoting a
    // right hand side object, but now takes a <code>DualFunctionalBase</code>
    // object that will assemble the right hand side vector of the dual
    // problem. The rest of the class is rather trivial.
    //
    // Since both primal and dual solver will use the same triangulation, but
    // different discretizations, it now becomes clear why we have made the
    // <code>Base</code> class a virtual one: since the final class will be
    // derived from both <code>PrimalSolver</code> as well as
    // <code>DualSolver</code>, it would have two <code>Base</code> instances,
    // would we not have marked the inheritance as virtual. Since in many
    // applications the base class would store much more information than just
    // the triangulation which needs to be shared between primal and dual
    // solvers, we do not usually want to use two such base classes.
    template <int dim>
    class DualSolver : public Solver<dim>
    {
    public:
      DualSolver (Triangulation<dim>       &triangulation,
                  const FiniteElement<dim> &fe,
                  const Quadrature<dim>    &quadrature,
                  const Quadrature<dim-1>  &face_quadrature,
                  const DualFunctional::DualFunctionalBase<dim> &dual_functional);

      virtual
      void
      solve_problem ();

      virtual
      unsigned int
      n_dofs () const;

      virtual
      void
      postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const;

    protected:
      const SmartPointer<const DualFunctional::DualFunctionalBase<dim> > dual_functional;
      virtual void assemble_rhs (Vector<double> &rhs) const;

      static const ZeroFunction<dim> boundary_values;

      // Same as above -- make a derived class a friend of this one:
      friend class WeightedResidual<dim>;
    };

    template <int dim>
    const ZeroFunction<dim> DualSolver<dim>::boundary_values;

    template <int dim>
    DualSolver<dim>::
    DualSolver (Triangulation<dim>       &triangulation,
                const FiniteElement<dim> &fe,
                const Quadrature<dim>    &quadrature,
                const Quadrature<dim-1>  &face_quadrature,
                const DualFunctional::DualFunctionalBase<dim> &dual_functional)
      :
      Base<dim> (triangulation),
      Solver<dim> (triangulation, fe,
                   quadrature, face_quadrature,
                   boundary_values),
      dual_functional (&dual_functional)
    {}


    template <int dim>
    void
    DualSolver<dim>::solve_problem ()
    {
      Solver<dim>::solve_problem ();
    }



    template <int dim>
    unsigned int
    DualSolver<dim>::n_dofs() const
    {
      return Solver<dim>::n_dofs();
    }


    template <int dim>
    void
    DualSolver<dim>::
    postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const
    {
      Solver<dim>::postprocess(postprocessor);
    }



    template <int dim>
    void
    DualSolver<dim>::
    assemble_rhs (Vector<double> &rhs) const
    {
      dual_functional->assemble_rhs (this->dof_handler, rhs);
    }


    // @sect4{The WeightedResidual class}

    // Here finally comes the main class of this program, the one that
    // implements the dual weighted residual error estimator. It joins the
    // primal and dual solver classes to use them for the computation of
    // primal and dual solutions, and implements the error representation
    // formula for use as error estimate and mesh refinement.
    //
    // The first few of the functions of this class are mostly overriders of
    // the respective functions of the base class:
    template <int dim>
    class WeightedResidual : public PrimalSolver<dim>,
      public DualSolver<dim>
    {
    public:
      WeightedResidual (Triangulation<dim>       &coarse_grid,
                        const FiniteElement<dim> &primal_fe,
                        const FiniteElement<dim> &dual_fe,
                        const Quadrature<dim>    &quadrature,
                        const Quadrature<dim-1>  &face_quadrature,
                        const Function<dim>      &rhs_function,
                        const Function<dim>      &boundary_values,
                        const DualFunctional::DualFunctionalBase<dim> &dual_functional);

      virtual
      void
      solve_problem ();

      virtual
      void
      postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const;

      virtual
      unsigned int
      n_dofs () const;

      virtual void refine_grid ();

      virtual
      void
      output_solution () const;

    private:
      // In the private section, we have two functions that are used to call
      // the <code>solve_problem</code> functions of the primal and dual base
      // classes. These two functions will be called in parallel by the
      // <code>solve_problem</code> function of this class.
      void solve_primal_problem ();
      void solve_dual_problem ();
      // Then declare abbreviations for active cell iterators, to avoid that
      // we have to write this lengthy name over and over again:

      typedef
      typename DoFHandler<dim>::active_cell_iterator
      active_cell_iterator;

      // Next, declare a data type that we will us to store the contribution
      // of faces to the error estimator. The idea is that we can compute the
      // face terms from each of the two cells to this face, as they are the
      // same when viewed from both sides. What we will do is to compute them
      // only once, based on some rules explained below which of the two
      // adjacent cells will be in charge to do so. We then store the
      // contribution of each face in a map mapping faces to their values, and
      // only collect the contributions for each cell by looping over the
      // cells a second time and grabbing the values from the map.
      //
      // The data type of this map is declared here:
      typedef
      typename std::map<typename DoFHandler<dim>::face_iterator,double>
      FaceIntegrals;

      // In the computation of the error estimates on cells and faces, we need
      // a number of helper objects, such as <code>FEValues</code> and
      // <code>FEFaceValues</code> functions, but also temporary objects
      // storing the values and gradients of primal and dual solutions, for
      // example. These fields are needed in the three functions that do the
      // integration on cells, and regular and irregular faces, respectively.
      //
      // There are three reasonable ways to provide these fields: first, as
      // local variables in the function that needs them; second, as member
      // variables of this class; third, as arguments passed to that function.
      //
      // These three alternatives all have drawbacks: the third that their
      // number is not negligible and would make calling these functions a
      // lengthy enterprise. The second has the drawback that it disallows
      // parallelization, since the threads that will compute the error
      // estimate have to have their own copies of these variables each, so
      // member variables of the enclosing class will not work. The first
      // approach, although straightforward, has a subtle but important
      // drawback: we will call these functions over and over again, many
      // thousands of times maybe; it now turns out that allocating
      // vectors and other objects that need memory from the heap is an
      // expensive business in terms of run-time, since memory allocation is
      // expensive when several threads are involved. It is thus
      // significantly better to allocate the memory only once, and recycle
      // the objects as often as possible.
      //
      // What to do? Our answer is to use a variant of the third strategy.
      // In fact, this is exactly what the WorkStream concept is supposed to
      // do (we have already introduced it above, but see also @ref threads).
      // To avoid that we have to give these functions a dozen or so
      // arguments, we pack all these variables into two structures, one which
      // is used for the computations on cells, the other doing them on the
      // faces. Both are then joined into the WeightedResidualScratchData class
      // that will serve as the "scratch data" class of the WorkStream concept:
      struct CellData
      {
        FEValues<dim>    fe_values;
        const SmartPointer<const Function<dim> > right_hand_side;

        std::vector<double> cell_residual;
        std::vector<double> rhs_values;
        std::vector<double> dual_weights;
        std::vector<double> cell_laplacians;
        CellData (const FiniteElement<dim> &fe,
                  const Quadrature<dim>    &quadrature,
                  const Function<dim>      &right_hand_side);
        CellData (const CellData &cell_data);
      };

      struct FaceData
      {
        FEFaceValues<dim>    fe_face_values_cell;
        FEFaceValues<dim>    fe_face_values_neighbor;
        FESubfaceValues<dim> fe_subface_values_cell;

        std::vector<double> jump_residual;
        std::vector<double> dual_weights;
        typename std::vector<Tensor<1,dim> > cell_grads;
        typename std::vector<Tensor<1,dim> > neighbor_grads;
        FaceData (const FiniteElement<dim> &fe,
                  const Quadrature<dim-1>  &face_quadrature);
        FaceData (const FaceData &face_data);
      };

      struct WeightedResidualScratchData
      {
        WeightedResidualScratchData(const PrimalSolver<dim> &primal_solver,
                                    const DualSolver<dim>   &dual_solver,
                                    const Vector<double>    &primal_solution,
                                    const Vector<double>    &dual_weights);

        WeightedResidualScratchData(const WeightedResidualScratchData &scratch_data);

        CellData       cell_data;
        FaceData       face_data;
        Vector<double> primal_solution;
        Vector<double> dual_weights;
      };


      // WorkStream::run generally wants both a scratch object and a copy object.
      // Here, for reasons similar to what we had in step-9 when discussing the
      // computation of an approximation of the gradient, we don't actually
      // need a "copy data" structure. Since WorkStream insists on having one of
      // these, we just declare an empty structure that does nothing other than
      // being there.
      struct WeightedResidualCopyData
      {};



      // Regarding the evaluation of the error estimator, we have one driver
      // function that uses WorkStream::run to call the second function on every
      // cell. The concept of using SynchronousIterators was already explained
      // in step-9:
      void estimate_error (Vector<float> &error_indicators) const;

      void estimate_on_one_cell (const SynchronousIterators<std_cxx11::tuple<
                                 active_cell_iterator,Vector<float>::iterator> > &cell_and_error,
                                 WeightedResidualScratchData                     &scratch_data,
                                 WeightedResidualCopyData                        &copy_data,
                                 FaceIntegrals                                   &face_integrals) const;

      // Then we have functions that do the actual integration of the error
      // representation formula. They will treat the terms on the cell
      // interiors, on those faces that have no hanging nodes, and on those
      // faces with hanging nodes, respectively:
      void
      integrate_over_cell (const SynchronousIterators<std_cxx11::tuple<
                           active_cell_iterator,Vector<float>::iterator> > &cell_and_error,
                           const Vector<double>                            &primal_solution,
                           const Vector<double>                            &dual_weights,
                           CellData                                        &cell_data) const;

      void
      integrate_over_regular_face (const active_cell_iterator &cell,
                                   const unsigned int          face_no,
                                   const Vector<double>       &primal_solution,
                                   const Vector<double>       &dual_weights,
                                   FaceData                   &face_data,
                                   FaceIntegrals              &face_integrals) const;
      void
      integrate_over_irregular_face (const active_cell_iterator &cell,
                                     const unsigned int          face_no,
                                     const Vector<double>       &primal_solution,
                                     const Vector<double>       &dual_weights,
                                     FaceData                   &face_data,
                                     FaceIntegrals              &face_integrals) const;
    };



    // In the implementation of this class, we first have the constructors of
    // the <code>CellData</code> and <code>FaceData</code> member classes, and
    // the <code>WeightedResidual</code> constructor. They only initialize
    // fields to their correct lengths, so we do not have to discuss them in
    // too much detail:
    template <int dim>
    WeightedResidual<dim>::CellData::
    CellData (const FiniteElement<dim> &fe,
              const Quadrature<dim>    &quadrature,
              const Function<dim>      &right_hand_side)
      :
      fe_values (fe, quadrature,
                 update_values   |
                 update_hessians |
                 update_quadrature_points |
                 update_JxW_values),
      right_hand_side (&right_hand_side),
      cell_residual (quadrature.size()),
      rhs_values (quadrature.size()),
      dual_weights (quadrature.size()),
      cell_laplacians (quadrature.size())
    {}



    template <int dim>
    WeightedResidual<dim>::CellData::
    CellData (const CellData &cell_data)
      :
      fe_values (cell_data.fe_values.get_fe(),
                 cell_data.fe_values.get_quadrature(),
                 update_values   |
                 update_hessians |
                 update_quadrature_points |
                 update_JxW_values),
      right_hand_side (cell_data.right_hand_side),
      cell_residual (cell_data.cell_residual),
      rhs_values (cell_data.rhs_values),
      dual_weights (cell_data.dual_weights),
      cell_laplacians (cell_data.cell_laplacians)
    {}



    template <int dim>
    WeightedResidual<dim>::FaceData::
    FaceData (const FiniteElement<dim> &fe,
              const Quadrature<dim-1>  &face_quadrature)
      :
      fe_face_values_cell (fe, face_quadrature,
                           update_values        |
                           update_gradients     |
                           update_JxW_values    |
                           update_normal_vectors),
      fe_face_values_neighbor (fe, face_quadrature,
                               update_values     |
                               update_gradients  |
                               update_JxW_values |
                               update_normal_vectors),
      fe_subface_values_cell (fe, face_quadrature,
                              update_gradients)
    {
      const unsigned int n_face_q_points
        = face_quadrature.size();

      jump_residual.resize(n_face_q_points);
      dual_weights.resize(n_face_q_points);
      cell_grads.resize(n_face_q_points);
      neighbor_grads.resize(n_face_q_points);
    }



    template <int dim>
    WeightedResidual<dim>::FaceData::
    FaceData (const FaceData &face_data)
      :
      fe_face_values_cell (face_data.fe_face_values_cell.get_fe(),
                           face_data.fe_face_values_cell.get_quadrature(),
                           update_values        |
                           update_gradients     |
                           update_JxW_values    |
                           update_normal_vectors),
      fe_face_values_neighbor (face_data.fe_face_values_neighbor.get_fe(),
                               face_data.fe_face_values_neighbor.get_quadrature(),
                               update_values     |
                               update_gradients  |
                               update_JxW_values |
                               update_normal_vectors),
      fe_subface_values_cell (face_data.fe_subface_values_cell.get_fe(),
                              face_data.fe_subface_values_cell.get_quadrature(),
                              update_gradients),
      jump_residual (face_data.jump_residual),
      dual_weights (face_data.dual_weights),
      cell_grads (face_data.cell_grads),
      neighbor_grads (face_data.neighbor_grads)
    {}



    template <int dim>
    WeightedResidual<dim>::WeightedResidualScratchData::
    WeightedResidualScratchData (const PrimalSolver<dim> &primal_solver,
                                 const DualSolver<dim>   &dual_solver,
                                 const Vector<double>    &primal_solution,
                                 const Vector<double>    &dual_weights)
      :
      cell_data (*dual_solver.fe,
                 *dual_solver.quadrature,
                 *primal_solver.rhs_function),
      face_data (*dual_solver.fe,
                 *dual_solver.face_quadrature),
      primal_solution(primal_solution),
      dual_weights(dual_weights)
    {}

    template <int dim>
    WeightedResidual<dim>::WeightedResidualScratchData::
    WeightedResidualScratchData (const WeightedResidualScratchData &scratch_data)
      :
      cell_data(scratch_data.cell_data),
      face_data(scratch_data.face_data),
      primal_solution(scratch_data.primal_solution),
      dual_weights(scratch_data.dual_weights)
    {}



    template <int dim>
    WeightedResidual<dim>::
    WeightedResidual (Triangulation<dim>       &coarse_grid,
                      const FiniteElement<dim> &primal_fe,
                      const FiniteElement<dim> &dual_fe,
                      const Quadrature<dim>    &quadrature,
                      const Quadrature<dim-1>  &face_quadrature,
                      const Function<dim>      &rhs_function,
                      const Function<dim>      &bv,
                      const DualFunctional::DualFunctionalBase<dim> &dual_functional)
      :
      Base<dim> (coarse_grid),
      PrimalSolver<dim> (coarse_grid, primal_fe,
                         quadrature, face_quadrature,
                         rhs_function, bv),
      DualSolver<dim> (coarse_grid, dual_fe,
                       quadrature, face_quadrature,
                       dual_functional)
    {}


    // The next five functions are boring, as they simply relay their work to
    // the base classes. The first calls the primal and dual solvers in
    // parallel, while postprocessing the solution and retrieving the number
    // of degrees of freedom is done by the primal class.
    template <int dim>
    void
    WeightedResidual<dim>::solve_problem ()
    {
      Threads::TaskGroup<> tasks;
      tasks += Threads::new_task (&WeightedResidual<dim>::solve_primal_problem,
                                  *this);
      tasks += Threads::new_task (&WeightedResidual<dim>::solve_dual_problem,
                                  *this);
      tasks.join_all();
    }


    template <int dim>
    void
    WeightedResidual<dim>::solve_primal_problem ()
    {
      PrimalSolver<dim>::solve_problem ();
    }

    template <int dim>
    void
    WeightedResidual<dim>::solve_dual_problem ()
    {
      DualSolver<dim>::solve_problem ();
    }


    template <int dim>
    void
    WeightedResidual<dim>::
    postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const
    {
      PrimalSolver<dim>::postprocess (postprocessor);
    }


    template <int dim>
    unsigned int
    WeightedResidual<dim>::n_dofs () const
    {
      return PrimalSolver<dim>::n_dofs();
    }



    // Now, it is becoming more interesting: the <code>refine_grid()</code>
    // function asks the error estimator to compute the cell-wise error
    // indicators, then uses their absolute values for mesh refinement.
    template <int dim>
    void
    WeightedResidual<dim>::refine_grid ()
    {
      // First call the function that computes the cell-wise and global error:
      Vector<float> error_indicators (this->triangulation->n_active_cells());
      estimate_error (error_indicators);

      // Then note that marking cells for refinement or coarsening only works
      // if all indicators are positive, to allow their comparison. Thus, drop
      // the signs on all these indicators:
      for (Vector<float>::iterator i=error_indicators.begin();
           i != error_indicators.end(); ++i)
        *i = std::fabs (*i);

      // Finally, we can select between different strategies for
      // refinement. The default here is to refine those cells with the
      // largest error indicators that make up for a total of 80 per cent of
      // the error, while we coarsen those with the smallest indicators that
      // make up for the bottom 2 per cent of the error.
      GridRefinement::refine_and_coarsen_fixed_fraction (*this->triangulation,
                                                         error_indicators,
                                                         0.8, 0.02);
      this->triangulation->execute_coarsening_and_refinement ();
    }


    // Since we want to output both the primal and the dual solution, we
    // overload the <code>output_solution</code> function. The only
    // interesting feature of this function is that the primal and dual
    // solutions are defined on different finite element spaces, which is not
    // the format the <code>DataOut</code> class expects. Thus, we have to
    // transfer them to a common finite element space. Since we want the
    // solutions only to see them qualitatively, we contend ourselves with
    // interpolating the dual solution to the (smaller) primal space. For the
    // interpolation, there is a library function, that takes a
    // <code>ConstraintMatrix</code> object including the hanging node
    // constraints. The rest is standard.
    //
    // There is, however, one work-around worth mentioning: in this function,
    // as in a couple of following ones, we have to access the
    // <code>DoFHandler</code> objects and solutions of both the primal as
    // well as of the dual solver. Since these are members of the
    // <code>Solver</code> base class which exists twice in the class
    // hierarchy leading to the present class (once as base class of the
    // <code>PrimalSolver</code> class, once as base class of the
    // <code>DualSolver</code> class), we have to disambiguate accesses to
    // them by telling the compiler a member of which of these two instances
    // we want to access. The way to do this would be identify the member by
    // pointing a path through the class hierarchy which disambiguates the
    // base class, for example writing <code>PrimalSolver::dof_handler</code>
    // to denote the member variable <code>dof_handler</code> from the
    // <code>Solver</code> base class of the <code>PrimalSolver</code>
    // class. Unfortunately, this confuses gcc's version 2.96 (a version that
    // was intended as a development snapshot, but delivered as system
    // compiler by Red Hat in their 7.x releases) so much that it bails out
    // and refuses to compile the code.
    //
    // Thus, we have to work around this problem. We do this by introducing
    // references to the <code>PrimalSolver</code> and <code>DualSolver</code>
    // components of the <code>WeightedResidual</code> object at the beginning
    // of the function. Since each of these has an unambiguous base class
    // <code>Solver</code>, we can access the member variables we want through
    // these references. However, we are now accessing protected member
    // variables of these classes through a pointer other than the
    // <code>this</code> pointer (in fact, this is of course the
    // <code>this</code> pointer, but not explicitly). This finally is the
    // reason why we had to declare the present class a friend of the classes
    // we so access.
    template <int dim>
    void
    WeightedResidual<dim>::output_solution () const
    {
      const PrimalSolver<dim> &primal_solver = *this;
      const DualSolver<dim>   &dual_solver   = *this;

      ConstraintMatrix primal_hanging_node_constraints;
      DoFTools::make_hanging_node_constraints (primal_solver.dof_handler,
                                               primal_hanging_node_constraints);
      primal_hanging_node_constraints.close();
      Vector<double> dual_solution (primal_solver.dof_handler.n_dofs());
      FETools::interpolate (dual_solver.dof_handler,
                            dual_solver.solution,
                            primal_solver.dof_handler,
                            primal_hanging_node_constraints,
                            dual_solution);

      DataOut<dim> data_out;
      data_out.attach_dof_handler (primal_solver.dof_handler);

      // Add the data vectors for which we want output. Add them both, the
      // <code>DataOut</code> functions can handle as many data vectors as you
      // wish to write to output:
      data_out.add_data_vector (primal_solver.solution,
                                "primal_solution");
      data_out.add_data_vector (dual_solution,
                                "dual_solution");

      data_out.build_patches ();

      std::ostringstream filename;
      filename << "solution-"
               << this->refinement_cycle
               << ".gnuplot"
               << std::ends;

      std::ofstream out (filename.str().c_str());
      data_out.write (out, DataOutBase::gnuplot);
    }


    // @sect3{Estimating errors}

    // @sect4{Error estimation driver functions}
    //
    // As for the actual computation of error estimates, let's start with the
    // function that drives all this, i.e. calls those functions that actually
    // do the work, and finally collects the results.
    template <int dim>
    void
    WeightedResidual<dim>::
    estimate_error (Vector<float> &error_indicators) const
    {
      const PrimalSolver<dim> &primal_solver = *this;
      const DualSolver<dim>   &dual_solver   = *this;

      // The first task in computing the error is to set up vectors that
      // denote the primal solution, and the weights (z-z_h)=(z-I_hz), both in
      // the finite element space for which we have computed the dual
      // solution. For this, we have to interpolate the primal solution to the
      // dual finite element space, and to subtract the interpolation of the
      // computed dual solution to the primal finite element
      // space. Fortunately, the library provides functions for the
      // interpolation into larger or smaller finite element spaces, so this
      // is mostly obvious.
      //
      // First, let's do that for the primal solution: it is cell-wise
      // interpolated into the finite element space in which we have solved
      // the dual problem: But, again as in the
      // <code>WeightedResidual::output_solution</code> function we first need
      // to create a ConstraintMatrix including the hanging node constraints,
      // but this time of the dual finite element space.
      ConstraintMatrix dual_hanging_node_constraints;
      DoFTools::make_hanging_node_constraints (dual_solver.dof_handler,
                                               dual_hanging_node_constraints);
      dual_hanging_node_constraints.close();
      Vector<double> primal_solution (dual_solver.dof_handler.n_dofs());
      FETools::interpolate (primal_solver.dof_handler,
                            primal_solver.solution,
                            dual_solver.dof_handler,
                            dual_hanging_node_constraints,
                            primal_solution);

      // Then for computing the interpolation of the numerically approximated
      // dual solution z into the finite element space of the primal solution
      // and subtracting it from z: use the
      // <code>interpolate_difference</code> function, that gives (z-I_hz) in
      // the element space of the dual solution.
      ConstraintMatrix primal_hanging_node_constraints;
      DoFTools::make_hanging_node_constraints (primal_solver.dof_handler,
                                               primal_hanging_node_constraints);
      primal_hanging_node_constraints.close();
      Vector<double> dual_weights (dual_solver.dof_handler.n_dofs());
      FETools::interpolation_difference (dual_solver.dof_handler,
                                         dual_hanging_node_constraints,
                                         dual_solver.solution,
                                         primal_solver.dof_handler,
                                         primal_hanging_node_constraints,
                                         dual_weights);

      // Note that this could probably have been more efficient since those
      // constraints have been used previously when assembling matrix and
      // right hand side for the primal problem and writing out the dual
      // solution. We leave the optimization of the program in this respect as
      // an exercise.

      // Having computed the dual weights we now proceed with computing the
      // cell and face residuals of the primal solution. First we set up a map
      // between face iterators and their jump term contributions of faces to
      // the error estimator. The reason is that we compute the jump terms
      // only once, from one side of the face, and want to collect them only
      // afterwards when looping over all cells a second time.
      //
      // We initialize this map already with a value of -1e20 for all faces,
      // since this value will stand out in the results if something should go
      // wrong and we fail to compute the value for a face for some
      // reason. Secondly, this initialization already makes the std::map
      // object allocate all objects it may possibly need. This is important
      // since we will write into this structure from parallel threads,
      // and doing so would not be thread-safe if the map needed to allocate
      // memory and thereby reshape its data structures. In other words, the
      // initial initialization relieves us from the necessity to synchronize the
      // threads through a mutex each time they write to (and modify the
      // structure of) this map.
      FaceIntegrals face_integrals;
      for (active_cell_iterator cell=dual_solver.dof_handler.begin_active();
           cell!=dual_solver.dof_handler.end();
           ++cell)
        for (unsigned int face_no=0;
             face_no<GeometryInfo<dim>::faces_per_cell;
             ++face_no)
          face_integrals[cell->face(face_no)] = -1e20;

      // Then set up a vector with error indicators and reserve one slot for
      // each cell and set it to zero. With this, we can then set up the
      // parallel iterator range just as we did in step-9, and hand it
      // all off to WorkStream::run to compute the estimators for all
      // cells in parallel:
      error_indicators.reinit (dual_solver.dof_handler
                               .get_tria().n_active_cells());

      typedef
      std_cxx11::tuple<active_cell_iterator,Vector<float>::iterator>
      IteratorTuple;

      SynchronousIterators<IteratorTuple>
      cell_and_error_begin(IteratorTuple (dual_solver.dof_handler.begin_active(),
                                          error_indicators.begin()));
      SynchronousIterators<IteratorTuple>
      cell_and_error_end  (IteratorTuple (dual_solver.dof_handler.end(),
                                          error_indicators.begin()));

      WorkStream::run(cell_and_error_begin,
                      cell_and_error_end,
                      std_cxx11::bind(&WeightedResidual<dim>::estimate_on_one_cell,
                                      this,
                                      std_cxx11::_1,
                                      std_cxx11::_2,
                                      std_cxx11::_3,
                                      std_cxx11::ref(face_integrals)),
                      std_cxx11::function<void (const WeightedResidualCopyData &)>(),
                      WeightedResidualScratchData (primal_solver,
                                                   dual_solver,
                                                   primal_solution,
                                                   dual_weights),
                      WeightedResidualCopyData());

      // Once the error contributions are computed, sum them up. For this,
      // note that the cell terms are already set, and that only the edge
      // terms need to be collected. Thus, loop over all cells and their
      // faces, make sure that the contributions of each of the faces are
      // there, and add them up. Only take minus one half of the jump term,
      // since the other half will be taken by the neighboring cell.
      unsigned int present_cell=0;
      for (active_cell_iterator cell=dual_solver.dof_handler.begin_active();
           cell!=dual_solver.dof_handler.end();
           ++cell, ++present_cell)
        for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell;
             ++face_no)
          {
            Assert(face_integrals.find(cell->face(face_no)) !=
                   face_integrals.end(),
                   ExcInternalError());
            error_indicators(present_cell)
            -= 0.5*face_integrals[cell->face(face_no)];
          }
      std::cout << "   Estimated error="
                << std::accumulate (error_indicators.begin(),
                                    error_indicators.end(), 0.)
                << std::endl;
    }


    // @sect4{Estimating on a single cell}

    // Next we have the function that is called to estimate the error on a
    // single cell. The function may be called multiple times if the library was
    // configured to use multithreading. Here it goes:
    template <int dim>
    void
    WeightedResidual<dim>::
    estimate_on_one_cell (const SynchronousIterators<std_cxx11::tuple<
                          active_cell_iterator,Vector<float>::iterator> > &cell_and_error,
                          WeightedResidualScratchData                       &scratch_data,
                          WeightedResidualCopyData                          &copy_data,
                          FaceIntegrals                                     &face_integrals) const
    {
      // First task on each cell is to compute the cell residual
      // contributions of this cell, and put them into the
      // <code>error_indicators</code> variable:
      active_cell_iterator cell = std_cxx11::get<0>(cell_and_error.iterators);

      integrate_over_cell (cell_and_error,
                           scratch_data.primal_solution,
                           scratch_data.dual_weights,
                           scratch_data.cell_data);

      // After computing the cell terms, turn to the face terms. For this,
      // loop over all faces of the present cell, and see whether
      // something needs to be computed on it:
      for (unsigned int face_no=0;
           face_no<GeometryInfo<dim>::faces_per_cell;
           ++face_no)
        {
          // First, if this face is part of the boundary, then there is
          // nothing to do. However, to make things easier when summing up
          // the contributions of the faces of cells, we enter this face
          // into the list of faces with a zero contribution to the error.
          if (cell->face(face_no)->at_boundary())
            {
              face_integrals[cell->face(face_no)] = 0;
              continue;
            }

          // Next, note that since we want to compute the jump terms on
          // each face only once although we access it twice (if it is not
          // at the boundary), we have to define some rules who is
          // responsible for computing on a face:
          //
          // First, if the neighboring cell is on the same level as this
          // one, i.e. neither further refined not coarser, then the one
          // with the lower index within this level does the work. In
          // other words: if the other one has a lower index, then skip
          // work on this face:
          if ((cell->neighbor(face_no)->has_children() == false) &&
              (cell->neighbor(face_no)->level() == cell->level()) &&
              (cell->neighbor(face_no)->index() < cell->index()))
            continue;

          // Likewise, we always work from the coarser cell if this and
          // its neighbor differ in refinement. Thus, if the neighboring
          // cell is less refined than the present one, then do nothing
          // since we integrate over the subfaces when we visit the coarse
          // cell.
          if (cell->at_boundary(face_no) == false)
            if (cell->neighbor(face_no)->level() < cell->level())
              continue;


          // Now we know that we are in charge here, so actually compute
          // the face jump terms. If the face is a regular one, i.e.  the
          // other side's cell is neither coarser not finer than this
          // cell, then call one function, and if the cell on the other
          // side is further refined, then use another function. Note that
          // the case that the cell on the other side is coarser cannot
          // happen since we have decided above that we handle this case
          // when we pass over that other cell.
          if (cell->face(face_no)->has_children() == false)
            integrate_over_regular_face (cell, face_no,
                                         scratch_data.primal_solution,
                                         scratch_data.dual_weights,
                                         scratch_data.face_data,
                                         face_integrals);
          else
            integrate_over_irregular_face (cell, face_no,
                                           scratch_data.primal_solution,
                                           scratch_data.dual_weights,
                                           scratch_data.face_data,
                                           face_integrals);
        }
    }


    // @sect4{Computing cell term error contributions}

    // As for the actual computation of the error contributions, first turn to
    // the cell terms:
    template <int dim>
    void WeightedResidual<dim>::
    integrate_over_cell (const SynchronousIterators<std_cxx11::tuple<
                         active_cell_iterator,Vector<float>::iterator> >   &cell_and_error,
                         const Vector<double>                              &primal_solution,
                         const Vector<double>                              &dual_weights,
                         CellData                                          &cell_data) const
    {
      // The tasks to be done are what appears natural from looking at the
      // error estimation formula: first get the right hand side and Laplacian
      // of the numerical solution at the quadrature points for the cell
      // residual,
      cell_data.fe_values.reinit (std_cxx11::get<0>(cell_and_error.iterators));
      cell_data.right_hand_side
      ->value_list (cell_data.fe_values.get_quadrature_points(),
                    cell_data.rhs_values);
      cell_data.fe_values.get_function_laplacians (primal_solution,
                                                   cell_data.cell_laplacians);

      // ...then get the dual weights...
      cell_data.fe_values.get_function_values (dual_weights,
                                               cell_data.dual_weights);

      // ...and finally build the sum over all quadrature points and store it
      // with the present cell:
      double sum = 0;
      for (unsigned int p=0; p<cell_data.fe_values.n_quadrature_points; ++p)
        sum += ((cell_data.rhs_values[p]+cell_data.cell_laplacians[p]) *
                cell_data.dual_weights[p] *
                cell_data.fe_values.JxW (p));
      *(std_cxx11::get<1>(cell_and_error.iterators)) += sum;
    }


    // @sect4{Computing edge term error contributions -- 1}

    // On the other hand, computation of the edge terms for the error estimate
    // is not so simple. First, we have to distinguish between faces with and
    // without hanging nodes. Because it is the simple case, we first consider
    // the case without hanging nodes on a face (let's call this the `regular'
    // case):
    template <int dim>
    void WeightedResidual<dim>::
    integrate_over_regular_face (const active_cell_iterator &cell,
                                 const unsigned int          face_no,
                                 const Vector<double>       &primal_solution,
                                 const Vector<double>       &dual_weights,
                                 FaceData                   &face_data,
                                 FaceIntegrals              &face_integrals) const
    {
      const unsigned int
      n_q_points = face_data.fe_face_values_cell.n_quadrature_points;

      // The first step is to get the values of the gradients at the
      // quadrature points of the finite element field on the present
      // cell. For this, initialize the <code>FEFaceValues</code> object
      // corresponding to this side of the face, and extract the gradients
      // using that object.
      face_data.fe_face_values_cell.reinit (cell, face_no);
      face_data.fe_face_values_cell.get_function_grads (primal_solution,
                                                        face_data.cell_grads);

      // The second step is then to extract the gradients of the finite
      // element solution at the quadrature points on the other side of the
      // face, i.e. from the neighboring cell.
      //
      // For this, do a sanity check before: make sure that the neighbor
      // actually exists (yes, we should not have come here if the neighbor
      // did not exist, but in complicated software there are bugs, so better
      // check this), and if this is not the case throw an error.
      Assert (cell->neighbor(face_no).state() == IteratorState::valid,
              ExcInternalError());
      // If we have that, then we need to find out with which face of the
      // neighboring cell we have to work, i.e. the <code>how-many'th</code> the
      // neighbor the present cell is of the cell behind the present face. For
      // this, there is a function, and we put the result into a variable with
      // the name <code>neighbor_neighbor</code>:
      const unsigned int
      neighbor_neighbor = cell->neighbor_of_neighbor (face_no);
      // Then define an abbreviation for the neighbor cell, initialize the
      // <code>FEFaceValues</code> object on that cell, and extract the
      // gradients on that cell:
      const active_cell_iterator neighbor = cell->neighbor(face_no);
      face_data.fe_face_values_neighbor.reinit (neighbor, neighbor_neighbor);
      face_data.fe_face_values_neighbor.get_function_grads (primal_solution,
                                                            face_data.neighbor_grads);

      // Now that we have the gradients on this and the neighboring cell,
      // compute the jump residual by multiplying the jump in the gradient
      // with the normal vector:
      for (unsigned int p=0; p<n_q_points; ++p)
        face_data.jump_residual[p]
          = ((face_data.cell_grads[p] - face_data.neighbor_grads[p]) *
             face_data.fe_face_values_cell.normal_vector(p));

      // Next get the dual weights for this face:
      face_data.fe_face_values_cell.get_function_values (dual_weights,
                                                         face_data.dual_weights);

      // Finally, we have to compute the sum over jump residuals, dual
      // weights, and quadrature weights, to get the result for this face:
      double face_integral = 0;
      for (unsigned int p=0; p<n_q_points; ++p)
        face_integral += (face_data.jump_residual[p] *
                          face_data.dual_weights[p]  *
                          face_data.fe_face_values_cell.JxW(p));

      // Double check that the element already exists and that it was not
      // already written to...
      Assert (face_integrals.find (cell->face(face_no)) != face_integrals.end(),
              ExcInternalError());
      Assert (face_integrals[cell->face(face_no)] == -1e20,
              ExcInternalError());

      // ...then store computed value at assigned location. Note that the
      // stored value does not contain the factor 1/2 that appears in the
      // error representation. The reason is that the term actually does not
      // have this factor if we loop over all faces in the triangulation, but
      // only appears if we write it as a sum over all cells and all faces of
      // each cell; we thus visit the same face twice. We take account of this
      // by using this factor -1/2 later, when we sum up the contributions for
      // each cell individually.
      face_integrals[cell->face(face_no)] = face_integral;
    }


    // @sect4{Computing edge term error contributions -- 2}

    // We are still missing the case of faces with hanging nodes. This is what
    // is covered in this function:
    template <int dim>
    void WeightedResidual<dim>::
    integrate_over_irregular_face (const active_cell_iterator &cell,
                                   const unsigned int          face_no,
                                   const Vector<double>       &primal_solution,
                                   const Vector<double>       &dual_weights,
                                   FaceData                   &face_data,
                                   FaceIntegrals              &face_integrals) const
    {
      // First again two abbreviations, and some consistency checks whether
      // the function is called only on faces for which it is supposed to be
      // called:
      const unsigned int
      n_q_points = face_data.fe_face_values_cell.n_quadrature_points;

      const typename DoFHandler<dim>::face_iterator
      face = cell->face(face_no);
      const typename DoFHandler<dim>::cell_iterator
      neighbor = cell->neighbor(face_no);
      Assert (neighbor.state() == IteratorState::valid,
              ExcInternalError());
      Assert (neighbor->has_children(),
              ExcInternalError());

      // Then find out which neighbor the present cell is of the adjacent
      // cell. Note that we will operate on the children of this adjacent
      // cell, but that their orientation is the same as that of their mother,
      // i.e. the neighbor direction is the same.
      const unsigned int
      neighbor_neighbor = cell->neighbor_of_neighbor (face_no);

      // Then simply do everything we did in the previous function for one
      // face for all the sub-faces now:
      for (unsigned int subface_no=0;
           subface_no<face->n_children(); ++subface_no)
        {
          // Start with some checks again: get an iterator pointing to the
          // cell behind the present subface and check whether its face is a
          // subface of the one we are considering. If that were not the case,
          // then there would be either a bug in the
          // <code>neighbor_neighbor</code> function called above, or -- worse
          // -- some function in the library did not keep to some underlying
          // assumptions about cells, their children, and their faces. In any
          // case, even though this assertion should not be triggered, it does
          // not harm to be cautious, and in optimized mode computations the
          // assertion will be removed anyway.
          const active_cell_iterator neighbor_child
            = cell->neighbor_child_on_subface (face_no, subface_no);
          Assert (neighbor_child->face(neighbor_neighbor) ==
                  cell->face(face_no)->child(subface_no),
                  ExcInternalError());

          // Now start the work by again getting the gradient of the solution
          // first at this side of the interface,
          face_data.fe_subface_values_cell.reinit (cell, face_no, subface_no);
          face_data.fe_subface_values_cell.get_function_grads (primal_solution,
                                                               face_data.cell_grads);
          // then at the other side,
          face_data.fe_face_values_neighbor.reinit (neighbor_child,
                                                    neighbor_neighbor);
          face_data.fe_face_values_neighbor.get_function_grads (primal_solution,
                                                                face_data.neighbor_grads);

          // and finally building the jump residuals. Since we take the normal
          // vector from the other cell this time, revert the sign of the
          // first term compared to the other function:
          for (unsigned int p=0; p<n_q_points; ++p)
            face_data.jump_residual[p]
              = ((face_data.neighbor_grads[p] - face_data.cell_grads[p]) *
                 face_data.fe_face_values_neighbor.normal_vector(p));

          // Then get dual weights:
          face_data.fe_face_values_neighbor.get_function_values (dual_weights,
                                                                 face_data.dual_weights);

          // At last, sum up the contribution of this sub-face, and set it in
          // the global map:
          double face_integral = 0;
          for (unsigned int p=0; p<n_q_points; ++p)
            face_integral += (face_data.jump_residual[p] *
                              face_data.dual_weights[p] *
                              face_data.fe_face_values_neighbor.JxW(p));
          face_integrals[neighbor_child->face(neighbor_neighbor)]
            = face_integral;
        }

      // Once the contributions of all sub-faces are computed, loop over all
      // sub-faces to collect and store them with the mother face for simple
      // use when later collecting the error terms of cells. Again make safety
      // checks that the entries for the sub-faces have been computed and do
      // not carry an invalid value.
      double sum = 0;
      for (unsigned int subface_no=0;
           subface_no<face->n_children(); ++subface_no)
        {
          Assert (face_integrals.find(face->child(subface_no)) !=
                  face_integrals.end(),
                  ExcInternalError());
          Assert (face_integrals[face->child(subface_no)] != -1e20,
                  ExcInternalError());

          sum += face_integrals[face->child(subface_no)];
        }
      // Finally store the value with the parent face.
      face_integrals[face] = sum;
    }

  }


  // @sect3{A simulation framework}

  // In the previous example program, we have had two functions that were used
  // to drive the process of solving on subsequently finer grids. We extend
  // this here to allow for a number of parameters to be passed to these
  // functions, and put all of that into framework class.
  //
  // You will have noted that this program is built up of a number of small
  // parts (evaluation functions, solver classes implementing various
  // refinement methods, different dual functionals, different problem and
  // data descriptions), which makes the program relatively simple to extend,
  // but also allows to solve a large number of different problems by
  // replacing one part by another. We reflect this flexibility by declaring a
  // structure in the following framework class that holds a number of
  // parameters that may be set to test various combinations of the parts of
  // this program, and which can be used to test it at various problems and
  // discretizations in a simple way.
  template <int dim>
  struct Framework
  {
  public:
    // First, we declare two abbreviations for simple use of the respective
    // data types:
    typedef Evaluation::EvaluationBase<dim> Evaluator;
    typedef std::list<Evaluator *>           EvaluatorList;


    // Then we have the structure which declares all the parameters that may
    // be set. In the default constructor of the structure, these values are
    // all set to default values, for simple use.
    struct ProblemDescription
    {
      // First allow for the degrees of the piecewise polynomials by which the
      // primal and dual problems will be discretized. They default to (bi-,
      // tri-)linear ansatz functions for the primal, and (bi-, tri-)quadratic
      // ones for the dual problem. If a refinement criterion is chosen that
      // does not need the solution of a dual problem, the value of the dual
      // finite element degree is of course ignored.
      unsigned int primal_fe_degree;
      unsigned int dual_fe_degree;

      // Then have an object that describes the problem type, i.e. right hand
      // side, domain, boundary values, etc. The pointer needed here defaults
      // to the Null pointer, i.e. you will have to set it in actual instances
      // of this object to make it useful.
      SmartPointer<const Data::SetUpBase<dim> > data;

      // Since we allow to use different refinement criteria (global
      // refinement, refinement by the Kelly error indicator, possibly with a
      // weight, and using the dual estimator), define a number of enumeration
      // values, and subsequently a variable of that type. It will default to
      // <code>dual_weighted_error_estimator</code>.
      enum RefinementCriterion
      {
        dual_weighted_error_estimator,
        global_refinement,
        kelly_indicator,
        weighted_kelly_indicator
      };

      RefinementCriterion refinement_criterion;

      // Next, an object that describes the dual functional. It is only needed
      // if the dual weighted residual refinement is chosen, and also defaults
      // to a Null pointer.
      SmartPointer<const DualFunctional::DualFunctionalBase<dim> > dual_functional;

      // Then a list of evaluation objects. Its default value is empty,
      // i.e. no evaluation objects.
      EvaluatorList evaluator_list;

      // Next to last, a function that is used as a weight to the
      // <code>RefinementWeightedKelly</code> class. The default value of this
      // pointer is zero, but you have to set it to some other value if you
      // want to use the <code>weighted_kelly_indicator</code> refinement
      // criterion.
      SmartPointer<const Function<dim> > kelly_weight;

      // Finally, we have a variable that denotes the maximum number of
      // degrees of freedom we allow for the (primal) discretization. If it is
      // exceeded, we stop the process of solving and intermittent mesh
      // refinement. Its default value is 20,000.
      unsigned int max_degrees_of_freedom;

      // Finally the default constructor of this class:
      ProblemDescription ();
    };

    // The driver framework class only has one method which calls solver and
    // mesh refinement intermittently, and does some other small tasks in
    // between. Since it does not need data besides the parameters given to
    // it, we make it static:
    static void run (const ProblemDescription &descriptor);
  };


  // As for the implementation, first the constructor of the parameter object,
  // setting all values to their defaults:
  template <int dim>
  Framework<dim>::ProblemDescription::ProblemDescription ()
    :
    primal_fe_degree (1),
    dual_fe_degree (2),
    refinement_criterion (dual_weighted_error_estimator),
    max_degrees_of_freedom (20000)
  {}



  // Then the function which drives the whole process:
  template <int dim>
  void Framework<dim>::run (const ProblemDescription &descriptor)
  {
    // First create a triangulation from the given data object,
    Triangulation<dim>
    triangulation (Triangulation<dim>::smoothing_on_refinement);
    descriptor.data->create_coarse_grid (triangulation);

    // then a set of finite elements and appropriate quadrature formula:
    const FE_Q<dim>     primal_fe(descriptor.primal_fe_degree);
    const FE_Q<dim>     dual_fe(descriptor.dual_fe_degree);
    const QGauss<dim>   quadrature(descriptor.dual_fe_degree+1);
    const QGauss<dim-1> face_quadrature(descriptor.dual_fe_degree+1);

    // Next, select one of the classes implementing different refinement
    // criteria.
    LaplaceSolver::Base<dim> *solver = 0;
    switch (descriptor.refinement_criterion)
      {
      case ProblemDescription::dual_weighted_error_estimator:
      {
        solver
          = new LaplaceSolver::WeightedResidual<dim> (triangulation,
                                                      primal_fe,
                                                      dual_fe,
                                                      quadrature,
                                                      face_quadrature,
                                                      descriptor.data->get_right_hand_side(),
                                                      descriptor.data->get_boundary_values(),
                                                      *descriptor.dual_functional);
        break;
      }

      case ProblemDescription::global_refinement:
      {
        solver
          = new LaplaceSolver::RefinementGlobal<dim> (triangulation,
                                                      primal_fe,
                                                      quadrature,
                                                      face_quadrature,
                                                      descriptor.data->get_right_hand_side(),
                                                      descriptor.data->get_boundary_values());
        break;
      }

      case ProblemDescription::kelly_indicator:
      {
        solver
          = new LaplaceSolver::RefinementKelly<dim> (triangulation,
                                                     primal_fe,
                                                     quadrature,
                                                     face_quadrature,
                                                     descriptor.data->get_right_hand_side(),
                                                     descriptor.data->get_boundary_values());
        break;
      }

      case ProblemDescription::weighted_kelly_indicator:
      {
        solver
          = new LaplaceSolver::RefinementWeightedKelly<dim> (triangulation,
                                                             primal_fe,
                                                             quadrature,
                                                             face_quadrature,
                                                             descriptor.data->get_right_hand_side(),
                                                             descriptor.data->get_boundary_values(),
                                                             *descriptor.kelly_weight);
        break;
      }

      default:
        AssertThrow (false, ExcInternalError());
      }

    // Now that all objects are in place, run the main loop. The stopping
    // criterion is implemented at the bottom of the loop.
    //
    // In the loop, first set the new cycle number, then solve the problem,
    // output its solution(s), apply the evaluation objects to it, then decide
    // whether we want to refine the mesh further and solve again on this
    // mesh, or jump out of the loop.
    for (unsigned int step=0; true; ++step)
      {
        std::cout << "Refinement cycle: "       << step
                  << std::endl;

        solver->set_refinement_cycle (step);
        solver->solve_problem ();
        solver->output_solution ();

        std::cout << "   Number of degrees of freedom="
                  << solver->n_dofs() << std::endl;

        for (typename EvaluatorList::const_iterator
             e = descriptor.evaluator_list.begin();
             e != descriptor.evaluator_list.end(); ++e)
          {
            (*e)->set_refinement_cycle (step);
            solver->postprocess (**e);
          }


        if (solver->n_dofs() < descriptor.max_degrees_of_freedom)
          solver->refine_grid ();
        else
          break;
      }

    // After the loop has run, clean up the screen, and delete objects no more
    // needed:
    std::cout << std::endl;
    delete solver;
    solver = 0;
  }

}



// @sect3{The main function}

// Here finally comes the main function. It drives the whole process by
// specifying a set of parameters to be used for the simulation (polynomial
// degrees, evaluation and dual functionals, etc), and passes them packed into
// a structure to the frame work class above.
int main ()
{
  try
    {
      using namespace dealii;
      using namespace Step14;

      deallog.depth_console (0);
      // Describe the problem we want to solve here by passing a descriptor
      // object to the function doing the rest of the work:
      const unsigned int dim = 2;
      Framework<dim>::ProblemDescription descriptor;

      // First set the refinement criterion we wish to use:
      descriptor.refinement_criterion
        = Framework<dim>::ProblemDescription::dual_weighted_error_estimator;
      // Here, we could as well have used <code>global_refinement</code> or
      // <code>weighted_kelly_indicator</code>. Note that the information
      // given about dual finite elements, dual functional, etc is only
      // important for the given choice of refinement criterion, and is
      // ignored otherwise.

      // Then set the polynomial degrees of primal and dual problem. We choose
      // here bi-linear and bi-quadratic ones:
      descriptor.primal_fe_degree = 1;
      descriptor.dual_fe_degree   = 2;

      // Then set the description of the test case, i.e. domain, boundary
      // values, and right hand side. These are prepackaged in classes. We
      // take here the description of <code>Exercise_2_3</code>, but you can
      // also use <code>CurvedRidges@<dim@></code>:
      descriptor.data = new Data::SetUp<Data::Exercise_2_3<dim>,dim> ();

      // Next set first a dual functional, then a list of evaluation
      // objects. We choose as default the evaluation of the value at an
      // evaluation point, represented by the classes
      // <code>PointValueEvaluation</code> in the namespaces of evaluation and
      // dual functional classes. You can also set the
      // <code>PointXDerivativeEvaluation</code> classes for the x-derivative
      // instead of the value at the evaluation point.
      //
      // Note that dual functional and evaluation objects should
      // match. However, you can give as many evaluation functionals as you
      // want, so you can have both point value and derivative evaluated after
      // each step.  One such additional evaluation is to output the grid in
      // each step.
      const Point<dim> evaluation_point (0.75, 0.75);
      descriptor.dual_functional
        = new DualFunctional::PointValueEvaluation<dim> (evaluation_point);

      Evaluation::PointValueEvaluation<dim>
      postprocessor1 (evaluation_point);
      Evaluation::GridOutput<dim>
      postprocessor2 ("grid");

      descriptor.evaluator_list.push_back (&postprocessor1);
      descriptor.evaluator_list.push_back (&postprocessor2);

      // Set the maximal number of degrees of freedom after which we want the
      // program to stop refining the mesh further:
      descriptor.max_degrees_of_freedom = 20000;

      // Finally pass the descriptor object to a function that runs the entire
      // solution with it:
      Framework<dim>::run (descriptor);
    }

  // Catch exceptions to give information about things that failed:
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

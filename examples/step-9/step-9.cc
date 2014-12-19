/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2000 - 2014 by the deal.II authors
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
 * Author: Wolfgang Bangerth, University of Heidelberg, 2000
 */


// Just as in previous examples, we have to include several files of which the
// meaning has already been discussed:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_out.h>

// The following two files provide classes and information for multithreaded
// programs. In the first one, the classes and functions are declared which we
// need to do assembly in parallel (i.e. the
// <code>WorkStream</code> namespace). The
// second file has a class <code>MultithreadInfo</code> (and a global object
// <code>multithread_info</code> of that type) which can be used to query the
// number of processors in your system, which is often useful when deciding
// how many threads to start in parallel.
#include <deal.II/base/work_stream.h>
#include <deal.II/base/multithread_info.h>

// The next new include file declares a base class <code>TensorFunction</code>
// not unlike the <code>Function</code> class, but with the difference that
// the return value is tensor-valued rather than scalar of vector-valued.
#include <deal.II/base/tensor_function.h>

#include <deal.II/numerics/error_estimator.h>

// This is C++, as we want to write some output to disk:
#include <fstream>
#include <iostream>


// The last step is as in previous programs:
namespace Step9
{
  using namespace dealii;

  // @sect3{AdvectionProblem class declaration}

  // Following we declare the main class of this program. It is very much
  // like the main classes of previous examples, so we again only comment on
  // the differences.
  template <int dim>
  class AdvectionProblem
  {
  public:
    AdvectionProblem ();
    ~AdvectionProblem ();
    void run ();

  private:
    void setup_system ();

    // The next set of functions will be used to assemble the
    // matrix. However, unlike in the previous examples, the
    // <code>assemble_system()</code> function will not do the work
    // itself, but rather will delegate the actual assembly to helper
    // functions <code>assemble_local_system()</code> and
    // <code>copy_local_to_global()</code>. The rationale is that
    // matrix assembly can be parallelized quite well, as the
    // computation of the local contributions on each cell is entirely
    // independent of other cells, and we only have to synchronize
    // when we add the contribution of a cell to the global
    // matrix.
    //
    // The strategy for parallelization we choose here is one of the
    // possibilities mentioned in detail in the @ref threads module in
    // the documentation. Specifically, we will use the WorkStream
    // approach discussed there. Since there is so much documentation
    // in this module, we will not repeat the rationale for the design
    // choices here (for example, if you read through the module
    // mentioned above, you will understand what the purpose of the
    // <code>AssemblyScratchData</code> and
    // <code>AssemblyCopyData</code> structures is). Rather, we will
    // only discuss the specific implementation.
    //
    // If you read the page mentioned above, you will find that in
    // order to parallelize assembly, we need two data structures --
    // one that corresponds to data that we need during local
    // integration ("scratch data", i.e., things we only need as
    // temporary storage), and one that carries information from the
    // local integration to the function that then adds the local
    // contributions to the corresponding elements of the global
    // matrix. The former of these typically contains the FEValues and
    // FEFaceValues objects, whereas the latter has the local matrix,
    // local right hand side, and information about which degrees of
    // freedom live on the cell for which we are assembling a local
    // contribution. With this information, the following should be
    // relatively self-explanatory:
    struct AssemblyScratchData
    {
      AssemblyScratchData (const FiniteElement<dim> &fe);
      AssemblyScratchData (const AssemblyScratchData &scratch_data);

      FEValues<dim>     fe_values;
      FEFaceValues<dim> fe_face_values;
    };

    struct AssemblyCopyData
    {
      FullMatrix<double>                   cell_matrix;
      Vector<double>                       cell_rhs;
      std::vector<types::global_dof_index> local_dof_indices;
    };

    void assemble_system ();
    void local_assemble_system (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                AssemblyScratchData                                  &scratch,
                                AssemblyCopyData                                     &copy_data);
    void copy_local_to_global (const AssemblyCopyData &copy_data);


    // The following functions again are as in previous examples, as are the
    // subsequent variables.
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;
    DoFHandler<dim>      dof_handler;

    FE_Q<dim>            fe;

    ConstraintMatrix     hanging_node_constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;
  };



  // @sect3{Equation data declaration}

  // Next we declare a class that describes the advection field. This, of
  // course, is a vector field with as many components as there are space
  // dimensions. One could now use a class derived from the
  // <code>Function</code> base class, as we have done for boundary values and
  // coefficients in previous examples, but there is another possibility in
  // the library, namely a base class that describes tensor valued
  // functions. In contrast to the usual <code>Function</code> objects, we
  // provide the compiler with knowledge on the size of the objects of the
  // return type. This enables the compiler to generate efficient code, which
  // is not so simple for usual vector-valued functions where memory has to be
  // allocated on the heap (thus, the <code>Function::vector_value</code>
  // function has to be given the address of an object into which the result
  // is to be written, in order to avoid copying and memory allocation and
  // deallocation on the heap). In addition to the known size, it is possible
  // not only to return vectors, but also tensors of higher rank; however,
  // this is not very often requested by applications, to be honest...
  //
  // The interface of the <code>TensorFunction</code> class is relatively
  // close to that of the <code>Function</code> class, so there is probably no
  // need to comment in detail the following declaration:
  template <int dim>
  class AdvectionField : public TensorFunction<1,dim>
  {
  public:
    AdvectionField () : TensorFunction<1,dim> () {}

    virtual Tensor<1,dim> value (const Point<dim> &p) const;

    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<Tensor<1,dim> >    &values) const;

    // In previous examples, we have used assertions that throw exceptions in
    // several places. However, we have never seen how such exceptions are
    // declared. This can be done as follows:
    DeclException2 (ExcDimensionMismatch,
                    unsigned int, unsigned int,
                    << "The vector has size " << arg1 << " but should have "
                    << arg2 << " elements.");
    // The syntax may look a little strange, but is reasonable. The format is
    // basically as follows: use the name of one of the macros
    // <code>DeclExceptionN</code>, where <code>N</code> denotes the number of
    // additional parameters which the exception object shall take. In this
    // case, as we want to throw the exception when the sizes of two vectors
    // differ, we need two arguments, so we use
    // <code>DeclException2</code>. The first parameter then describes the
    // name of the exception, while the following declare the data types of
    // the parameters. The last argument is a sequence of output directives
    // that will be piped into the <code>std::cerr</code> object, thus the
    // strange format with the leading <code>@<@<</code> operator and the
    // like. Note that we can access the parameters which are passed to the
    // exception upon construction (i.e. within the <code>Assert</code> call)
    // by using the names <code>arg1</code> through <code>argN</code>, where
    // <code>N</code> is the number of arguments as defined by the use of the
    // respective macro <code>DeclExceptionN</code>.
    //
    // To learn how the preprocessor expands this macro into actual code,
    // please refer to the documentation of the exception classes in the base
    // library. Suffice it to say that by this macro call, the respective
    // exception class is declared, which also has error output functions
    // already implemented.
  };



  // The following two functions implement the interface described above. The
  // first simply implements the function as described in the introduction,
  // while the second uses the same trick to avoid calling a virtual function
  // as has already been introduced in the previous example program. Note the
  // check for the right sizes of the arguments in the second function, which
  // should always be present in such functions; it is our experience that
  // many if not most programming errors result from incorrectly initialized
  // arrays, incompatible parameters to functions and the like; using
  // assertion as in this case can eliminate many of these problems.
  template <int dim>
  Tensor<1,dim>
  AdvectionField<dim>::value (const Point<dim> &p) const
  {
    Point<dim> value;
    value[0] = 2;
    for (unsigned int i=1; i<dim; ++i)
      value[i] = 1+0.8*std::sin(8*numbers::PI*p[0]);

    return value;
  }



  template <int dim>
  void
  AdvectionField<dim>::value_list (const std::vector<Point<dim> > &points,
                                   std::vector<Tensor<1,dim> >    &values) const
  {
    Assert (values.size() == points.size(),
            ExcDimensionMismatch (values.size(), points.size()));

    for (unsigned int i=0; i<points.size(); ++i)
      values[i] = AdvectionField<dim>::value (points[i]);
  }




  // Besides the advection field, we need two functions describing the source
  // terms (<code>right hand side</code>) and the boundary values. First for
  // the right hand side, which follows the same pattern as in previous
  // examples. As described in the introduction, the source is a constant
  // function in the vicinity of a source point, which we denote by the
  // constant static variable <code>center_point</code>. We set the values of
  // this center using the same template tricks as we have shown in the step-7
  // example program. The rest is simple and has been shown previously,
  // including the way to avoid virtual function calls in the
  // <code>value_list</code> function.
  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const;

  private:
    static const Point<dim> center_point;
  };


  template <>
  const Point<1> RightHandSide<1>::center_point = Point<1> (-0.75);

  template <>
  const Point<2> RightHandSide<2>::center_point = Point<2> (-0.75, -0.75);

  template <>
  const Point<3> RightHandSide<3>::center_point = Point<3> (-0.75, -0.75, -0.75);



  // The only new thing here is that we check for the value of the
  // <code>component</code> parameter. As this is a scalar function, it is
  // obvious that it only makes sense if the desired component has the index
  // zero, so we assert that this is indeed the
  // case. <code>ExcIndexRange</code> is a global predefined exception
  // (probably the one most often used, we therefore made it global instead of
  // local to some class), that takes three parameters: the index that is
  // outside the allowed range, the first element of the valid range and the
  // one past the last (i.e. again the half-open interval so often used in the
  // C++ standard library):
  template <int dim>
  double
  RightHandSide<dim>::value (const Point<dim>   &p,
                             const unsigned int  component) const
  {
    Assert (component == 0, ExcIndexRange (component, 0, 1));
    const double diameter = 0.1;
    return ( (p-center_point).square() < diameter*diameter ?
             .1/std::pow(diameter,dim) :
             0);
  }



  template <int dim>
  void
  RightHandSide<dim>::value_list (const std::vector<Point<dim> > &points,
                                  std::vector<double>            &values,
                                  const unsigned int              component) const
  {
    Assert (values.size() == points.size(),
            ExcDimensionMismatch (values.size(), points.size()));

    for (unsigned int i=0; i<points.size(); ++i)
      values[i] = RightHandSide<dim>::value (points[i], component);
  }



  // Finally for the boundary values, which is just another class derived from
  // the <code>Function</code> base class:
  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    BoundaryValues () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const;
  };



  template <int dim>
  double
  BoundaryValues<dim>::value (const Point<dim>   &p,
                              const unsigned int  component) const
  {
    Assert (component == 0, ExcIndexRange (component, 0, 1));

    const double sine_term = std::sin(16*numbers::PI*std::sqrt(p.square()));
    const double weight    = std::exp(-5*p.square()) / std::exp(-5.);
    return sine_term * weight;
  }



  template <int dim>
  void
  BoundaryValues<dim>::value_list (const std::vector<Point<dim> > &points,
                                   std::vector<double>            &values,
                                   const unsigned int              component) const
  {
    Assert (values.size() == points.size(),
            ExcDimensionMismatch (values.size(), points.size()));

    for (unsigned int i=0; i<points.size(); ++i)
      values[i] = BoundaryValues<dim>::value (points[i], component);
  }



  // @sect3{GradientEstimation class declaration}

  // Now, finally, here comes the class that will compute the difference
  // approximation of the gradient on each cell and weighs that with a power
  // of the mesh size, as described in the introduction.  This class is a
  // simple version of the <code>DerivativeApproximation</code> class in the
  // library, that uses similar techniques to obtain finite difference
  // approximations of the gradient of a finite element field, or of higher
  // derivatives.
  //
  // The class has one public static function <code>estimate</code> that is
  // called to compute a vector of error indicators, and a few private functions
  // that do the actual work on all active cells. As in other parts of the
  // library, we follow an informal convention to use vectors of floats for
  // error indicators rather than the common vectors of doubles, as the
  // additional accuracy is not necessary for estimated values.
  //
  // In addition to these two functions, the class declares two exceptions
  // which are raised when a cell has no neighbors in each of the space
  // directions (in which case the matrix described in the introduction would
  // be singular and can't be inverted), while the other one is used in the
  // more common case of invalid parameters to a function, namely a vector of
  // wrong size.
  //
  // Two other comments: first, the class has no non-static member functions
  // or variables, so this is not really a class, but rather serves the
  // purpose of a <code>namespace</code> in C++. The reason that we chose a
  // class over a namespace is that this way we can declare functions that are
  // private. This can be done with namespaces as well, if one declares some
  // functions in header files in the namespace and implements these and other
  // functions in the implementation file. The functions not declared in the
  // header file are still in the namespace but are not callable from
  // outside. However, as we have only one file here, it is not possible to
  // hide functions in the present case.
  //
  // The second comment is that the dimension template parameter is attached
  // to the function rather than to the class itself. This way, you don't have
  // to specify the template parameter yourself as in most other cases, but
  // the compiler can figure its value out itself from the dimension of the
  // DoF handler object that one passes as first argument.
  //
  // Before jumping into the fray with the implementation, let us also comment
  // on the parallelization strategy. We have already introduced the necessary
  // framework for using the WorkStream concept in the declaration of the main
  // class of this program above. We will use it again here. In the current
  // context, this means that we have to define (i) classes for scratch and
  // copy objects, (ii) a function that does the local computation on one
  // cell, and (iii) a function that copies the local result into a global
  // object. Given this general framework, we will, however, deviate from it a
  // bit. In particular, WorkStream was generally invented for cases where
  // each local computation on a cell <i>adds</i> to a global object -- for
  // example, when assembling linear systems where we add local contributions
  // into a global matrix and right hand side. Here, however, the situation is
  // slightly different: we compute contributions from every cell
  // individually, but then all we need to do is put them into an element of
  // an output vector that is unique to each cell. Consequently, there is no
  // risk that the write operations from two cells might conflict, and the
  // elaborate machinery of WorkStream to avoid conflicting writes is not
  // necessary. Consequently, what we will do is this: We still need a scratch
  // object that holds, for example, the FEValues object.  However, we only
  // create an fake, empty copy data structure. Likewise, we do need the
  // function that computes local contributions, but since it can already put
  // the result into its final location, we do not need a copy-local-to-global
  // function and will instead give the WorkStream::run function an empty
  // function object -- the equivalent to a NULL function pointer.
  //
  // The second idea to make this approach work is this: If we want to write
  // the result into its final destination right away, then the local worker
  // function needs to already know where this destination is. Here, this is
  // an element of a vector -- but which element is something that the local
  // worker function (or, if we wanted to use one, a copy-local-to-global
  // function) can not determine easily just knowing an iterator to a cell it
  // is supposed to work on. Consequently, in addition to a cell, we need to
  // pass a second piece of identifying information along: the element of the
  // output vector to write into. What this means is that the work items are
  // identified by two iterators: to a cell, and to an output vector
  // element. Moving from one work item to the next requires incrementing both
  // iterators. deal.II has a class for this, called SynchronousIterators,
  // that takes a tuple of iterator types as arguments and stores an iterator
  // of each type. Whenever the SynchronousIterators object is incremented, it
  // increments the stored iterators in turn. Thus, this class is exactly what
  // we need to do our work, and we consequently use it as the first argument
  // of the worker function. We will further down below show how to create
  // such an object.
  class GradientEstimation
  {
  public:
    template <int dim>
    static void estimate (const DoFHandler<dim> &dof,
                          const Vector<double>  &solution,
                          Vector<float>         &error_per_cell);

    DeclException2 (ExcInvalidVectorLength,
                    int, int,
                    << "Vector has length " << arg1 << ", but should have "
                    << arg2);
    DeclException0 (ExcInsufficientDirections);

  private:
    template <int dim>
    struct EstimateScratchData
    {
      EstimateScratchData (const FiniteElement<dim> &fe,
                           const Vector<double>     &solution);
      EstimateScratchData (const EstimateScratchData &data);

      FEValues<dim> fe_midpoint_value;
      Vector<double> solution;
    };

    struct EstimateCopyData
    {};

    template <int dim>
    static
    void estimate_cell (const SynchronousIterators<std_cxx11::tuple<typename DoFHandler<dim>::active_cell_iterator,
                        Vector<float>::iterator> >     &cell,
                        EstimateScratchData<dim>       &scratch_data,
                        const EstimateCopyData         &copy_data);
  };



  // @sect3{AdvectionProblem class implementation}


  // Now for the implementation of the main class. Constructor, destructor and
  // the function <code>setup_system</code> follow the same pattern that was
  // used previously, so we need not comment on these three function:
  template <int dim>
  AdvectionProblem<dim>::AdvectionProblem ()
    :
    dof_handler (triangulation),
    fe(1)
  {}



  template <int dim>
  AdvectionProblem<dim>::~AdvectionProblem ()
  {
    dof_handler.clear ();
  }



  template <int dim>
  void AdvectionProblem<dim>::setup_system ()
  {
    dof_handler.distribute_dofs (fe);

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

    system_matrix.reinit (sparsity_pattern);

    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());
  }



  // In the following function, the matrix and right hand side are
  // assembled. As stated in the documentation of the main class above, it
  // does not do this itself, but rather delegates to the function following
  // next, utilizing the WorkStream concept discussed in @ref threads .
  //
  // If you have looked through the @ref threads module, you will have
  // seen that assembling in parallel does not take an incredible
  // amount of extra code as long as you diligently describe what the
  // scratch and copy data objects are, and if you define suitable
  // functios for the local assembly and the copy operation from local
  // contributions to global objects. This done, the following will do
  // all the heavy lifting to get these operations done on multiple
  // threads on as many cores as you have in your system:
  template <int dim>
  void AdvectionProblem<dim>::assemble_system ()
  {
    WorkStream::run(dof_handler.begin_active(),
                    dof_handler.end(),
                    *this,
                    &AdvectionProblem::local_assemble_system,
                    &AdvectionProblem::copy_local_to_global,
                    AssemblyScratchData(fe),
                    AssemblyCopyData());


    // After the matrix has been assembled in parallel, we still have to
    // eliminate hanging node constraints. This is something that can't be
    // done on each of the threads separately, so we have to do it now.
    // Note also, that unlike in previous examples, there are no boundary
    // conditions to be applied to the system of equations. This, of course,
    // is due to the fact that we have included them into the weak formulation
    // of the problem.
    hanging_node_constraints.condense (system_matrix);
    hanging_node_constraints.condense (system_rhs);
  }



  // As already mentioned above, we need to have scratch objects for
  // the parallel computation of local contributions. These objects
  // contain FEValues and FEFaceValues objects, and so we will need to
  // have constructors and copy constructors that allow us to create
  // them. In initializing them, note first that we use bilinear
  // elements, soGauss formulae with two points in each space
  // direction are sufficient.  For the cell terms we need the values
  // and gradients of the shape functions, the quadrature points in
  // order to determine the source density and the advection field at
  // a given point, and the weights of the quadrature points times the
  // determinant of the Jacobian at these points. In contrast, for the
  // boundary integrals, we don't need the gradients, but rather the
  // normal vectors to the cells. This determines which update flags
  // we will have to pass to the constructors of the members of the
  // class:
  template <int dim>
  AdvectionProblem<dim>::AssemblyScratchData::
  AssemblyScratchData (const FiniteElement<dim> &fe)
    :
    fe_values (fe,
               QGauss<dim>(2),
               update_values   | update_gradients |
               update_quadrature_points | update_JxW_values),
    fe_face_values (fe,
                    QGauss<dim-1>(2),
                    update_values     | update_quadrature_points   |
                    update_JxW_values | update_normal_vectors)
  {}



  template <int dim>
  AdvectionProblem<dim>::AssemblyScratchData::
  AssemblyScratchData (const AssemblyScratchData &scratch_data)
    :
    fe_values (scratch_data.fe_values.get_fe(),
               scratch_data.fe_values.get_quadrature(),
               update_values   | update_gradients |
               update_quadrature_points | update_JxW_values),
    fe_face_values (scratch_data.fe_face_values.get_fe(),
                    scratch_data.fe_face_values.get_quadrature(),
                    update_values     | update_quadrature_points   |
                    update_JxW_values | update_normal_vectors)
  {}




  // Now, this is the function that does the actual work. It is not very
  // different from the <code>assemble_system</code> functions of previous
  // example programs, so we will again only comment on the differences. The
  // mathematical stuff follows closely what we have said in the introduction.
  //
  // There are a number of points worth mentioning here, though. The
  // first one is that we have moved the FEValues and FEFaceValues
  // objects into the ScratchData object. We have done so because the
  // alternative would have been to simply create one every time we
  // get into this function -- i.e., on every cell. It now turns out
  // that the FEValues classes were written with the explicit goal of
  // moving everything that remains the same from cell to cell into
  // the construction of the object, and only do as little work as
  // possible in FEValues::reinit() whenever we move to a new
  // cell. What this means is that it would be very expensive to
  // create a new object of this kind in this function as we would
  // have to do it for every cell -- exactly the thing we wanted to
  // avoid with the FEValues class. Instead, what we do is create it
  // only once (or a small number of times) in the scratch objects and
  // then re-use it as often as we can.
  //
  // This begs the question of whether there are other objects we
  // create in this function whose creation is expensive compared to
  // its use. Indeed, at the top of the function, we declare all sorts
  // of objects. The <code>AdvectionField</code>,
  // <code>RightHandSide</code> and <code>BoundaryValues</code> do not
  // cost much to create, so there is no harm here. However,
  // allocating memory in creating the <code>rhs_values</code> and
  // similar variables below typically costs a significant amount of
  // time, compared to just accessing the (temporary) values we store
  // in them. Consequently, these would be candidates for moving into
  // the <code>AssemblyScratchData</code> class. We will leave this as
  // an exercise.
  template <int dim>
  void
  AdvectionProblem<dim>::
  local_assemble_system (const typename DoFHandler<dim>::active_cell_iterator &cell,
                         AssemblyScratchData                                  &scratch_data,
                         AssemblyCopyData                                     &copy_data)
  {
    // First of all, we will need some objects that describe boundary values,
    // right hand side function and the advection field. As we will only
    // perform actions on these objects that do not change them, we declare
    // them as constant, which can enable the compiler in some cases to
    // perform additional optimizations.
    const AdvectionField<dim> advection_field;
    const RightHandSide<dim>  right_hand_side;
    const BoundaryValues<dim> boundary_values;

    // Then we define some abbreviations to avoid unnecessarily long lines:
    const unsigned int dofs_per_cell   = fe.dofs_per_cell;
    const unsigned int n_q_points      = scratch_data.fe_values.get_quadrature().size();
    const unsigned int n_face_q_points = scratch_data.fe_face_values.get_quadrature().size();

    // We declare cell matrix and cell right hand side...
    copy_data.cell_matrix.reinit (dofs_per_cell, dofs_per_cell);
    copy_data.cell_rhs.reinit (dofs_per_cell);

    // ... an array to hold the global indices of the degrees of freedom of
    // the cell on which we are presently working...
    copy_data.local_dof_indices.resize(dofs_per_cell);

    // ... and array in which the values of right hand side, advection
    // direction, and boundary values will be stored, for cell and face
    // integrals respectively:
    std::vector<double>         rhs_values (n_q_points);
    std::vector<Tensor<1,dim> > advection_directions (n_q_points);
    std::vector<double>         face_boundary_values (n_face_q_points);
    std::vector<Tensor<1,dim> > face_advection_directions (n_face_q_points);


    // ... then initialize the <code>FEValues</code> object...
    scratch_data.fe_values.reinit (cell);

    // ... obtain the values of right hand side and advection directions
    // at the quadrature points...
    advection_field.value_list (scratch_data.fe_values.get_quadrature_points(),
                                advection_directions);
    right_hand_side.value_list (scratch_data.fe_values.get_quadrature_points(),
                                rhs_values);

    // ... set the value of the streamline diffusion parameter as
    // described in the introduction...
    const double delta = 0.1 * cell->diameter ();

    // ... and assemble the local contributions to the system matrix and
    // right hand side as also discussed above:
    for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            copy_data.cell_matrix(i,j) += ((advection_directions[q_point] *
                                            scratch_data.fe_values.shape_grad(j,q_point)   *
                                            (scratch_data.fe_values.shape_value(i,q_point) +
                                             delta *
                                             (advection_directions[q_point] *
                                              scratch_data.fe_values.shape_grad(i,q_point)))) *
                                           scratch_data.fe_values.JxW(q_point));

          copy_data.cell_rhs(i) += ((scratch_data.fe_values.shape_value(i,q_point) +
                                     delta *
                                     (advection_directions[q_point] *
                                      scratch_data.fe_values.shape_grad(i,q_point))        ) *
                                    rhs_values[q_point] *
                                    scratch_data.fe_values.JxW (q_point));
        }

    // Besides the cell terms which we have built up now, the bilinear
    // form of the present problem also contains terms on the boundary of
    // the domain. Therefore, we have to check whether any of the faces of
    // this cell are on the boundary of the domain, and if so assemble the
    // contributions of this face as well. Of course, the bilinear form
    // only contains contributions from the <code>inflow</code> part of
    // the boundary, but to find out whether a certain part of a face of
    // the present cell is part of the inflow boundary, we have to have
    // information on the exact location of the quadrature points and on
    // the direction of flow at this point; we obtain this information
    // using the FEFaceValues object and only decide within the main loop
    // whether a quadrature point is on the inflow boundary.
    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
      if (cell->face(face)->at_boundary())
        {
          // Ok, this face of the present cell is on the boundary of the
          // domain. Just as for the usual FEValues object which we have
          // used in previous examples and also above, we have to
          // reinitialize the FEFaceValues object for the present face:
          scratch_data.fe_face_values.reinit (cell, face);

          // For the quadrature points at hand, we ask for the values of
          // the inflow function and for the direction of flow:
          boundary_values.value_list (scratch_data.fe_face_values.get_quadrature_points(),
                                      face_boundary_values);
          advection_field.value_list (scratch_data.fe_face_values.get_quadrature_points(),
                                      face_advection_directions);

          // Now loop over all quadrature points and see whether it is on
          // the inflow or outflow part of the boundary. This is
          // determined by a test whether the advection direction points
          // inwards or outwards of the domain (note that the normal
          // vector points outwards of the cell, and since the cell is at
          // the boundary, the normal vector points outward of the domain,
          // so if the advection direction points into the domain, its
          // scalar product with the normal vector must be negative):
          for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
            if (scratch_data.fe_face_values.normal_vector(q_point) *
                face_advection_directions[q_point]
                < 0)
              // If the is part of the inflow boundary, then compute the
              // contributions of this face to the global matrix and right
              // hand side, using the values obtained from the
              // FEFaceValues object and the formulae discussed in the
              // introduction:
              for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                  for (unsigned int j=0; j<dofs_per_cell; ++j)
                    copy_data.cell_matrix(i,j) -= (face_advection_directions[q_point] *
                                                   scratch_data.fe_face_values.normal_vector(q_point) *
                                                   scratch_data.fe_face_values.shape_value(i,q_point) *
                                                   scratch_data.fe_face_values.shape_value(j,q_point) *
                                                   scratch_data.fe_face_values.JxW(q_point));

                  copy_data.cell_rhs(i) -= (face_advection_directions[q_point] *
                                            scratch_data.fe_face_values.normal_vector(q_point) *
                                            face_boundary_values[q_point]         *
                                            scratch_data.fe_face_values.shape_value(i,q_point) *
                                            scratch_data.fe_face_values.JxW(q_point));
                }
        }


    // Now go on by transferring the local contributions to the system of
    // equations into the global objects. The first step was to obtain the
    // global indices of the degrees of freedom on this cell.
    cell->get_dof_indices (copy_data.local_dof_indices);
  }



  // The second function we needed to write was the one that copies
  // the local contributions the previous function has computed and
  // put into the copy data object, into the global matrix and right
  // hand side vector objects. This is essentially what we always had
  // as the last block of code when assembling something on every
  // cell. The following should therefore be pretty obvious:
  template <int dim>
  void
  AdvectionProblem<dim>::copy_local_to_global (const AssemblyCopyData &copy_data)
  {
    for (unsigned int i=0; i<copy_data.local_dof_indices.size(); ++i)
      {
        for (unsigned int j=0; j<copy_data.local_dof_indices.size(); ++j)
          system_matrix.add (copy_data.local_dof_indices[i],
                             copy_data.local_dof_indices[j],
                             copy_data.cell_matrix(i,j));

        system_rhs(copy_data.local_dof_indices[i]) += copy_data.cell_rhs(i);
      }
  }




  // Following is the function that solves the linear system of equations. As
  // the system is no more symmetric positive definite as in all the previous
  // examples, we can't use the Conjugate Gradients method anymore. Rather, we
  // use a solver that is tailored to nonsymmetric systems like the one at
  // hand, the BiCGStab method. As preconditioner, we use the Jacobi method.
  template <int dim>
  void AdvectionProblem<dim>::solve ()
  {
    SolverControl           solver_control (1000, 1e-12);
    SolverBicgstab<>        bicgstab (solver_control);

    PreconditionJacobi<> preconditioner;
    preconditioner.initialize(system_matrix, 1.0);

    bicgstab.solve (system_matrix, solution, system_rhs,
                    preconditioner);

    hanging_node_constraints.distribute (solution);
  }


  // The following function refines the grid according to the quantity
  // described in the introduction. The respective computations are made in
  // the class <code>GradientEstimation</code>. The only difference to
  // previous examples is that we refine a little more aggressively (0.5
  // instead of 0.3 of the number of cells).
  template <int dim>
  void AdvectionProblem<dim>::refine_grid ()
  {
    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

    GradientEstimation::estimate (dof_handler,
                                  solution,
                                  estimated_error_per_cell);

    GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                     estimated_error_per_cell,
                                                     0.5, 0.03);

    triangulation.execute_coarsening_and_refinement ();
  }



  // Writing output to disk is done in the same way as in the previous
  // examples...
  template <int dim>
  void AdvectionProblem<dim>::output_results (const unsigned int cycle) const
  {
    std::string filename = "grid-";
    filename += ('0' + cycle);
    Assert (cycle < 10, ExcInternalError());

    filename += ".eps";
    std::ofstream output (filename.c_str());

    GridOut grid_out;
    grid_out.write_eps (triangulation, output);
  }


  // ... as is the main loop (setup -- solve -- refine)
  template <int dim>
  void AdvectionProblem<dim>::run ()
  {
    for (unsigned int cycle=0; cycle<6; ++cycle)
      {
        std::cout << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          {
            GridGenerator::hyper_cube (triangulation, -1, 1);
            triangulation.refine_global (4);
          }
        else
          {
            refine_grid ();
          }


        std::cout << "   Number of active cells:       "
                  << triangulation.n_active_cells()
                  << std::endl;

        setup_system ();

        std::cout << "   Number of degrees of freedom: "
                  << dof_handler.n_dofs()
                  << std::endl;

        assemble_system ();
        solve ();
        output_results (cycle);
      }

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "solution");
    data_out.build_patches ();

    std::ofstream output ("final-solution.vtk");
    data_out.write_vtk (output);
  }



  // @sect3{GradientEstimation class implementation}

  // Now for the implementation of the <code>GradientEstimation</code> class.
  // Let us start by defining constructors for the
  // <code>EstimateScratchData</code> class used by the
  // <code>estimate_cell()</code> function:
  template <int dim>
  GradientEstimation::EstimateScratchData<dim>::
  EstimateScratchData (const FiniteElement<dim> &fe,
                       const Vector<double>     &solution)
    :
    fe_midpoint_value(fe,
                      QMidpoint<dim> (),
                      update_values | update_quadrature_points),
    solution(solution)
  {}


  template <int dim>
  GradientEstimation::EstimateScratchData<dim>::
  EstimateScratchData(const EstimateScratchData &scratch_data)
    :
    fe_midpoint_value(scratch_data.fe_midpoint_value.get_fe(),
                      scratch_data.fe_midpoint_value.get_quadrature(),
                      update_values | update_quadrature_points),
    solution(scratch_data.solution)
  {}


  // Next for the implementation of the <code>GradientEstimation</code>
  // class. The first function does not much except for delegating work to the
  // other function, but there is a bit of setup at the top.
  //
  // Before starting with the work, we check that the vector into which the
  // results are written has the right size. It is a common error that such
  // parameters have the wrong size, but the resulting damage by not
  // catching these errors are very subtle as they are usually corruption of
  // data somewhere in memory. Often, the problems emerging from this are
  // not reproducible, and it is well worth the effort to
  // check for such things.
  //
  // The second piece is to set up the iterator that goes in lockstep over the
  // cells of the domain and the corresponding elements of the output vector
  // (see above where we introduced the <code>SynchronousIterators</code>
  // class). We can abbreviate the process slightly by introducing a
  // <code>typedef</code> that denotes a pair of iterators. This being set up,
  // we can hand the whole thing off to WorkStream::run, keeping in mind that
  // we do not need a copy-local-to-global function here but can get away by
  // simply using a default-constructed function object (the equivalent to a
  // NULL function pointer).
  template <int dim>
  void
  GradientEstimation::estimate (const DoFHandler<dim> &dof_handler,
                                const Vector<double>  &solution,
                                Vector<float>         &error_per_cell)
  {
    Assert (error_per_cell.size() == dof_handler.get_tria().n_active_cells(),
            ExcInvalidVectorLength (error_per_cell.size(),
                                    dof_handler.get_tria().n_active_cells()));

    typedef std_cxx11::tuple<typename DoFHandler<dim>::active_cell_iterator,Vector<float>::iterator>
    IteratorTuple;

    SynchronousIterators<IteratorTuple>
    begin_sync_it (IteratorTuple (dof_handler.begin_active(),
                                  error_per_cell.begin())),
                                                       end_sync_it (IteratorTuple (dof_handler.end(),
                                                           error_per_cell.end()));

    WorkStream::run (begin_sync_it,
                     end_sync_it,
                     &GradientEstimation::template estimate_cell<dim>,
                     std_cxx11::function<void (const EstimateCopyData &)> (),
                     EstimateScratchData<dim> (dof_handler.get_fe(),
                                               solution),
                     EstimateCopyData ());
  }


  // Following now the function that actually computes the finite difference
  // approximation to the gradient. The general outline of the function is to
  // first compute the list of active neighbors of the present cell and then
  // compute the quantities described in the introduction for each of the
  // neighbors. The reason for this order is that it is not a one-liner to
  // find a given neighbor with locally refined meshes. In principle, an
  // optimized implementation would find neighbors and the quantities
  // depending on them in one step, rather than first building a list of
  // neighbors and in a second step their contributions but we will gladly
  // leave this as an exercise. As discussed before, the worker function
  // passed to WorkStream::run works on "scratch" objects that keep all
  // temporary objects. This way, we do not need to create and initialize
  // objects that are expensive to initialize within the function that does
  // the work, every time it is called for a given cell. Such an argument is
  // passed as the second argument. The third argument would be a "copy-data"
  // object (see @ref threads for more information) but we do not actually use
  // any of these here. Because WorkStream::run insists on passing three
  // arguments, we declare this function with three arguments, but simply
  // ignore the last one.
  //
  // (This is unsatisfactory from an esthetic perspective. It can be avoided,
  // at the cost of some other trickery. If you allow, let us here show
  // how. First, assume that we had declared this function to only take two
  // arguments by omitting the unused last one. Now, WorkStream::run still
  // wants to call this function with three arguments, so we need to find a
  // way to "forget" the third argument in the call. Simply passing
  // WorkStream::run the pointer to the function as we do above will not do
  // this -- the compiler will complain that a function declared to have two
  // arguments is called with three arguments.  However, we can do this by
  // passing the following as the third argument when calling WorkStream::run
  // above:
  // @code
  //    std_cxx11::function<void (const SynchronousIterators<IteratorTuple> &,
  //                              EstimateScratchData<dim>                  &,
  //                              EstimateCopyData                          &)>
  //      (std_cxx11::bind (&GradientEstimation::template estimate_cell<dim>,
  //                        std_cxx11::_1,
  //                        std_cxx11::_2))
  // @endcode
  // This creates a function object taking three arguments, but when it calls
  // the underlying function object, it simply only uses the first and second
  // argument -- we simply "forget" to use the third argument :-)
  // In the end, this isn't completely obvious either, and so we didn't implement
  // it, but hey -- it can be done!)
  //
  // Now for the details:
  template <int dim>
  void
  GradientEstimation::estimate_cell (const SynchronousIterators<std_cxx11::tuple<typename DoFHandler<dim>::active_cell_iterator,
                                     Vector<float>::iterator> > &cell,
                                     EstimateScratchData<dim>                                               &scratch_data,
                                     const EstimateCopyData &)
  {
    // We need space for the tensor <code>Y</code>, which is the sum of
    // outer products of the y-vectors.
    Tensor<2,dim> Y;


    // Then we allocate a vector to hold iterators to all active neighbors of
    // a cell. We reserve the maximal number of active neighbors in order to
    // avoid later reallocations. Note how this maximal number of active
    // neighbors is computed here.
    std::vector<typename DoFHandler<dim>::active_cell_iterator> active_neighbors;
    active_neighbors.reserve (GeometryInfo<dim>::faces_per_cell *
                              GeometryInfo<dim>::max_children_per_face);

    typename DoFHandler<dim>::active_cell_iterator cell_it(std_cxx11::get<0>(cell.iterators));

    // First initialize the <code>FEValues</code> object, as well as the
    // <code>Y</code> tensor:
    scratch_data.fe_midpoint_value.reinit (cell_it);

    // Then allocate the vector that will be the sum over the y-vectors
    // times the approximate directional derivative:
    Tensor<1,dim> projected_gradient;


    // Now before going on first compute a list of all active neighbors of
    // the present cell. We do so by first looping over all faces and see
    // whether the neighbor there is active, which would be the case if it
    // is on the same level as the present cell or one level coarser (note
    // that a neighbor can only be once coarser than the present cell, as
    // we only allow a maximal difference of one refinement over a face in
    // deal.II). Alternatively, the neighbor could be on the same level
    // and be further refined; then we have to find which of its children
    // are next to the present cell and select these (note that if a child
    // of of neighbor of an active cell that is next to this active cell,
    // needs necessarily be active itself, due to the one-refinement rule
    // cited above).
    //
    // Things are slightly different in one space dimension, as there the
    // one-refinement rule does not exist: neighboring active cells may
    // differ in as many refinement levels as they like. In this case, the
    // computation becomes a little more difficult, but we will explain
    // this below.
    //
    // Before starting the loop over all neighbors of the present cell, we
    // have to clear the array storing the iterators to the active
    // neighbors, of course.
    active_neighbors.clear ();
    for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
      if (! std_cxx11::get<0>(cell.iterators)->at_boundary(face_no))
        {
          // First define an abbreviation for the iterator to the face and
          // the neighbor
          const typename DoFHandler<dim>::face_iterator
          face = std_cxx11::get<0>(cell.iterators)->face(face_no);
          const typename DoFHandler<dim>::cell_iterator
          neighbor = std_cxx11::get<0>(cell.iterators)->neighbor(face_no);

          // Then check whether the neighbor is active. If it is, then it
          // is on the same level or one level coarser (if we are not in
          // 1D), and we are interested in it in any case.
          if (neighbor->active())
            active_neighbors.push_back (neighbor);
          else
            {
              // If the neighbor is not active, then check its children.
              if (dim == 1)
                {
                  // To find the child of the neighbor which bounds to the
                  // present cell, successively go to its right child if
                  // we are left of the present cell (n==0), or go to the
                  // left child if we are on the right (n==1), until we
                  // find an active cell.
                  typename DoFHandler<dim>::cell_iterator
                  neighbor_child = neighbor;
                  while (neighbor_child->has_children())
                    neighbor_child = neighbor_child->child (face_no==0 ? 1 : 0);

                  // As this used some non-trivial geometrical intuition,
                  // we might want to check whether we did it right,
                  // i.e. check whether the neighbor of the cell we found
                  // is indeed the cell we are presently working
                  // on. Checks like this are often useful and have
                  // frequently uncovered errors both in algorithms like
                  // the line above (where it is simple to involuntarily
                  // exchange <code>n==1</code> for <code>n==0</code> or
                  // the like) and in the library (the assumptions
                  // underlying the algorithm above could either be wrong,
                  // wrongly documented, or are violated due to an error
                  // in the library). One could in principle remove such
                  // checks after the program works for some time, but it
                  // might be a good things to leave it in anyway to check
                  // for changes in the library or in the algorithm above.
                  //
                  // Note that if this check fails, then this is certainly
                  // an error that is irrecoverable and probably qualifies
                  // as an internal error. We therefore use a predefined
                  // exception class to throw here.
                  Assert (neighbor_child->neighbor(face_no==0 ? 1 : 0)
                          ==std_cxx11::get<0>(cell.iterators),ExcInternalError());

                  // If the check succeeded, we push the active neighbor
                  // we just found to the stack we keep:
                  active_neighbors.push_back (neighbor_child);
                }
              else
                // If we are not in 1d, we collect all neighbor children
                // `behind' the subfaces of the current face
                for (unsigned int subface_no=0; subface_no<face->n_children(); ++subface_no)
                  active_neighbors.push_back (
                    std_cxx11::get<0>(cell.iterators)->neighbor_child_on_subface(face_no,subface_no));
            }
        }

    // OK, now that we have all the neighbors, lets start the computation
    // on each of them. First we do some preliminaries: find out about the
    // center of the present cell and the solution at this point. The
    // latter is obtained as a vector of function values at the quadrature
    // points, of which there are only one, of course. Likewise, the
    // position of the center is the position of the first (and only)
    // quadrature point in real space.
    const Point<dim> this_center = scratch_data.fe_midpoint_value.quadrature_point(0);

    std::vector<double> this_midpoint_value(1);
    scratch_data.fe_midpoint_value.get_function_values (scratch_data.solution, this_midpoint_value);


    // Now loop over all active neighbors and collect the data we
    // need. Allocate a vector just like <code>this_midpoint_value</code>
    // which we will use to store the value of the solution in the
    // midpoint of the neighbor cell. We allocate it here already, since
    // that way we don't have to allocate memory repeatedly in each
    // iteration of this inner loop (memory allocation is a rather
    // expensive operation):
    std::vector<double> neighbor_midpoint_value(1);
    typename std::vector<typename DoFHandler<dim>::active_cell_iterator>::const_iterator
    neighbor_ptr = active_neighbors.begin();
    for (; neighbor_ptr!=active_neighbors.end(); ++neighbor_ptr)
      {
        // First define an abbreviation for the iterator to the active
        // neighbor cell:
        const typename DoFHandler<dim>::active_cell_iterator
        neighbor = *neighbor_ptr;

        // Then get the center of the neighbor cell and the value of the
        // finite element function thereon. Note that for this information
        // we have to reinitialize the <code>FEValues</code> object for
        // the neighbor cell.
        scratch_data.fe_midpoint_value.reinit (neighbor);
        const Point<dim> neighbor_center = scratch_data.fe_midpoint_value.quadrature_point(0);

        scratch_data.fe_midpoint_value.get_function_values (scratch_data.solution,
                                                            neighbor_midpoint_value);

        // Compute the vector <code>y</code> connecting the centers of the
        // two cells. Note that as opposed to the introduction, we denote
        // by <code>y</code> the normalized difference vector, as this is
        // the quantity used everywhere in the computations.
        Point<dim>   y        = neighbor_center - this_center;
        const double distance = std::sqrt(y.square());
        y /= distance;

        // Then add up the contribution of this cell to the Y matrix...
        for (unsigned int i=0; i<dim; ++i)
          for (unsigned int j=0; j<dim; ++j)
            Y[i][j] += y[i] * y[j];

        // ... and update the sum of difference quotients:
        projected_gradient += (neighbor_midpoint_value[0] -
                               this_midpoint_value[0]) /
                              distance *
                              y;
      }

    // If now, after collecting all the information from the neighbors, we
    // can determine an approximation of the gradient for the present
    // cell, then we need to have passed over vectors <code>y</code> which
    // span the whole space, otherwise we would not have all components of
    // the gradient. This is indicated by the invertibility of the matrix.
    //
    // If the matrix should not be invertible, this means that the present
    // cell had an insufficient number of active neighbors. In contrast to
    // all previous cases, where we raised exceptions, this is, however,
    // not a programming error: it is a runtime error that can happen in
    // optimized mode even if it ran well in debug mode, so it is
    // reasonable to try to catch this error also in optimized mode. For
    // this case, there is the <code>AssertThrow</code> macro: it checks
    // the condition like the <code>Assert</code> macro, but not only in
    // debug mode; it then outputs an error message, but instead of
    // terminating the program as in the case of the <code>Assert</code>
    // macro, the exception is thrown using the <code>throw</code> command
    // of C++. This way, one has the possibility to catch this error and
    // take reasonable counter actions. One such measure would be to
    // refine the grid globally, as the case of insufficient directions
    // can not occur if every cell of the initial grid has been refined at
    // least once.
    AssertThrow (determinant(Y) != 0,
                 ExcInsufficientDirections());

    // If, on the other hand the matrix is invertible, then invert it,
    // multiply the other quantity with it and compute the estimated error
    // using this quantity and the right powers of the mesh width:
    const Tensor<2,dim> Y_inverse = invert(Y);

    Point<dim> gradient;
    contract (gradient, Y_inverse, projected_gradient);

    // The last part of this function is the one where we
    // write into the element of the output vector what
    // we have just computed. As above, we need to get
    // at the second element of the pair of iterators, which requires
    // slightly awkward syntax but is not otherwise particularly
    // difficult:
    *(std_cxx11::get<1>(cell.iterators)) = (std::pow(std_cxx11::get<0>(cell.iterators)->diameter(),
                                                     1+1.0*dim/2) *
                                            std::sqrt(gradient.square()));

  }
}


// @sect3{Main function}

// The <code>main</code> function is exactly like in previous examples, with
// the only difference in the name of the main class that actually does the
// computation.
int main ()
{
  try
    {
      dealii::deallog.depth_console (0);

      Step9::AdvectionProblem<2> advection_problem_2d;
      advection_problem_2d.run ();
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

/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2007 - 2016 by the deal.II authors
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
 * Author: Moritz Allmaras, Texas A&M University, 2007
 */


// @sect3{Include files}

// The following header files are unchanged from step-7 and have been
// discussed before:

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <fstream>

// This header file contains the necessary declarations for the
// ParameterHandler class that we will use to read our parameters from a
// configuration file:
#include <deal.II/base/parameter_handler.h>

// For solving the linear system, we'll use the sparse LU-decomposition
// provided by UMFPACK (see the SparseDirectUMFPACK class), for which the
// following header file is needed.  Note that in order to compile this
// tutorial program, the deal.II-library needs to be built with UMFPACK
// support, which is enabled by default:
#include <deal.II/lac/sparse_direct.h>

// The FESystem class allows us to stack several FE-objects to one compound,
// vector-valued finite element field. The necessary declarations for this
// class are provided in this header file:
#include <deal.II/fe/fe_system.h>

// Finally, include the header file that declares the Timer class that we will
// use to determine how much time each of the operations of our program takes:
#include <deal.II/base/timer.h>

// As the last step at the beginning of this program, we put everything that
// is in this program into its namespace and, within it, make everything that
// is in the deal.II namespace globally available, without the need to prefix
// everything with <code>dealii</code><code>::</code>:
namespace Step29
{
  using namespace dealii;


  // @sect3{The <code>DirichletBoundaryValues</code> class}

  // First we define a class for the function representing the Dirichlet
  // boundary values. This has been done many times before and therefore does
  // not need much explanation.
  //
  // Since there are two values $v$ and $w$ that need to be prescribed at the
  // boundary, we have to tell the base class that this is a vector-valued
  // function with two components, and the <code>vector_value</code> function
  // and its cousin <code>vector_value_list</code> must return vectors with
  // two entries. In our case the function is very simple, it just returns 1
  // for the real part $v$ and 0 for the imaginary part $w$ regardless of the
  // point where it is evaluated.
  template <int dim>
  class DirichletBoundaryValues : public Function<dim>
  {
  public:
    DirichletBoundaryValues() : Function<dim> (2) {};

    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &values) const;

    virtual void vector_value_list (const std::vector<Point<dim> > &points,
                                    std::vector<Vector<double> >   &value_list) const;
  };


  template <int dim>
  inline
  void DirichletBoundaryValues<dim>::vector_value (const Point<dim> &/*p*/,
                                                   Vector<double>   &values) const
  {
    Assert (values.size() == 2, ExcDimensionMismatch (values.size(), 2));

    values(0) = 1;
    values(1) = 0;
  }


  template <int dim>
  void DirichletBoundaryValues<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                                        std::vector<Vector<double> >   &value_list) const
  {
    Assert (value_list.size() == points.size(),
            ExcDimensionMismatch (value_list.size(), points.size()));

    for (unsigned int p=0; p<points.size(); ++p)
      DirichletBoundaryValues<dim>::vector_value (points[p], value_list[p]);
  }

  // @sect3{The <code>ParameterReader</code> class}

  // The next class is responsible for preparing the ParameterHandler object
  // and reading parameters from an input file.  It includes a function
  // <code>declare_parameters</code> that declares all the necessary
  // parameters and a <code>read_parameters</code> function that is called
  // from outside to initiate the parameter reading process.
  class ParameterReader : public Subscriptor
  {
  public:
    ParameterReader(ParameterHandler &);
    void read_parameters(const std::string);

  private:
    void declare_parameters();
    ParameterHandler &prm;
  };

  // The constructor stores a reference to the ParameterHandler object that is
  // passed to it:
  ParameterReader::ParameterReader(ParameterHandler &paramhandler)
    :
    prm(paramhandler)
  {}

  // @sect4{<code>ParameterReader::declare_parameters</code>}

  // The <code>declare_parameters</code> function declares all the parameters
  // that our ParameterHandler object will be able to read from input files,
  // along with their types, range conditions and the subsections they appear
  // in. We will wrap all the entries that go into a section in a pair of
  // braces to force the editor to indent them by one level, making it simpler
  // to read which entries together form a section:
  void ParameterReader::declare_parameters()
  {
    // Parameters for mesh and geometry include the number of global
    // refinement steps that are applied to the initial coarse mesh and the
    // focal distance $d$ of the transducer lens. For the number of refinement
    // steps, we allow integer values in the range $[0,\infty)$, where the
    // omitted second argument to the Patterns::Integer object denotes the
    // half-open interval.  For the focal distance any number greater than
    // zero is accepted:
    prm.enter_subsection ("Mesh & geometry parameters");
    {
      prm.declare_entry("Number of refinements", "6",
                        Patterns::Integer(0),
                        "Number of global mesh refinement steps "
                        "applied to initial coarse grid");

      prm.declare_entry("Focal distance", "0.3",
                        Patterns::Double(0),
                        "Distance of the focal point of the lens "
                        "to the x-axis");
    }
    prm.leave_subsection ();

    // The next subsection is devoted to the physical parameters appearing in
    // the equation, which are the frequency $\omega$ and wave speed
    // $c$. Again, both need to lie in the half-open interval $[0,\infty)$
    // represented by calling the Patterns::Double class with only the left
    // end-point as argument:
    prm.enter_subsection ("Physical constants");
    {
      prm.declare_entry("c", "1.5e5",
                        Patterns::Double(0),
                        "Wave speed");

      prm.declare_entry("omega", "5.0e7",
                        Patterns::Double(0),
                        "Frequency");
    }
    prm.leave_subsection ();


    // Last but not least we would like to be able to change some properties
    // of the output, like filename and format, through entries in the
    // configuration file, which is the purpose of the last subsection:
    prm.enter_subsection ("Output parameters");
    {
      prm.declare_entry("Output file", "solution",
                        Patterns::Anything(),
                        "Name of the output file (without extension)");

      // Since different output formats may require different parameters for
      // generating output (like for example, postscript output needs
      // viewpoint angles, line widths, colors etc), it would be cumbersome if
      // we had to declare all these parameters by hand for every possible
      // output format supported in the library. Instead, each output format
      // has a <code>FormatFlags::declare_parameters</code> function, which
      // declares all the parameters specific to that format in an own
      // subsection. The following call of
      // DataOutInterface<1>::declare_parameters executes
      // <code>declare_parameters</code> for all available output formats, so
      // that for each format an own subsection will be created with
      // parameters declared for that particular output format. (The actual
      // value of the template parameter in the call, <code>@<1@></code>
      // above, does not matter here: the function does the same work
      // independent of the dimension, but happens to be in a
      // template-parameter-dependent class.)  To find out what parameters
      // there are for which output format, you can either consult the
      // documentation of the DataOutBase class, or simply run this program
      // without a parameter file present. It will then create a file with all
      // declared parameters set to their default values, which can
      // conveniently serve as a starting point for setting the parameters to
      // the values you desire.
      DataOutInterface<1>::declare_parameters (prm);
    }
    prm.leave_subsection ();
  }

  // @sect4{<code>ParameterReader::read_parameters</code>}

  // This is the main function in the ParameterReader class.  It gets called
  // from outside, first declares all the parameters, and then reads them from
  // the input file whose filename is provided by the caller. After the call
  // to this function is complete, the <code>prm</code> object can be used to
  // retrieve the values of the parameters read in from the file:
  void ParameterReader::read_parameters (const std::string parameter_file)
  {
    declare_parameters();

    prm.read_input (parameter_file);
  }



  // @sect3{The <code>ComputeIntensity</code> class}

  // As mentioned in the introduction, the quantity that we are really after
  // is the spatial distribution of the intensity of the ultrasound wave,
  // which corresponds to $|u|=\sqrt{v^2+w^2}$. Now we could just be content
  // with having $v$ and $w$ in our output, and use a suitable visualization
  // or postprocessing tool to derive $|u|$ from the solution we
  // computed. However, there is also a way to output data derived from the
  // solution in deal.II, and we are going to make use of this mechanism here.

  // So far we have always used the DataOut::add_data_vector function to add
  // vectors containing output data to a DataOut object.  There is a special
  // version of this function that in addition to the data vector has an
  // additional argument of type DataPostprocessor. What happens when this
  // function is used for output is that at each point where output data is to
  // be generated, the DataPostprocessor::compute_derived_quantities_scalar or
  // DataPostprocessor::compute_derived_quantities_vector function of the
  // specified DataPostprocessor object is invoked to compute the output
  // quantities from the values, the gradients and the second derivatives of
  // the finite element function represented by the data vector (in the case
  // of face related data, normal vectors are available as well). Hence, this
  // allows us to output any quantity that can locally be derived from the
  // values of the solution and its derivatives.  Of course, the ultrasound
  // intensity $|u|$ is such a quantity and its computation doesn't even
  // involve any derivatives of $v$ or $w$.

  // In practice, the DataPostprocessor class only provides an interface to
  // this functionality, and we need to derive our own class from it in order
  // to implement the functions specified by the interface. In the most
  // general case one has to implement several member functions but if the
  // output quantity is a single scalar then some of this boilerplate code can
  // be handled by a more specialized class, DataPostprocessorScalar and we
  // can derive from that one instead. This is what the
  // <code>ComputeIntensity</code> class does:
  template <int dim>
  class ComputeIntensity : public DataPostprocessorScalar<dim>
  {
  public:
    ComputeIntensity ();

    virtual
    void
    compute_derived_quantities_vector (const std::vector<Vector<double> >               &uh,
                                       const std::vector<std::vector<Tensor<1, dim> > > &duh,
                                       const std::vector<std::vector<Tensor<2, dim> > > &dduh,
                                       const std::vector<Point<dim> >                   &normals,
                                       const std::vector<Point<dim> >                   &evaluation_points,
                                       std::vector<Vector<double> >                     &computed_quantities) const;
  };

  // In the constructor, we need to call the constructor of the base class
  // with two arguments. The first denotes the name by which the single scalar
  // quantity computed by this class should be represented in output files. In
  // our case, the postprocessor has $|u|$ as output, so we use "Intensity".
  //
  // The second argument is a set of flags that indicate which data is needed
  // by the postprocessor in order to compute the output quantities.  This can
  // be any subset of update_values, update_gradients and update_hessians
  // (and, in the case of face data, also update_normal_vectors), which are
  // documented in UpdateFlags.  Of course, computation of the derivatives
  // requires additional resources, so only the flags for data that are really
  // needed should be given here, just as we do when we use FEValues objects.
  // In our case, only the function values of $v$ and $w$ are needed to
  // compute $|u|$, so we're good with the update_values flag.
  template <int dim>
  ComputeIntensity<dim>::ComputeIntensity ()
    :
    DataPostprocessorScalar<dim> ("Intensity",
                                  update_values)
  {}


  // The actual postprocessing happens in the following function.  Its inputs
  // are a vector representing values of the function (which is here
  // vector-valued) representing the data vector given to
  // DataOut::add_data_vector, evaluated at all evaluation points where we
  // generate output, and some tensor objects representing derivatives (that
  // we don't use here since $|u|$ is computed from just $v$ and $w$, and for
  // which we assign no name to the corresponding function argument).  The
  // derived quantities are returned in the <code>computed_quantities</code>
  // vector.  Remember that this function may only use data for which the
  // respective update flag is specified by
  // <code>get_needed_update_flags</code>. For example, we may not use the
  // derivatives here, since our implementation of
  // <code>get_needed_update_flags</code> requests that only function values
  // are provided.
  template <int dim>
  void
  ComputeIntensity<dim>::compute_derived_quantities_vector (
    const std::vector<Vector<double> >                 &uh,
    const std::vector<std::vector<Tensor<1, dim> > >   & /*duh*/,
    const std::vector<std::vector<Tensor<2, dim> > >   & /*dduh*/,
    const std::vector<Point<dim> >                     & /*normals*/,
    const std::vector<Point<dim> >                     & /*evaluation_points*/,
    std::vector<Vector<double> >                       &computed_quantities
  ) const
  {
    Assert(computed_quantities.size() == uh.size(),
           ExcDimensionMismatch (computed_quantities.size(), uh.size()));

    // The computation itself is straightforward: We iterate over each entry
    // in the output vector and compute $|u|$ from the corresponding values of
    // $v$ and $w$:
    for (unsigned int i=0; i<computed_quantities.size(); i++)
      {
        Assert(computed_quantities[i].size() == 1,
               ExcDimensionMismatch (computed_quantities[i].size(), 1));
        Assert(uh[i].size() == 2, ExcDimensionMismatch (uh[i].size(), 2));

        computed_quantities[i](0) = std::sqrt(uh[i](0)*uh[i](0) + uh[i](1)*uh[i](1));
      }
  }


  // @sect3{The <code>UltrasoundProblem</code> class}

  // Finally here is the main class of this program.  It's member functions
  // are very similar to the previous examples, in particular step-4, and the
  // list of member variables does not contain any major surprises either.
  // The ParameterHandler object that is passed to the constructor is stored
  // as a reference to allow easy access to the parameters from all functions
  // of the class.  Since we are working with vector valued finite elements,
  // the FE object we are using is of type FESystem.
  template <int dim>
  class UltrasoundProblem
  {
  public:
    UltrasoundProblem (ParameterHandler &);
    ~UltrasoundProblem ();
    void run ();

  private:
    void make_grid ();
    void setup_system ();
    void assemble_system ();
    void solve ();
    void output_results () const;

    ParameterHandler      &prm;

    Triangulation<dim>     triangulation;
    DoFHandler<dim>        dof_handler;
    FESystem<dim>          fe;

    SparsityPattern        sparsity_pattern;
    SparseMatrix<double>   system_matrix;
    Vector<double>         solution, system_rhs;
  };



  // The constructor takes the ParameterHandler object and stores it in a
  // reference. It also initializes the DoF-Handler and the finite element
  // system, which consists of two copies of the scalar Q1 field, one for $v$
  // and one for $w$:
  template <int dim>
  UltrasoundProblem<dim>::UltrasoundProblem (ParameterHandler  &param)
    :
    prm(param),
    dof_handler(triangulation),
    fe(FE_Q<dim>(1), 2)
  {}


  template <int dim>
  UltrasoundProblem<dim>::~UltrasoundProblem ()
  {
    dof_handler.clear();
  }

  // @sect4{<code>UltrasoundProblem::make_grid</code>}

  // Here we setup the grid for our domain.  As mentioned in the exposition,
  // the geometry is just a unit square (in 2d) with the part of the boundary
  // that represents the transducer lens replaced by a sector of a circle.
  template <int dim>
  void UltrasoundProblem<dim>::make_grid ()
  {
    // First we generate some logging output and start a timer so we can
    // compute execution time when this function is done:
    deallog << "Generating grid... ";
    Timer timer;
    timer.start ();

    // Then we query the values for the focal distance of the transducer lens
    // and the number of mesh refinement steps from our ParameterHandler
    // object:
    prm.enter_subsection ("Mesh & geometry parameters");

    const double                focal_distance = prm.get_double("Focal distance");
    const unsigned int  n_refinements  = prm.get_integer("Number of refinements");

    prm.leave_subsection ();

    // Next, two points are defined for position and focal point of the
    // transducer lens, which is the center of the circle whose segment will
    // form the transducer part of the boundary. Notice that this is the only
    // point in the program where things are slightly different in 2D and 3D.
    // Even though this tutorial only deals with the 2D case, the necessary
    // additions to make this program functional in 3D are so minimal that we
    // opt for including them:
    const Point<dim>    transducer = (dim == 2) ?
                                     Point<dim> (0.5, 0.0) :
                                     Point<dim> (0.5, 0.5, 0.0);
    const Point<dim>   focal_point = (dim == 2) ?
                                     Point<dim> (0.5, focal_distance) :
                                     Point<dim> (0.5, 0.5, focal_distance);


    // As initial coarse grid we take a simple unit square with 5 subdivisions
    // in each direction. The number of subdivisions is chosen so that the
    // line segment $[0.4,0.6]$ that we want to designate as the transducer
    // boundary is spanned by a single face. Then we step through all cells to
    // find the faces where the transducer is to be located, which in fact is
    // just the single edge from 0.4 to 0.6 on the x-axis. This is where we
    // want the refinements to be made according to a circle shaped boundary,
    // so we mark this edge with a different manifold indicator. Since we will
    // Dirichlet boundary conditions on the transducer, we also change its
    // boundary indicator.
    GridGenerator::subdivided_hyper_cube (triangulation, 5, 0, 1);

    typename Triangulation<dim>::cell_iterator
    cell = triangulation.begin (),
    endc = triangulation.end();

    for (; cell!=endc; ++cell)
      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
        if ( cell->face(face)->at_boundary() &&
             ((cell->face(face)->center() - transducer).norm_square() < 0.01) )
          {

            cell->face(face)->set_boundary_id (1);
            cell->face(face)->set_manifold_id (1);
          }
    // For the circle part of the transducer lens, a SphericalManifold object
    // is used (which, of course, in 2D just represents a circle), with center
    // computed as above. By marking this object as <code>static</code>, we
    // ensure that it lives until the end of the program and thereby longer
    // than the triangulation object we will associate with it. We then assign
    // this boundary-object to the part of the boundary with boundary indicator 1:
    static const SphericalManifold<dim> boundary(focal_point);
    triangulation.set_manifold(1, boundary);

    // Now global refinement is executed. Cells near the transducer location
    // will be automatically refined according to the circle shaped boundary
    // of the transducer lens:
    triangulation.refine_global (n_refinements);

    // Lastly, we generate some more logging output. We stop the timer and
    // query the number of CPU seconds elapsed since the beginning of the
    // function:
    timer.stop ();
    deallog << "done ("
            << timer()
            << "s)"
            << std::endl;

    deallog << "  Number of active cells:  "
            << triangulation.n_active_cells()
            << std::endl;
  }


  // @sect4{<code>UltrasoundProblem::setup_system</code>}
  //
  // Initialization of the system matrix, sparsity patterns and vectors are
  // the same as in previous examples and therefore do not need further
  // comment. As in the previous function, we also output the run time of what
  // we do here:
  template <int dim>
  void UltrasoundProblem<dim>::setup_system ()
  {
    deallog << "Setting up system... ";
    Timer timer;
    timer.start();

    dof_handler.distribute_dofs (fe);

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, dsp);
    sparsity_pattern.copy_from (dsp);

    system_matrix.reinit (sparsity_pattern);
    system_rhs.reinit (dof_handler.n_dofs());
    solution.reinit (dof_handler.n_dofs());

    timer.stop ();
    deallog << "done ("
            << timer()
            << "s)"
            << std::endl;

    deallog << "  Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << std::endl;
  }


  // @sect4{<code>UltrasoundProblem::assemble_system</code>}

  // As before, this function takes care of assembling the system matrix and
  // right hand side vector:
  template <int dim>
  void UltrasoundProblem<dim>::assemble_system ()
  {
    deallog << "Assembling system matrix... ";
    Timer timer;
    timer.start ();

    // First we query wavespeed and frequency from the ParameterHandler object
    // and store them in local variables, as they will be used frequently
    // throughout this function.

    prm.enter_subsection ("Physical constants");

    const double omega = prm.get_double("omega"),
                 c     = prm.get_double("c");

    prm.leave_subsection ();

    // As usual, for computing integrals ordinary Gauss quadrature rule is
    // used. Since our bilinear form involves boundary integrals on
    // $\Gamma_2$, we also need a quadrature rule for surface integration on
    // the faces, which are $dim-1$ dimensional:
    QGauss<dim>    quadrature_formula(2);
    QGauss<dim-1>  face_quadrature_formula(2);

    const unsigned int n_q_points       = quadrature_formula.size(),
                       n_face_q_points  = face_quadrature_formula.size(),
                       dofs_per_cell    = fe.dofs_per_cell;

    // The FEValues objects will evaluate the shape functions for us.  For the
    // part of the bilinear form that involves integration on $\Omega$, we'll
    // need the values and gradients of the shape functions, and of course the
    // quadrature weights.  For the terms involving the boundary integrals,
    // only shape function values and the quadrature weights are necessary.
    FEValues<dim>  fe_values (fe, quadrature_formula,
                              update_values | update_gradients |
                              update_JxW_values);

    FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula,
                                      update_values | update_JxW_values);

    // As usual, the system matrix is assembled cell by cell, and we need a
    // matrix for storing the local cell contributions as well as an index
    // vector to transfer the cell contributions to the appropriate location
    // in the global system matrix after.
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    for (; cell!=endc; ++cell)
      {

        // On each cell, we first need to reset the local contribution matrix
        // and request the FEValues object to compute the shape functions for
        // the current cell:
        cell_matrix = 0;
        fe_values.reinit (cell);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {

                // At this point, it is important to keep in mind that we are
                // dealing with a finite element system with two
                // components. Due to the way we constructed this FESystem,
                // namely as the Cartesian product of two scalar finite
                // element fields, each shape function has only a single
                // nonzero component (they are, in deal.II lingo, @ref
                // GlossPrimitive "primitive").  Hence, each shape function
                // can be viewed as one of the $\phi$'s or $\psi$'s from the
                // introduction, and similarly the corresponding degrees of
                // freedom can be attributed to either $\alpha$ or $\beta$.
                // As we iterate through all the degrees of freedom on the
                // current cell however, they do not come in any particular
                // order, and so we cannot decide right away whether the DoFs
                // with index $i$ and $j$ belong to the real or imaginary part
                // of our solution.  On the other hand, if you look at the
                // form of the system matrix in the introduction, this
                // distinction is crucial since it will determine to which
                // block in the system matrix the contribution of the current
                // pair of DoFs will go and hence which quantity we need to
                // compute from the given two shape functions.  Fortunately,
                // the FESystem object can provide us with this information,
                // namely it has a function
                // FESystem::system_to_component_index, that for each local
                // DoF index returns a pair of integers of which the first
                // indicates to which component of the system the DoF
                // belongs. The second integer of the pair indicates which
                // index the DoF has in the scalar base finite element field,
                // but this information is not relevant here. If you want to
                // know more about this function and the underlying scheme
                // behind primitive vector valued elements, take a look at
                // step-8 or the @ref vector_valued module, where these topics
                // are explained in depth.
                if (fe.system_to_component_index(i).first ==
                    fe.system_to_component_index(j).first)
                  {

                    // If both DoFs $i$ and $j$ belong to same component,
                    // i.e. their shape functions are both $\phi$'s or both
                    // $\psi$'s, the contribution will end up in one of the
                    // diagonal blocks in our system matrix, and since the
                    // corresponding entries are computed by the same formula,
                    // we do not bother if they actually are $\phi$ or $\psi$
                    // shape functions. We can simply compute the entry by
                    // iterating over all quadrature points and adding up
                    // their contributions, where values and gradients of the
                    // shape functions are supplied by our FEValues object.

                    for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
                      cell_matrix(i,j) += (((fe_values.shape_value(i,q_point) *
                                             fe_values.shape_value(j,q_point)) *
                                            (- omega * omega)
                                            +
                                            (fe_values.shape_grad(i,q_point) *
                                             fe_values.shape_grad(j,q_point)) *
                                            c * c) *
                                           fe_values.JxW(q_point));

                    // You might think that we would have to specify which
                    // component of the shape function we'd like to evaluate
                    // when requesting shape function values or gradients from
                    // the FEValues object. However, as the shape functions
                    // are primitive, they have only one nonzero component,
                    // and the FEValues class is smart enough to figure out
                    // that we are definitely interested in this one nonzero
                    // component.
                  }
              }
          }


        // We also have to add contributions due to boundary terms. To this
        // end, we loop over all faces of the current cell and see if first it
        // is at the boundary, and second has the correct boundary indicator
        // associated with $\Gamma_2$, the part of the boundary where we have
        // absorbing boundary conditions:
        for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
          if (cell->face(face)->at_boundary() &&
              (cell->face(face)->boundary_id() == 0) )
            {


              // These faces will certainly contribute to the off-diagonal
              // blocks of the system matrix, so we ask the FEFaceValues
              // object to provide us with the shape function values on this
              // face:
              fe_face_values.reinit (cell, face);


              // Next, we loop through all DoFs of the current cell to find
              // pairs that belong to different components and both have
              // support on the current face:
              for (unsigned int i=0; i<dofs_per_cell; ++i)
                for (unsigned int j=0; j<dofs_per_cell; ++j)
                  if ((fe.system_to_component_index(i).first !=
                       fe.system_to_component_index(j).first) &&
                      fe.has_support_on_face(i, face) &&
                      fe.has_support_on_face(j, face))
                    // The check whether shape functions have support on a
                    // face is not strictly necessary: if we don't check for
                    // it we would simply add up terms to the local cell
                    // matrix that happen to be zero because at least one of
                    // the shape functions happens to be zero. However, we can
                    // save that work by adding the checks above.

                    // In either case, these DoFs will contribute to the
                    // boundary integrals in the off-diagonal blocks of the
                    // system matrix. To compute the integral, we loop over
                    // all the quadrature points on the face and sum up the
                    // contribution weighted with the quadrature weights that
                    // the face quadrature rule provides.  In contrast to the
                    // entries on the diagonal blocks, here it does matter
                    // which one of the shape functions is a $\psi$ and which
                    // one is a $\phi$, since that will determine the sign of
                    // the entry.  We account for this by a simple conditional
                    // statement that determines the correct sign. Since we
                    // already checked that DoF $i$ and $j$ belong to
                    // different components, it suffices here to test for one
                    // of them to which component it belongs.
                    for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
                      cell_matrix(i,j) += ((fe.system_to_component_index(i).first == 0) ? -1 : 1) *
                                          fe_face_values.shape_value(i,q_point) *
                                          fe_face_values.shape_value(j,q_point) *
                                          c *
                                          omega *
                                          fe_face_values.JxW(q_point);
            }

        // Now we are done with this cell and have to transfer its
        // contributions from the local to the global system matrix. To this
        // end, we first get a list of the global indices of the this cells
        // DoFs...
        cell->get_dof_indices (local_dof_indices);


        // ...and then add the entries to the system matrix one by one:
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            system_matrix.add (local_dof_indices[i],
                               local_dof_indices[j],
                               cell_matrix(i,j));
      }


    // The only thing left are the Dirichlet boundary values on $\Gamma_1$,
    // which is characterized by the boundary indicator 1. The Dirichlet
    // values are provided by the <code>DirichletBoundaryValues</code> class
    // we defined above:
    std::map<types::global_dof_index,double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler,
                                              1,
                                              DirichletBoundaryValues<dim>(),
                                              boundary_values);

    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        solution,
                                        system_rhs);

    timer.stop ();
    deallog << "done ("
            << timer()
            << "s)"
            << std::endl;
  }



  // @sect4{<code>UltrasoundProblem::solve</code>}

  // As already mentioned in the introduction, the system matrix is neither
  // symmetric nor definite, and so it is not quite obvious how to come up
  // with an iterative solver and a preconditioner that do a good job on this
  // matrix.  We chose instead to go a different way and solve the linear
  // system with the sparse LU decomposition provided by UMFPACK. This is
  // often a good first choice for 2D problems and works reasonably well even
  // for a large number of DoFs.  The deal.II interface to UMFPACK is given by
  // the SparseDirectUMFPACK class, which is very easy to use and allows us to
  // solve our linear system with just 3 lines of code.

  // Note again that for compiling this example program, you need to have the
  // deal.II library built with UMFPACK support.
  template <int dim>
  void UltrasoundProblem<dim>::solve ()
  {
    deallog << "Solving linear system... ";
    Timer timer;
    timer.start ();

    // The code to solve the linear system is short: First, we allocate an
    // object of the right type. The following <code>initialize</code> call
    // provides the matrix that we would like to invert to the
    // SparseDirectUMFPACK object, and at the same time kicks off the
    // LU-decomposition. Hence, this is also the point where most of the
    // computational work in this program happens.
    SparseDirectUMFPACK  A_direct;
    A_direct.initialize(system_matrix);

    // After the decomposition, we can use <code>A_direct</code> like a matrix
    // representing the inverse of our system matrix, so to compute the
    // solution we just have to multiply with the right hand side vector:
    A_direct.vmult (solution, system_rhs);

    timer.stop ();
    deallog << "done ("
            << timer ()
            << "s)"
            << std::endl;
  }



  // @sect4{<code>UltrasoundProblem::output_results</code>}

  // Here we output our solution $v$ and $w$ as well as the derived quantity
  // $|u|$ in the format specified in the parameter file. Most of the work for
  // deriving $|u|$ from $v$ and $w$ was already done in the implementation of
  // the <code>ComputeIntensity</code> class, so that the output routine is
  // rather straightforward and very similar to what is done in the previous
  // tutorials.
  template <int dim>
  void UltrasoundProblem<dim>::output_results () const
  {
    deallog << "Generating output... ";
    Timer timer;
    timer.start ();

    // Define objects of our <code>ComputeIntensity</code> class and a DataOut
    // object:
    ComputeIntensity<dim> intensities;
    DataOut<dim> data_out;

    data_out.attach_dof_handler (dof_handler);

    // Next we query the output-related parameters from the ParameterHandler.
    // The DataOut::parse_parameters call acts as a counterpart to the
    // DataOutInterface<1>::declare_parameters call in
    // <code>ParameterReader::declare_parameters</code>. It collects all the
    // output format related parameters from the ParameterHandler and sets the
    // corresponding properties of the DataOut object accordingly.
    prm.enter_subsection("Output parameters");

    const std::string output_file    = prm.get("Output file");
    data_out.parse_parameters(prm);

    prm.leave_subsection ();

    // Now we put together the filename from the base name provided by the
    // ParameterHandler and the suffix which is provided by the DataOut class
    // (the default suffix is set to the right type that matches the one set
    // in the .prm file through parse_parameters()):
    const std::string filename = output_file +
                                 data_out.default_suffix();

    std::ofstream output (filename.c_str());

    // The solution vectors $v$ and $w$ are added to the DataOut object in the
    // usual way:
    std::vector<std::string> solution_names;
    solution_names.push_back ("Re_u");
    solution_names.push_back ("Im_u");

    data_out.add_data_vector (solution, solution_names);

    // For the intensity, we just call <code>add_data_vector</code> again, but
    // this with our <code>ComputeIntensity</code> object as the second
    // argument, which effectively adds $|u|$ to the output data:
    data_out.add_data_vector (solution, intensities);

    // The last steps are as before. Note that the actual output format is now
    // determined by what is stated in the input file, i.e. one can change the
    // output format without having to re-compile this program:
    data_out.build_patches ();
    data_out.write (output);

    timer.stop ();
    deallog << "done ("
            << timer()
            << "s)"
            << std::endl;
  }



  // @sect4{<code>UltrasoundProblem::run</code>}

  // Here we simply execute our functions one after the other:
  template <int dim>
  void UltrasoundProblem<dim>::run ()
  {
    make_grid ();
    setup_system ();
    assemble_system ();
    solve ();
    output_results ();
  }
}


// @sect4{The <code>main</code> function}

// Finally the <code>main</code> function of the program. It has the same
// structure as in almost all of the other tutorial programs. The only
// exception is that we define ParameterHandler and
// <code>ParameterReader</code> objects, and let the latter read in the
// parameter values from a textfile called <code>step-29.prm</code>. The
// values so read are then handed over to an instance of the UltrasoundProblem
// class:
int main ()
{
  try
    {
      using namespace dealii;
      using namespace Step29;

      deallog.depth_console(5);

      ParameterHandler  prm;
      ParameterReader   param(prm);
      param.read_parameters("step-29.prm");

      UltrasoundProblem<2>  ultrasound_problem (prm);
      ultrasound_problem.run ();
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

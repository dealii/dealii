/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2009 - 2013 by the deal.II authors
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
 * Author: Guido Kanschat, Texas A&M University, 2009
 */


// The first few files have already been covered in previous examples and will
// thus not be further commented on:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/mapping_q1.h>
// Here the discontinuous finite elements are defined. They are used in the
// same way as all other finite elements, though -- as you have seen in
// previous tutorial programs -- there isn't much user interaction with finite
// element classes at all: they are passed to <code>DoFHandler</code> and
// <code>FEValues</code> objects, and that is about it.
#include <deal.II/fe/fe_dgq.h>
// We are going to use the simplest possible solver, called Richardson
// iteration, that represents a simple defect correction. This, in combination
// with a block SSOR preconditioner (defined in precondition_block.h), that
// uses the special block matrix structure of system matrices arising from DG
// discretizations.
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/precondition_block.h>
// We are going to use gradients as refinement indicator.
#include <deal.II/numerics/derivative_approximation.h>

// Here come the new include files for using the MeshWorker framework. The
// first contains the class MeshWorker::DoFInfo, which provides local
// integrators with a mapping between local and global degrees of freedom. It
// stores the results of local integrals as well in its base class
// Meshworker::LocalResults.  In the second of these files, we find an object
// of type MeshWorker::IntegrationInfo, which is mostly a wrapper around a
// group of FEValues objects. The file <tt>meshworker/simple.h</tt> contains
// classes assembling locally integrated data into a global system containing
// only a single matrix. Finally, we will need the file that runs the loop
// over all mesh cells and faces.
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/simple.h>
#include <deal.II/meshworker/loop.h>

// Like in all programs, we finish this section by including the needed C++
// headers and declaring we want to use objects in the dealii namespace
// without prefix.
#include <iostream>
#include <fstream>


namespace Step12
{
  using namespace dealii;

  // @sect3{Equation data}
  //
  // First, we define a class describing the inhomogeneous boundary
  // data. Since only its values are used, we implement value_list(), but
  // leave all other functions of Function undefined.
  template <int dim>
  class BoundaryValues:  public Function<dim>
  {
  public:
    BoundaryValues () {};
    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double> &values,
                             const unsigned int component=0) const;
  };

  // Given the flow direction, the inflow boundary of the unit square
  // $[0,1]^2$ are the right and the lower boundaries. We prescribe
  // discontinuous boundary values 1 and 0 on the x-axis and value 0 on the
  // right boundary. The values of this function on the outflow boundaries
  // will not be used within the DG scheme.
  template <int dim>
  void BoundaryValues<dim>::value_list(const std::vector<Point<dim> > &points,
                                       std::vector<double> &values,
                                       const unsigned int) const
  {
    Assert(values.size()==points.size(),
           ExcDimensionMismatch(values.size(),points.size()));

    for (unsigned int i=0; i<values.size(); ++i)
      {
        if (points[i](0)<0.5)
          values[i]=1.;
        else
          values[i]=0.;
      }
  }
  // @sect3{The AdvectionProblem class}
  //
  // After this preparations, we proceed with the main class of this program,
  // called AdvectionProblem. It is basically the main class of step-6. We do
  // not have a ConstraintMatrix, because there are no hanging node
  // constraints in DG discretizations.

  // Major differences will only come up in the implementation of the assemble
  // functions, since here, we not only need to cover the flux integrals over
  // faces, we also use the MeshWorker interface to simplify the loops
  // involved.
  template <int dim>
  class AdvectionProblem
  {
  public:
    AdvectionProblem ();
    void run ();

  private:
    void setup_system ();
    void assemble_system ();
    void solve (Vector<double> &solution);
    void refine_grid ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;
    const MappingQ1<dim> mapping;

    // Furthermore we want to use DG elements of degree 1 (but this is only
    // specified in the constructor). If you want to use a DG method of a
    // different degree the whole program stays the same, only replace 1 in
    // the constructor by the desired polynomial degree.
    FE_DGQ<dim>          fe;
    DoFHandler<dim>      dof_handler;

    // The next four members represent the linear system to be
    // solved. <code>system_matrix</code> and <code>right_hand_side</code> are
    // generated by <code>assemble_system()</code>, the <code>solution</code>
    // is computed in <code>solve()</code>. The <code>sparsity_pattern</code>
    // is used to determine the location of nonzero elements in
    // <code>system_matrix</code>.
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       right_hand_side;

    // Finally, we have to provide functions that assemble the cell, boundary,
    // and inner face terms. Within the MeshWorker framework, the loop over
    // all cells and much of the setup of operations will be done outside this
    // class, so all we have to provide are these three operations. They will
    // then work on intermediate objects for which first, we here define
    // typedefs to the info objects handed to the local integration functions
    // in order to make our life easier below.
    typedef MeshWorker::DoFInfo<dim> DoFInfo;
    typedef MeshWorker::IntegrationInfo<dim> CellInfo;

    // The following three functions are then the ones that get called inside
    // the generic loop over all cells and faces. They are the ones doing the
    // actual integration.
    //
    // In our code below, these functions do not access member variables of
    // the current class, so we can mark them as <code>static</code> and
    // simply pass pointers to these functions to the MeshWorker
    // framework. If, however, these functions would want to access member
    // variables (or needed additional arguments beyond the ones specified
    // below), we could use the facilities of boost::bind (or std::bind,
    // respectively) to provide the MeshWorker framework with objects that act
    // as if they had the required number and types of arguments, but have in
    // fact other arguments already bound.
    static void integrate_cell_term (DoFInfo &dinfo,
                                     CellInfo &info);
    static void integrate_boundary_term (DoFInfo &dinfo,
                                         CellInfo &info);
    static void integrate_face_term (DoFInfo &dinfo1,
                                     DoFInfo &dinfo2,
                                     CellInfo &info1,
                                     CellInfo &info2);
  };


  // We start with the constructor. The 1 in the constructor call of
  // <code>fe</code> is the polynomial degree.
  template <int dim>
  AdvectionProblem<dim>::AdvectionProblem ()
    :
    mapping (),
    fe (1),
    dof_handler (triangulation)
  {}


  template <int dim>
  void AdvectionProblem<dim>::setup_system ()
  {
    // In the function that sets up the usual finite element data structures,
    // we first need to distribute the DoFs.
    dof_handler.distribute_dofs (fe);

    // We start by generating the sparsity pattern. To this end, we first fill
    // an intermediate object of type CompressedSparsityPattern with the
    // couplings appearing in the system. After building the pattern, this
    // object is copied to <code>sparsity_pattern</code> and can be discarded.

    // To build the sparsity pattern for DG discretizations, we can call the
    // function analogue to DoFTools::make_sparsity_pattern, which is called
    // DoFTools::make_flux_sparsity_pattern:
    CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
    DoFTools::make_flux_sparsity_pattern (dof_handler, c_sparsity);
    sparsity_pattern.copy_from(c_sparsity);

    // Finally, we set up the structure of all components of the linear
    // system.
    system_matrix.reinit (sparsity_pattern);
    solution.reinit (dof_handler.n_dofs());
    right_hand_side.reinit (dof_handler.n_dofs());
  }

  // @sect4{The assemble_system function}

  // Here we see the major difference to assembling by hand. Instead of
  // writing loops over cells and faces, we leave all this to the MeshWorker
  // framework. In order to do so, we just have to define local integration
  // functions and use one of the classes in namespace MeshWorker::Assembler
  // to build the global system.
  template <int dim>
  void AdvectionProblem<dim>::assemble_system ()
  {
    // This is the magic object, which knows everything about the data
    // structures and local integration.  This is the object doing the work in
    // the function MeshWorker::loop(), which is implicitly called by
    // MeshWorker::integration_loop() below. After the functions to which we
    // provide pointers did the local integration, the
    // MeshWorker::Assembler::SystemSimple object distributes these into the
    // global sparse matrix and the right hand side vector.
    MeshWorker::IntegrationInfoBox<dim> info_box;

    // First, we initialize the quadrature formulae and the update flags in
    // the worker base class. For quadrature, we play safe and use a QGauss
    // formula with number of points one higher than the polynomial degree
    // used. Since the quadratures for cells, boundary and interior faces can
    // be selected independently, we have to hand over this value three times.
    const unsigned int n_gauss_points = dof_handler.get_fe().degree+1;
    info_box.initialize_gauss_quadrature(n_gauss_points,
                                         n_gauss_points,
                                         n_gauss_points);

    // These are the types of values we need for integrating our system. They
    // are added to the flags used on cells, boundary and interior faces, as
    // well as interior neighbor faces, which is forced by the four @p true
    // values.
    info_box.initialize_update_flags();
    UpdateFlags update_flags = update_quadrature_points |
                               update_values            |
                               update_gradients;
    info_box.add_update_flags(update_flags, true, true, true, true);

    // After preparing all data in <tt>info_box</tt>, we initialize the
    // FEValues objects in there.
    info_box.initialize(fe, mapping);

    // The object created so far helps us do the local integration on each
    // cell and face. Now, we need an object which receives the integrated
    // (local) data and forwards them to the assembler.
    MeshWorker::DoFInfo<dim> dof_info(dof_handler);

    // Now, we have to create the assembler object and tell it, where to put
    // the local data. These will be our system matrix and the right hand
    // side.
    MeshWorker::Assembler::SystemSimple<SparseMatrix<double>, Vector<double> >
    assembler;
    assembler.initialize(system_matrix, right_hand_side);

    // Finally, the integration loop over all active cells (determined by the
    // first argument, which is an active iterator).
    //
    // As noted in the discussion when declaring the local integration
    // functions in the class declaration, the arguments expected by the
    // assembling integrator class are not actually function pointers. Rather,
    // they are objects that can be called like functions with a certain
    // number of arguments. Consequently, we could also pass objects with
    // appropriate operator() implementations here, or the result of std::bind
    // if the local integrators were, for example, non-static member
    // functions.
    MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>, MeshWorker::IntegrationInfoBox<dim> >
    (dof_handler.begin_active(), dof_handler.end(),
     dof_info, info_box,
     &AdvectionProblem<dim>::integrate_cell_term,
     &AdvectionProblem<dim>::integrate_boundary_term,
     &AdvectionProblem<dim>::integrate_face_term,
     assembler);
  }


  // @sect4{The local integrators}

  // These are the functions given to the MeshWorker::integration_loop()
  // called just above. They compute the local contributions to the system
  // matrix and right hand side on cells and faces.
  template <int dim>
  void AdvectionProblem<dim>::integrate_cell_term (DoFInfo &dinfo,
                                                   CellInfo &info)
  {
    // First, let us retrieve some of the objects used here from @p info. Note
    // that these objects can handle much more complex structures, thus the
    // access here looks more complicated than might seem necessary.
    const FEValuesBase<dim> &fe_v = info.fe_values();
    FullMatrix<double> &local_matrix = dinfo.matrix(0).matrix;
    const std::vector<double> &JxW = fe_v.get_JxW_values ();

    // With these objects, we continue local integration like always. First,
    // we loop over the quadrature points and compute the advection vector in
    // the current point.
    for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
      {
        Point<dim> beta;
        beta(0) = -fe_v.quadrature_point(point)(1);
        beta(1) = fe_v.quadrature_point(point)(0);
        beta /= beta.norm();

        // We solve a homogeneous equation, thus no right hand side shows up
        // in the cell term.  What's left is integrating the matrix entries.
        for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
          for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
            local_matrix(i,j) -= beta*fe_v.shape_grad(i,point)*
                                 fe_v.shape_value(j,point) *
                                 JxW[point];
      }
  }

  // Now the same for the boundary terms. Note that now we use FEValuesBase,
  // the base class for both FEFaceValues and FESubfaceValues, in order to get
  // access to normal vectors.
  template <int dim>
  void AdvectionProblem<dim>::integrate_boundary_term (DoFInfo &dinfo,
                                                       CellInfo &info)
  {
    const FEValuesBase<dim> &fe_v = info.fe_values();
    FullMatrix<double> &local_matrix = dinfo.matrix(0).matrix;
    Vector<double> &local_vector = dinfo.vector(0).block(0);

    const std::vector<double> &JxW = fe_v.get_JxW_values ();
    const std::vector<Point<dim> > &normals = fe_v.get_normal_vectors ();

    std::vector<double> g(fe_v.n_quadrature_points);

    static BoundaryValues<dim> boundary_function;
    boundary_function.value_list (fe_v.get_quadrature_points(), g);

    for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
      {
        Point<dim> beta;
        beta(0) = -fe_v.quadrature_point(point)(1);
        beta(1) = fe_v.quadrature_point(point)(0);
        beta /= beta.norm();

        const double beta_n=beta * normals[point];
        if (beta_n>0)
          for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
            for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
              local_matrix(i,j) += beta_n *
                                   fe_v.shape_value(j,point) *
                                   fe_v.shape_value(i,point) *
                                   JxW[point];
        else
          for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
            local_vector(i) -= beta_n *
                               g[point] *
                               fe_v.shape_value(i,point) *
                               JxW[point];
      }
  }

  // Finally, the interior face terms. The difference here is that we receive
  // two info objects, one for each cell adjacent to the face and we assemble
  // four matrices, one for each cell and two for coupling back and forth.
  template <int dim>
  void AdvectionProblem<dim>::integrate_face_term (DoFInfo &dinfo1,
                                                   DoFInfo &dinfo2,
                                                   CellInfo &info1,
                                                   CellInfo &info2)
  {
    // For quadrature points, weights, etc., we use the FEValuesBase object of
    // the first argument.
    const FEValuesBase<dim> &fe_v = info1.fe_values();

    // For additional shape functions, we have to ask the neighbors
    // FEValuesBase.
    const FEValuesBase<dim> &fe_v_neighbor = info2.fe_values();

    // Then we get references to the four local matrices. The letters u and v
    // refer to trial and test functions, respectively. The %numbers indicate
    // the cells provided by info1 and info2. By convention, the two matrices
    // in each info object refer to the test functions on the respective
    // cell. The first matrix contains the interior couplings of that cell,
    // while the second contains the couplings between cells.
    FullMatrix<double> &u1_v1_matrix = dinfo1.matrix(0,false).matrix;
    FullMatrix<double> &u2_v1_matrix = dinfo1.matrix(0,true).matrix;
    FullMatrix<double> &u1_v2_matrix = dinfo2.matrix(0,true).matrix;
    FullMatrix<double> &u2_v2_matrix = dinfo2.matrix(0,false).matrix;

    // Here, following the previous functions, we would have the local right
    // hand side vectors. Fortunately, the interface terms only involve the
    // solution and the right hand side does not receive any contributions.

    const std::vector<double> &JxW = fe_v.get_JxW_values ();
    const std::vector<Point<dim> > &normals = fe_v.get_normal_vectors ();

    for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
      {
        Point<dim> beta;
        beta(0) = -fe_v.quadrature_point(point)(1);
        beta(1) = fe_v.quadrature_point(point)(0);
        beta /= beta.norm();

        const double beta_n=beta * normals[point];
        if (beta_n>0)
          {
            // This term we've already seen:
            for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
              for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
                u1_v1_matrix(i,j) += beta_n *
                                     fe_v.shape_value(j,point) *
                                     fe_v.shape_value(i,point) *
                                     JxW[point];

            // We additionally assemble the term $(\beta\cdot n u,\hat
            // v)_{\partial \kappa_+}$,
            for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
              for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
                u1_v2_matrix(k,j) -= beta_n *
                                     fe_v.shape_value(j,point) *
                                     fe_v_neighbor.shape_value(k,point) *
                                     JxW[point];
          }
        else
          {
            // This one we've already seen, too:
            for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
              for (unsigned int l=0; l<fe_v_neighbor.dofs_per_cell; ++l)
                u2_v1_matrix(i,l) += beta_n *
                                     fe_v_neighbor.shape_value(l,point) *
                                     fe_v.shape_value(i,point) *
                                     JxW[point];

            // And this is another new one: $(\beta\cdot n \hat u,\hat
            // v)_{\partial \kappa_-}$:
            for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
              for (unsigned int l=0; l<fe_v_neighbor.dofs_per_cell; ++l)
                u2_v2_matrix(k,l) -= beta_n *
                                     fe_v_neighbor.shape_value(l,point) *
                                     fe_v_neighbor.shape_value(k,point) *
                                     JxW[point];
          }
      }
  }


  // @sect3{All the rest}
  //
  // For this simple problem we use the simplest possible solver, called
  // Richardson iteration, that represents a simple defect correction. This,
  // in combination with a block SSOR preconditioner, that uses the special
  // block matrix structure of system matrices arising from DG
  // discretizations. The size of these blocks are the number of DoFs per
  // cell. Here, we use a SSOR preconditioning as we have not renumbered the
  // DoFs according to the flow field. If the DoFs are renumbered in the
  // downstream direction of the flow, then a block Gauss-Seidel
  // preconditioner (see the PreconditionBlockSOR class with relaxation=1)
  // does a much better job.
  template <int dim>
  void AdvectionProblem<dim>::solve (Vector<double> &solution)
  {
    SolverControl           solver_control (1000, 1e-12);
    SolverRichardson<>      solver (solver_control);

    // Here we create the preconditioner,
    PreconditionBlockSSOR<SparseMatrix<double> > preconditioner;

    // then assign the matrix to it and set the right block size:
    preconditioner.initialize(system_matrix, fe.dofs_per_cell);

    // After these preparations we are ready to start the linear solver.
    solver.solve (system_matrix, solution, right_hand_side,
                  preconditioner);
  }


  // We refine the grid according to a very simple refinement criterion,
  // namely an approximation to the gradient of the solution. As here we
  // consider the DG(1) method (i.e. we use piecewise bilinear shape
  // functions) we could simply compute the gradients on each cell. But we do
  // not want to base our refinement indicator on the gradients on each cell
  // only, but want to base them also on jumps of the discontinuous solution
  // function over faces between neighboring cells. The simplest way of doing
  // that is to compute approximative gradients by difference quotients
  // including the cell under consideration and its neighbors. This is done by
  // the <code>DerivativeApproximation</code> class that computes the
  // approximate gradients in a way similar to the
  // <code>GradientEstimation</code> described in step-9 of this tutorial. In
  // fact, the <code>DerivativeApproximation</code> class was developed
  // following the <code>GradientEstimation</code> class of step-9. Relating
  // to the discussion in step-9, here we consider $h^{1+d/2}|\nabla_h
  // u_h|$. Furthermore we note that we do not consider approximate second
  // derivatives because solutions to the linear advection equation are in
  // general not in $H^2$ but in $H^1$ (to be more precise, in $H^1_\beta$)
  // only.
  template <int dim>
  void AdvectionProblem<dim>::refine_grid ()
  {
    // The <code>DerivativeApproximation</code> class computes the gradients
    // to float precision. This is sufficient as they are approximate and
    // serve as refinement indicators only.
    Vector<float> gradient_indicator (triangulation.n_active_cells());

    // Now the approximate gradients are computed
    DerivativeApproximation::approximate_gradient (mapping,
                                                   dof_handler,
                                                   solution,
                                                   gradient_indicator);

    // and they are cell-wise scaled by the factor $h^{1+d/2}$
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (unsigned int cell_no=0; cell!=endc; ++cell, ++cell_no)
      gradient_indicator(cell_no)*=std::pow(cell->diameter(), 1+1.0*dim/2);

    // Finally they serve as refinement indicator.
    GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                     gradient_indicator,
                                                     0.3, 0.1);

    triangulation.execute_coarsening_and_refinement ();
  }


  // The output of this program consists of eps-files of the adaptively
  // refined grids and the numerical solutions given in gnuplot format. This
  // was covered in previous examples and will not be further commented on.
  template <int dim>
  void AdvectionProblem<dim>::output_results (const unsigned int cycle) const
  {
    // Write the grid in eps format.
    std::string filename = "grid-";
    filename += ('0' + cycle);
    Assert (cycle < 10, ExcInternalError());

    filename += ".eps";
    deallog << "Writing grid to <" << filename << ">" << std::endl;
    std::ofstream eps_output (filename.c_str());

    GridOut grid_out;
    grid_out.write_eps (triangulation, eps_output);

    // Output of the solution in gnuplot format.
    filename = "sol-";
    filename += ('0' + cycle);
    Assert (cycle < 10, ExcInternalError());

    filename += ".gnuplot";
    deallog << "Writing solution to <" << filename << ">" << std::endl;
    std::ofstream gnuplot_output (filename.c_str());

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "u");

    data_out.build_patches ();

    data_out.write_gnuplot(gnuplot_output);
  }


  // The following <code>run</code> function is similar to previous examples.
  template <int dim>
  void AdvectionProblem<dim>::run ()
  {
    for (unsigned int cycle=0; cycle<6; ++cycle)
      {
        deallog << "Cycle " << cycle << std::endl;

        if (cycle == 0)
          {
            GridGenerator::hyper_cube (triangulation);

            triangulation.refine_global (3);
          }
        else
          refine_grid ();


        deallog << "Number of active cells:       "
                << triangulation.n_active_cells()
                << std::endl;

        setup_system ();

        deallog << "Number of degrees of freedom: "
                << dof_handler.n_dofs()
                << std::endl;

        assemble_system ();
        solve (solution);

        output_results (cycle);
      }
  }
}


// The following <code>main</code> function is similar to previous examples as
// well, and need not be commented on.
int main ()
{
  try
    {
      Step12::AdvectionProblem<2> dgmethod;
      dgmethod.run ();
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

/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2009 - 2015 by the deal.II authors
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
 * Author: Luca Heltai, Cataldo Manigrasso, 2009
 */


// @sect3{Include files}

// The program starts with including a bunch of include files that we will use
// in the various parts of the program. Most of them have been discussed in
// previous tutorials already:
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_selector.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

// And here are a few C++ standard header files that we will need:
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

// The last part of this preamble is to import everything in the dealii
// namespace into the one into which everything in this program will go:
namespace Step34
{
  using namespace dealii;


  // @sect3{Single and double layer operator kernels}

  // First, let us define a bit of the boundary integral equation machinery.

  // The following two functions are the actual calculations of the single and
  // double layer potential kernels, that is $G$ and $\nabla G$. They are well
  // defined only if the vector $R = \mathbf{y}-\mathbf{x}$ is different from
  // zero.
  namespace LaplaceKernel
  {
    template <int dim>
    double single_layer(const Tensor<1,dim> &R)
    {
      switch (dim)
        {
        case 2:
          return (-std::log(R.norm()) / (2*numbers::PI) );

        case 3:
          return (1./( R.norm()*4*numbers::PI ) );

        default:
          Assert(false, ExcInternalError());
          return 0.;
        }
    }



    template <int dim>
    Tensor<1,dim> double_layer(const Tensor<1,dim> &R)
    {
      switch (dim)
        {
        case 2:
          return R / ( -2*numbers::PI * R.norm_square());
        case 3:
          return R / ( -4*numbers::PI * R.norm_square() * R.norm() );

        default:
          Assert(false, ExcInternalError());
          return Tensor<1,dim>();
        }
    }
  }


  // @sect3{The BEMProblem class}

  // The structure of a boundary element method code is very similar to the
  // structure of a finite element code, and so the member functions of this
  // class are like those of most of the other tutorial programs. In
  // particular, by now you should be familiar with reading parameters from an
  // external file, and with the splitting of the different tasks into
  // different modules. The same applies to boundary element methods, and we
  // won't comment too much on them, except on the differences.
  template <int dim>
  class BEMProblem
  {
  public:
    BEMProblem(const unsigned int fe_degree = 1,
               const unsigned int mapping_degree = 1);

    void run();

  private:

    void read_parameters (const std::string &filename);

    void read_domain();

    void refine_and_resize();

    // The only really different function that we find here is the assembly
    // routine. We wrote this function in the most possible general way, in
    // order to allow for easy generalization to higher order methods and to
    // different fundamental solutions (e.g., Stokes or Maxwell).
    //
    // The most noticeable difference is the fact that the final matrix is
    // full, and that we have a nested loop inside the usual loop on cells
    // that visits all support points of the degrees of freedom.  Moreover,
    // when the support point lies inside the cell which we are visiting, then
    // the integral we perform becomes singular.
    //
    // The practical consequence is that we have two sets of quadrature
    // formulas, finite element values and temporary storage, one for standard
    // integration and one for the singular integration, which are used where
    // necessary.
    void assemble_system();

    // There are two options for the solution of this problem. The first is to
    // use a direct solver, and the second is to use an iterative solver. We
    // opt for the second option.
    //
    // The matrix that we assemble is not symmetric, and we opt to use the
    // GMRES method; however the construction of an efficient preconditioner
    // for boundary element methods is not a trivial issue. Here we use a non
    // preconditioned GMRES solver. The options for the iterative solver, such
    // as the tolerance, the maximum number of iterations, are selected
    // through the parameter file.
    void solve_system();

    // Once we obtained the solution, we compute the $L^2$ error of the
    // computed potential as well as the $L^\infty$ error of the approximation
    // of the solid angle. The mesh we are using is an approximation of a
    // smooth curve, therefore the computed diagonal matrix of fraction of
    // angles or solid angles $\alpha(\mathbf{x})$ should be constantly equal
    // to $\frac 12$. In this routine we output the error on the potential and
    // the error in the approximation of the computed angle. Notice that the
    // latter error is actually not the error in the computation of the angle,
    // but a measure of how well we are approximating the sphere and the
    // circle.
    //
    // Experimenting a little with the computation of the angles gives very
    // accurate results for simpler geometries. To verify this you can comment
    // out, in the read_domain() method, the tria.set_manifold(1, manifold)
    // line, and check the alpha that is generated by the program. By removing
    // this call, whenever the mesh is refined new nodes will be placed along
    // the straight lines that made up the coarse mesh, rather than be pulled
    // onto the surface that we really want to approximate. In the three
    // dimensional case, the coarse grid of the sphere is obtained starting
    // from a cube, and the obtained values of alphas are exactly $\frac 12$
    // on the nodes of the faces, $\frac 34$ on the nodes of the edges and
    // $\frac 78$ on the 8 nodes of the vertices.
    void compute_errors(const unsigned int cycle);

    // Once we obtained a solution on the codimension one domain, we want to
    // interpolate it to the rest of the space. This is done by performing
    // again the convolution of the solution with the kernel in the
    // compute_exterior_solution() function.
    //
    // We would like to plot the velocity variable which is the gradient of
    // the potential solution. The potential solution is only known on the
    // boundary, but we use the convolution with the fundamental solution to
    // interpolate it on a standard dim dimensional continuous finite element
    // space. The plot of the gradient of the extrapolated solution will give
    // us the velocity we want.
    //
    // In addition to the solution on the exterior domain, we also output the
    // solution on the domain's boundary in the output_results() function, of
    // course.
    void compute_exterior_solution();

    void output_results(const unsigned int cycle);

    // To allow for dimension independent programming, we specialize this
    // single function to extract the singular quadrature formula needed to
    // integrate the singular kernels in the interior of the cells.
    const Quadrature<dim-1> & get_singular_quadrature(
      const typename DoFHandler<dim-1, dim>::active_cell_iterator &cell,
      const unsigned int index) const;


    // The usual deal.II classes can be used for boundary element methods by
    // specifying the "codimension" of the problem. This is done by setting
    // the optional second template arguments to Triangulation, FiniteElement
    // and DoFHandler to the dimension of the embedding space. In our case we
    // generate either 1 or 2 dimensional meshes embedded in 2 or 3
    // dimensional spaces.
    //
    // The optional argument by default is equal to the first argument, and
    // produces the usual finite element classes that we saw in all previous
    // examples.
    //
    // The class is constructed in a way to allow for arbitrary order of
    // approximation of both the domain (through high order mapping) and the
    // finite element space. The order of the finite element space and of the
    // mapping can be selected in the constructor of the class.

    Triangulation<dim-1, dim>   tria;
    FE_Q<dim-1,dim>             fe;
    DoFHandler<dim-1,dim>       dh;
    MappingQ<dim-1, dim>      mapping;

    // In BEM methods, the matrix that is generated is dense. Depending on the
    // size of the problem, the final system might be solved by direct LU
    // decomposition, or by iterative methods. In this example we use an
    // unpreconditioned GMRES method. Building a preconditioner for BEM method
    // is non trivial, and we don't treat this subject here.

    FullMatrix<double>    system_matrix;
    Vector<double>        system_rhs;

    // The next two variables will denote the solution $\phi$ as well as a
    // vector that will hold the values of $\alpha(\mathbf x)$ (the fraction
    // of $\Omega$ visible from a point $\mathbf x$) at the support points of
    // our shape functions.

    Vector<double>              phi;
    Vector<double>              alpha;

    // The convergence table is used to output errors in the exact solution
    // and in the computed alphas.

    ConvergenceTable  convergence_table;

    // The following variables are the ones that we fill through a parameter
    // file.  The new objects that we use in this example are the
    // Functions::ParsedFunction object and the QuadratureSelector object.
    //
    // The Functions::ParsedFunction class allows us to easily and quickly
    // define new function objects via parameter files, with custom
    // definitions which can be very complex (see the documentation of that
    // class for all the available options).
    //
    // We will allocate the quadrature object using the QuadratureSelector
    // class that allows us to generate quadrature formulas based on an
    // identifying string and on the possible degree of the formula itself. We
    // used this to allow custom selection of the quadrature formulas for the
    // standard integration, and to define the order of the singular
    // quadrature rule.
    //
    // We also define a couple of parameters which are used in case we wanted
    // to extend the solution to the entire domain.

    Functions::ParsedFunction<dim> wind;
    Functions::ParsedFunction<dim> exact_solution;

    unsigned int singular_quadrature_order;
    std_cxx11::shared_ptr<Quadrature<dim-1> > quadrature;

    SolverControl solver_control;

    unsigned int n_cycles;
    unsigned int external_refinement;

    bool run_in_this_dimension;
    bool extend_solution;
  };


  // @sect4{BEMProblem::BEMProblem and BEMProblem::read_parameters}

  // The constructor initializes the various object in much the same way as
  // done in the finite element programs such as step-4 or step-6. The only
  // new ingredient here is the ParsedFunction object, which needs, at
  // construction time, the specification of the number of components.
  //
  // For the exact solution the number of vector components is one, and no
  // action is required since one is the default value for a ParsedFunction
  // object. The wind, however, requires dim components to be
  // specified. Notice that when declaring entries in a parameter file for the
  // expression of the Functions::ParsedFunction, we need to specify the
  // number of components explicitly, since the function
  // Functions::ParsedFunction::declare_parameters is static, and has no
  // knowledge of the number of components.
  template <int dim>
  BEMProblem<dim>::BEMProblem(const unsigned int fe_degree,
                              const unsigned int mapping_degree)
    :
    fe(fe_degree),
    dh(tria),
    mapping(mapping_degree, true),
    wind(dim)
  {}


  template <int dim>
  void BEMProblem<dim>::read_parameters (const std::string &filename)
  {
    deallog << std::endl << "Parsing parameter file " << filename << std::endl
            << "for a " << dim << " dimensional simulation. " << std::endl;

    ParameterHandler prm;

    prm.declare_entry("Number of cycles", "4",
                      Patterns::Integer());
    prm.declare_entry("External refinement", "5",
                      Patterns::Integer());
    prm.declare_entry("Extend solution on the -2,2 box", "true",
                      Patterns::Bool());
    prm.declare_entry("Run 2d simulation", "true",
                      Patterns::Bool());
    prm.declare_entry("Run 3d simulation", "true",
                      Patterns::Bool());

    prm.enter_subsection("Quadrature rules");
    {
      prm.declare_entry("Quadrature type", "gauss",
                        Patterns::Selection(QuadratureSelector<(dim-1)>::get_quadrature_names()));
      prm.declare_entry("Quadrature order", "4", Patterns::Integer());
      prm.declare_entry("Singular quadrature order", "5", Patterns::Integer());
    }
    prm.leave_subsection();

    // For both two and three dimensions, we set the default input data to be
    // such that the solution is $x+y$ or $x+y+z$. The actually computed
    // solution will have value zero at infinity. In this case, this coincide
    // with the exact solution, and no additional corrections are needed, but
    // you should be aware of the fact that we arbitrarily set $\phi_\infty$,
    // and the exact solution we pass to the program needs to have the same
    // value at infinity for the error to be computed correctly.
    //
    // The use of the Functions::ParsedFunction object is pretty straight
    // forward. The Functions::ParsedFunction::declare_parameters function
    // takes an additional integer argument that specifies the number of
    // components of the given function. Its default value is one. When the
    // corresponding Functions::ParsedFunction::parse_parameters method is
    // called, the calling object has to have the same number of components
    // defined here, otherwise an exception is thrown.
    //
    // When declaring entries, we declare both 2 and three dimensional
    // functions. However only the dim-dimensional one is ultimately
    // parsed. This allows us to have only one parameter file for both 2 and 3
    // dimensional problems.
    //
    // Notice that from a mathematical point of view, the wind function on the
    // boundary should satisfy the condition $\int_{\partial\Omega}
    // \mathbf{v}\cdot \mathbf{n} d \Gamma = 0$, for the problem to have a
    // solution. If this condition is not satisfied, then no solution can be
    // found, and the solver will not converge.
    prm.enter_subsection("Wind function 2d");
    {
      Functions::ParsedFunction<2>::declare_parameters(prm, 2);
      prm.set("Function expression", "1; 1");
    }
    prm.leave_subsection();

    prm.enter_subsection("Wind function 3d");
    {
      Functions::ParsedFunction<3>::declare_parameters(prm, 3);
      prm.set("Function expression", "1; 1; 1");
    }
    prm.leave_subsection();

    prm.enter_subsection("Exact solution 2d");
    {
      Functions::ParsedFunction<2>::declare_parameters(prm);
      prm.set("Function expression", "x+y");
    }
    prm.leave_subsection();

    prm.enter_subsection("Exact solution 3d");
    {
      Functions::ParsedFunction<3>::declare_parameters(prm);
      prm.set("Function expression", "x+y+z");
    }
    prm.leave_subsection();


    // In the solver section, we set all SolverControl parameters. The object
    // will then be fed to the GMRES solver in the solve_system() function.
    prm.enter_subsection("Solver");
    SolverControl::declare_parameters(prm);
    prm.leave_subsection();

    // After declaring all these parameters to the ParameterHandler object,
    // let's read an input file that will give the parameters their values. We
    // then proceed to extract these values from the ParameterHandler object:
    prm.read_input(filename);

    n_cycles = prm.get_integer("Number of cycles");
    external_refinement = prm.get_integer("External refinement");
    extend_solution = prm.get_bool("Extend solution on the -2,2 box");

    prm.enter_subsection("Quadrature rules");
    {
      quadrature =
        std_cxx11::shared_ptr<Quadrature<dim-1> >
        (new QuadratureSelector<dim-1> (prm.get("Quadrature type"),
                                        prm.get_integer("Quadrature order")));
      singular_quadrature_order = prm.get_integer("Singular quadrature order");
    }
    prm.leave_subsection();

    prm.enter_subsection(std::string("Wind function ")+
                         Utilities::int_to_string(dim)+std::string("d"));
    {
      wind.parse_parameters(prm);
    }
    prm.leave_subsection();

    prm.enter_subsection(std::string("Exact solution ")+
                         Utilities::int_to_string(dim)+std::string("d"));
    {
      exact_solution.parse_parameters(prm);
    }
    prm.leave_subsection();

    prm.enter_subsection("Solver");
    solver_control.parse_parameters(prm);
    prm.leave_subsection();


    // Finally, here's another example of how to use parameter files in
    // dimension independent programming.  If we wanted to switch off one of
    // the two simulations, we could do this by setting the corresponding "Run
    // 2d simulation" or "Run 3d simulation" flag to false:
    run_in_this_dimension = prm.get_bool("Run " +
                                         Utilities::int_to_string(dim) +
                                         "d simulation");
  }


  // @sect4{BEMProblem::read_domain}

  // A boundary element method triangulation is basically the same as a
  // (dim-1) dimensional triangulation, with the difference that the vertices
  // belong to a (dim) dimensional space.
  //
  // Some of the mesh formats supported in deal.II use by default three
  // dimensional points to describe meshes. These are the formats which are
  // compatible with the boundary element method capabilities of deal.II. In
  // particular we can use either UCD or GMSH formats. In both cases, we have
  // to be particularly careful with the orientation of the mesh, because,
  // unlike in the standard finite element case, no reordering or
  // compatibility check is performed here.  All meshes are considered as
  // oriented, because they are embedded in a higher dimensional space. (See
  // the documentation of the GridIn and of the Triangulation for further
  // details on orientation of cells in a triangulation.) In our case, the
  // normals to the mesh are external to both the circle in 2d or the sphere
  // in 3d.
  //
  // The other detail that is required for appropriate refinement of
  // the boundary element mesh, is an accurate description of the
  // manifold that the mesh is approximating. We already saw this
  // several times for the boundary of standard finite element meshes
  // (for example in step-5 and step-6), and here the principle and
  // usage is the same, except that the SphericalManifold class takes
  // an additional template parameter that specifies the embedding
  // space dimension. The function object still has to be static to
  // live at least as long as the triangulation object to which it is
  // attached.

  template <int dim>
  void BEMProblem<dim>::read_domain()
  {
    static const Point<dim> center = Point<dim>();
    static const SphericalManifold<dim-1, dim> manifold(center);

    std::ifstream in;
    switch (dim)
      {
      case 2:
        in.open ("coarse_circle.inp");
        break;

      case 3:
        in.open ("coarse_sphere.inp");
        break;

      default:
        Assert (false, ExcNotImplemented());
      }

    GridIn<dim-1, dim> gi;
    gi.attach_triangulation (tria);
    gi.read_ucd (in);

    tria.set_all_manifold_ids(1);
    tria.set_manifold(1, manifold);
  }


  // @sect4{BEMProblem::refine_and_resize}

  // This function globally refines the mesh, distributes degrees of freedom,
  // and resizes matrices and vectors.

  template <int dim>
  void BEMProblem<dim>::refine_and_resize()
  {
    tria.refine_global(1);

    dh.distribute_dofs(fe);

    const unsigned int n_dofs =  dh.n_dofs();

    system_matrix.reinit(n_dofs, n_dofs);

    system_rhs.reinit(n_dofs);
    phi.reinit(n_dofs);
    alpha.reinit(n_dofs);
  }


  // @sect4{BEMProblem::assemble_system}

  // The following is the main function of this program, assembling the matrix
  // that corresponds to the boundary integral equation.
  template <int dim>
  void BEMProblem<dim>::assemble_system()
  {

    // First we initialize an FEValues object with the quadrature formula for
    // the integration of the kernel in non singular cells. This quadrature is
    // selected with the parameter file, and needs to be quite precise, since
    // the functions we are integrating are not polynomial functions.
    FEValues<dim-1,dim> fe_v(mapping, fe, *quadrature,
                             update_values |
                             update_cell_normal_vectors |
                             update_quadrature_points |
                             update_JxW_values);

    const unsigned int n_q_points = fe_v.n_quadrature_points;

    std::vector<types::global_dof_index> local_dof_indices(fe.dofs_per_cell);

    std::vector<Vector<double> > cell_wind(n_q_points, Vector<double>(dim) );
    double normal_wind;

    // Unlike in finite element methods, if we use a collocation boundary
    // element method, then in each assembly loop we only assemble the
    // information that refers to the coupling between one degree of freedom
    // (the degree associated with support point $i$) and the current
    // cell. This is done using a vector of fe.dofs_per_cell elements, which
    // will then be distributed to the matrix in the global row $i$. The
    // following object will hold this information:
    Vector<double>      local_matrix_row_i(fe.dofs_per_cell);

    // The index $i$ runs on the collocation points, which are the support
    // points of the $i$th basis function, while $j$ runs on inner integration
    // points.

    // We construct a vector of support points which will be used in the local
    // integrations:
    std::vector<Point<dim> > support_points(dh.n_dofs());
    DoFTools::map_dofs_to_support_points<dim-1, dim>( mapping, dh, support_points);


    // After doing so, we can start the integration loop over all cells, where
    // we first initialize the FEValues object and get the values of
    // $\mathbf{\tilde v}$ at the quadrature points (this vector field should
    // be constant, but it doesn't hurt to be more general):
    typename DoFHandler<dim-1,dim>::active_cell_iterator
    cell = dh.begin_active(),
    endc = dh.end();

    for (cell = dh.begin_active(); cell != endc; ++cell)
      {
        fe_v.reinit(cell);
        cell->get_dof_indices(local_dof_indices);

        const std::vector<Point<dim> >    &q_points = fe_v.get_quadrature_points();
        const std::vector<Tensor<1,dim> > &normals  = fe_v.get_all_normal_vectors();
        wind.vector_value_list(q_points, cell_wind);

        // We then form the integral over the current cell for all degrees of
        // freedom (note that this includes degrees of freedom not located on
        // the current cell, a deviation from the usual finite element
        // integrals). The integral that we need to perform is singular if one
        // of the local degrees of freedom is the same as the support point
        // $i$. A the beginning of the loop we therefore check whether this is
        // the case, and we store which one is the singular index:
        for (unsigned int i=0; i<dh.n_dofs() ; ++i)
          {

            local_matrix_row_i = 0;

            bool is_singular = false;
            unsigned int singular_index = numbers::invalid_unsigned_int;

            for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
              if (local_dof_indices[j] == i)
                {
                  singular_index = j;
                  is_singular = true;
                  break;
                }

            // We then perform the integral. If the index $i$ is not one of
            // the local degrees of freedom, we simply have to add the single
            // layer terms to the right hand side, and the double layer terms
            // to the matrix:
            if (is_singular == false)
              {
                for (unsigned int q=0; q<n_q_points; ++q)
                  {
                    normal_wind = 0;
                    for (unsigned int d=0; d<dim; ++d)
                      normal_wind += normals[q][d]*cell_wind[q](d);

                    const Tensor<1,dim> R = q_points[q] - support_points[i];

                    system_rhs(i) += ( LaplaceKernel::single_layer(R)   *
                                       normal_wind                      *
                                       fe_v.JxW(q) );

                    for (unsigned int j=0; j<fe.dofs_per_cell; ++j)

                      local_matrix_row_i(j) -= ( ( LaplaceKernel::double_layer(R)     *
                                                   normals[q] )            *
                                                 fe_v.shape_value(j,q)     *
                                                 fe_v.JxW(q)       );
                  }
              }
            else
              {
                // Now we treat the more delicate case. If we are here, this
                // means that the cell that runs on the $j$ index contains
                // support_point[i]. In this case both the single and the
                // double layer potential are singular, and they require
                // special treatment.
                //
                // Whenever the integration is performed with the singularity
                // inside the given cell, then a special quadrature formula is
                // used that allows one to integrate arbitrary functions
                // against a singular weight on the reference cell.
                //
                // The correct quadrature formula is selected by the
                // get_singular_quadrature function, which is explained in
                // detail below.
                Assert(singular_index != numbers::invalid_unsigned_int,
                       ExcInternalError());

                const Quadrature<dim-1> & singular_quadrature =
                  get_singular_quadrature(cell, singular_index);

                FEValues<dim-1,dim> fe_v_singular (mapping, fe, singular_quadrature,
                                                   update_jacobians |
                                                   update_values |
                                                   update_cell_normal_vectors |
                                                   update_quadrature_points );

                fe_v_singular.reinit(cell);

                std::vector<Vector<double> > singular_cell_wind( singular_quadrature.size(),
                                                                 Vector<double>(dim) );

                const std::vector<Tensor<1,dim> > &singular_normals  = fe_v_singular.get_all_normal_vectors();
                const std::vector<Point<dim> >    &singular_q_points = fe_v_singular.get_quadrature_points();

                wind.vector_value_list(singular_q_points, singular_cell_wind);

                for (unsigned int q=0; q<singular_quadrature.size(); ++q)
                  {
                    const Tensor<1,dim> R = singular_q_points[q] - support_points[i];
                    double normal_wind = 0;
                    for (unsigned int d=0; d<dim; ++d)
                      normal_wind += (singular_cell_wind[q](d)*
                                      singular_normals[q][d]);

                    system_rhs(i) += ( LaplaceKernel::single_layer(R) *
                                       normal_wind                         *
                                       fe_v_singular.JxW(q) );

                    for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
                      {
                        local_matrix_row_i(j) -= (( LaplaceKernel::double_layer(R) *
                                                    singular_normals[q])                *
                                                  fe_v_singular.shape_value(j,q)        *
                                                  fe_v_singular.JxW(q)       );
                      }
                  }
              }

            // Finally, we need to add the contributions of the current cell
            // to the global matrix.
            for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
              system_matrix(i,local_dof_indices[j])
              += local_matrix_row_i(j);
          }
      }

    // The second part of the integral operator is the term
    // $\alpha(\mathbf{x}_i) \phi_j(\mathbf{x}_i)$. Since we use a collocation
    // scheme, $\phi_j(\mathbf{x}_i)=\delta_{ij}$ and the corresponding matrix
    // is a diagonal one with entries equal to $\alpha(\mathbf{x}_i)$.

    // One quick way to compute this diagonal matrix of the solid angles, is
    // to use the Neumann matrix itself. It is enough to multiply the matrix
    // with a vector of elements all equal to -1, to get the diagonal matrix
    // of the alpha angles, or solid angles (see the formula in the
    // introduction for this). The result is then added back onto the system
    // matrix object to yield the final form of the matrix:
    Vector<double> ones(dh.n_dofs());
    ones.add(-1.);

    system_matrix.vmult(alpha, ones);
    alpha.add(1);
    for (unsigned int i = 0; i<dh.n_dofs(); ++i)
      system_matrix(i,i) +=  alpha(i);
  }


  // @sect4{BEMProblem::solve_system}

  // The next function simply solves the linear system.
  template <int dim>
  void BEMProblem<dim>::solve_system()
  {
    SolverGMRES<Vector<double> > solver (solver_control);
    solver.solve (system_matrix, phi, system_rhs, PreconditionIdentity());
  }


  // @sect4{BEMProblem::compute_errors}

  // The computation of the errors is exactly the same in all other example
  // programs, and we won't comment too much. Notice how the same methods that
  // are used in the finite element methods can be used here.
  template <int dim>
  void BEMProblem<dim>::compute_errors(const unsigned int cycle)
  {
    Vector<float> difference_per_cell (tria.n_active_cells());
    VectorTools::integrate_difference (mapping, dh, phi,
                                       exact_solution,
                                       difference_per_cell,
                                       QGauss<(dim-1)>(2*fe.degree+1),
                                       VectorTools::L2_norm);
    const double L2_error = VectorTools::compute_global_error(tria,
                                                              difference_per_cell,
                                                              VectorTools::L2_norm);

    // The error in the alpha vector can be computed directly using the
    // Vector::linfty_norm() function, since on each node, the value should be
    // $\frac 12$. All errors are then output and appended to our
    // ConvergenceTable object for later computation of convergence rates:
    Vector<double> difference_per_node(alpha);
    difference_per_node.add(-.5);

    const double alpha_error = difference_per_node.linfty_norm();
    const unsigned int n_active_cells=tria.n_active_cells();
    const unsigned int n_dofs=dh.n_dofs();

    deallog << "Cycle " << cycle << ':'
            << std::endl
            << "   Number of active cells:       "
            << n_active_cells
            << std::endl
            << "   Number of degrees of freedom: "
            << n_dofs
            << std::endl;

    convergence_table.add_value("cycle", cycle);
    convergence_table.add_value("cells", n_active_cells);
    convergence_table.add_value("dofs", n_dofs);
    convergence_table.add_value("L2(phi)", L2_error);
    convergence_table.add_value("Linfty(alpha)", alpha_error);
  }


  // Singular integration requires a careful selection of the quadrature
  // rules. In particular the deal.II library provides quadrature rules which
  // are tailored for logarithmic singularities (QGaussLog, QGaussLogR), as
  // well as for 1/R singularities (QGaussOneOverR).
  //
  // Singular integration is typically obtained by constructing weighted
  // quadrature formulas with singular weights, so that it is possible to
  // write
  //
  // \f[ \int_K f(x) s(x) dx = \sum_{i=1}^N w_i f(q_i) \f]
  //
  // where $s(x)$ is a given singularity, and the weights and quadrature
  // points $w_i,q_i$ are carefully selected to make the formula above an
  // equality for a certain class of functions $f(x)$.
  //
  // In all the finite element examples we have seen so far, the weight of the
  // quadrature itself (namely, the function $s(x)$), was always constantly
  // equal to 1.  For singular integration, we have two choices: we can use
  // the definition above, factoring out the singularity from the integrand
  // (i.e., integrating $f(x)$ with the special quadrature rule), or we can
  // ask the quadrature rule to "normalize" the weights $w_i$ with $s(q_i)$:
  //
  // \f[ \int_K f(x) s(x) dx = \int_K g(x) dx = \sum_{i=1}^N
  //   \frac{w_i}{s(q_i)} g(q_i) \f]
  //
  // We use this second option, through the @p factor_out_singularity
  // parameter of both QGaussLogR and QGaussOneOverR.
  //
  // These integrals are somewhat delicate, especially in two dimensions, due
  // to the transformation from the real to the reference cell, where the
  // variable of integration is scaled with the determinant of the
  // transformation.
  //
  // In two dimensions this process does not result only in a factor appearing
  // as a constant factor on the entire integral, but also on an additional
  // integral altogether that needs to be evaluated:
  //
  // \f[ \int_0^1 f(x)\ln(x/\alpha) dx = \int_0^1 f(x)\ln(x) dx - \int_0^1
  //  f(x) \ln(\alpha) dx.  \f]
  //
  // This process is taken care of by the constructor of the QGaussLogR class,
  // which adds additional quadrature points and weights to take into
  // consideration also the second part of the integral.
  //
  // A similar reasoning should be done in the three dimensional case, since
  // the singular quadrature is tailored on the inverse of the radius $r$ in
  // the reference cell, while our singular function lives in real space,
  // however in the three dimensional case everything is simpler because the
  // singularity scales linearly with the determinant of the
  // transformation. This allows us to build the singular two dimensional
  // quadrature rules only once and, reuse them over all cells.
  //
  // In the one dimensional singular integration this is not possible, since
  // we need to know the scaling parameter for the quadrature, which is not
  // known a priori. Here, the quadrature rule itself depends also on the size
  // of the current cell. For this reason, it is necessary to create a new
  // quadrature for each singular integration.
  //
  // The different quadrature rules are built inside the
  // get_singular_quadrature, which is specialized for dim=2 and dim=3, and
  // they are retrieved inside the assemble_system function. The index given
  // as an argument is the index of the unit support point where the
  // singularity is located.

  template<>
  const Quadrature<2> &BEMProblem<3>::get_singular_quadrature(
    const DoFHandler<2,3>::active_cell_iterator &,
    const unsigned int index) const
  {
    Assert(index < fe.dofs_per_cell,
           ExcIndexRange(0, fe.dofs_per_cell, index));

    static std::vector<QGaussOneOverR<2> > quadratures;
    if (quadratures.size() == 0)
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        quadratures.push_back(QGaussOneOverR<2>(singular_quadrature_order,
                                                fe.get_unit_support_points()[i],
                                                true));
    return quadratures[index];
  }


  template<>
  const Quadrature<1> &BEMProblem<2>::get_singular_quadrature(
    const DoFHandler<1,2>::active_cell_iterator &cell,
    const unsigned int index) const
  {
    Assert(index < fe.dofs_per_cell,
           ExcIndexRange(0, fe.dofs_per_cell, index));

    static Quadrature<1> *q_pointer = NULL;
    if (q_pointer) delete q_pointer;

    q_pointer = new QGaussLogR<1>(singular_quadrature_order,
                                  fe.get_unit_support_points()[index],
                                  1./cell->measure(), true);
    return (*q_pointer);
  }



  // @sect4{BEMProblem::compute_exterior_solution}

  // We'd like to also know something about the value of the potential $\phi$
  // in the exterior domain: after all our motivation to consider the boundary
  // integral problem was that we wanted to know the velocity in the exterior
  // domain!
  //
  // To this end, let us assume here that the boundary element domain is
  // contained in the box $[-2,2]^{\text{dim}}$, and we extrapolate the actual
  // solution inside this box using the convolution with the fundamental
  // solution. The formula for this is given in the introduction.
  //
  // The reconstruction of the solution in the entire space is done on a
  // continuous finite element grid of dimension dim. These are the usual
  // ones, and we don't comment any further on them. At the end of the
  // function, we output this exterior solution in, again, much the usual way.
  template <int dim>
  void BEMProblem<dim>::compute_exterior_solution()
  {
    Triangulation<dim>  external_tria;
    GridGenerator::hyper_cube(external_tria, -2, 2);

    FE_Q<dim>           external_fe(1);
    DoFHandler<dim>     external_dh (external_tria);
    Vector<double>      external_phi;

    external_tria.refine_global(external_refinement);
    external_dh.distribute_dofs(external_fe);
    external_phi.reinit(external_dh.n_dofs());

    typename DoFHandler<dim-1,dim>::active_cell_iterator
    cell = dh.begin_active(),
    endc = dh.end();


    FEValues<dim-1,dim> fe_v(mapping, fe, *quadrature,
                             update_values |
                             update_cell_normal_vectors |
                             update_quadrature_points |
                             update_JxW_values);

    const unsigned int n_q_points = fe_v.n_quadrature_points;

    std::vector<types::global_dof_index> dofs(fe.dofs_per_cell);

    std::vector<double> local_phi(n_q_points);
    std::vector<double> normal_wind(n_q_points);
    std::vector<Vector<double> > local_wind(n_q_points, Vector<double>(dim) );

    std::vector<Point<dim> > external_support_points(external_dh.n_dofs());
    DoFTools::map_dofs_to_support_points<dim>(StaticMappingQ1<dim>::mapping,
                                              external_dh, external_support_points);

    for (cell = dh.begin_active(); cell != endc; ++cell)
      {
        fe_v.reinit(cell);

        const std::vector<Point<dim> >    &q_points = fe_v.get_quadrature_points();
        const std::vector<Tensor<1,dim> > &normals  = fe_v.get_all_normal_vectors();

        cell->get_dof_indices(dofs);
        fe_v.get_function_values(phi, local_phi);

        wind.vector_value_list(q_points, local_wind);

        for (unsigned int q=0; q<n_q_points; ++q)
          {
            normal_wind[q] = 0;
            for (unsigned int d=0; d<dim; ++d)
              normal_wind[q] += normals[q][d]*local_wind[q](d);
          }

        for (unsigned int i=0; i<external_dh.n_dofs(); ++i)
          for (unsigned int q=0; q<n_q_points; ++q)
            {

              const Tensor<1,dim> R = q_points[q] - external_support_points[i];

              external_phi(i) += ( ( LaplaceKernel::single_layer(R) *
                                     normal_wind[q]
                                     +
                                     (LaplaceKernel::double_layer(R) *
                                      normals[q] )            *
                                     local_phi[q] )           *
                                   fe_v.JxW(q) );
            }
      }

    DataOut<dim> data_out;

    data_out.attach_dof_handler(external_dh);
    data_out.add_data_vector(external_phi, "external_phi");
    data_out.build_patches();

    const std::string
    filename = Utilities::int_to_string(dim) + "d_external.vtk";
    std::ofstream file(filename.c_str());

    data_out.write_vtk(file);
  }


  // @sect4{BEMProblem::output_results}

  // Outputting the results of our computations is a rather mechanical
  // tasks. All the components of this function have been discussed before.
  template <int dim>
  void BEMProblem<dim>::output_results(const unsigned int cycle)
  {
    DataOut<dim-1, DoFHandler<dim-1, dim> > dataout;

    dataout.attach_dof_handler(dh);
    dataout.add_data_vector(phi, "phi",
                            DataOut<dim-1, DoFHandler<dim-1, dim> >::type_dof_data);
    dataout.add_data_vector(alpha, "alpha",
                            DataOut<dim-1, DoFHandler<dim-1, dim> >::type_dof_data);
    dataout.build_patches(mapping,
                          mapping.get_degree(),
                          DataOut<dim-1, DoFHandler<dim-1, dim> >::curved_inner_cells);

    std::string filename = ( Utilities::int_to_string(dim) +
                             "d_boundary_solution_" +
                             Utilities::int_to_string(cycle) +
                             ".vtk" );
    std::ofstream file(filename.c_str());

    dataout.write_vtk(file);

    if (cycle == n_cycles-1)
      {
        convergence_table.set_precision("L2(phi)", 3);
        convergence_table.set_precision("Linfty(alpha)", 3);

        convergence_table.set_scientific("L2(phi)", true);
        convergence_table.set_scientific("Linfty(alpha)", true);

        convergence_table
        .evaluate_convergence_rates("L2(phi)", ConvergenceTable::reduction_rate_log2);
        convergence_table
        .evaluate_convergence_rates("Linfty(alpha)", ConvergenceTable::reduction_rate_log2);
        deallog << std::endl;
        convergence_table.write_text(std::cout);
      }
  }


  // @sect4{BEMProblem::run}

  // This is the main function. It should be self explanatory in its
  // briefness:
  template <int dim>
  void BEMProblem<dim>::run()
  {

    read_parameters("parameters.prm");

    if (run_in_this_dimension == false)
      {
        deallog << "Run in dimension " << dim
                << " explicitly disabled in parameter file. "
                << std::endl;
        return;
      }

    read_domain();

    for (unsigned int cycle=0; cycle<n_cycles; ++cycle)
      {
        refine_and_resize();
        assemble_system();
        solve_system();
        compute_errors(cycle);
        output_results(cycle);
      }

    if (extend_solution == true)
      compute_exterior_solution();
  }
}


// @sect3{The main() function}

// This is the main function of this program. It is exactly like all previous
// tutorial programs:
int main ()
{
  try
    {
      using namespace dealii;
      using namespace Step34;

      const unsigned int degree = 1;
      const unsigned int mapping_degree = 1;

      deallog.depth_console (3);
      BEMProblem<2> laplace_problem_2d(degree, mapping_degree);
      laplace_problem_2d.run();

      BEMProblem<3> laplace_problem_3d(degree, mapping_degree);
      laplace_problem_3d.run();
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

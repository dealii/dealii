/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2000 - 2025 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */


// @sect3{Include files}

// First we include header files that we have already used in previous
// example programs:
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/numerics/data_out.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

// This again is C++:
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>

// The last step is as in all previous programs:
namespace Step101
{
  using namespace dealii;


  // @sect3{Parameter struct for CMS}

  // We start with the parameter struct for CMS. We need to have it at the beginning
  // of the code because we call it in other parts of the code and else, we would
  // get an error message.
  struct CmsParameters
  {
    // Material parameters (Lamé parameters lambda and mu + density)
    double lambda = 1.0;
    double mu = 1.0;
    double rho = 1.0;

    // We use the parameter thickness to compute the mass for the assessment of the
    // impact of the artificially added mass to the model. We still have a 2D problem,
    // this is a non-existing thickness only used for the assessment, therefore, we
    // set it to 1.0.
    double thickness = 1.0; 

    // To select a bigger critical step, we define the parameter dt_factor.
    // The newly selected critical step, the target critical step,
    // <code>dt_target</code> will then be calculated through
    // <code>dt_target = dt_factor * dt_crit_global</code>, which is based on the selected
    // critical cells.
    double dt_factor = 1.1;

    // Acceptability criteria

    // Safety cap used in the assessment of the impact of the artificially added
    // mass to the system's response. If the mass increase of the model is greater
    // than <code>max_added_mass_fraction</code>, the target critical step
    // <code>dt_target</code> will be reduced iteratively until the mass of the model
    // equals or is smaller than <code>max_added_mass_fraction</code> or until
    // the maximum number of iterations is reached.
    double max_added_mass_fraction = 0.10;
    double max_scaled_cell_fraction = 0.10;

    bool auto_reduce_dt_target = true;

    // When the mass increase is too large, shrink <code>dt_target</code> by
    // the factor <code>dt_shrink_factor</code> and retry.
    double dt_shrink_factor = 0.95;

    // Maximum number of iterations in the loop that reduces the target critical
    // step so that the artificially added mass to the system does not have
    // a major impact on the system's response. If this number of iterations
    // is reached and the artificially added mass is bigger than
    // <code>max_added_mass_fraction</code> multiplied by the original mass,
    // there are two possible options. The first option is that you increase
    // the factor <code>max_added_mass_fraction</code>. However, you should do
    // this carefully considering your specific system. Will an increase of this
    // factor result in too big inertia forces? Then you should not do this.
    // The second option is that CMS is not suitable for your model. You may consider
    // remeshing it with a coarser mesh.
    unsigned int max_dt_adjustment_steps = 25;

    // Output a VTK field with the resulting density scaling factor.
    bool write_density_scaling_output = true;
  };

  // @sect3{The <code>ElasticProblem</code> class template}

  // The first part of the program is the declaration of the main class. All of
  // this has been copied verbatim from step-8.

  // TODO: If something is changed/added compared to step-8, comment it
  template <int dim>
  // TODO: Check if this class is needed
   class ElasticProblem
   {
   public:
     //ElasticProblem();
     ElasticProblem(const CmsParameters &parameters); // TODO: Check if this is correct
     void run();

   private:
   struct CmsResult
     {
      double dt_crit_global_before = 0.0;
      double dt_target = 0.0;
      double dt_crit_global_after = 0.0;

      double total_mass = 0.0;
      double added_mass = 0.0;
      double added_mass_ratio = 0.0;

      unsigned int n_scaled = 0;
      double scaled_fraction = 0.0;

      std::vector<double> h_char;

      // Per-cell density scaling alpha = rho_new / rho_old
      std::vector<double> density_scaling;
     };

     //void setup_system();
     // void assemble_system();
     // void solve();
     // void refine_grid();
     // void output_results(const unsigned int cycle) const;
     //void assemble_stiffness_and_rhs();
     //void assemble_consistent_mass_matrix();

     // Mesh and geometry
     void make_conforming_graded_mesh();

     static double compute_min_edge_length_of_cell(
      const typename Triangulation<dim>::active_cell_iterator &cell
     );

     void compute_characteristic_length(std::vector<double> &h_char) const;

     // CFL / dt computations
     double compute_dt_cell(const double h, const double rho_local) const;
     double compute_dt_crit_global(const std::vector<double> &h_char, const std::vector<double> &density_scaling) const;

     // CMS
     CmsResult apply_cms_for_target_dt(const std::vector<double> &h_char, const double dt_target) const;

     CmsResult choose_dt_and_apply_cms(const std::vector<double> &h_char) const;

     void verify_mesh_is_quad_and_non_degenerate() const;

     void print_cell_debug_table(const std::vector<double> &h_char, const std::vector<double> &density_scaling, const double dt_target) const;

     // Output
     void output_density_scaling(const std::vector<double> &density_scaling) const;

     const CmsParameters parameters;

     Triangulation<dim> triangulation;
     DoFHandler<dim> dof_handler;
     const FESystem<dim> fe;

  //    const CmsParameters parameters;

  //    Triangulation<dim> triangulation{Triangulation<dim>::limit_level_difference_at_vertices};
  //    DoFHandler<dim>    dof_handler;

  //    const FESystem<dim> fe;

  //    AffineConstraints<double> constraints;

  //    SparsityPattern      sparsity_pattern;
  //    SparseMatrix<double> system_matrix;
  //    SparseMatrix<double> stiffness_matrix;
  //    SparseMatrix<double> mass_matrix;

  // //    Vector<double> solution;
  //    Vector<double> system_rhs;
  //    static double compute_min_edge_length_of_cell(
  //     const typename DoFHandler<dim>::active_cell_iterator &cell
  //    );

  //    std::pair<double, std::vector<unsigned int>>
  //    select_critical_cells_by_percentile(const std::vector<double> &h_char) const;

  //    double compute_pressure_wave_speed() const;

  //    CmsResult apply_classical_mass_scaling(
  //     const std::vector<double> &h_char, const std::vector<unsigned int> &critical_cell_ids, const double dt_target
  //    ) const;

  //    CmsResult choose_dt_and_apply_cms(const std::vector<double> &h_char) const;

  //    void output_density_scaling(const std::vector<double> &density_scaling) const;


    };

    template <int dim>
    ElasticProblem<dim>::ElasticProblem(const CmsParameters &parameters)
      : parameters(parameters)
      , dof_handler(triangulation)
      , fe(FE_Q<dim>(1) ^ dim)
    {}

    // Mesh generation without hanging nodes
    // Build a structured quad mesh on [-1,1]² but non-uniform spacing in x/y (clustering).
    // This yields different cell sizes AND is conforming (no hanging modes / no T-junctions).

    template <int dim>
    void ElasticProblem<dim>::make_conforming_graded_mesh()
    {
      Assert(dim == 2, ExcMessage("This demo is implemented for dim=2 only."));

      // Coordinates with clustering in several regions:
      //  - coarse on the left/bottom
      //  - finer near x ~ 0.6 and y ~ 0.6
      //  - also some refinement near x ~ -0.4, y ~ 0.2 (second region)
      const std::vector<double> x = {
        -1.0, -0.7, -0.4, -0.2, 0.0, 0.2, 0.4,
        0.52, 0.58, 0.62, 0.8, 0.9, 1.0
      };

      const std::vector<double> y = {
        -1.0, -0.7, -0.4, -0.2, 0.0, 0.15, 0.22, 0.3,
        0.45, 0.55, 0.6, 0.63, 0.7, 0.85, 1.0
      };

      std::vector<Point<2>> vertices;
      vertices.reserve(x.size() * y.size());

      auto vid = [&](unsigned int i, unsigned int j){
        return j * x.size() + i;
      };

      for (unsigned int j = 0; j < y.size(); ++j)
        for (unsigned int i = 0; i < x.size(); ++i)
          vertices.emplace_back(x[i],y[j]);

      std::vector<CellData<2>> cells;
      cells.reserve((x.size()-1) * (y.size()-1));

      for (unsigned int j = 0; j < y.size()-1; ++j)
        for (unsigned int i = 0; i < x.size()-1; ++i)
          {
            CellData<2> c;
            c.vertices[0] = vid(i,j); // bottom-left
            c.vertices[1] = vid(i+1,j); // bottom-right
            c.vertices[2] = vid(i,j+1); // top-left
            c.vertices[3] = vid(i+1,j+1); //top-right
            c.material_id = 0;
            cells.push_back(c);
          }

      SubCellData subcell_data;
      triangulation.create_triangulation(vertices,cells,subcell_data);

      // No global / local refine here -> mesh stays conforming.
      // No hanging nodes can appear because we never do adaptive refinement

    }


  // @sect3{Right hand side values}

  // TODO: DELETE THIS LARGE COMMENT, IT'S FROM STEP-8, JUST BRIEFLY EXPLAIN
  // Note to self: However, I can probably leave the right hand side values section as it is

  // Before going over to the implementation of the main class, we declare and
  // define the function which describes the right hand side. This time, the
  // right hand side is vector-valued, as is the solution, so we will describe
  // the changes required for this in some more detail.
  //
  // To prevent cases where the return vector has not previously been set to
  // the right size we test for this case and otherwise throw an exception at
  // the beginning of the function. This could be done by writing
  // `Assert (values.size() == points.size(), some exception text)`, but
  // because checking for the equality in the sizes of two objects is
  // such a common operation, there is a short-cut: `AssertDimension`.
  // The operation behind this command is that it compares the two given
  // sizes and, if they are not equal, aborts the program with a suitable
  // error message that we don't have to write from scratch in all of the
  // places where we want to have this kind of check. (As for the other
  // `Assert` variations, the check is removed in optimized mode.)
  // Note that enforcing that output arguments
  // already have the correct size is a convention in deal.II, and enforced
  // almost everywhere. The reason is that we would otherwise have to check at
  // the beginning of the function and possibly change the size of the output
  // vector. This is expensive, and would almost always be unnecessary (the
  // first call to the function would set the vector to the right size, and
  // subsequent calls would only have to do redundant checks). In addition,
  // checking and possibly resizing the vector is an operation that can not be
  // removed if we can't rely on the assumption that the vector already has
  // the correct size; this is in contrast to the call to `Assert` that is
  // completely removed if the program is compiled in optimized mode.
  //
  // Likewise, if by some accident someone tried to compile and run the
  // program in only one space dimension (in which the elastic equations do
  // not make much sense since they reduce to the ordinary Laplace equation),
  // we terminate the program in the second assertion. The program will work
  // just fine in 3d, however.
  // template <int dim>
  // void right_hand_side(const std::vector<Point<dim>> &points,
  //                      std::vector<Tensor<1, dim>>   &values)
  // {
  //   AssertDimension(values.size(), points.size());
  //   Assert(dim >= 2, ExcNotImplemented());

    // The rest of the function implements computing force values. We will use
    // a constant (unit) force in x-direction located in two little circles
    // (or spheres, in 3d) around points (0.5,0) and (-0.5,0), and y-force in
    // an area around the origin; in 3d, the z-component of these centers is
    // zero as well.
    //
    // For this, let us first define two objects that denote the centers of
    // these areas. Note that upon construction of the Point objects, all
    // components are set to zero.
  //   Point<dim> point_1, point_2;
  //   point_1[0] = 0.5;
  //   point_2[0] = -0.5;

  //   for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
  //     {
  //       // If <code>points[point_n]</code> is in a circle (sphere) of radius
  //       // 0.2 around one of these points, then set the force in x-direction
  //       // to one, otherwise to zero:
  //       if (((points[point_n] - point_1).norm_square() < 0.2 * 0.2) ||
  //           ((points[point_n] - point_2).norm_square() < 0.2 * 0.2))
  //         values[point_n][0] = 1.0;
  //       else
  //         values[point_n][0] = 0.0;

  //       // Likewise, if <code>points[point_n]</code> is in the vicinity of the
  //       // origin, then set the y-force to one, otherwise to zero:
  //       if (points[point_n].norm_square() < 0.2 * 0.2)
  //         values[point_n][1] = 1.0;
  //       else
  //         values[point_n][1] = 0.0;
  //     }
  // }


  // @sect3{The <code>ElasticProblem</code> class implementation}

  // @sect4{ElasticProblem::ElasticProblem constructor}

  // Following is the constructor of the main class. As said before, we would
  // like to construct a vector-valued finite element that is composed of
  // several scalar finite elements (i.e., we want to build the vector-valued
  // element so that each of its vector components consists of the shape
  // functions of a scalar element). Of course, the number of scalar finite
  // elements we would like to stack together equals the number of components
  // the solution function has, which is <code>dim</code> since we consider
  // displacement in each space direction. The FESystem class can handle this:
  // we pass it the finite element of which we would like to compose the
  // system of, and how often to repeat it. There are different ways to
  // tell the FESystem constructor how to do this, but the one that is
  // closest to mathematical notation is to write out what we want to do
  // mathematically: We want to construct the finite element space
  // $Q_1^d$ where the index 1 corresponds to the polynomial degree and
  // the exponent $d$ to the space dimension -- because the *displacement*
  // we try to simulate here is a vector with exactly $d$ components. The
  // FESystem class then lets us create this space by initialization with
  // `FE_Q<dim>(1)^dim`, emulating the mathematical notation.
  //
  // (We could also have written `fe(FE_Q<dim>(1), dim)`, which would simply
  // have called a different constructor of the FESystem class that first
  // takes the "base element" and then a "multiplicity", i.e., a number that
  // indicates how many times the base element is to be repeated. The two
  // ways of writing things are entirely equivalent; we choose the one that
  // is closer to mathematical notation.)

  // TODO: CHECK IF THIS IS REALLY OK

  // @sect4{Constructor}

  // template <int dim>
  // // ElasticProblem<dim>::ElasticProblem()
  // ElasticProblem<dim>::ElasticProblem(const CmsParameters &parameters) //TODO: Check
  //   : parameters(parameters)
  //   , dof_handler(triangulation)
  //   , fe(FE_Q<dim>(1) ^ dim)
  // {}

  // In fact, the FESystem class has several more constructors which can
  // perform more complex operations than just stacking together several
  // scalar finite elements of the same type into one; we will get to know
  // these possibilities in later examples.


  // @sect4{ElasticProblem::setup_system}

  // Setting up the system of equations is identical to the function used in
  // the step-6 example. The DoFHandler class and all other classes used here
  // are fully aware that the finite element we want to use is vector-valued,
  // and take care of the vector-valuedness of the finite element
  // themselves. (In fact, they do not, but this does not need to bother you:
  // since they only need to know how many degrees of freedom there are per
  // vertex, line and cell, and they do not ask what they represent,
  // i.e. whether the finite element under consideration is vector-valued or
  // whether it is, for example, a scalar Hermite element with several degrees
  // of freedom on each vertex).
  // template <int dim>
  // void ElasticProblem<dim>::setup_system()
  // {
  //   dof_handler.distribute_dofs(fe);
  //   // solution.reinit(dof_handler.n_dofs());
  //   system_rhs.reinit(dof_handler.n_dofs());

  //   constraints.clear();
  //   DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    
  //   // Same boundary condition as step-8: u=0 on boundary id 0.
  //   VectorTools::interpolate_boundary_values(dof_handler,
  //                                            types::boundary_id(0),
  //                                            Functions::ZeroFunction<dim>(dim),
  //                                            constraints);
  //   constraints.close();

  //   DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  //   DoFTools::make_sparsity_pattern(dof_handler,
  //                                   dsp,
  //                                   constraints,
  //                                   /*keep_constrained_dofs = */ false);
  //   sparsity_pattern.copy_from(dsp);

  //   system_matrix.reinit(sparsity_pattern);
  //   mass_matrix.reinit(sparsity_pattern);
  //   stiffness_matrix.reinit(sparsity_pattern);
  // }


  // @sect4{ElasticProblem::assemble_system}

  // The big changes in this program are in the creation of matrix and right
  // hand side, since they are problem-dependent. We will go through that
  // process step-by-step, since it is a bit more complicated than in previous
  // examples.
  //
  // The first parts of this function are the same as before, however: setting
  // up a suitable quadrature formula, initializing an FEValues object for the
  // (vector-valued) finite element we use as well as the quadrature object,
  // and declaring a number of auxiliary arrays. In addition, we declare the
  // ever same two abbreviations: <code>n_q_points</code> and
  // <code>dofs_per_cell</code>. The number of degrees of freedom per cell we
  // now obviously ask from the composed finite element rather than from the
  // underlying scalar Q1 element. Here, it is <code>dim</code> times the
  // number of degrees of freedom per cell of the Q1 element, though this is
  // not explicit knowledge we need to care about:
  // template <int dim>
  // void ElasticProblem<dim>::assemble_stiffness_and_rhs() // TODO: CHECK IF NEW NAME MAKES SENSE
  // {
  //   const QGauss<dim> quadrature_formula(fe.degree + 1);

  //   FEValues<dim> fe_values(fe,
  //                           quadrature_formula,
  //                           update_values | update_gradients |
  //                             update_quadrature_points | update_JxW_values);

  //   const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  //   const unsigned int n_q_points    = quadrature_formula.size();

  //   FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  //   Vector<double>     cell_rhs(dofs_per_cell);

  //   std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  //   // As was shown in previous examples as well, we need a place where to
  //   // store the values of the coefficients at all the quadrature points on a
  //   // cell. In the present situation, we have two coefficients, lambda and
  //   // mu.

  //   const Functions::ConstantFunction<dim> lambda(parameters.lambda);
  //   const Functions::ConstantFunction<dim> mu(parameters.mu);

  //   std::vector<double> lambda_values(n_q_points);
  //   std::vector<double> mu_values(n_q_points);

    // Well, we could as well have omitted the above two arrays since we will
    // use constant coefficients for both lambda and mu, which can be declared
    // like this. They both represent functions always returning the constant
    // value 1.0. Although we could omit the respective factors in the
    // assemblage of the matrix, we use them here for purpose of
    // demonstration.
    // Functions::ConstantFunction<dim> lambda(1.), mu(1.);

    // Like the two constant functions above, we will call the function
    // right_hand_side just once per cell to make things simpler.
    // std::vector<Tensor<1, dim>> rhs_values(n_q_points);

    // stiffness_matrix = 0;
    // system_rhs = 0;

    // // Now we can begin with the loop over all cells:
    // for (const auto &cell : dof_handler.active_cell_iterators())
    //   {
    //     fe_values.reinit(cell);

    //     cell_matrix = 0;
    //     cell_rhs    = 0;

    //     // Next we get the values of the coefficients at the quadrature
    //     // points. Likewise for the right hand side:
    //     lambda.value_list(fe_values.get_quadrature_points(), lambda_values);
    //     mu.value_list(fe_values.get_quadrature_points(), mu_values);
    //     right_hand_side(fe_values.get_quadrature_points(), rhs_values);

        // Then assemble the entries of the local @ref GlossStiffnessMatrix "stiffness matrix" and right
        // hand side vector. This follows almost one-to-one the pattern
        // described in the introduction of this example.  One of the few
        // comments in place is that we can compute the number
        // <code>comp(i)</code>, i.e. the index of the only nonzero vector
        // component of shape function <code>i</code> using the
        // <code>fe.system_to_component_index(i).first</code> function call
        // below.
        //
        // (By accessing the <code>first</code> variable of the return value
        // of the <code>system_to_component_index</code> function, you might
        // already have guessed that there is more in it. In fact, the
        // function returns a <code>std::pair@<unsigned int, unsigned
        // int@></code>, of which the first element is <code>comp(i)</code>
        // and the second is the value <code>base(i)</code> also noted in the
        // introduction, i.e.  the index of this shape function within all the
        // shape functions that are nonzero in this component,
        // i.e. <code>base(i)</code> in the diction of the introduction. This
        // is not a number that we are usually interested in, however.)
        //
        // With this knowledge, we can assemble the local matrix
        // contributions:
        // for (const unsigned int i : fe_values.dof_indices())
        //   {
        //     const unsigned int component_i =
        //       fe.system_to_component_index(i).first;

        //     for (const unsigned int j : fe_values.dof_indices())
        //       {
        //         const unsigned int component_j =
        //           fe.system_to_component_index(j).first;

        //         for (const unsigned int q_point :
        //              fe_values.quadrature_point_indices())
        //           {
        //             cell_matrix(i, j) +=
                      // The first term is $(\lambda \partial_i u_i, \partial_j
                      // v_j) + (\mu \partial_i u_j, \partial_j v_i)$. Note
                      // that <code>shape_grad(i,q_point)</code> returns the
                      // gradient of the only nonzero component of the i-th
                      // shape function at quadrature point q_point. The
                      // component <code>comp(i)</code> of the gradient, which
                      // is the derivative of this only nonzero vector
                      // component of the i-th shape function with respect to
                      // the comp(i)th coordinate is accessed by the appended
                      // brackets.
                      // (                                                  //
                      //   (fe_values.shape_grad(i, q_point)[component_i] * //
                      //    fe_values.shape_grad(j, q_point)[component_j] * //
                      //    lambda_values[q_point])                         //
                      //   +                                                //
                      //   (fe_values.shape_grad(i, q_point)[component_j] * //
                      //    fe_values.shape_grad(j, q_point)[component_i] * //
                      //    mu_values[q_point])                             //
                      //   +                                                //
                        // The second term is $(\mu \nabla u_i, \nabla
                        // v_j)$. We need not access a specific component of
                        // the gradient, since we only have to compute the
                        // scalar product of the two gradients, of which an
                        // overloaded version of <tt>operator*</tt> takes
                        // care, as in previous examples.
                        //
                        // Note that by using the <tt>?:</tt> operator, we only
                        // do this if <tt>component_i</tt> equals
                        // <tt>component_j</tt>, otherwise a zero is added
                        // (which will be optimized away by the compiler).
  //                       ((component_i == component_j) ?        //
  //                          (fe_values.shape_grad(i, q_point) * //
  //                           fe_values.shape_grad(j, q_point) * //
  //                           mu_values[q_point]) :              //
  //                          0)                                  //
  //                       ) *                                    //
  //                     fe_values.JxW(q_point);                  //
  //                 }
  //             }
  //         }

  //       // Assembling the right hand side is also just as discussed in the
  //       // introduction:
  //       for (const unsigned int i : fe_values.dof_indices())
  //         {
  //           const unsigned int component_i =
  //             fe.system_to_component_index(i).first;

  //           for (const unsigned int q_point :
  //                fe_values.quadrature_point_indices())
  //             cell_rhs(i) += fe_values.shape_value(i, q_point) *
  //                            rhs_values[q_point][component_i] *
  //                            fe_values.JxW(q_point);
  //         }

  //       // The transfer from local degrees of freedom into the global matrix
  //       // and right hand side vector does not depend on the equation under
  //       // consideration, and is thus the same as in all previous
  //       // examples.
  //       cell->get_dof_indices(local_dof_indices);
  //       constraints.distribute_local_to_global(
  //         cell_matrix, cell_rhs, local_dof_indices, stiffness_matrix, system_rhs);
  //     }
  // }

  // @sect4{assemble_consistent_mass_matrix}

  // template <int dim>
  // void ElasticProblem<dim>::assemble_consistent_mass_matrix()
  // {
  //   const QGauss<dim> quadrature_formula(fe.degree + 1);

  //   const Functions::ConstantFunction<dim> density(parameters.rho);

  //   mass_matrix = 0;
  //   MatrixCreator::create_mass_matrix(
  //     dof_handler, quadrature_formula, mass_matrix, &density, constraints
  //   );
  // }

  // @sect4{Characteristic element length: minimum edge length}

  template <int dim>
  double ElasticProblem<dim>::compute_min_edge_length_of_cell(
    //const typename DoFHandler<dim>::active_cell_iterator &cell
    const typename Triangulation<dim>::active_cell_iterator &cell
  )
  {
    double min_length = std::numeric_limits<double>::max();

    for (unsigned int line = 0; line < GeometryInfo<dim>::lines_per_cell; ++line)
    {
      const unsigned int v0_index = GeometryInfo<dim>::line_to_cell_vertices(line,0);
      const unsigned int v1_index = GeometryInfo<dim>::line_to_cell_vertices(line,1);

      const Point<dim> &v0 = cell->vertex(v0_index);
      const Point<dim> &v1 = cell->vertex(v1_index);

      min_length = std::min(min_length, v0.distance(v1));
    }

    Assert(min_length > 0.0, ExcInternalError());
    return min_length;
  }

  template <int dim>
  void ElasticProblem<dim>::compute_characteristic_length(
    std::vector<double> &h_char
  ) const
  {
    h_char.clear();
    h_char.reserve(triangulation.n_active_cells());

    for (const auto &cell : dof_handler.active_cell_iterators())
      h_char.push_back(compute_min_edge_length_of_cell(cell));

    Assert(h_char.size() == triangulation.n_active_cells(), ExcInternalError());
  }

  // @sect4{Select critical cells by percentile}

  // template <int dim>
  // std::pair<double, std::vector<unsigned int>>
  // ElasticProblem<dim>::select_critical_cells_by_percentile(
  //   const std::vector<double> &h_char
  // ) const
  // {
  //   Assert(!h_char.empty(), ExcInternalError());
  //   Assert(parameters.critical_fraction > 0.0 && parameters.critical_fraction <= 1.0, ExcMessage("critical_fraction must be in (0,1]."));

  //   std::vector<unsigned int> ids(h_char.size());
  //   std::iota(ids.begin(), ids.end(), 0);

  //   std::sort(ids.begin(), ids.end(), [&](const unsigned int a, const unsigned int b){
  //     return h_char[a] < h_char[b];
  //   });

  //   const unsigned int n_select = std::max<unsigned int>(1, static_cast<unsigned int>(std::ceil(parameters.critical_fraction * static_cast<double>(ids.size()))));

  //   ids.resize(n_select);

  //   double h_reference = 0.0;
  //   for (const unsigned int id: ids)
  //     h_reference = std::max(h_reference, h_char[id]);

  //   Assert(h_reference > 0.0, ExcInternalError());
  //   return {h_reference, ids};
  // }


  // @sect4{Wave speed}

  // template <int dim>
  // double ElasticProblem<dim>::compute_pressure_wave_speed() const
  // {
  //   Assert(parameters.rho > 0.0, ExcMessage("Density must be positive."));
  //   Assert(parameters.lambda + 2.0 * parameters.mu > 0.0, ExcMessage("lambda + 2 mu must be positive."));

  //   return std::sqrt((parameters.lambda + 2.0 * parameters.mu) / parameters.rho);
  // }

  // @sect{Compute global critical time step}
  template <int dim>
  double ElasticProblem<dim>::compute_dt_cell(
    const double h,
    const double rho_local
    //const std::vector<double> &h_char,
    //const std::vector<double> &density_scaling
  ) const
  {
    // Assert(h_char.size() == triangulation.n_active_cells(), ExcInternalError());
    // Assert(density_scaling.size() == triangulation.n_active_cells(), ExcInternalError());

    const double denom = parameters.lambda + 2.0 * parameters.mu;
    Assert(denom > 0.0, ExcMessage("lambda + 2 must be positive."));
    Assert(rho_local > 0.0, ExcMessage("rho_local must be positive."));
    Assert(h > 0.0, ExcMessage("h must be positive."));

    // double dt_min = std::numeric_limits<double>::max();

    // for (unsigned int i = 0; i < h_char.size(); ++i)
    // {
    //   const double rho_i = parameters.rho * density_scaling[i];
    //   Assert(rho_i > 0.0, ExcMessage("Per-cell rho must be positive."));

    //   const double dt_i = h_char[i] * std:::sqrt(rho_i / denom);
    //   dt_min = std::min(dt_min, dt_i);
    // }

    return h * std::sqrt(rho_local / denom);
  }

  template <int dim>
  double ElasticProblem<dim>::compute_dt_crit_global(
    const std::vector<double> &h_char,
    const std::vector<double> &density_scaling
  ) const
  {
    Assert(h_char.size() == triangulation.n_active_cells(), ExcInternalError());
    Assert(density_scaling.size() == triangulation.n_active_cells(), ExcInternalError());

    double dt_min = std::numeric_limits<double>::max();

    for (unsigned int i = 0; i < h_char.size(); ++i)
    {
      const double rho_i = parameters.rho * density_scaling[i];
      const double dt_i = compute_dt_cell(h_char[i],rho_i);
      dt_min = std::min(dt_min, dt_i);
    }

    return dt_min;
  }

  // @sect4{Apply classical mass scaling (CMS)}

  template <int dim>
  typename ElasticProblem<dim>::CmsResult
  ElasticProblem<dim>::apply_cms_for_target_dt(
    const std::vector<double> &h_char, const double dt_target
  ) const
  {
    Assert(dt_target > 0.0, ExcMessage("dt_target must be positive."));
    Assert(h_char.size() == triangulation.n_active_cells(), ExcInternalError());

    CmsResult result;
    result.h_char = h_char;
    result.dt_target = dt_target;

    const double rho0 = parameters.rho;
    const double denom = parameters.lambda + 2.0 * parameters.mu;

    result.density_scaling.assign(h_char.size(), 1.0);

    // First compute dt_crit_global_before (all scaling=1)
    std::vector<double> ones(h_char.size(),1.0);
    result.dt_crit_global_before = compute_dt_crit_global(h_char, ones);

    // Mass accounting
    double total_mass = 0.0;
    double added_mass = 0.0;

    unsigned int cell_index = 0;
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      const double cell_volume = cell->measure() * parameters.thickness;
      total_mass += rho0 * cell_volume;

      // Per-cell dt before scaling
      const double dt_i = compute_dt_cell(h_char[cell_index],rho0);

      // Scale exactly those cells that violate dt_target:
      if (dt_i < dt_target)
      {
        // Want: dt_target = h * sqrt(rho_required / denom)
        // => rho_required = denom * (dt_target / h)^2
        const double h = h_char[cell_index];
        const double rho_required = denom * (dt_target * dt_target) / (h * h);

        if (rho_required > rho0)
        {
          const double alpha = rho_required / rho0;
          result.density_scaling[cell_index] = alpha;
          added_mass += (rho_required - rho0) * cell_volume;
        }
      }

      // const double h = h_char[cell_index];
      // const double rho_required = (dt_target * dt_target) * (parameters.lambda + 2.0 * parameters.mu) / (h * h);

      // Only apply mass scaling on the selected critical cells:
      // if (is_critical[cell_index] && rho_required > rho0)
      // {
      //   result.density_scaling[cell_index] = rho_required / rho0;
      //   added_mass += (rho_required - rho0) * cell_volume;
      // }

      ++cell_index;
    }

    result.total_mass = total_mass;
    result.added_mass = added_mass;
    result.added_mass_ratio = (total_mass > 0.0 ? added_mass / total_mass : 0.0);

    // Count scaled cells
    for (const double a : result.density_scaling)
      if (a > 1.0 + 1e-12)
        ++result.n_scaled;

    result.scaled_fraction = (
      result.density_scaling.empty() ? 0.0 : static_cast<double>(result.n_scaled) / static_cast<double>(result.density_scaling.size()));

    // dt_crit_global_after
    result.dt_crit_global_after = compute_dt_crit_global(h_char, result.density_scaling);

    // double dt_ref = std::numeric_limits<double>::max();
    // for (const unsigned int id : critical_cell_ids)
    // {
    //   const double h = h_char[id];
    //   const double dt_cell = h * std::sqrt(parameters.rho / (parameters.lambda + 2.0 * parameters.mu));
    //   dt_ref = std::min(dt_ref, dt_cell);
    // }

    // result.dt_crit_reference = dt_ref;

    // unsigned int n_scaled = 0;
    // for (const double a : result.density_scaling)
    //   if (a > 1.0 + 1e-12)
    //     ++n_scaled;

    // // Print
    // std::cout << "Scaled cells: " << n_scaled
    //           << " / " << result.density_scaling.size() << std::endl;

    return result;

  }


  //@sect4{Choose dt_target and ensure mass increase stays acceptable}

  template <int dim>
  typename ElasticProblem<dim>::CmsResult
  ElasticProblem<dim>::choose_dt_and_apply_cms(const std::vector<double> &h_char) const
  {
    // global dt before CMS
    std::vector<double> ones(h_char.size(), 1.0);
    const double dt_global_before = compute_dt_crit_global(h_char, ones);

    double dt_target = parameters.dt_factor * dt_global_before;

    CmsResult best = apply_cms_for_target_dt(h_char, dt_target);

    if (!parameters.auto_reduce_dt_target)
      return best;

    // If constraints violated, reduce dt_target iteratively.
    for (unsigned int k = 0; k < parameters.max_dt_adjustment_steps; ++k)
    {
      const bool too_much_mass = (best.added_mass_ratio > parameters.max_added_mass_fraction);

      const bool too_many_cells = (best.scaled_fraction > parameters.max_scaled_cell_fraction);

      if(!too_much_mass && !too_many_cells)
        return best;

      dt_target *= parameters.dt_shrink_factor;
      best = apply_cms_for_target_dt(h_char, dt_target);
    }

    // Return best effort after max iterations
    return best;

    // const auto [h_ref, critical_ids] = select_critical_cells_by_percentile(h_char);

    // const double dt_crit_reference = h_ref * std::sqrt(parameters.rho / (parameters.lambda + 2.0 * parameters.mu));

    // double dt_target = parameters.dt_factor * dt_crit_reference;

    // for (unsigned int attempt = 0; attempt < parameters.max_dt_adjustment_steps; ++attempt)
    // {
    //   CmsResult cms = apply_classical_mass_scaling(h_char, critical_ids, dt_target);

    //   cms.dt_crit_reference = dt_crit_reference;

    //   if (cms.added_mass_ratio <= parameters.max_added_mass_fraction)
    //     return cms;

    //   dt_target *= parameters.dt_shrink_factor;
    // }
    // return apply_classical_mass_scaling(h_char,critical_ids,dt_target);
  }

  // @sect4{}

  template <int dim>
  void ElasticProblem<dim>::verify_mesh_is_quad_and_non_degenerate() const
  {
    Assert(dim == 2, ExcNotImplemented());

    // 1) Check number of vertices per cell (quad should have 4)
    for (const auto &cell : triangulation.active_cell_iterators())
    {
      Assert(cell->n_vertices() == 4, ExcInternalError());
      Assert(cell->measure() > 0.0, ExcMessage("Found a cell with non-positive area."));
    }

    // 2) Check total area
    double area_sum = 0.0;
    for (const auto &cell : triangulation.active_cell_iterators())
      area_sum += cell->measure();

    const double expected_area = 4.0; // [-1,1]x[-1,1]
    std::cout << "Mesh check: sum(cell areas) = " << area_sum << " (expected ~ " << expected_area << ")\n";

    AssertThrow(std::abs(area_sum - expected_area) < 1e-12, ExcMessage("Total mesh area does not match expected domain area."));
  }

  // @sect4{}
  template <int dim>
  void ElasticProblem<dim>::print_cell_debug_table
  (
    const std::vector<double> &h_char, const std::vector<double> &density_scaling, const double dt_target
  ) const
  {
    Assert(h_char.size() == triangulation.n_active_cells(), ExcInternalError());
    Assert(density_scaling.size() == triangulation.n_active_cells(), ExcInternalError());

    const double denom = parameters.lambda + 2.0 * parameters.mu;
    const double rho0 = parameters.rho;

    std::cout << "\n--- Per-cell debug (first all cells) ---\n";
    std::cout << "idx   area     e0       e1       e2     e3      h_min      dt_cell      alpha       scaled\n";

    unsigned int idx = 0;
    for (const auto &cell : triangulation.active_cell_iterators())
    {
      // edge lengths for quad
      double e[4];
      for (unsigned int l = 0; l < 4; ++l)
      {
        const unsigned int v0 = GeometryInfo<dim>::line_to_cell_vertices(l, 0);
        const unsigned int v1 = GeometryInfo<dim>::line_to_cell_vertices(l, 1);
        e[l] = cell->vertex(v0).distance(cell->vertex(v1));
      }

      const double hmin = *std::min_element(std::begin(e), std::end(e));
      const double dt_i = hmin * std::sqrt(rho0 / denom);

      const double alpha = density_scaling[idx];
      const bool scaled = (alpha > 1.0 + 1e-12);

      std::cout << idx << "   "
                << cell->measure() << "    "
                << e[0] << "    " << e[1] << "    " << e[2] << "    " << e[3] << "   "
                << hmin << "    "
                << dt_i << "    "
                << alpha << "    "
                << (scaled ? "YES" : "NO")
                << "\n";

      // check h_char matches what we just computed
      AssertThrow(std::abs(h_char[idx] - hmin) < 1e-14, ExcMessage("h_char mismatch: stored h_char != recomputed min edge length"));

      ++idx;
    }

    std::cout << "--- End per-cell debug ---\n";
  }

  // @sect4{Output density scaling to VTK}

  template <int dim>
  void ElasticProblem<dim>::output_density_scaling(
    const std::vector<double> &density_scaling
  ) const
  {
    Assert(density_scaling.size() == triangulation.n_active_cells(), ExcInternalError());

    Vector<float> scaling_field(triangulation.n_active_cells());

    for (unsigned int i = 0; i < scaling_field.size(); ++i)
      scaling_field[i] = static_cast<float>(density_scaling[i]);

    DataOut<dim> data_out;
    data_out.attach_triangulation(triangulation);
    data_out.add_data_vector(scaling_field, "rho_scaling");
    data_out.build_patches();

    std::ofstream output("cms_rho_scaling.vtk");
    data_out.write_vtk(output);
  }

  // @sect4{ElasticProblem::solve}

  // The solver does not care about where the system of equations comes from, as
  // long as it is positive definite and symmetric (which are the
  // requirements for the use of the CG solver), which the system indeed
  // is. Therefore, we need not change anything.
  // template <int dim>
  // void ElasticProblem<dim>::solve()
  // {
  //   SolverControl            solver_control(1000, 1e-6 * system_rhs.l2_norm());
  //   SolverCG<Vector<double>> cg(solver_control);

  //   PreconditionSSOR<SparseMatrix<double>> preconditioner;
  //   preconditioner.initialize(system_matrix, 1.2);

  //   cg.solve(system_matrix, solution, system_rhs, preconditioner);

  //   constraints.distribute(solution);
  // }


  // @sect4{ElasticProblem::refine_grid}

  // The function that does the refinement of the grid is the same as in the
  // step-6 example. The quadrature formula is adapted to the linear elements
  // again. Note that the error estimator by default adds up the estimated
  // obtained from all components of the finite element solution, i.e., it
  // uses the displacement in all directions with the same weight. If we would
  // like the grid to be adapted to the x-displacement only, we could pass the
  // function an additional parameter which tells it to do so and do not
  // consider the displacements in all other directions for the error
  // indicators. However, for the current problem, it seems appropriate to
  // consider all displacement components with equal weight.
  // template <int dim>
  // void ElasticProblem<dim>::refine_grid()
  // {
  //   Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

  //   KellyErrorEstimator<dim>::estimate(dof_handler,
  //                                      QGauss<dim - 1>(fe.degree + 1),
  //                                      {},
  //                                      solution,
  //                                      estimated_error_per_cell);

  //   GridRefinement::refine_and_coarsen_fixed_number(triangulation,
  //                                                   estimated_error_per_cell,
  //                                                   0.3,
  //                                                   0.03);

  //   triangulation.execute_coarsening_and_refinement();
  // }


  // @sect4{ElasticProblem::output_results}

  // The output happens mostly as has been shown in previous examples
  // already. The only difference is that the solution function is vector
  // valued. The DataOut class takes care of this automatically, but we have
  // to give each component of the solution vector a different name.
  //
  // To do this, the DataOut::add_vector() function wants a vector of
  // strings. Since the number of components is the same as the number
  // of dimensions we are working in, we use the <code>switch</code>
  // statement below.
  //
  // We note that some graphics programs have restriction on what
  // characters are allowed in the names of variables. deal.II therefore
  // supports only the minimal subset of these characters that is supported
  // by all programs. Basically, these are letters, numbers, underscores,
  // and some other characters, but in particular no whitespace and
  // minus/hyphen. The library will throw an exception otherwise, at least
  // if in debug mode.
  //
  // After listing the 1d, 2d, and 3d case, it is good style to let the
  // program die if we run into a case which we did not consider. You have
  // previously already seen the use of the `Assert` macro that generates
  // aborts the program with an error message if a condition is not satisfied
  // (see step-5, for example). We could use this in the `default` case
  // below, in the form `Assert(false, ExcNotImplemented())` -- in other words,
  // the "condition" here is always `false`, and so the assertion always fails
  // and always aborts the program whenever it gets to the default statement.
  // This is perhaps more difficult to read than necessary, and consequently
  // there is a short-cut: `DEAL_II_NOT_IMPLEMENTED()`. It does the same
  // as the form above (with the minor difference that it also aborts the
  // program in release mode). It is written in all-caps because that makes
  // it stand out visually (and also because it is not actually a function,
  // but a macro).
  // template <int dim>
  // void ElasticProblem<dim>::output_results(const unsigned int cycle) const
  // {
  //   DataOut<dim> data_out;
  //   data_out.attach_dof_handler(dof_handler);

  //   std::vector<std::string> solution_names;
  //   switch (dim)
  //     {
  //       case 1:
  //         solution_names.emplace_back("displacement");
  //         break;
  //       case 2:
  //         solution_names.emplace_back("x_displacement");
  //         solution_names.emplace_back("y_displacement");
  //         break;
  //       case 3:
  //         solution_names.emplace_back("x_displacement");
  //         solution_names.emplace_back("y_displacement");
  //         solution_names.emplace_back("z_displacement");
  //         break;
  //       default:
  //         DEAL_II_NOT_IMPLEMENTED();
  //     }

    // After setting up the names for the different components of the
    // solution vector, we can add the solution vector to the list of
    // data vectors scheduled for output. Note that the following
    // function takes a vector of strings as second argument, whereas
    // the one which we have used in all previous examples took a
    // single string there (which was the right choice because
    // we had only a single solution variable in all previous examples).
  //   data_out.add_data_vector(solution, solution_names);
  //   data_out.build_patches();

  //   std::ofstream output("solution-" + std::to_string(cycle) + ".vtk");
  //   data_out.write_vtk(output);
  // }



  // @sect4{ElasticProblem::run}

  // The <code>run</code> function does the same things as in step-6, for
  // example. This time, we use the square [-1,1]^d as domain, and we refine
  // it globally four times before starting the first iteration.
  //
  // The reason for refining is a bit accidental: we use the QGauss
  // quadrature formula with two points in each direction for integration of the
  // right hand side; that means that there are four quadrature points on each
  // cell (in 2d). If we only refine the initial grid once globally, then there
  // will be only four quadrature points in each direction on the
  // domain. However, the right hand side function was chosen to be rather
  // localized and in that case, by pure chance, it happens that all quadrature
  // points lie at points where the right hand side function is zero (in
  // mathematical terms, the quadrature points happen to be at points outside
  // the <i>support</i> of the right hand side function). The right hand side
  // vector computed with quadrature will then contain only zeroes (even though
  // it would of course be nonzero if we had computed the right hand side vector
  // exactly using the integral) and the solution of the system of
  // equations is the zero vector, i.e., a finite element function that is zero
  // everywhere. In a sense, we
  // should not be surprised that this is happening since we have chosen
  // an initial grid that is totally unsuitable for the problem at hand.
  //
  // The unfortunate thing is that if the discrete solution is constant, then
  // the error indicators computed by the KellyErrorEstimator class are zero
  // for each cell as well, and the call to
  // Triangulation::refine_and_coarsen_fixed_number() will not flag any cells
  // for refinement (why should it if the indicated error is zero for each
  // cell?). The grid in the next iteration will therefore consist of four
  // cells only as well, and the same problem occurs again.
  //
  // The conclusion needs to be: while of course we will not choose the
  // initial grid to be well-suited for the accurate solution of the problem,
  // we must at least choose it such that it has the chance to capture the
  // important features of the solution. In this case, it needs to be able to
  // see the right hand side. Thus, we refine globally four times. (Any larger
  // number of global refinement steps would of course also work.)
  template <int dim>
  void ElasticProblem<dim>::run()
  {
    Assert(dim == 2, ExcMessage("This step-101 is written for dim=2."));

    make_conforming_graded_mesh();

    dof_handler.distribute_dofs(fe);

    verify_mesh_is_quad_and_non_degenerate();

    // GridGenerator::hyper_cube(triangulation,-1,1);
    
    // // Base mesh
    // triangulation.refine_global(3);

    // // Locally refine a small region to create a few very small cells.
    // // Here: refine cells whose center is close to (0.6, 0.6).
    // for (unsigned int step = 0; step < 4; ++step)
    // {
    //   for (const auto &cell : triangulation.active_cell_iterators())
    //     if (cell->center().distance(Point<2>(0.6,0.6))<0.35)
    //       cell->set_refine_flag();

    //   triangulation.execute_coarsening_and_refinement();
    // }

    // setup_system();

    std::cout << "Step-101 (CMS) (CMS timestep targetting)\n";
    std::cout << "Active cells: " << triangulation.n_active_cells() << '\n';
    std::cout << "DoFs: " << dof_handler.n_dofs() << '\n';

    // assemble_stiffness_and_rhs();
    // assemble_consistent_mass_matrix();

    // Characteristic sizes
    std::vector<double> h_char;
    compute_characteristic_length(h_char);

    // const double c_p = compute_pressure_wave_speed();

    // const std::vector<double> density_scaling_before(triangulation.n_active_cells(),1.0);
    // const double dt_crit_before = compute_dt_crit_global(h_char, density_scaling_before);

    // const auto cms = choose_dt_and_apply_cms(h_char);

    // const double dt_crit_after = compute_dt_crit_global(h_char, cms.density_scaling);

    const auto cms = choose_dt_and_apply_cms(h_char);

    print_cell_debug_table(h_char, cms.density_scaling, cms.dt_target);

    double total_area = 0.0;
    for (const auto &cell : triangulation.active_cell_iterators())
      total_area += cell->measure();

    std::cout << "\nTotal model area = " << total_area << "\n";

    std::cout << "\nTargetting:\n";
    std::cout << "dt_factor = " << parameters.dt_factor << '\n';
    std::cout << "dt_crit_global (before) = " << cms.dt_crit_global_before << '\n';
    std::cout << "dt_target = " << cms.dt_target << '\n';
    std::cout << "dt_crit_global (after) = " << cms.dt_crit_global_after << "\n";
    std::cout << "improvement factor = " << (cms.dt_crit_global_after / cms.dt_crit_global_before) << "\n";

    std::cout << "\nCMS impact:\n";
    std::cout << "Scaled cells: " << cms.n_scaled << " / " << cms.density_scaling.size() << " (" << 100.0 * cms.scaled_fraction << "%)\n";

    const bool too_many_cells = (cms.scaled_fraction > parameters.max_scaled_cell_fraction);

    const bool too_much_mass = (cms.added_mass_ratio > parameters.max_added_mass_fraction);

    if (too_many_cells || too_much_mass)
    {
      std::cout << "\nWARNING: CMS may not be reasonable for this setup.\n";
      
      if (too_many_cells)
      {
        std::cout << " - Scaled-cell fraction exceeds limit: " << 100.0 * cms.scaled_fraction << " %\n";
      }

      if (too_much_mass)
      {
        std::cout << " - Added-mass ratio exceeds limit: " << 100.0 * cms.added_mass_ratio << " % > " << 100.0 * parameters.max_added_mass_fraction << " %\n";
      }

      std::cout << "Recommendations:\n";
      std::cout << " - Reduce dt_factor (smaller timestep increase), or \n";
      std::cout << " - Remesh with fewer extremely small cells / smoother grading.\n";
    } 
    
    if (parameters.write_density_scaling_output)
    {
      output_density_scaling(cms.density_scaling);
      std::cout << "\nWrote: cms_rho_scaling.vtk (cell data: rho_scaling)\n";
    }
  }

    // std::cout << "\nGlobal critical time step check:\n";
    // std::cout << "dt_crit_before = " << dt_crit_before << '\n';
    // std::cout << "dt_crit_after = " << dt_crit_after << '\n';
    // std::cout << "ratio (after/before) = " << (dt_crit_after / dt_crit_before) << '\n';

    // std::cout << "\nMaterial parameters: \n";
    // std::cout << "lambda = " << parameters.lambda << '\n';
    // std::cout << "mu = " << parameters.mu << '\n';
    // std::cout << "rho = " << parameters.rho << '\n';
    // std::cout << "c_p = " << c_p << " (pressure wave speed)\n";

    // std::cout << "\nCFL / CMS diagnostics:\n";
    // std::cout << "critical_fraction = " << parameters.critical_fraction << '\n';
    // std::cout << "dt_crit_reference = " << cms.dt_crit_reference << '\n';
    // std::cout << "dt_target = " << cms.dt_target << '\n';

    // std::cout << "\nMass scaling assessment:\n";
    // std::cout << "total_mass = " << cms.total_mass << '\n';
    // std::cout << "added_mass = " << cms.added_mass << '\n';
    // std::cout << "added_mass_ratio = " << 100.0 * cms.added_mass_ratio << " %\n";
    // std::cout << "max_allowed = " << 100.0 * parameters.max_added_mass_fraction << " %\n";

  //   if (parameters.write_density_scaling_output)
  //   {
  //     output_density_scaling(cms.density_scaling);
  //     std::cout << "\nWrote: cms_rho_scaling.vtk (cell data: rho_scaling)\n";
  //   }
  // }
    // for (unsigned int cycle = 0; cycle < 8; ++cycle)
    //   {
    //     std::cout << "Cycle " << cycle << ':' << std::endl;

    //     if (cycle == 0)
    //       {
    //         GridGenerator::hyper_cube(triangulation, -1, 1);
    //         triangulation.refine_global(4);
    //       }
    //     else
    //       refine_grid();

    //     std::cout << "   Number of active cells:       "
    //               << triangulation.n_active_cells() << std::endl;

    //     setup_system();

    //     std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
    //               << std::endl;

    //     assemble_system();
    //     solve();
    //     output_results(cycle);
    //   }
 // }
} // namespace Step101

// @sect3{The <code>main</code> function}

// After closing the <code>Step101</code> namespace in the last line above, the
// following is the main function of the program and is again exactly like in
// step-6 (apart from the changed class names, of course).
int main()
{
  try
    {
      //Step101::ElasticProblem<2> elastic_problem_2d;
      Step101::CmsParameters prm;

      prm.lambda = 1.0;
      prm.mu = 1.0;
      prm.rho = 1.0;

      prm.thickness = 1.0;

      prm.dt_factor = 1.05;

      prm.max_added_mass_fraction = 0.10;
      prm.max_scaled_cell_fraction = 0.10;

      prm.auto_reduce_dt_target = true;
      prm.dt_shrink_factor = 0.95;
      prm.max_dt_adjustment_steps = 25;

      prm.write_density_scaling_output = true;

      Step101::ElasticProblem<2> problem(prm);
      problem.run();
      //elastic_problem_2d.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
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
      std::cerr << std::endl
                << std::endl
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

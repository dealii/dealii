/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2021 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth, Colorado State University, 2021.
 * Based on step-15 by Sven Wetterauer, University of Heidelberg, 2012.
 */


// @sect3{Include files}

// The majority of the include files used in this program are
// well known from step-6 and similar programs:

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/fe/mapping_q.h>
#include <deal.II/matrix_free/fe_point_evaluation.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>


// The ones that are new are only the following three: The first declares the
// DiscreteTime class that helps us keep track of time in a time-dependent
// simulation. The latter two provide all of the particle functionality,
// namely a way to keep track of particles located on a mesh (the
// Particles::ParticleHandler class) and the ability to output these
// particles' locations and their properties for the purposes of
// visualization (the Particles::DataOut class).
#include <deal.II/base/discrete_time.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/data_out.h>

#include <fstream>

using namespace dealii;


// @sect3{Global definitions}

// As is customary, we put everything that corresponds to the details of the
// program into a namespace of its own. At the top, we define a few constants
// for which we would rather use symbolic names than hard-coded numbers.
//
// Specifically, we define numbers for
// @ref GlossBoundaryIndicator "boundary indicators"
// for the various parts of the geometry, as well as the physical properties
// of electrons and other specifics of the setup we use here.
//
// For the boundary indicators, let us start enumerating at some
// random value 101. The principle here is to use numbers that are
// *uncommon*. If there are pre-defined boundary indicators previously
// set by the `GridGenerator` functions, they will likely be small
// integers starting from zero, but not in this rather randomly chosen
// range. Using numbers such as those below avoids the possibility for
// conflicts, and also reduces the temptation to just spell these
// numbers out in the program (because you will probably never
// remember which is which, whereas you might have been tempted if
// they had started at 0).
namespace Step83
{
  namespace BoundaryIds
  {
    constexpr types::boundary_id open          = 101;
    constexpr types::boundary_id cathode       = 102;
    constexpr types::boundary_id focus_element = 103;
    constexpr types::boundary_id anode         = 104;
  } // namespace BoundaryIds

  namespace Constants
  {
    constexpr double electron_mass   = 9.1093837015e-31;
    constexpr double electron_charge = 1.602176634e-19;

    constexpr double V0 = 1;

    constexpr double E_threshold = 0.05;

    constexpr double electrons_per_particle = 3e15;
  } // namespace Constants


  // @sect3{The main class}

  // The following is then the main class of this program. It has,
  // fundamentally, the same structure as step-6 and many other
  // tutorial programs. This includes the majority of the member
  // functions (with the purpose of the rest probably self-explanatory
  // from their names) as well as only a small number of member
  // variables beyond those of step-6, all of which are related to
  // dealing with particles.
  template <int dim>
  class CathodeRaySimulator
  {
  public:
    CathodeRaySimulator();

    void run();

  private:
    void make_grid();
    void setup_system();
    void assemble_system();
    void solve_field();
    void refine_grid();

    void create_particles();
    void move_particles();
    void track_lost_particle(
      const typename Particles::ParticleIterator<dim> &        particle,
      const typename Triangulation<dim>::active_cell_iterator &cell);


    void update_timestep_size();
    void output_results() const;

    Triangulation<dim>        triangulation;
    MappingQGeneric<dim>      mapping;
    FE_Q<dim>                 fe;
    DoFHandler<dim>           dof_handler;
    AffineConstraints<double> constraints;

    SparseMatrix<double> system_matrix;
    SparsityPattern      sparsity_pattern;

    Vector<double> solution;
    Vector<double> system_rhs;

    Particles::ParticleHandler<dim> particle_handler;
    types::particle_index           next_unused_particle_id;
    types::particle_index           n_recently_lost_particles;
    types::particle_index           n_total_lost_particles;
    types::particle_index           n_particles_lost_through_anode;

    DiscreteTime time;
  };



  // @sect3{The <code>CathodeRaySimulator</code> class implementation}

  // @sect4{The <code>CathodeRaySimulator</code> constructor}

  // So then let us get started on the implementation. What the constructor
  // does is really only a straight-forward initialization of all of the member
  // variables at the top. The only two worth mentioning are the
  // `particle_handler`, which is handed a reference to the triangulation
  // on which the particles will live (currently of course still empty,
  // but the particle handler stores the reference and will use it once
  // particles are added -- which happens after the triangulation is built).
  // The other piece of information it gets is how many "properties"
  // each particle needs to store. Here, all we need each particle to
  // remember is its current velocity, i.e., a vector with `dim`
  // components. There are, however, other intrinsic properties that
  // each particle has and that the Particles::ParticleHandler class
  // automatically and always makes sure are available; in particular,
  // these are the current location of a particle, the cell it is on,
  // it's reference location within that cell, and the particle's ID.
  //
  // The only other variable of interest is `time`, an object of type
  // DiscreteTime. It keeps track of the current time we are in a
  // time-dependent simulation, and is initialized with the start time
  // (zero) and end time ($10^{-4}$). We will later set the time step
  // size in `update_timestep_size()`.
  //
  // The body of the constructor consists of a piece of code we have
  // already discussed in the introduction. Namely, we make sure that the
  // `track_lost_particle()` function is called by the `particle_handler`
  // object every time a particle leaves the domain.
  template <int dim>
  CathodeRaySimulator<dim>::CathodeRaySimulator()
    : mapping(1)
    , fe(2)
    , dof_handler(triangulation)
    , particle_handler(triangulation, mapping, /*n_properties=*/dim)
    , next_unused_particle_id(0)
    , n_recently_lost_particles(0)
    , n_total_lost_particles(0)
    , n_particles_lost_through_anode(0)
    , time(0, 1e-4)
  {
    particle_handler.signals.particle_lost.connect(
      [this](const typename Particles::ParticleIterator<dim> &        particle,
             const typename Triangulation<dim>::active_cell_iterator &cell) {
        this->track_lost_particle(particle, cell);
      });
  }



  // @sect4{The <code>CathodeRaySimulator::make_grid</code> function}

  // The next function is then responsible for generating the mesh on which
  // we want to solve. Recall how the domain looks like:
  //   <p align="center">
  //     <img
  //     src="https://www.dealii.org/images/steps/developer/step-19.geometry.png"
  //          alt="The geometry used in this program"
  //          width="600">
  //   </p>
  // We subdivide this geometry into a mesh of $4\times 2$ cells that looks
  // like this:
  // @code
  //   *---*---*---*---*
  //   \   |   |   |   |
  //    *--*---*---*---*
  //   /   |   |   |   |
  //   *---*---*---*---*
  // @endcode
  // The way this is done is by first defining where the $15=5\times 3$
  // vertices are located -- here, we say that they are on integer points
  // with the middle one on the left side moved to the right by a value of
  // `delta=0.5`.
  //
  // In the following, we then have to say which vertices together form
  // the 8 cells. The following code is then entirely equivalent to what
  // we also do in step-14:
  template <int dim>
  void CathodeRaySimulator<dim>::make_grid()
  {
    static_assert(dim == 2,
                  "This function is currently only implemented for 2d.");

    const double       delta = 0.5;
    const unsigned int nx    = 5;
    const unsigned int ny    = 3;

    const std::vector<Point<dim>> vertices //
      = {{0, 0},
         {1, 0},
         {2, 0},
         {3, 0},
         {4, 0},
         {delta, 1},
         {1, 1},
         {2, 1},
         {3, 1},
         {4, 1},
         {0, 2},
         {1, 2},
         {2, 2},
         {3, 2},
         {4, 2}};
    AssertDimension(vertices.size(), nx * ny);

    const std::vector<unsigned int> cell_vertices[(nx - 1) * (ny - 1)] = {
      {0, 1, nx + 0, nx + 1},
      {1, 2, nx + 1, nx + 2},
      {2, 3, nx + 2, nx + 3},
      {3, 4, nx + 3, nx + 4},

      {5, nx + 1, 2 * nx + 0, 2 * nx + 1},
      {nx + 1, nx + 2, 2 * nx + 1, 2 * nx + 2},
      {nx + 2, nx + 3, 2 * nx + 2, 2 * nx + 3},
      {nx + 3, nx + 4, 2 * nx + 3, 2 * nx + 4}};

    // With these arrays out of the way, we can move to slightly higher
    // higher-level data structures. We create a vector of CellData
    // objects that store for each cell to be created the vertices in
    // question as well as the @ref GlossMaterialId "material id" (which
    // we will here simply set to zero since we don't use it in the program).
    //
    // This information is then handed to the
    // Triangulation::create_triangulation() function, and the mesh is twice
    // globally refined.
    std::vector<CellData<dim>> cells((nx - 1) * (ny - 1), CellData<dim>());
    for (unsigned int i = 0; i < cells.size(); ++i)
      {
        cells[i].vertices    = cell_vertices[i];
        cells[i].material_id = 0;
      }

    triangulation.create_triangulation(
      vertices,
      cells,
      SubCellData()); // No boundary information

    triangulation.refine_global(2);

    // The remaining part of the function loops over all cells and their faces,
    // and if a face is at the boundary determines which boundary indicator
    // should be applied to it. The various conditions should make sense if
    // you compare the code with the picture of the geometry above.
    //
    // Once done with this step, we refine the mesh once more globally.
    for (auto &cell : triangulation.active_cell_iterators())
      for (auto &face : cell->face_iterators())
        if (face->at_boundary())
          {
            if ((face->center()[0] > 0) && (face->center()[0] < 0.5) &&
                (face->center()[1] > 0) && (face->center()[1] < 2))
              face->set_boundary_id(BoundaryIds::cathode);
            else if ((face->center()[0] > 0) && (face->center()[0] < 2))
              face->set_boundary_id(BoundaryIds::focus_element);
            else if ((face->center()[0] > 4 - 1e-12) &&
                     ((face->center()[1] > 1.5) || (face->center()[1] < 0.5)))
              face->set_boundary_id(BoundaryIds::anode);
            else
              face->set_boundary_id(BoundaryIds::open);
          }

    triangulation.refine_global(1);
  }


  // @sect4{The <code>CathodeRaySimulator::setup_system</code> function}

  // The next function in this program deals with setting up the various
  // objects related to solving the partial differential equations. It is
  // in essence a copy of the corresponding function in step-6 and requires
  // no further discussion.
  template <int dim>
  void CathodeRaySimulator<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    VectorTools::interpolate_boundary_values(dof_handler,
                                             BoundaryIds::cathode,
                                             Functions::ConstantFunction<dim>(
                                               -Constants::V0),
                                             constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             BoundaryIds::focus_element,
                                             Functions::ConstantFunction<dim>(
                                               -Constants::V0),
                                             constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             BoundaryIds::anode,
                                             Functions::ConstantFunction<dim>(
                                               +Constants::V0),
                                             constraints);
    constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ false);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
  }


  // @sect4{The <code>CathodeRaySimulator::assemble_system</code> function}

  // The function that computes
  // the matrix entries is again in essence a copy of the
  // corresponding function in step-6:
  template <int dim>
  void CathodeRaySimulator<dim>::assemble_system()
  {
    system_matrix = 0;
    system_rhs    = 0;

    const QGauss<dim> quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_matrix = 0;
        cell_rhs    = 0;

        fe_values.reinit(cell);

        for (const unsigned int q_index : fe_values.quadrature_point_indices())
          for (const unsigned int i : fe_values.dof_indices())
            {
              for (const unsigned int j : fe_values.dof_indices())
                cell_matrix(i, j) +=
                  (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                   fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                   fe_values.JxW(q_index));           // dx
            }

        // The only interesting part of this function is how it forms the right
        // hand side of the linear system. Recall that the right hand side
        // of the PDE is
        // @f[
        //   \sum_p (N e)\delta(\mathbf x-\mathbf x_p),
        // @f]
        // where we have used $p$ to index the particles here to avoid
        // confusion with the shape function $\varphi_i$; $\mathbf x_p$
        // is the position of the $p$th particle.
        //
        // When multiplied by a test function $\varphi_i$ and integrated over
        // the domain results in a right hand side vector
        // @f{align*}{
        //   F_i &= \int_\Omega \varphi_i (\mathbf x)\left[
        //                \sum_p (N e)\delta(\mathbf x-\mathbf x_p) \right] dx
        //   \\  &=  \sum_p (N e) \varphi_i(\mathbf x_p).
        // @f}
        // Note that the final line no longer contains an integral, and
        // consequently also no occurrence of $dx$ which would require the
        // appearance of the `JxW` symbol in our code.
        //
        // For a given cell $K$, this cell's contribution to the right hand
        // side is then
        // @f{align*}{
        //   F_i^K &= \sum_{p, \mathbf x_p\in K} (N e) \varphi_i(\mathbf x_p),
        // @f}
        // i.e., we only have to worry about those particles that are actually
        // located on the current cell $K$.
        //
        // In practice, what we do here is the following: If there are any
        // particles on the current cell, then we first obtain an iterator range
        // pointing to the first particle of that cell as well as the particle
        // past the last one on this cell (or the end iterator) -- i.e., a
        // half-open range as is common for C++ functions. Knowing now the list
        // of particles, we query their reference locations (with respect to
        // the reference cell), evaluate the shape functions in these reference
        // locations, and compute the force according to the formula above
        // (without any FEValues::JxW).
        //
        // @note It is worth pointing out that calling the
        //   Particles::ParticleHandler::particles_in_cell() and
        //   Particles::ParticleHandler::n_particles_in_cell() functions is not
        //   very efficient on problems with a large number of particles. But it
        //   illustrates the easiest way to write this algorithm, and so we are
        //   willing to incur this cost for the moment for expository purposes.
        //   We discuss the issue in more detail in the
        //   <a href="#extensions">"possibilities for extensions" section</a>
        //   below, and use a better approach in step-70, for example.
        if (particle_handler.n_particles_in_cell(cell) > 0)
          for (const auto &particle : particle_handler.particles_in_cell(cell))
            {
              const Point<dim> &reference_location =
                particle.get_reference_location();
              for (const unsigned int i : fe_values.dof_indices())
                cell_rhs(i) +=
                  (fe.shape_value(i, reference_location) * // phi_i(x_p)
                   (-Constants::electrons_per_particle *   // N
                    Constants::electron_charge));          // e
            }

        // Finally, we can copy the contributions of this cell into
        // the global matrix and right hand side vector:
        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
  }


  // @sect4{CathodeRaySimulator::solve}

  // The function that solves the linear system is then again exactly as in
  // step-6:
  template <int dim>
  void CathodeRaySimulator<dim>::solve_field()
  {
    SolverControl            solver_control(1000, 1e-12);
    SolverCG<Vector<double>> solver(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    solver.solve(system_matrix, solution, system_rhs, preconditioner);

    constraints.distribute(solution);
  }


  // @sect4{CathodeRaySimulator::refine_grid}

  // The final field-related function is the one that refines the grid. We will
  // call it a number of times in the first time step to obtain a mesh that
  // is well-adapted to the structure of the solution and, in particular,
  // resolves the various singularities in the solution that are due to
  // re-entrant corners and places where the boundary condition type
  // changes. You might want to refer to step-6 again for more details:
  template <int dim>
  void CathodeRaySimulator<dim>::refine_grid()
  {
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

    KellyErrorEstimator<dim>::estimate(dof_handler,
                                       QGauss<dim - 1>(fe.degree + 1),
                                       {},
                                       solution,
                                       estimated_error_per_cell);

    GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                    estimated_error_per_cell,
                                                    0.1,
                                                    0.03);

    triangulation.execute_coarsening_and_refinement();
  }


  // @sect4{CathodeRaySimulator::create_particles}

  // Let us now turn to the functions that deal with particles. The first one
  // is about the creation of particles. As mentioned in the introduction,
  // we want to create a particle at points of the cathode if the the electric
  // field $\mathbf E=\nabla V$ exceeds a certain threshold, i.e., if
  // $|\mathbf E| \ge E_\text{threshold}$, and if furthermore the electric field
  // points into the domain (i.e., if $\mathbf E \cdot \mathbf n < 0$). As is
  // common in the finite element method, we evaluate fields (and their
  // derivatives) at specific evaluation points; typically, these are
  // "quadrature points", and so we create a "quadrature formula" that we will
  // use to designate the points at which we want to evaluate the solution.
  // Here, we will simply take QMidpoint implying that we will only check the
  // threshold condition at the midpoints of faces. We then use this to
  // initialize an object of type FEFaceValues to evaluate the solution at these
  // points.
  //
  // All of this will then be used in a loop over all cells, their faces, and
  // specifically those faces that are at the boundary and, moreover, the
  // cathode part of the boundary.
  template <int dim>
  void CathodeRaySimulator<dim>::create_particles()
  {
    FEFaceValues<dim> fe_face_values(fe,
                                     QMidpoint<dim - 1>(),
                                     update_quadrature_points |
                                       update_gradients |
                                       update_normal_vectors);

    std::vector<Tensor<1, dim>> solution_gradients(
      fe_face_values.n_quadrature_points);

    for (const auto &cell : dof_handler.active_cell_iterators())
      for (const auto &face : cell->face_iterators())
        if (face->at_boundary() &&
            (face->boundary_id() == BoundaryIds::cathode))
          {
            fe_face_values.reinit(cell, face);

            // So we have found a face on the cathode. Next, we let the
            // FEFaceValues object compute the gradient of the solution at each
            // "quadrature" point, and extract the electric field vector from
            // the gradient in the form of a Tensor variable through the methods
            // discussed in the
            // @ref vector_valued "vector-valued problems" documentation module.
            const FEValuesExtractors::Scalar electric_potential(0);
            fe_face_values[electric_potential].get_function_gradients(
              solution, solution_gradients);
            for (const unsigned int q_point :
                 fe_face_values.quadrature_point_indices())
              {
                const Tensor<1, dim> E = solution_gradients[q_point];

                // Electrons can only escape the cathode if the electric field
                // strength exceeds a threshold and,
                // crucially, if the electric field points *into* the domain.
                // Once we have that checked, we create a new
                // Particles::Particle object at this location and insert it
                // into the Particles::ParticleHandler object with a unique ID.
                //
                // The only thing that may be not obvious here is that we also
                // associate with this particle the location in the reference
                // coordinates of the cell we are currently on. This is done
                // because we will in downstream functions compute quantities
                // such as the electric field at the location of the particle
                // (e.g., to compute the forces that act on it when updating its
                // position in each time step). Evaluating a finite element
                // field at arbitrary coordinates is quite an expensive
                // operation because shape functions are really only defined on
                // the reference cell, and so when asking for the electric field
                // at an arbitrary point requires us first to determine what the
                // reference coordinates of that point are. To avoid having to
                // do this over and over, we determine these coordinates once
                // and for all and then store these reference coordinates
                // directly with the particle.
                if ((E * fe_face_values.normal_vector(q_point) < 0) &&
                    (E.norm() > Constants::E_threshold))
                  {
                    const Point<dim> &location =
                      fe_face_values.quadrature_point(q_point);

                    Particles::Particle<dim> new_particle;
                    new_particle.set_location(location);
                    new_particle.set_reference_location(
                      mapping.transform_real_to_unit_cell(cell, location));
                    new_particle.set_id(next_unused_particle_id);
                    particle_handler.insert_particle(new_particle, cell);

                    ++next_unused_particle_id;
                  }
              }
          }

    // At the end of all of these insertions, we let the `particle_handler`
    // update some internal statistics about the particles it stores.
    particle_handler.update_cached_numbers();
  }


  // @sect4{CathodeRaySimulator::move_particles}

  // The second particle-related function is the one that moves the particles
  // in each time step. To do this, we have to loop over all cells, the
  // particles in each cell, and evaluate the electric field at each of the
  // particles' positions.
  //
  // The approach used here is conceptually the same used in the
  // `assemble_system()` function: We loop over all cells, find the particles
  // located there (with the same caveat about the inefficiency of the algorithm
  // used here to find these particles), and use FEPointEvaluation object to
  // evaluate the gradient at these positions:
  template <int dim>
  void CathodeRaySimulator<dim>::move_particles()
  {
    const double dt = time.get_next_step_size();

    Vector<double>            solution_values(fe.n_dofs_per_cell());
    FEPointEvaluation<1, dim> evaluator(mapping, fe, update_gradients);

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (particle_handler.n_particles_in_cell(cell) > 0)
        {
          const typename Particles::ParticleHandler<
            dim>::particle_iterator_range particles_in_cell =
            particle_handler.particles_in_cell(cell);

          std::vector<Point<dim>> particle_positions;
          for (const auto &particle : particles_in_cell)
            particle_positions.push_back(particle.get_reference_location());

          cell->get_dof_values(solution, solution_values);

          // Then we can ask the FEPointEvaluation object for the gradients of
          // the solution (i.e., the electric field $\mathbf E$) at these
          // locations and loop over the individual particles:
          evaluator.reinit(cell, particle_positions);
          evaluator.evaluate(make_array_view(solution_values),
                             EvaluationFlags::gradients);

          {
            typename Particles::ParticleHandler<dim>::particle_iterator
              particle = particles_in_cell.begin();
            for (unsigned int particle_index = 0;
                 particle != particles_in_cell.end();
                 ++particle, ++particle_index)
              {
                const Tensor<1, dim> &E =
                  evaluator.get_gradient(particle_index);

                // Having now obtained the electric field at the location of one
                // of the particles, we use this to update first the velocity
                // and then the position. To do so, let us first get the old
                // velocity out of the properties stored with the particle,
                // compute the acceleration, update the velocity, and store this
                // new velocity again in the properties of the particle. Recall
                // that this corresponds to the first of the following set of
                // update equations discussed in the introduction:
                // @f{align*}{
                //     \frac{{\mathbf v}_i^{(n)}
                //           -{\mathbf v}_i^{(n-1)}}{\Delta t}
                //     &= \frac{e\nabla V^{(n)}}{m}
                //  \\ \frac{{\mathbf x}_i^{(n)}-{\mathbf x}_i^{(n-1)}}
                //          {\Delta t} &= {\mathbf v}_i^{(n)}.
                // @f}
                const Tensor<1, dim> old_velocity(particle->get_properties());

                const Tensor<1, dim> acceleration =
                  Constants::electron_charge / Constants::electron_mass * E;

                const Tensor<1, dim> new_velocity =
                  old_velocity + acceleration * dt;

                particle->set_properties(make_array_view(new_velocity));

                // With the new velocity, we can then also update the location
                // of the particle and tell the particle about it.
                const Point<dim> new_location =
                  particle->get_location() + dt * new_velocity;
                particle->set_location(new_location);
              }
          }
        }

    // Having updated the locations and properties (i.e., velocities) of all
    // particles, we need to make sure that the `particle_handler` again knows
    // which cells they are in, and what their locations in the coordinate
    // system of the reference cell are. The following function does that. (It
    // also makes sure that, in parallel computations, particles are moved from
    // one processor to another processor if a particle moves from the subdomain
    // owned by the former to the subdomain owned by the latter.)
    particle_handler.sort_particles_into_subdomains_and_cells();
  }


  // @sect4{CathodeRaySimulator::track_lost_particle}

  // The final particle-related function is the one that is called whenever a
  // particle is lost from the simulation. This typically happens if it leaves
  // the domain. If that happens, this function is called both the cell (which
  // we can ask for its new location) and the cell it was previously on. The
  // function then keeps track of updating the number of particles lost in this
  // time step, the total number of lost particles, and then estimates whether
  // the particle left through the hole in the middle of the anode. We do so by
  // first checking whether the cell it was in last had an $x$ coordinate to the
  // left of the right boundary (located at $x=4$) and the particle now has a
  // position to the right of the right boundary. If that is so, we compute a
  // direction vector of its motion that is normalized so that the $x$ component
  // of the direction vector is equal to $1$. With this direction vector, we can
  // compute where it would have intersected the line $x=4$. If this intersect
  // is between $0.5$ and $1.5$, then we claim that the particle left through
  // the hole and increment a counter.
  template <int dim>
  void CathodeRaySimulator<dim>::track_lost_particle(
    const typename Particles::ParticleIterator<dim> &        particle,
    const typename Triangulation<dim>::active_cell_iterator &cell)
  {
    ++n_recently_lost_particles;
    ++n_total_lost_particles;

    const Point<dim> current_location              = particle->get_location();
    const Point<dim> approximate_previous_location = cell->center();

    if ((approximate_previous_location[0] < 4) && (current_location[0] > 4))
      {
        const Tensor<1, dim> direction =
          (current_location - approximate_previous_location) /
          (current_location[0] - approximate_previous_location[0]);

        const double right_boundary_intercept =
          approximate_previous_location[1] +
          (4 - approximate_previous_location[0]) * direction[1];
        if ((right_boundary_intercept > 0.5) &&
            (right_boundary_intercept < 1.5))
          ++n_particles_lost_through_anode;
      }
  }



  // @sect4{CathodeRaySimulator::update_timestep_size}

  // As discussed at length in the introduction, we need to respect a time step
  // condition whereby particles can not move further than one cell in one time
  // step. To ensure that this is the case, we again first compute the maximal
  // speed of all particles on each cell, and divide the cell size by that
  // speed. We then compute the next time step size as the minimum of this
  // quantity over all cells, using the safety factor discussed in the
  // introduction, and set this as the desired time step size using the
  // DiscreteTime::set_desired_time_step_size() function.
  template <int dim>
  void CathodeRaySimulator<dim>::update_timestep_size()
  {
    if (time.get_step_number() > 0)
      {
        double min_cell_size_over_velocity = std::numeric_limits<double>::max();

        for (const auto &cell : dof_handler.active_cell_iterators())
          if (particle_handler.n_particles_in_cell(cell) > 0)
            {
              const double cell_size = cell->minimum_vertex_distance();

              double max_particle_velocity(0.0);

              for (const auto &particle :
                   particle_handler.particles_in_cell(cell))
                {
                  const Tensor<1, dim> velocity(particle.get_properties());
                  max_particle_velocity =
                    std::max(max_particle_velocity, velocity.norm());
                }

              if (max_particle_velocity > 0)
                min_cell_size_over_velocity =
                  std::min(min_cell_size_over_velocity,
                           cell_size / max_particle_velocity);
            }

        constexpr double c_safety = 0.5;
        time.set_desired_next_step_size(c_safety * 0.5 *
                                        min_cell_size_over_velocity);
      }
    // As mentioned in the introduction, we have to treat the very first
    // time step differently since there, particles are not available yet or
    // do not yet have the information associated that we need for the
    // computation of a reasonable step length. The formulas below follow the
    // discussion in the introduction.
    else
      {
        const QTrapezoid<dim> vertex_quadrature;
        FEValues<dim> fe_values(fe, vertex_quadrature, update_gradients);

        std::vector<Tensor<1, dim>> field_gradients(vertex_quadrature.size());

        double min_timestep = std::numeric_limits<double>::max();

        for (const auto &cell : dof_handler.active_cell_iterators())
          if (particle_handler.n_particles_in_cell(cell) > 0)
            {
              const double cell_size = cell->minimum_vertex_distance();

              fe_values.reinit(cell);
              fe_values.get_function_gradients(solution, field_gradients);

              double max_E = 0;
              for (const auto q_point : fe_values.quadrature_point_indices())
                max_E = std::max(max_E, field_gradients[q_point].norm());

              if (max_E > 0)
                min_timestep =
                  std::min(min_timestep,
                           std::sqrt(0.5 * cell_size *
                                     Constants::electron_mass /
                                     Constants::electron_charge / max_E));
            }

        time.set_desired_next_step_size(min_timestep);
      }
  }



  // @sect4{The <code>CathodeRaySimulator::output_results()</code> function}

  // The final function implementing pieces of the overall algorithm is the one
  // that generates graphical output. In the current context, we want to output
  // both the electric potential field as well as the particle locations and
  // velocities. But we also want to output the electric field, i.e., the
  // gradient of the solution.
  //
  // deal.II has a general way how one can compute derived quantities from
  // the solution and output those as well. Here, this is the electric
  // field, but it could also be some other quantity -- say, the norm of the
  // electric field, or in fact anything else one could want to compute from
  // the solution $V_h(\mathbf x)$ or its derivatives. This general solution
  // uses the DataPostprocessor class and, in cases like the one here where we
  // want to output a quantity that represents a vector field, the
  // DataPostprocessorVector class.
  //
  // Rather than try and explain how this class works, let us simply refer to
  // the documentation of the DataPostprocessorVector class that has essentially
  // this case as a well-documented example.
  template <int dim>
  class ElectricFieldPostprocessor : public DataPostprocessorVector<dim>
  {
  public:
    ElectricFieldPostprocessor()
      : DataPostprocessorVector<dim>("electric_field", update_gradients)
    {}

    virtual void evaluate_scalar_field(
      const DataPostprocessorInputs::Scalar<dim> &input_data,
      std::vector<Vector<double>> &computed_quantities) const override
    {
      AssertDimension(input_data.solution_gradients.size(),
                      computed_quantities.size());

      for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
        {
          AssertDimension(computed_quantities[p].size(), dim);
          for (unsigned int d = 0; d < dim; ++d)
            computed_quantities[p][d] = input_data.solution_gradients[p][d];
        }
    }
  };



  // With this, the `output_results()` function becomes relatively
  // straightforward: We use the DataOut class as we have in almost every one of
  // the previous tutorial programs to output the solution (the "electric
  // potential") and we use the postprocessor defined above to also output its
  // gradient (the "electric field"). This all is then written into a file in
  // VTU format after also associating the current time and time step number
  // with this file.
  template <int dim>
  void CathodeRaySimulator<dim>::output_results() const
  {
    {
      ElectricFieldPostprocessor<dim> electric_field;
      DataOut<dim>                    data_out;
      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(solution, "electric_potential");
      data_out.add_data_vector(solution, electric_field);
      data_out.build_patches();

      data_out.set_flags(
        DataOutBase::VtkFlags(time.get_current_time(), time.get_step_number()));

      std::ofstream output("solution-" +
                           Utilities::int_to_string(time.get_step_number(), 4) +
                           ".vtu");
      data_out.write_vtu(output);
    }

    // Output the particle positions and properties is not more complicated. The
    // Particles::DataOut class plays the role of the DataOut class for
    // particles, and all we have to do is tell that class where to take
    // particles from and how to interpret the `dim` components of the
    // properties -- namely, as a single vector indicating the velocity, rather
    // than as `dim` scalar properties. The rest is then the same as above:
    {
      Particles::DataOut<dim, dim> particle_out;
      particle_out.build_patches(
        particle_handler,
        std::vector<std::string>(dim, "velocity"),
        std::vector<DataComponentInterpretation::DataComponentInterpretation>(
          dim, DataComponentInterpretation::component_is_part_of_vector));

      particle_out.set_flags(
        DataOutBase::VtkFlags(time.get_current_time(), time.get_step_number()));

      std::ofstream output("particles-" +
                           Utilities::int_to_string(time.get_step_number(), 4) +
                           ".vtu");
      particle_out.write_vtu(output);
    }
  }


  // @sect4{CathodeRaySimulator::run}

  // The last member function of the principal class of this program is then the
  // driver. At the top, it refines the mesh a number of times by solving the
  // problem (with not particles yet created) on a sequence of finer and finer
  // meshes.
  template <int dim>
  void CathodeRaySimulator<dim>::run()
  {
    make_grid();

    // do a few refinement cycles up front
    const unsigned int n_pre_refinement_cycles = 3;
    for (unsigned int refinement_cycle = 0;
         refinement_cycle < n_pre_refinement_cycles;
         ++refinement_cycle)
      {
        setup_system();
        assemble_system();
        solve_field();
        refine_grid();
      }


    // Now do the loop over time. The sequence of steps follows closely the
    // outline of the algorithm discussed in the introduction. As discussed in
    // great detail in the documentation of the DiscreteTime class, while we
    // move the field and particle information forward by one time step, the
    // time stored in the `time` variable is not consistent with where (some of)
    // these quantities are (in the diction of DiscreteTime, this is the "update
    // stage"). The call to `time.advance_time()` makes everything consistent
    // again by setting the `time` variable to the time at which the field and
    // particles already are, and once we are in this "consistent stage", we can
    // generate graphical output and write information about the current state
    // of the simulation to screen.
    setup_system();
    do
      {
        std::cout << "Timestep " << time.get_step_number() + 1 << std::endl;
        std::cout << "  Field degrees of freedom:                 "
                  << dof_handler.n_dofs() << std::endl;

        assemble_system();
        solve_field();

        create_particles();
        std::cout << "  Total number of particles in simulation:  "
                  << particle_handler.n_global_particles() << std::endl;

        n_recently_lost_particles = 0;
        update_timestep_size();
        move_particles();

        time.advance_time();

        output_results();

        std::cout << "  Number of particles lost this time step:  "
                  << n_recently_lost_particles << std::endl;
        if (n_total_lost_particles > 0)
          std::cout << "  Fraction of particles lost through anode: "
                    << 1. * n_particles_lost_through_anode /
                         n_total_lost_particles
                    << std::endl;

        std::cout << std::endl
                  << "  Now at t=" << time.get_current_time()
                  << ", dt=" << time.get_previous_step_size() << '.'
                  << std::endl
                  << std::endl;
      }
    while (time.is_at_end() == false);
  }
} // namespace Step83



// @sect3{The <code>main</code> function}

// The final function of the program is then again the `main()` function. It is
// unchanged in all tutorial programs since step-6 and so there is nothing new
// to discuss:
int main()
{
  try
    {
      Step83::CathodeRaySimulator<2> cathode_ray_simulator_2d;
      cathode_ray_simulator_2d.run();
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

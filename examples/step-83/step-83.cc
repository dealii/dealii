/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2023 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * Author: Pasquale Africa, SISSA, 2024,
 *         Wolfgang Bangerth, Colorado State University, 2024,
 *         Bruno Blais, Polytechnique Montreal, 2024.
 */


// @sect3{Include files}

// This program, with the exception of the checkpointing component
// is identical to step-19, and so the following include files are
// all the same:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/discrete_time.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/fe/mapping_q.h>
#include <deal.II/matrix_free/fe_point_evaluation.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/data_out.h>

// The only thing new are the following two include files. They are the ones
// that declare the classes we use as archives for reading (`iarchive` = input
// archive) and writing (`oarchive` = output archive) serialized data:
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <filesystem>
#include <fstream>
#include <string>

// @sect3{Global definitions}

// As is customary, we put everything that corresponds to the details of the
// program into a namespace of its own.
namespace Step83
{
  using namespace dealii;

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

  // The following is then the main class of this program. It is,
  // fundamentally, identical to step-19 with the exception of
  // the `checkpoint()` and `restart()` functions, along with the
  // `serialize()` function we use to serialize and deserialize the
  // data this class stores. The `serialize()` function is called
  // by the BOOST serialization framework, and consequently has to
  // have exactly the set of arguments used here. Furthermore, because
  // it is called by BOOST functions, it has to be `public`; the other
  // two new functions are as always made `private`.
  //
  // The `run()` function has also been modified to enable simulation restart
  // via its new argument `do_restart` that indicates whether or not to
  // start the simulation from a checkpoint.
  template <int dim>
  class CathodeRaySimulator
  {
  public:
    CathodeRaySimulator();

    void run(const bool do_restart);

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version);

  private:
    void make_grid();
    void setup_system();
    void assemble_system();
    void solve_field();
    void refine_grid();

    void create_particles();
    void move_particles();
    void track_lost_particle(
      const Particles::ParticleIterator<dim>                  &particle,
      const typename Triangulation<dim>::active_cell_iterator &cell);

    void update_timestep_size();
    void output_results() const;

    void checkpoint();
    void restart();

    Triangulation<dim>        triangulation;
    const MappingQ<dim>       mapping;
    const FE_Q<dim>           fe;
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

  // @sect4{The unchanged parts of the class}

  // Let us start with those parts of the class that are all unchanged
  // from step-19 and about which you can learn there. We will
  // then pick up with commentary again when we get to the two new
  // functions, `checkpoint()` and `restart()`, along with how the
  // `run()` function needs to be modified:
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
      [this](const typename Particles::ParticleIterator<dim>         &particle,
             const typename Triangulation<dim>::active_cell_iterator &cell) {
        this->track_lost_particle(particle, cell);
      });
  }



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

    std::vector<CellData<dim>> cells((nx - 1) * (ny - 1), CellData<dim>());
    for (unsigned int i = 0; i < cells.size(); ++i)
      {
        cells[i].vertices    = cell_vertices[i];
        cells[i].material_id = 0;
      }

    GridTools::consistently_order_cells(cells);
    triangulation.create_triangulation(
      vertices,
      cells,
      SubCellData()); // No boundary information

    triangulation.refine_global(2);

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
                                    /* keep_constrained_dofs = */ false);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
  }



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

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
  }



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

            const FEValuesExtractors::Scalar electric_potential(0);
            fe_face_values[electric_potential].get_function_gradients(
              solution, solution_gradients);
            for (const unsigned int q_point :
                 fe_face_values.quadrature_point_indices())
              {
                const Tensor<1, dim> E = solution_gradients[q_point];

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

    particle_handler.update_cached_numbers();
  }



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

                const Tensor<1, dim> old_velocity(particle->get_properties());

                const Tensor<1, dim> acceleration =
                  Constants::electron_charge / Constants::electron_mass * E;

                const Tensor<1, dim> new_velocity =
                  old_velocity + acceleration * dt;

                particle->set_properties(new_velocity);

                const Point<dim> new_location =
                  particle->get_location() + dt * new_velocity;
                particle->set_location(new_location);
              }
          }
        }

    particle_handler.sort_particles_into_subdomains_and_cells();
  }



  template <int dim>
  void CathodeRaySimulator<dim>::track_lost_particle(
    const typename Particles::ParticleIterator<dim>         &particle,
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

      DataOutBase::VtkFlags output_flags;
      output_flags.time  = time.get_current_time();
      output_flags.cycle = time.get_step_number();
      output_flags.physical_units["electric_potential"] = "V";
      output_flags.physical_units["electric_field"]     = "V/m";

      data_out.set_flags(output_flags);

      std::ofstream output("solution-" +
                           Utilities::int_to_string(time.get_step_number(), 4) +
                           ".vtu");
      data_out.write_vtu(output);
    }

    {
      Particles::DataOut<dim> particle_out;
      particle_out.build_patches(
        particle_handler,
        std::vector<std::string>(dim, "velocity"),
        std::vector<DataComponentInterpretation::DataComponentInterpretation>(
          dim, DataComponentInterpretation::component_is_part_of_vector));

      DataOutBase::VtkFlags output_flags;
      output_flags.time                       = time.get_current_time();
      output_flags.cycle                      = time.get_step_number();
      output_flags.physical_units["velocity"] = "m/s";

      particle_out.set_flags(output_flags);

      std::ofstream output("particles-" +
                           Utilities::int_to_string(time.get_step_number(), 4) +
                           ".vtu");
      particle_out.write_vtu(output);
    }
  }



  // @sect4{CathodeRaySimulator::serialize()}

  // The first of the new function is the one that is called by the
  // BOOST Serialization framework to serialize and deserialize the
  // data of this class. It has already been discussed in the introduction
  // to this program and so does not provide any surprises. All it does is
  // write those member variables of the current class that cannot be
  // re-created easily into an archive, or read these members from it.
  // (Whether `operator &` facilitates a write or read operation depends on
  // whether the `Archive` type is an output or input archive.)
  //
  // The function takes a second argument, `version`, that can be used to
  // create checkpoints that have version numbers. This is useful if one
  // evolves programs by adding more member variables, but still wants to
  // retain the ability to read checkpoint files created with earlier
  // versions of the program. The `version` variable would, in that case,
  // be used to represent which version of the program wrote the file,
  // and if necessary to read only those variables that were written with
  // that past version, finding a different way to initialize the new member
  // variables that have been added since then. We will not make use of this
  // ability here.
  //
  // Finally, while the program that indents all deal.II source files format
  // the following code as
  // @code
  //   ar &solution;
  // @endcode
  // as if we are taking the address of the `triangulation` variable, the
  // way you *should* read the code is as
  // @code
  //   ar & solution;
  // @endcode
  // where `operator &` is a binary operator that could either be interpreted
  // as `operator <<` for output or `operator >>` for input.
  //
  // As discussed in the introduction, we do not serialize the `triangulation`
  // member variable, instead leaving that to separate calls in the
  // `checkpoint()` and `restart()` functions below.
  template <int dim>
  template <class Archive>
  void CathodeRaySimulator<dim>::serialize(Archive &ar,
                                           const unsigned int /* version */)
  {
    ar &solution;
    ar &particle_handler;
    ar &next_unused_particle_id;
    ar &n_recently_lost_particles;
    ar &n_total_lost_particles;
    ar &n_particles_lost_through_anode;
    ar &time;
  }



  // @sect4{CathodeRaySimulator::checkpoint()}

  // The checkpoint function of the principal class of this program is then
  // quite straightforward: We create an output file (and check that it is
  // writable), create an output archive, and then move the serialized
  // contents of the current object (i.e., the `*this` object) into the
  // archive. The use of `operator<<` here calls the `serialize()` function
  // above with an output archive as argument. When the destructor of the
  // `archive` variable is called at the end of the code block within which
  // it lives, the entire archive is written into the output file stream it
  // is associated with.
  //
  // As mentioned in the introduction, "real" applications would not use
  // text-based archives as provided by the `boost::archive::text_oarchive`
  // class, but use binary and potentially compressed versions. This can
  // easily be achieved by using differently named classes, and the BOOST
  // documentation explains how to do that.
  template <int dim>
  void CathodeRaySimulator<dim>::checkpoint()
  {
    std::cout << "--- Writing checkpoint... ---" << std::endl << std::endl;

    {
      std::ofstream checkpoint_file("tmp.checkpoint_step_83");
      AssertThrow(checkpoint_file,
                  ExcMessage(
                    "Could not write to the <tmp.checkpoint_step_83> file."));

      boost::archive::text_oarchive archive(checkpoint_file);

      archive << *this;
    }

    // The second part of the serialization is all of the data that we can
    // attach to cells -- see the discussion about this in the introduction.
    // Here, the only data we attach to cells are the particles. We then
    // let the triangulation save these into a set of files that all start
    // with the same prefix as we chose above, namely "tmp.checkpoint":
    particle_handler.prepare_for_serialization();
    triangulation.save("tmp.checkpoint");


    // At this point, the serialized data of this file has ended up in a number
    // of files that all start with `tmp.checkpoint` file. As mentioned in the
    // introduction, we do not want to directly overwrite the checkpointing
    // files from the previous checkpoint operation, for fear that the program
    // may be interrupted *while writing the checkpoint files*. This would
    // result in corrupted files, and defeat the whole purpose of checkpointing
    // because one cannot restart from such a file. On the other hand, if we got
    // here, we know that the "tmp.checkpoint*" files were successfully written,
    // and we can rename it to "checkpoint*", in the process replacing the old
    // file.
    //
    // We do this move operation by calling the
    // [C++ function that does the renaming of
    // files](https://en.cppreference.com/w/cpp/filesystem/rename).
    // Note that it is documented as being for all practical purposes
    // "atomic", i.e., we do not need to worry that the program may be
    // interrupted somewhere within the renaming operation itself. Of
    // course, it is possible that we get interrupted somewhere between
    // renaming one file and the next, and that presents problems in
    // itself -- in essence, we want the entire renaming operation of all of
    // these files to be atomic. With a couple dozen lines of extra code, one
    // could address this issue (using strategies that databases use frequently:
    // if one operation fails, we need to
    // [rollback](https://en.wikipedia.org/wiki/Rollback_(data_management))
    // the entire transaction). For the purposes of this program, this is
    // perhaps too much, and we will simply hope that that doesn't happen,
    // perhaps based on the belief that renaming files is much faster than
    // writing them, and that unlike writing checkpoint files, renaming does not
    // require much memory or disk space and so does not risk running out of
    // either.
    //
    // As a consequence, the following code first loops over all files in
    // the current directory, picks out those that start with the string
    // "tmp.checkpoint", and puts them into a list. In a second step,
    // we loop over the list and rename each of these files to one
    // whose name consists of the "tmp.checkpoint*" file but stripped off
    // its first four characters (i.e., only the "checkpoint*" part). We use
    // this approach, rather than listing the files we want to rename,
    // because we do not actually know the names of the files written by
    // the Triangulation::save() function, though we know how each of these
    // file names starts.
    std::list<std::string> tmp_checkpoint_files;
    for (const auto &dir_entry : std::filesystem::directory_iterator("."))
      if (dir_entry.is_regular_file() &&
          (dir_entry.path().filename().string().find("tmp.checkpoint") == 0))
        tmp_checkpoint_files.push_back(dir_entry.path().filename().string());

    for (const std::string &filename : tmp_checkpoint_files)
      std::filesystem::rename(filename, filename.substr(4, std::string::npos));
  }


  // @sect4{CathodeRaySimulator::restart()}

  // The restart function of this class then simply does the opposite:
  // It opens an input file (and triggers an error if that file cannot be
  // opened), associates an input archive with it, and then reads the
  // contents of the current object from it, again using the
  // `serialize()` function from above. Clearly, since we have written
  // data into a text-based archive above, we need to use the corresponding
  // `boost::archive::text_iarchive` class for reading.
  //
  // In a second step, we ask the triangulation to read in cell-attached
  // data, and then tell the Particles::ParticleHandler object to re-create
  // its information about all of the particles from the data just read.
  //
  // The function ends by printing a status message about having restarted:
  template <int dim>
  void CathodeRaySimulator<dim>::restart()
  {
    {
      std::ifstream checkpoint_file("checkpoint_step_83");
      AssertThrow(checkpoint_file,
                  ExcMessage(
                    "Could not read from the <checkpoint_step_83> file."));

      boost::archive::text_iarchive archive(checkpoint_file);
      archive >> *this;
    }

    triangulation.load("checkpoint");
    particle_handler.deserialize();

    std::cout << "--- Restarting at t=" << time.get_current_time()
              << ", dt=" << time.get_next_step_size() << std::endl
              << std::endl;
  }


  // @sect4{CathodeRaySimulator::run()}

  // The last member function of the principal class of this program is then the
  // driver. The driver takes a single argument to indicate if the simulation
  // is a restart. If it is not a restart, the mesh is set up and the problem is
  // solved like in step-19. If it is a restart, then we read in everything that
  // is a history variable from the checkpoint file via the `restart()`
  // function. Recall that everything that is inside the `if` block at the top
  // of the function is exactly like in step-19, as is almost everything that
  // follows:
  template <int dim>
  void CathodeRaySimulator<dim>::run(const bool do_restart)
  {
    if (do_restart == false)
      {
        make_grid();

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
      }
    else
      {
        restart();
      }

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

        // The only other difference between this program and step-19 is that
        // we checkpoint the simulation every ten time steps:
        if (time.get_step_number() % 10 == 0)
          checkpoint();
      }
    while (time.is_at_end() == false);
  }
} // namespace Step83


// @sect3{The <code>main</code> function}

// The final function of the program is then again the `main()` function. Its
// overall structure is unchanged in all tutorial programs since step-6 and
// so there is nothing new to discuss about this aspect.
//
// The only difference is that we need to figure out whether a restart was
// requested, or whether the program should simply start from scratch when
// called. We do this using a command line argument: The `argc` argument to
// `main()` indicates how many command line arguments were provided when
// the program was called (counting the name of the program as the zeroth
// argument), and `argv` is an array of strings with as many elements as
// `argc` that contains these command line arguments. So if you call
// the program as
// @code
//   ./step-83
// @endcode
// then `argc` will be 1, and `argv` will be the array with one element
// and content `[ "./step-83" ]`. On the other hand, if you call the program
// as
// @code
//   ./step-83 restart
// @endcode
// then `argc` will be 2, and `argv` will be the array with two elements
// and content `[ "./step-83", "restart" ]`. Every other choice should be
// flagged as an error. The following try block then does exactly this:
int main(int argc, char *argv[])
{
  try
    {
      Step83::CathodeRaySimulator<2> cathode_ray_simulator;

      if (argc == 1)
        cathode_ray_simulator.run(false); // no restart
      else if ((argc == 2) && (std::string(argv[1]) == "restart"))
        cathode_ray_simulator.run(true); // yes restart
      else
        {
          std::cerr << "Error: The only allowed command line argument to this\n"
                       "       program is 'restart'."
                    << std::endl;
          return 1;
        }
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
